/* ==========================================================================
 * GRAMSCI OpenCL compute kernels  (3PCF, equilateral-3PCF, 4PCF, 4PCF-parity)
 *
 * Portable single-precision port of the OpenACC kernels in ../src_gpu.  Apple's
 * OpenCL on Apple Silicon has no fp64 (CL_DEVICE_DOUBLE_FP_CONFIG == 0) and no
 * float atomics, so the design differs from the nvfortran build in two ways:
 *
 *   1. All device arithmetic is float.  Host reduction of the partials is done
 *      in double, and the parity chirality sign is looked up from a host-built
 *      table (computed in double) so it is bit-exact with the CPU reference.
 *      Accumulation uses Kahan compensated summation (kadd below, one comp
 *      array per partial array): a plain fp32 += stops registering increments
 *      once a column partial passes ~2/eps times the increment (~1.7e7 tuples
 *      per column for unit weights), which silently and one-sidedly undercounts
 *      the monotone RRR/RRRR channels on production-size runs.  With the
 *      compensation term the error stays O(eps) independent of the count.
 *
 *   2. No atomics.  Each work-item IS one "gang": it owns column g of the
 *      partial accumulators and processes a strided set of hubs
 *      (i = istart+g, istart+g+ngang, ...).  Because no two work-items ever
 *      touch the same accumulator slot, the histogram needs no atomics; the
 *      host sums the ngang columns afterwards (in double).
 *
 * CSR layout (1-based offsets, matching the Fortran arrays byte-for-byte):
 *   ptr[i-1] .. ptr[i]-1   are the (1-based) edge indices of node i
 *   id[e-1], dist[e-1]     neighbor id and distance-bin of edge e (1..tot)
 * Distance bins are 1..nbins; bin 0 means "no such edge".
 * ========================================================================== */

/* Kahan compensated add into sum[idx] (comp[idx] carries the rounding error).
 * OpenCL C defaults to strict IEEE float add (no reassociation), so the
 * compensation cannot be optimised away.  The true total is sum - comp; the
 * host subtracts the comp array in its double reduction. */
inline void kadd(__global float *sum, __global float *comp, long idx, float val)
{
    float y = val - comp[idx];
    float t = sum[idx] + y;
    comp[idx] = (t - sum[idx]) - y;
    sum[idx]  = t;
}

/* CAS-emulated float atomic add (Apple OpenCL 1.2 has integer atomics only,
 * atomic_cmpxchg on 32-bit ints being core since 1.1).  Used ONLY for the
 * jackknife touching sums: their fp32 precision is protected by the host,
 * which commits and re-zeros the buffers after every bucketed window
 * (cl_run_bucketed), so a single window's per-slot sum stays far below the
 * ~2/eps saturation that motivates kadd for the full-run channels.  With
 * -njk the query always runs bucketed for exactly this reason. */
inline void atomic_fadd(volatile __global float *p, float v)
{
    if (v == 0.0f) return;
    volatile __global uint *up = (volatile __global uint *)p;
    uint old = *up, assumed;
    do {
        assumed = old;
        old = atomic_cmpxchg(up, assumed, as_uint(as_float(assumed) + v));
    } while (old != assumed);
}

/* Binary search node `from`'s sorted adjacency row for neighbor `to`.
 * Returns the distance-bin (1..nbins) or 0 if the edge is absent.
 * `long` indexing throughout: edge counts may exceed 2^31. */
inline char find_bin(__global const long *ptr,
                     __global const int  *id,
                     __global const char *dist,
                     int from, int to)
{
    long s  = ptr[from - 1];          /* 1-based first edge index of `from` */
    int  nn = (int)(ptr[from] - s);
    if (nn == 0) return 0;
    long b0 = s - 1;                  /* 0-based offset of first edge */
    if (id[b0] > to || id[b0 + nn - 1] < to) return 0;  /* outside sorted range */

    int lo = 0, hi = nn - 1;
    while (lo <= hi) {
        int mid = (lo + hi) >> 1;
        int v   = id[b0 + mid];
        if (v == to)      return dist[b0 + mid];
        else if (v < to)  lo = mid + 1;
        else              hi = mid - 1;
    }
    return 0;
}

/* --------------------------------------------------------------------------
 * 3PCF, all triangle configurations (binary-search over the CSR half-graph).
 * part_nnn[g*cb + bin-1] gets +w3; part_rrr likewise gets -w3 for RRR triples
 * (matching the OpenACC sign convention: N3 += sum(part_rrr)).
 * -------------------------------------------------------------------------- */
__kernel void k_3pcf_all(__global const long  *ptr,
                         __global const int   *id,
                         __global const char  *dist,
                         __global const float *w,
                         __global const int   *bt3,    /* nbins^3 config table */
                         __global float       *part_nnn,
                         __global float       *part_rrr,
                         const int nbins, const int cb,
                         const int num_data, const int istart,
                         const int iend, const int ngang,
                         const int nbuckets, const int blo, const int bhi,
                         __global char *done_flag,
                         __global float *comp_nnn,
                         __global float *comp_rrr,
                         __global const int   *region, /* jackknife labels, 1..njk */
                         __global float       *jkn,    /* touching sums, [cb x njk] */
                         __global float       *jkr,
                         const int njk)
{
    int g = (int)get_global_id(0);
    long gcol = (long)g * cb;

    for (int i = istart + g; i <= iend; i += ngang) {
        /* interleaved hub batching: process only hubs in the current bucket
         * window [blo,bhi).  Buckets are i mod nbuckets, so every launch's hubs
         * are spread across the whole catalog -> balanced load, full occupancy
         * (contiguous slices would isolate the expensive clustered hubs). */
        if (bhi - blo < nbuckets) {   /* partial window: skip out-of-window hubs */
            int bk = (i - istart) % nbuckets;
            if (bk < blo || bk >= bhi) continue;
        }
        long s  = ptr[i - 1];
        int  nn = (int)(ptr[i] - s);
        if (nn < 2) continue;
        long b0 = s - 1;
        int jr1 = (njk > 0) ? region[i - 1] : 0;

        for (int k1 = 1; k1 <= nn; k1++) {
            char c1 = dist[b0 + k1 - 1];
            int  n1 = id  [b0 + k1 - 1];
            float wi_w1 = w[i - 1] * w[n1 - 1];
            int jr2 = (njk > 0) ? region[n1 - 1] : 0;

            for (int k2 = k1 + 1; k2 <= nn; k2++) {
                char c2 = dist[b0 + k2 - 1];
                int  n2 = id  [b0 + k2 - 1];

                char c3 = find_bin(ptr, id, dist, n1, n2);
                if (c3 == 0) continue;

                int bin = bt3[(c1 - 1) + nbins * ((c2 - 1) + nbins * (c3 - 1))];
                float w3 = wi_w1 * w[n2 - 1];
                int rrr = (i > num_data && n1 > num_data && n2 > num_data);

                kadd(part_nnn, comp_nnn, gcol + bin - 1, w3);
                if (rrr)
                    kadd(part_rrr, comp_rrr, gcol + bin - 1, -w3);

                /* Jackknife: count this triplet once per DISTINCT region it
                 * touches (unrolled dedup, matching the CPU/OpenACC kernels).
                 * RRR carries the minus sign, like part_rrr above. */
                if (njk > 0) {
                    int jr3 = region[n2 - 1];
                    if (jr1 > 0) {
                        atomic_fadd(&jkn[(bin - 1) + (long)cb * (jr1 - 1)], w3);
                        if (rrr) atomic_fadd(&jkr[(bin - 1) + (long)cb * (jr1 - 1)], -w3);
                    }
                    if (jr2 > 0 && jr2 != jr1) {
                        atomic_fadd(&jkn[(bin - 1) + (long)cb * (jr2 - 1)], w3);
                        if (rrr) atomic_fadd(&jkr[(bin - 1) + (long)cb * (jr2 - 1)], -w3);
                    }
                    if (jr3 > 0 && jr3 != jr1 && jr3 != jr2) {
                        atomic_fadd(&jkn[(bin - 1) + (long)cb * (jr3 - 1)], w3);
                        if (rrr) atomic_fadd(&jkr[(bin - 1) + (long)cb * (jr3 - 1)], -w3);
                    }
                }
            }
        }
    }
    done_flag[g] = 1;   /* work-item completed its full hub sweep (watchdog check) */
}

/* --------------------------------------------------------------------------
 * Equilateral 3PCF: only triangles whose three sides share one radial bin.
 * Accumulate per radial bin (cb == nbins here, bin index = c1-1).
 * -------------------------------------------------------------------------- */
__kernel void k_3pcf_equi(__global const long  *ptr,
                          __global const int   *id,
                          __global const char  *dist,
                          __global const float *w,
                          __global float       *part_nnn,
                          __global float       *part_rrr,
                          const int nbins, const int cb,
                          const int num_data, const int istart,
                          const int iend, const int ngang,
                         const int nbuckets, const int blo, const int bhi,
                         __global char *done_flag,
                         __global float *comp_nnn,
                         __global float *comp_rrr)
{
    int g = (int)get_global_id(0);
    long gcol = (long)g * cb;

    for (int i = istart + g; i <= iend; i += ngang) {
        /* interleaved hub batching: process only hubs in the current bucket
         * window [blo,bhi).  Buckets are i mod nbuckets, so every launch's hubs
         * are spread across the whole catalog -> balanced load, full occupancy
         * (contiguous slices would isolate the expensive clustered hubs). */
        if (bhi - blo < nbuckets) {   /* partial window: skip out-of-window hubs */
            int bk = (i - istart) % nbuckets;
            if (bk < blo || bk >= bhi) continue;
        }
        long s  = ptr[i - 1];
        int  nn = (int)(ptr[i] - s);
        if (nn < 2) continue;
        long b0 = s - 1;

        for (int k1 = 1; k1 <= nn; k1++) {
            char c1 = dist[b0 + k1 - 1];
            int  n1 = id  [b0 + k1 - 1];
            float wi_w1 = w[i - 1] * w[n1 - 1];

            for (int k2 = k1 + 1; k2 <= nn; k2++) {
                char c2 = dist[b0 + k2 - 1];
                if (c2 != c1) continue;            /* equilateral filter */
                int  n2 = id[b0 + k2 - 1];
                char c3 = find_bin(ptr, id, dist, n1, n2);
                if (c3 != c1) continue;

                int bin = (int)c1;                 /* 1-based radial bin */
                float w3 = wi_w1 * w[n2 - 1];

                kadd(part_nnn, comp_nnn, gcol + bin - 1, w3);
                if (i > num_data && n1 > num_data && n2 > num_data)
                    kadd(part_rrr, comp_rrr, gcol + bin - 1, -w3);
            }
        }
    }
    done_flag[g] = 1;   /* work-item completed its full hub sweep (watchdog check) */
}

/* --------------------------------------------------------------------------
 * 4PCF, all configurations — pure binary-search strategy (matches the CPU
 * query_graph_4pcf_bsearch).  For each k1<k2<k3 the three closing edges are
 * found with find_bin (id1-id2, id1-id3, id2-id3).
 *
 * NOTE: unlike the OpenACC build, this does NOT precompute a per-hub
 * connectivity matrix (lmat).  The lmat optimisation needs a per-work-item
 * read-write GLOBAL scratch buffer; the bsearch kernel uses only read-only
 * inputs + a private accumulator column (structurally identical to the
 * reliable 3PCF kernel), so it has no scratch to mismanage and is
 * deterministic.  It does ~25x more searches per hub — a fine trade for
 * simplicity/correctness on a laptop GPU.  The host drives this kernel in
 * short hub-range batches so no single launch trips the GPU watchdog timer
 * (which would otherwise return partial results and silently drop counts).
 * -------------------------------------------------------------------------- */
__kernel void k_4pcf_all(__global const long  *ptr,
                         __global const int   *id,
                         __global const char  *dist,
                         __global const float *w,
                         __global const int   *bt6,    /* nbins^6 config table */
                         __global float       *part_n4,
                         __global float       *part_r4,
                         const int nbins, const int ncfg,
                         const int num_data,
                         const int istart, const int iend, const int ngang,
                         const int nbuckets, const int blo, const int bhi,
                         __global char *done_flag,
                         __global float *comp_n4,
                         __global float *comp_r4,
                         __global const int   *region, /* jackknife labels, 1..njk */
                         __global float       *jkn,    /* touching sums, [ncfg x njk] */
                         __global float       *jkr,
                         const int njk)
{
    int  g    = (int)get_global_id(0);
    long gcol = (long)g * ncfg;
    long n    = nbins;

    for (int i = istart + g; i <= iend; i += ngang) {
        /* interleaved hub batching: process only hubs in the current bucket
         * window [blo,bhi).  Buckets are i mod nbuckets, so every launch's hubs
         * are spread across the whole catalog -> balanced load, full occupancy
         * (contiguous slices would isolate the expensive clustered hubs). */
        if (bhi - blo < nbuckets) {   /* partial window: skip out-of-window hubs */
            int bk = (i - istart) % nbuckets;
            if (bk < blo || bk >= bhi) continue;
        }
        long s  = ptr[i - 1];
        int  nn = (int)(ptr[i] - s);
        if (nn <= 2) continue;
        long b0 = s - 1;
        int jr1 = (njk > 0) ? region[i - 1] : 0;

        for (int k1 = 1; k1 <= nn; k1++) {
            char c1 = dist[b0 + k1 - 1];
            int  n1 = id  [b0 + k1 - 1];
            int jr2 = (njk > 0) ? region[n1 - 1] : 0;

            for (int k2 = k1 + 1; k2 <= nn; k2++) {
                int  n2 = id[b0 + k2 - 1];
                char c4 = find_bin(ptr, id, dist, n1, n2);
                if (c4 == 0) continue;
                char c2 = dist[b0 + k2 - 1];
                float w12 = w[i - 1] * w[n1 - 1] * w[n2 - 1];
                int rand12 = (i > num_data && n1 > num_data && n2 > num_data);
                int jr3 = (njk > 0) ? region[n2 - 1] : 0;

                for (int k3 = k2 + 1; k3 <= nn; k3++) {
                    int  n3 = id[b0 + k3 - 1];
                    char c5 = find_bin(ptr, id, dist, n1, n3);
                    if (c5 == 0) continue;
                    char c6 = find_bin(ptr, id, dist, n2, n3);
                    if (c6 == 0) continue;
                    char c3 = dist[b0 + k3 - 1];

                    long off = (c1 - 1) + n * ((c2 - 1) + n * ((c3 - 1) +
                               n * ((c4 - 1) + n * ((c5 - 1) + n * (c6 - 1)))));
                    int cfgidx = abs(bt6[off]);
                    if (cfgidx == 0) continue;

                    float w4 = w12 * w[n3 - 1];
                    int isr = (rand12 && n3 > num_data);
                    kadd(part_n4, comp_n4, gcol + cfgidx - 1, w4);
                    if (isr)
                        kadd(part_r4, comp_r4, gcol + cfgidx - 1, w4);

                    /* Jackknife: once per DISTINCT region touched (RRRR is
                     * accumulated with PLUS, matching the CPU convention). */
                    if (njk > 0) {
                        int jr4 = region[n3 - 1];
                        long jb = (long)(cfgidx - 1);
                        if (jr1 > 0) {
                            atomic_fadd(&jkn[jb + (long)ncfg * (jr1 - 1)], w4);
                            if (isr) atomic_fadd(&jkr[jb + (long)ncfg * (jr1 - 1)], w4);
                        }
                        if (jr2 > 0 && jr2 != jr1) {
                            atomic_fadd(&jkn[jb + (long)ncfg * (jr2 - 1)], w4);
                            if (isr) atomic_fadd(&jkr[jb + (long)ncfg * (jr2 - 1)], w4);
                        }
                        if (jr3 > 0 && jr3 != jr1 && jr3 != jr2) {
                            atomic_fadd(&jkn[jb + (long)ncfg * (jr3 - 1)], w4);
                            if (isr) atomic_fadd(&jkr[jb + (long)ncfg * (jr3 - 1)], w4);
                        }
                        if (jr4 > 0 && jr4 != jr1 && jr4 != jr2 && jr4 != jr3) {
                            atomic_fadd(&jkn[jb + (long)ncfg * (jr4 - 1)], w4);
                            if (isr) atomic_fadd(&jkr[jb + (long)ncfg * (jr4 - 1)], w4);
                        }
                    }
                }
            }
        }
    }
    done_flag[g] = 1;   /* work-item completed its full hub sweep (watchdog check) */
}

/* --------------------------------------------------------------------------
 * 4PCF with parity decomposition (pure binary-search; see k_4pcf_all note).
 * Even channel  += w4;  odd channel += (parity_flip * sign_V) * w4.
 * parity_flip = sign(bt6 entry) * chiral[cfg]; sign_V is looked up from the
 * host-built signv table (computed in double, so the chirality matches the
 * CPU exactly).  chiral[cfg] (config_module chiral_4pcf) is 0 for binned
 * configurations invariant under an odd vertex permutation: their odd 4PCF is
 * identically zero, and any sign the estimator gave them would depend on the
 * vertex labelling, i.e. on catalogue row order.
 *   signv[(p1-1)*ndir*ndir + (p2-1)*ndir + (p3-1)] in {-1,0,+1}
 * -------------------------------------------------------------------------- */
__kernel void k_4pcf_parity(__global const long  *ptr,
                            __global const int   *id,
                            __global const char  *dist,
                            __global const short *phi,   /* int16: pixel index reaches ntheta*nphi=256+ */
                            __global const float *w,
                            __global const int   *bt6,
                            __global const char  *signv,  /* ndir^3 sign table */
                            __global float       *pn_even,
                            __global float       *pn_odd,
                            __global float       *pr_even,
                            __global float       *pr_odd,
                            const int nbins, const int ncfg,
                            const int num_data, const int ndir,
                            const int istart, const int iend, const int ngang,
                         const int nbuckets, const int blo, const int bhi,
                         __global char *done_flag,
                         __global float *cn_even,
                         __global float *cn_odd,
                         __global float *cr_even,
                         __global float *cr_odd,
                         __global const int   *region, /* jackknife labels, 1..njk */
                         __global float       *jkn,    /* [ncfg x 2ch x njk], Fortran N4jk layout */
                         __global float       *jkr,
                         const int njk,
                         __global const char  *chiral) /* [ncfg]: 1 chiral, 0 achiral */
{
    int  g    = (int)get_global_id(0);
    long gcol = (long)g * ncfg;
    long n    = nbins;

    for (int i = istart + g; i <= iend; i += ngang) {
        /* interleaved hub batching: process only hubs in the current bucket
         * window [blo,bhi).  Buckets are i mod nbuckets, so every launch's hubs
         * are spread across the whole catalog -> balanced load, full occupancy
         * (contiguous slices would isolate the expensive clustered hubs). */
        if (bhi - blo < nbuckets) {   /* partial window: skip out-of-window hubs */
            int bk = (i - istart) % nbuckets;
            if (bk < blo || bk >= bhi) continue;
        }
        long s  = ptr[i - 1];
        int  nn = (int)(ptr[i] - s);
        if (nn <= 2) continue;
        long b0 = s - 1;
        int jr1 = (njk > 0) ? region[i - 1] : 0;

        for (int k1 = 1; k1 <= nn; k1++) {
            char c1 = dist[b0 + k1 - 1];
            int  n1 = id  [b0 + k1 - 1];
            int  p1 = (int)phi[b0 + k1 - 1];
            int jr2 = (njk > 0) ? region[n1 - 1] : 0;

            for (int k2 = k1 + 1; k2 <= nn; k2++) {
                int  n2 = id[b0 + k2 - 1];
                char c4 = find_bin(ptr, id, dist, n1, n2);
                if (c4 == 0) continue;
                char c2 = dist[b0 + k2 - 1];
                int  p2 = (int)phi[b0 + k2 - 1];
                float w12 = w[i - 1] * w[n1 - 1] * w[n2 - 1];
                int rand12 = (i > num_data && n1 > num_data && n2 > num_data);
                int jr3 = (njk > 0) ? region[n2 - 1] : 0;

                for (int k3 = k2 + 1; k3 <= nn; k3++) {
                    int  n3 = id[b0 + k3 - 1];
                    char c5 = find_bin(ptr, id, dist, n1, n3);
                    if (c5 == 0) continue;
                    char c6 = find_bin(ptr, id, dist, n2, n3);
                    if (c6 == 0) continue;
                    char c3 = dist[b0 + k3 - 1];

                    long off = (c1 - 1) + n * ((c2 - 1) + n * ((c3 - 1) +
                               n * ((c4 - 1) + n * ((c5 - 1) + n * (c6 - 1)))));
                    int raw = bt6[off];
                    int cfgidx = abs(raw);
                    if (cfgidx == 0) continue;
                    int pflip = ((raw > 0) ? 1 : -1) * (int)chiral[cfgidx - 1];

                    int p3 = (int)phi[b0 + k3 - 1];
                    int sgn = (int)signv[(long)(p1 - 1) * ndir * ndir
                                         + (p2 - 1) * ndir + (p3 - 1)];

                    float w4 = w12 * w[n3 - 1];
                    float w4o = (float)(pflip * sgn) * w4;
                    int isr = (rand12 && n3 > num_data);

                    kadd(pn_even, cn_even, gcol + cfgidx - 1, w4);
                    kadd(pn_odd,  cn_odd,  gcol + cfgidx - 1, w4o);
                    if (isr) {
                        kadd(pr_even, cr_even, gcol + cfgidx - 1, w4);
                        kadd(pr_odd,  cr_odd,  gcol + cfgidx - 1, w4o);
                    }

                    /* Jackknife: once per DISTINCT region touched; even and
                     * odd channels laid out as Fortran N4jk(cfg, ch, region):
                     * index = (cfg-1) + ncfg*((ch-1) + 2*(region-1)). */
                    if (njk > 0) {
                        int jr4 = region[n3 - 1];
                        long jb = (long)(cfgidx - 1);
                        if (jr1 > 0) {
                            long je = jb + (long)ncfg * 2 * (jr1 - 1);
                            atomic_fadd(&jkn[je], w4);
                            atomic_fadd(&jkn[je + ncfg], w4o);
                            if (isr) { atomic_fadd(&jkr[je], w4);
                                       atomic_fadd(&jkr[je + ncfg], w4o); }
                        }
                        if (jr2 > 0 && jr2 != jr1) {
                            long je = jb + (long)ncfg * 2 * (jr2 - 1);
                            atomic_fadd(&jkn[je], w4);
                            atomic_fadd(&jkn[je + ncfg], w4o);
                            if (isr) { atomic_fadd(&jkr[je], w4);
                                       atomic_fadd(&jkr[je + ncfg], w4o); }
                        }
                        if (jr3 > 0 && jr3 != jr1 && jr3 != jr2) {
                            long je = jb + (long)ncfg * 2 * (jr3 - 1);
                            atomic_fadd(&jkn[je], w4);
                            atomic_fadd(&jkn[je + ncfg], w4o);
                            if (isr) { atomic_fadd(&jkr[je], w4);
                                       atomic_fadd(&jkr[je + ncfg], w4o); }
                        }
                        if (jr4 > 0 && jr4 != jr1 && jr4 != jr2 && jr4 != jr3) {
                            long je = jb + (long)ncfg * 2 * (jr4 - 1);
                            atomic_fadd(&jkn[je], w4);
                            atomic_fadd(&jkn[je + ncfg], w4o);
                            if (isr) { atomic_fadd(&jkr[je], w4);
                                       atomic_fadd(&jkr[je + ncfg], w4o); }
                        }
                    }
                }
            }
        }
    }
    done_flag[g] = 1;   /* work-item completed its full hub sweep (watchdog check) */
}
