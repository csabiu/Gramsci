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
                         __global const int   *buffer,
                         __global const int   *bt3,    /* nbins^3 config table */
                         __global float       *part_nnn,
                         __global float       *part_rrr,
                         const int nbins, const int cb,
                         const int num_data, const int istart,
                         const int iend, const int ngang,
                         const int nbuckets, const int blo, const int bhi,
                         __global char *done_flag)
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
        if (buffer[i - 1] == 1) continue;
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
                int  n2 = id  [b0 + k2 - 1];

                char c3 = find_bin(ptr, id, dist, n1, n2);
                if (c3 == 0) continue;

                int bin = bt3[(c1 - 1) + nbins * ((c2 - 1) + nbins * (c3 - 1))];
                float w3 = wi_w1 * w[n2 - 1];

                part_nnn[gcol + bin - 1] += w3;
                if (i > num_data && n1 > num_data && n2 > num_data)
                    part_rrr[gcol + bin - 1] -= w3;
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
                          __global const int   *buffer,
                          __global float       *part_nnn,
                          __global float       *part_rrr,
                          const int nbins, const int cb,
                          const int num_data, const int istart,
                          const int iend, const int ngang,
                         const int nbuckets, const int blo, const int bhi,
                         __global char *done_flag)
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
        if (buffer[i - 1] == 1) continue;
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

                part_nnn[gcol + bin - 1] += w3;
                if (i > num_data && n1 > num_data && n2 > num_data)
                    part_rrr[gcol + bin - 1] -= w3;
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
                         __global const int   *buffer,
                         __global const int   *bt6,    /* nbins^6 config table */
                         __global float       *part_n4,
                         __global float       *part_r4,
                         const int nbins, const int ncfg,
                         const int num_data,
                         const int istart, const int iend, const int ngang,
                         const int nbuckets, const int blo, const int bhi,
                         __global char *done_flag)
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
        if (buffer[i - 1] == 1) continue;
        long s  = ptr[i - 1];
        int  nn = (int)(ptr[i] - s);
        if (nn <= 2) continue;
        long b0 = s - 1;

        for (int k1 = 1; k1 <= nn; k1++) {
            char c1 = dist[b0 + k1 - 1];
            int  n1 = id  [b0 + k1 - 1];

            for (int k2 = k1 + 1; k2 <= nn; k2++) {
                int  n2 = id[b0 + k2 - 1];
                char c4 = find_bin(ptr, id, dist, n1, n2);
                if (c4 == 0) continue;
                char c2 = dist[b0 + k2 - 1];
                float w12 = w[i - 1] * w[n1 - 1] * w[n2 - 1];
                int rand12 = (i > num_data && n1 > num_data && n2 > num_data);

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
                    part_n4[gcol + cfgidx - 1] += w4;
                    if (rand12 && n3 > num_data)
                        part_r4[gcol + cfgidx - 1] += w4;
                }
            }
        }
    }
    done_flag[g] = 1;   /* work-item completed its full hub sweep (watchdog check) */
}

/* --------------------------------------------------------------------------
 * 4PCF with parity decomposition (pure binary-search; see k_4pcf_all note).
 * Even channel  += w4;  odd channel += (parity_flip * sign_V) * w4.
 * parity_flip = sign(bt6 entry); sign_V is looked up from the host-built
 * signv table (computed in double, so the chirality matches the CPU exactly).
 *   signv[(p1-1)*ndir*ndir + (p2-1)*ndir + (p3-1)] in {-1,0,+1}
 * -------------------------------------------------------------------------- */
__kernel void k_4pcf_parity(__global const long  *ptr,
                            __global const int   *id,
                            __global const char  *dist,
                            __global const char  *phi,
                            __global const float *w,
                            __global const int   *buffer,
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
                         __global char *done_flag)
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
        if (buffer[i - 1] == 1) continue;
        long s  = ptr[i - 1];
        int  nn = (int)(ptr[i] - s);
        if (nn <= 2) continue;
        long b0 = s - 1;

        for (int k1 = 1; k1 <= nn; k1++) {
            char c1 = dist[b0 + k1 - 1];
            int  n1 = id  [b0 + k1 - 1];
            int  p1 = (int)phi[b0 + k1 - 1];

            for (int k2 = k1 + 1; k2 <= nn; k2++) {
                int  n2 = id[b0 + k2 - 1];
                char c4 = find_bin(ptr, id, dist, n1, n2);
                if (c4 == 0) continue;
                char c2 = dist[b0 + k2 - 1];
                int  p2 = (int)phi[b0 + k2 - 1];
                float w12 = w[i - 1] * w[n1 - 1] * w[n2 - 1];
                int rand12 = (i > num_data && n1 > num_data && n2 > num_data);

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
                    int pflip = (raw > 0) ? 1 : -1;

                    int p3 = (int)phi[b0 + k3 - 1];
                    int sgn = (int)signv[(long)(p1 - 1) * ndir * ndir
                                         + (p2 - 1) * ndir + (p3 - 1)];

                    float w4 = w12 * w[n3 - 1];
                    float w4o = (float)(pflip * sgn) * w4;

                    pn_even[gcol + cfgidx - 1] += w4;
                    pn_odd [gcol + cfgidx - 1] += w4o;
                    if (rand12 && n3 > num_data) {
                        pr_even[gcol + cfgidx - 1] += w4;
                        pr_odd [gcol + cfgidx - 1] += w4o;
                    }
                }
            }
        }
    }
    done_flag[g] = 1;   /* work-item completed its full hub sweep (watchdog check) */
}
