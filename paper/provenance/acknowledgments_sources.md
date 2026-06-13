# Acknowledgments — source provenance and submission checklist

The mandated acknowledgment text was retrieved and verified against
official sources (multi-agent research, 2026-06-13). Keep for the record.

## DESI DR1 funding paragraph (MANDATORY — verbatim)
Verified element-by-element as ACCEPTABLE against:
- https://data.desi.lbl.gov/doc/acknowledgments/
- DR1 papers arXiv:2511.20757, arXiv:2512.03229 (Dec 2025) — identical funder list
Confirmed-correct details often gotten wrong:
- "CONAHCYT" (current), NOT the older "CONACYT"
- "Ministry of Science, Innovation and Universities of Spain
  (MICIU/AEI/10.13039/501100011033)", NOT the older "MICINN"
- Contract DE-AC02-05CH11231 + NERSC named.

### ⚠️ BEFORE SUBMISSION
1. Re-confirm the paragraph against the LIVE page
   https://data.desi.lbl.gov/doc/acknowledgments/ (DESI updates it periodically).
2. Insert author (CGS) personal funding/grant support — the one remaining
   \todo in gramsci_v2.tex.
3. If the final analysis uses any DESI Value-Added Catalog, Legacy Imaging
   Surveys, or NERSC SCGSR support, append those extra blocks from the same page.
4. Optional: Tohono O'odham Nation / Kitt Peak land-acknowledgment sentence
   (some DESI papers include it).
5. Add the DESI DR1 LSS-catalogue construction paper + systematics papers
   (fiber-assignment arXiv:2411.12025, spectroscopic-systematics arXiv:2405.17208)
   to the citations if the referee/policy requires the specific data-product refs.
6. EZmock DESI-DR1 production paper: currently cited only via the method
   (Chuang et al. 2015). Add Zhao et al. (DESI DR1 mocks) once it has an arXiv id.

## Data products / citations
- DESI DR1 release: DESI Collaboration 2026, AJ 171, 285 (arXiv:2503.14745) -> bib key `desidr1`
- Quijote: Villaescusa-Navarro et al. 2020 (arXiv:1909.05273) -> `villaescusa20` (no mandated ack text)
- EZmock method: Chuang et al. 2015 (arXiv:1409.1124) -> `chuang15`

## Software citations added to bib
astropy13/18/22, hunter07 (Matplotlib), virtanen20 (SciPy); reused harris20 (NumPy),
kennel04 (kdtree2), sabiu19 (GRAMSCI). NVIDIA HPC SDK + OpenACC: URL footnotes only
(no canonical paper). \facility omitted (no observational facility).
