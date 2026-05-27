# Project status

## Paper #1 (this repo): cross-sectional OD × activity in NHANES 2013–2014

| | |
|---|---|
| Status | Submitted to medRxiv and *Lancet Healthy Longevity* (May 2026) |
| Manuscript | `paper/manuscript.pdf` |
| Combined submission PDF | `paper/manuscript_complete.pdf` (text + tables + figures) |
| Supplementary appendix | `paper/supplementary/supplementary_appendix.pdf` |
| STROBE checklist | `paper/STROBE_checklist.pdf` |
| Word counts | Abstract 259, body 5,183 |
| Pre-registration | None for paper #1 (acknowledged in Methods); paper #2 pre-registered before data acquisition |

## Paper #2 (separate repo, deferred): mortality follow-up

| | |
|---|---|
| Status | Pre-registered, awaiting completion of paper #1 review cycle |
| Pre-registration | OSF [10.17605/OSF.IO/ZX8RN](https://doi.org/10.17605/OSF.IO/ZX8RN) |
| Decision threshold | 30% attenuation of OD–mortality HR when activity signature added to Cox model |
| Data source | NHANES Linked Mortality File (public-use, downloaded after pre-registration) |
| Lock-in date | Pre-registration locked May 2026 |

## Compute and reproducibility

- Code in `src/`, pure Python and R
- Stage 1 reconstructs the analytic dataset from raw NHANES XPTs
- Stage 8 takes ~25 min on a modern laptop; everything else completes in seconds
- All sensitivity analyses (MICE, IPSW bootstrap, split-half, frailty, severity-gradient,
  cut-point variants, self-report comparison, race-stratified) are in `src/stage30_analysis.R`
- Pinned dependency versions in `requirements.txt` and `r_requirements.txt`

## License summary

| Asset | License |
|-------|---------|
| Code (Python, R) | MIT |
| Manuscript, figures, supplementary appendix | CC BY 4.0 |
| NHANES source data | U.S. public domain (CDC/NCHS) |
| OSF pre-registration | CC0 |
