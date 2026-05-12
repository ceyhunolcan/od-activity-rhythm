# Pipeline scope

This directory contains the analytic pipeline scripts. Coverage and known
caveats:

## Coverage

- `stage1_build_analytic.py` — merges 13 NHANES XPT files into the analytic
  dataset (n = 2,327 after age/PST/exam/accelerometer inclusion). Produces
  `analytic_full.csv`, `analytic_seqn_list.csv`, `attrition_log.csv`.

- `stage25_extract_hourly.py` — per-participant per-hour mean MIMS from the
  PAXMIN_H minute file, using PAXFTIME for real clock-time hour-of-day.
  Produces `hourly_profile.csv` underlying Figure 2.

- `stage8_minute_level_fragmentation.py` — minute-level features from PAXMIN_H:
  ASTP/SATP, bout-length distributions, prolonged-sitting fractions, discrete-
  time hazards, and time-of-day decomposition. Produces
  `fragmentation_features.csv` (paper Table S14).

- `stage30_analysis.R` — primary regression analysis (M1–M4 progressive
  adjustment, survey-weighted with cluster-robust SE), BH-FDR within outcome
  family, and MICE-pooled sensitivity (via `mitools::MIcombine` for proper
  survey-design variance under multiple imputation). Produces Tables 2, 3,
  S1, S2, S5.

## Dependency on intermediate stages not in this directory

`stage30_analysis.R` joins `analytic_full.csv` (from stage 1) and
`fragmentation_features.csv` (from stage 8) and then expects an
`activity_summary.csv` containing the per-day rhythm metrics
(`mean_mims`, `mvpa_min`, `IS`, `IV`, `ASTP`, `RA`, `M10`, `L5`,
`total_sleep_min`, `WASO`, `sleep_efficiency`). That summary file is built
upstream from the PAXMIN_H per-day file using Karas et al. cut-points for
MVPA and van Hees et al. nonparametric rhythm metrics; the build step lives
outside this repository because the per-day file is intermediate output that
is not redistributed with this archive. `stage30_analysis.R` exits with a
clear error message if `activity_summary.csv` is missing.

## Coding notes

- All NHANES Yes/No covariates (DIQ010 diabetes, BPQ020 hypertension, CSQ240
  head injury, CSQ260 sinus history, CSQ010 self-reported smell) preserve
  NaN for refused/don't-know/missing rather than collapsing into 0.
- Smoking is coded `never` / `former` / `current` based on the SMQ020 +
  SMQ040 combination; sentinel codes (7 = Refused, 9 = Don't know) preserve
  NaN.
- `od_binary` (PST ≤ 5) preserves NaN when any of the 8 PST items is missing.
- `comorbidity_count` from MCQ160A–N + MCQ220 treats refused/DK as "not
  having" the condition, the standard NHANES analytic simplification.
- The 8 PST item names and per-item correct response codes are documented
  inline in `stage1_build_analytic.py` (`PST_ITEMS` constant) against the
  CSX_H codebook.
- Stage 8 computes bouts and transitions within contiguous wake blocks
  (PAXSSNMP gap ≤ 1 minute between consecutive wake minutes), preventing
  across-sleep concatenation from inflating bout lengths or transition counts.
- Stage 30 uses `mitools::MIcombine` rather than `mice::pool` for MICE-pooled
  inference, because the standard `pool()` does not propagate complex-survey
  variance through Rubin's rules. Cohen's *d* is standardized by the
  survey-weighted SD of the outcome (`svyvar` on the full design) to keep
  the standardizer on the same population scale as the regression coefficient.
