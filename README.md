# od-activity-rhythm

Analysis code, manuscript, and supplementary materials for:

**Olfactory dysfunction, daytime activity reduction, and 24-hour rhythm fragmentation - but not food-odor recognition - in NHANES 2013–2014**

Ceyhun Olcan, Dartmouth College.

- **medRxiv preprint**: https://doi.org/10.1101/2026.MM.DD.XXXXXXXX *(DOI pending screening; typically 2–4 business days after submission)*
- **OSF pre-registration** (paper #2 follow-up): https://doi.org/10.17605/OSF.IO/ZX8RN
- **Zenodo archive of this repo** (citable concept DOI): https://doi.org/10.5281/zenodo.20132927

---

## What's here

```
paper/
  manuscript.md              Markdown source
  manuscript.pdf             Text-only PDF (14 pages)
  manuscript_complete.pdf    Text + tables + 5 figures (22 pages)
  STROBE_checklist.md/.pdf   STROBE Statement v4 compliance
  tables/                    Main tables (1–4) as md, pdf, and CSV
  figures/                   Main figures (1–5) as pdf and png
  supplementary/             Supplementary appendix + 21 supp tables + 2 supp figures
src/
  stage1_build_analytic.py             construct analytic dataset from 18 NHANES XPTs
  stage25_extract_hourly.py            hourly mean MIMS per participant
  stage8_minute_level_fragmentation.py bouts/transitions/hazards from PAXMIN_H
  stage30_analysis.R                   primary regression, FDR, MICE pooled
data/
  attrition_log.csv          STROBE attrition counts
  analytic_seqn_list.csv     2,327 SEQNs in the final analytic sample
docs/
  variable_dictionary.md     data dictionary for analytic_full.csv
```

The actual NHANES XPT files are not redistributed - they live at
https://wwwn.cdc.gov/nchs/nhanes/ and are downloaded into the working dir
before running stage 1.

## Reproducing the analysis

1. Download NHANES 2013–2014 cycle files: DEMO_H, BMX_H, BPX_H, CSX_H, CSQ_H,
   PAXMIN_H, PAXHD_H, PAQ_H, SMQ_H, DIQ_H, GHB_H, MCQ_H, BPQ_H, DPQ_H, INQ_H,
   HUQ_H, RXQ_RX_H, RDQ_H, SLQ_H.

2. From the directory containing the XPTs:
   ```bash
   python src/stage1_build_analytic.py
   python src/stage25_extract_hourly.py \
       --paxmin PAXMIN_H.xpt --paxhd PAXHD_H.xpt \
       --features paxmin_features.csv \
       --seqn analytic_seqn_list.csv \
       --out paxmin_output
   python src/stage8_minute_level_fragmentation.py
   Rscript src/stage30_analysis.R
   ```

   Stage 8 takes ~25 min on a 2022 MacBook Air; everything else is faster.

## Dependencies

- Python 3.10+ with pandas, numpy
- R 4.3+ with survey, mice, dplyr, readr

See `requirements.txt` and `r_requirements.txt`.

## Data dictionary

See `docs/variable_dictionary.md`. The most-used variables:

| name | what |
|------|------|
| `SEQN` | NHANES participant ID |
| `PST_correct` | Pocket Smell Test, # correctly identified (0–8) |
| `od_binary` | 1 if PST_correct ≤ 5, else 0 |
| `mean_mims` | per-day mean MIMS units |
| `mvpa_min` | min/day in MVPA (Karas: ≥ 37.5 MIMS) |
| `IS`, `IV` | interdaily stability, intradaily variability |
| `ASTP_minute` | active-to-sedentary transition probability (minute-level) |
| `SATP_minute` | sedentary-to-active transition probability |
| `WTMEC2YR`, `SDMVSTRA`, `SDMVPSU` | NHANES MEC weights and clustering |

## Citing

If you use this code or build on these findings, please cite:

> Olcan C. Olfactory dysfunction, daytime activity reduction, and 24-hour
> rhythm fragmentation - but not food-odor recognition - in NHANES 2013–2014.
> medRxiv 2026 (preprint, DOI pending).

For the code archive specifically: cite the Zenodo concept DOI
[10.5281/zenodo.20132927](https://doi.org/10.5281/zenodo.20132927).

## License

- **Code**: MIT (see [LICENSE](LICENSE))
- **Manuscript and figures**: CC BY 4.0 (matches Lancet Healthy Longevity policy)
- **NHANES data**: U.S. public domain (CDC/NCHS)

## Contact

ceyhun.olcan.27@dartmouth.edu | ORCID [0000-0002-6326-6071](https://orcid.org/0000-0002-6326-6071)
