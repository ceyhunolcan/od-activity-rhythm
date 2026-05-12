# od-activity-rhythm

Analysis code and pipeline for:

**Olfactory dysfunction, daytime activity reduction, and 24-hour rhythm fragmentation — but not food-odor recognition — in NHANES 2013–2014**

Ceyhun Olcan, Dartmouth College.
medRxiv preprint: [DOI pending].

---

## What's here

```
src/
  stage1_build_analytic.py       construct analytic dataset from 18 NHANES XPTs
  stage25_extract_hourly.py      hourly mean MIMS per participant
  stage8_minute_level_fragmentation.py    bouts/transitions/hazards from PAXMIN_H
  stage30_analysis.R             primary regression, FDR, MICE pooled
data/
  attrition_log.csv              STROBE attrition counts
  analytic_seqn_list.csv         2,327 SEQNs in the final analytic sample
docs/
  variable_dictionary.md         data dictionary for analytic_full.csv
```

The actual NHANES XPT files are not redistributed — they live at
https://wwwn.cdc.gov/nchs/nhanes/ and are downloaded into the working dir
before running stage 1.

## Reproducing

1. Download NHANES 2013–2014 cycle files: DEMO_H, BMX_H, BPX_H, CSX_H, CSQ_H,
   PAXMIN_H, PAXHD_H, PAQ_H, SMQ_H, DIQ_H, GHB_H, MCQ_H, BPQ_H, DPQ_H, INQ_H,
   HUQ_H, RXQ_RX_H, RDQ_H, SLQ_H.

2. From the directory containing the XPTs:
   ```
   python src/stage1_build_analytic.py
   python src/stage25_extract_hourly.py --paxmin PAXMIN_H.xpt --paxhd PAXHD_H.xpt \
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
| `od_binary` | 1 if PST_correct <= 5, else 0 |
| `mean_mims` | per-day mean MIMS units |
| `mvpa_min` | min/day in MVPA (Karas: >= 37.5 MIMS) |
| `IS`, `IV` | interdaily stability, intradaily variability |
| `ASTP_minute` | active-to-sedentary transition probability (minute-level) |
| `SATP_minute` | sedentary-to-active transition probability |
| `WTMEC2YR`, `SDMVSTRA`, `SDMVPSU` | NHANES MEC weights and clustering |

## License

Code: MIT (see LICENSE).
Data: NHANES is public domain (US CDC/NCHS).

## Citing

If you use this code please cite the manuscript above. The exact analytic
SEQN list and attrition counts are reproduced in `data/`.

## Contact

ceyhun.olcan.27@dartmouth.edu
