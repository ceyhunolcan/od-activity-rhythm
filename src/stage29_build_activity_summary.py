#!/usr/bin/env python3
"""
Stage 29: build activity_summary.csv for stage 30.

Stage 2 writes paxmin_features.csv with per-participant activity, rhythm and
sleep measures, but its column names are cut-point-set specific and differ from
the names stage30_analysis.R expects. This stage selects and renames them.

The --cutset argument chooses which cut-point set supplies the sedentary,
light and MVPA minute counts, so the same pipeline can be run under each of
the sets defined in CUT_SETS in stage 2.

Usage:
    python src/stage29_build_activity_summary.py \
        --features paxmin_output/paxmin_features.csv \
        --cutset karas_mixed \
        --out activity_summary.csv
"""

import argparse

import pandas as pd

# stage 2 name -> stage 30 name, for columns that do not depend on the cut-point set
RENAMES = {
    'mims_mean_per_day':    'mean_mims',
    'astp':                 'ASTP',
    'total_sleep_time_min': 'total_sleep_min',
    'waso_min':             'WASO',
}

# columns stage 30 reads under their own names, no rename needed
PASSTHROUGH = ['SEQN', 'IS', 'IV', 'RA', 'M10', 'L5', 'sleep_efficiency',
               'n_valid_days', 'meets_4day_inclusion']


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--features', required=True,
                    help='paxmin_features.csv written by stage 2')
    ap.add_argument('--cutset', default='karas_mixed',
                    help='cut-point set supplying the intensity minute counts')
    ap.add_argument('--out', default='activity_summary.csv')
    args = ap.parse_args()

    feats = pd.read_csv(args.features)

    mvpa_col = f'mvpa_min_day_{args.cutset}'
    sed_col = f'sedentary_min_day_{args.cutset}'
    light_col = f'light_min_day_{args.cutset}'
    for col in (mvpa_col, sed_col, light_col):
        if col not in feats.columns:
            raise SystemExit(
                f"{args.features} has no column {col!r}. Available cut-point "
                f"sets: "
                + ", ".join(sorted(
                    c.replace('mvpa_min_day_', '')
                    for c in feats.columns if c.startswith('mvpa_min_day_')
                ))
            )

    keep = [c for c in PASSTHROUGH if c in feats.columns]
    out = feats[keep].copy()

    for old, new in RENAMES.items():
        if old in feats.columns:
            out[new] = feats[old]

    out['mvpa_min'] = feats[mvpa_col]
    out['sedentary_min'] = feats[sed_col]
    out['light_min'] = feats[light_col]
    out['cutset'] = args.cutset

    out.to_csv(args.out, index=False)
    print(f'Wrote {args.out}: rows={len(out)}, cols={out.shape[1]}, '
          f'cutset={args.cutset}')
    missing = [n for n in ('mean_mims', 'mvpa_min', 'IS', 'IV', 'ASTP')
               if n not in out.columns]
    if missing:
        print(f'WARNING: stage 30 primary outcomes missing: {missing}')


if __name__ == '__main__':
    main()
