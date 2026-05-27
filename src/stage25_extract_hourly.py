"""
Compute per-participant per-hour mean MIMS from PAXMIN_H.

Uses PAXFTIME (start time of day from PAXHD_H) to convert PAXSSNMP
(samples since monitor start, 80Hz) into real clock-time hour-of-day.

Output: hourly_profile.csv -- one row per SEQN x hour pair.
"""

import argparse
import os
import time
import numpy as np
import pandas as pd


SAMPLES_PER_MINUTE = 4800


def parse_ftime(x):
    # PAXFTIME may come back as bytes, string, or numeric (seconds since midnight)
    if x is None or (isinstance(x, float) and np.isnan(x)):
        return np.nan
    if isinstance(x, bytes):
        try:
            x = x.decode('utf-8', errors='ignore')
        except Exception:
            return np.nan
    if isinstance(x, (int, float, np.integer, np.floating)):
        return float(x) / 60.0
    s = str(x).strip().strip("'\"")
    if not s:
        return np.nan
    if ':' in s:
        try:
            parts = s.split(':')
            return int(parts[0]) * 60 + int(parts[1]) + (int(parts[2]) if len(parts) > 2 else 0) / 60.0
        except Exception:
            return np.nan
    try:
        return float(s) / 60.0
    except ValueError:
        return np.nan


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--paxmin', required=True)
    ap.add_argument('--paxhd', required=True)
    ap.add_argument('--features', required=True,
                    help='paxmin_features.csv with meets_4day_inclusion flag')
    ap.add_argument('--seqn', required=True, help='analytic_seqn_list.csv')
    ap.add_argument('--out', required=True)
    ap.add_argument('--chunksize', type=int, default=500_000)
    args = ap.parse_args()

    os.makedirs(args.out, exist_ok=True)

    # who to keep
    seqns = set(pd.read_csv(args.seqn)['SEQN'].astype('int64'))
    feats = pd.read_csv(args.features)
    if 'meets_4day_inclusion' not in feats.columns:
        raise KeyError(
            f"{args.features} is missing 'meets_4day_inclusion' column "
            f"(have: {list(feats.columns)})"
        )
    incl = set(feats[feats['meets_4day_inclusion'] == 1]['SEQN'].astype('int64'))
    keep = seqns & incl
    print(f'Keeping n = {len(keep)} included participants')

    # PAXHD -> start time per participant
    print('Reading PAXHD...')
    if args.paxhd.endswith(('.xpt', '.xpt.txt')):
        ph = pd.read_sas(args.paxhd, format='xport')
    else:
        ph = pd.read_csv(args.paxhd)
    ph['SEQN'] = ph['SEQN'].astype('int64')
    ftime_col = 'PAXFTIME' if 'PAXFTIME' in ph.columns else 'PAXFTM'
    if ftime_col not in ph.columns:
        raise KeyError(f'No PAXFTIME/PAXFTM in PAXHD; columns: {list(ph.columns)}')
    ph['ftime_min'] = ph[ftime_col].apply(parse_ftime)
    seqn_to_ftime = dict(zip(ph['SEQN'], ph['ftime_min']))

    # stream PAXMIN
    print('Streaming PAXMIN_H...')
    chunks = []
    n_total = n_kept = 0
    t0 = time.time()
    for i, ch in enumerate(pd.read_sas(args.paxmin, format='xport', chunksize=args.chunksize), 1):
        n_total += len(ch)
        ch['SEQN'] = ch['SEQN'].astype('int64')
        ch = ch[ch['SEQN'].isin(keep)]
        if len(ch):
            ch = ch[['SEQN', 'PAXSSNMP', 'PAXMTSM']].copy()
            ch['PAXMTSM'] = ch['PAXMTSM'].astype('float32')
            chunks.append(ch)
            n_kept += len(ch)
        if i % 20 == 0:
            print(f'  chunk {i}: total={n_total:,} kept={n_kept:,} ({time.time()-t0:.0f}s)')
    df = pd.concat(chunks, ignore_index=True)
    del chunks

    # PAXSSNMP -> minute-of-day = (ftime_min + minutes_since_start) mod 1440
    df['ftime_min'] = df['SEQN'].map(seqn_to_ftime)

    # drop participants whose PAXFTIME couldn't be parsed (otherwise the int8
    # cast below silently turns NaN into hour=0)
    n_before = len(df)
    df = df.dropna(subset=['ftime_min']).copy()
    n_dropped = n_before - len(df)
    if n_dropped:
        bad_seqns = sorted({s for s, t in seqn_to_ftime.items() if pd.isna(t)})
        print(f'  dropped {n_dropped:,} rows from {len(bad_seqns)} participants '
              f'with unparseable PAXFTIME (first few: {bad_seqns[:5]})')

    df['min_since_start'] = df['PAXSSNMP'].astype('int64') // SAMPLES_PER_MINUTE
    df['minute_of_day'] = ((df['ftime_min'] + df['min_since_start']) % 1440).astype('float32')
    df['hour'] = (df['minute_of_day'] // 60).astype('int8')

    # aggregate
    hourly = (df.groupby(['SEQN', 'hour'])['PAXMTSM']
                .mean()
                .reset_index()
                .rename(columns={'PAXMTSM': 'mean_mims'}))

    out_path = os.path.join(args.out, 'hourly_profile.csv')
    hourly.to_csv(out_path, index=False)
    print(f'Wrote {out_path}: {len(hourly):,} rows ({hourly["SEQN"].nunique()} participants x 24h)')

    # quick check: peak hours should be morning/evening
    print('\nOverall mean MIMS by hour-of-day:')
    overall = hourly.groupby('hour')['mean_mims'].mean()
    for h in range(24):
        v = overall.get(h, np.nan)
        print(f'  {h:02d}:00  {v:5.2f}' if not np.isnan(v) else f'  {h:02d}:00    n/a')


if __name__ == '__main__':
    main()
