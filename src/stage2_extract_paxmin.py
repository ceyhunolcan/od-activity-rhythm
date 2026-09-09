#!/usr/bin/env python3
"""
Paper #1: Olfactory dysfunction and accelerometer-derived physical activity
Stage 2: PAXMIN_H feature extraction (run LOCALLY, patched sleep algorithm).
"""
import argparse
import os
import time
import numpy as np
import pandas as pd

CUT_SETS = {
    # Karas 2022 Table 4 maps activity-count cut-points onto MIMS: count 1853 ->
    # 10.558 (sedentary/active in older adults) and count 3940 -> 19.614
    # (light to moderate-vigorous). karas_table4 uses those published values.
    'karas_table4': {'sed_max': 10.558, 'light_max': 19.614},
    # Retained for comparison with the original analysis. The 37.5 upper bound
    # is not a Karas value and its provenance could not be established.
    'karas':       {'sed_max': 10.558, 'light_max': 37.5},
    'wolff_hughes':{'sed_max': 8.6,    'light_max': 24.6},
    'belcher':     {'sed_max': 12.0,   'light_max': 37.6},
}

VALID_WEAR_MIN_PER_DAY = 960
MIN_VALID_DAYS         = 4
SAMPLES_PER_MINUTE     = 4800
MINUTES_PER_DAY        = 1440
SLEEP_DAY_OFFSET_MIN   = 720
MIN_MAIN_SLEEP_MIN     = 240
SLEEP_MERGE_GAP_MIN    = 30


def log(msg, fh=None):
    line = f'[{time.strftime("%H:%M:%S")}] {msg}'
    print(line, flush=True)
    if fh is not None:
        fh.write(line + '\n'); fh.flush()


def derive_minute_index(df):
    samples = df['PAXSSNMP'].astype('int64')
    df['total_minute']  = (samples // SAMPLES_PER_MINUTE).astype('int32')
    df['minute_of_day'] = (df['total_minute'] % MINUTES_PER_DAY).astype('int16')
    return df


def find_main_sleep_period(states, sleep_state=2,
                           merge_gap=SLEEP_MERGE_GAP_MIN,
                           min_period=MIN_MAIN_SLEEP_MIN):
    """Longest run of mostly-sleep minutes, allowing intra-sleep awakenings
    up to merge_gap. Choi 2011 style. Returns (start, end_excl) or (None, None).
    """
    n = len(states)
    if n == 0:
        return None, None
    is_sleep = (states == sleep_state).astype(np.int8)
    if is_sleep.sum() == 0:
        return None, None
    diffs  = np.diff(is_sleep, prepend=0, append=0)
    starts = np.where(diffs == 1)[0]
    ends   = np.where(diffs == -1)[0]
    if len(starts) == 0:
        return None, None
    merged = [(int(starts[0]), int(ends[0]))]
    for s, e in zip(starts[1:], ends[1:]):
        last_s, last_e = merged[-1]
        if int(s) - last_e <= merge_gap:
            merged[-1] = (last_s, int(e))
        else:
            merged.append((int(s), int(e)))
    longest = max(merged, key=lambda x: x[1] - x[0])
    if (longest[1] - longest[0]) < min_period:
        return None, None
    return longest


def compute_features_for_participant(g):
    g = g.sort_values('total_minute').reset_index(drop=True)
    if len(g) == 0:
        return None

    day_rows = []
    for day, gd in g.groupby('PAXDAYM', sort=True):
        st = gd['PAXPREDM'].values.astype(np.int8)
        mm = gd['PAXMTSM'].astype(float).values
        wake_mask    = (st == 1)
        sleep_mask   = (st == 2)
        nonwear_mask = (st == 3)
        wake_min  = int(wake_mask.sum())
        sleep_min = int(sleep_mask.sum())
        wear_min  = wake_min + sleep_min

        cat = {}
        if wake_min > 0:
            wm = mm[wake_mask]
            for name, c in CUT_SETS.items():
                cat[f'sed_min_{name}']   = int((wm <  c['sed_max']).sum())
                cat[f'light_min_{name}'] = int(((wm >= c['sed_max']) &
                                                (wm <  c['light_max'])).sum())
                cat[f'mvpa_min_{name}']  = int((wm >= c['light_max']).sum())
            mims_mean = float(wm.mean())
            peak1m  = float(wm.max())
            peak30m = (float(pd.Series(wm).rolling(30, min_periods=30).mean().max())
                       if wake_min >= 30 else np.nan)
            active = wm >= CUT_SETS['karas']['sed_max']
            if active.sum() > 1:
                trans = int(((active[:-1]) & (~active[1:])).sum())
                astp  = trans / max(int(active[:-1].sum()), 1)
            else:
                astp = np.nan
        else:
            for name in CUT_SETS:
                cat[f'sed_min_{name}']   = 0
                cat[f'light_min_{name}'] = 0
                cat[f'mvpa_min_{name}']  = 0
            mims_mean = peak1m = peak30m = astp = np.nan

        valid_day = int(wear_min >= VALID_WEAR_MIN_PER_DAY)
        row = {'PAXDAYM': int(day), 'wear_min': wear_min,
               'wake_min': wake_min, 'sleep_min_calendar': sleep_min,
               'nonwear_min': int(nonwear_mask.sum()),
               'mims_mean_wake': mims_mean,
               'peak1m_mims': peak1m, 'peak30m_mims': peak30m,
               'astp': astp, 'valid_day': valid_day}
        row.update(cat)
        day_rows.append(row)
    daily = pd.DataFrame(day_rows)

    g['sleep_day'] = (g['total_minute'] - SLEEP_DAY_OFFSET_MIN) // MINUTES_PER_DAY
    sleep_rows = []
    for sd, gs in g.groupby('sleep_day', sort=True):
        st  = gs['PAXPREDM'].values.astype(np.int8)
        mod = gs['minute_of_day'].values.astype(np.int32)
        s, e = find_main_sleep_period(st)
        if s is None:
            continue
        if s == 0 or e == len(st):
            continue
        block_st = st[s:e]
        tib = e - s
        tst = int((block_st == 2).sum())
        waso = tib - tst
        eff  = tst / tib if tib > 0 else np.nan
        onset  = int(mod[s])
        offset = (int(mod[e-1]) + 1) % MINUTES_PER_DAY
        if offset <= onset:
            mid = ((onset + offset + MINUTES_PER_DAY) / 2) % MINUTES_PER_DAY
        else:
            mid = (onset + offset) / 2
        sleep_rows.append({'sleep_onset_min': onset,
                           'sleep_offset_min': offset,
                           'sleep_midpoint_min': float(mid),
                           'time_in_bed_min': int(tib),
                           'total_sleep_time_min': tst,
                           'waso_min': int(waso),
                           'sleep_efficiency': float(eff)})
    sleep_df = pd.DataFrame(sleep_rows)

    n_valid = int(daily['valid_day'].sum())
    if n_valid < MIN_VALID_DAYS:
        return {'SEQN': int(g['SEQN'].iloc[0]),
                'n_valid_days': n_valid,
                'meets_4day_inclusion': 0}

    valid = daily[daily['valid_day'] == 1]
    out = {
        'SEQN': int(g['SEQN'].iloc[0]),
        'n_valid_days': n_valid,
        'meets_4day_inclusion': 1,
        'mean_wear_hours_per_day': float(valid['wear_min'].mean()) / 60.0,
        'mims_mean_per_day': float(valid['mims_mean_wake'].mean()),
        'peak1m_mims':       float(valid['peak1m_mims'].mean()),
        'peak30m_mims':      float(valid['peak30m_mims'].mean()),
    }
    for k in CUT_SETS:
        out[f'sedentary_min_day_{k}'] = float(valid[f'sed_min_{k}'].mean())
        out[f'light_min_day_{k}']     = float(valid[f'light_min_{k}'].mean())
        out[f'mvpa_min_day_{k}']      = float(valid[f'mvpa_min_{k}'].mean())
    out['astp'] = float(valid['astp'].mean(skipna=True))
    if not sleep_df.empty:
        out.update({
            'sleep_n_nights':          int(len(sleep_df)),
            'sleep_onset_hour_mean':   float(sleep_df['sleep_onset_min'].mean()) / 60.0,
            'sleep_offset_hour_mean':  float(sleep_df['sleep_offset_min'].mean()) / 60.0,
            'sleep_midpoint_hour':     float(sleep_df['sleep_midpoint_min'].mean()) / 60.0,
            'time_in_bed_min':         float(sleep_df['time_in_bed_min'].mean()),
            'total_sleep_time_min':    float(sleep_df['total_sleep_time_min'].mean()),
            'waso_min':                float(sleep_df['waso_min'].mean()),
            'sleep_efficiency':        float(sleep_df['sleep_efficiency'].mean()),
        })
    else:
        for k in ['sleep_n_nights','sleep_onset_hour_mean','sleep_offset_hour_mean',
                  'sleep_midpoint_hour','time_in_bed_min','total_sleep_time_min',
                  'waso_min','sleep_efficiency']:
            out[k] = np.nan

    valid_days = set(int(x) for x in valid['PAXDAYM'])
    g_v = g[g['PAXDAYM'].astype(int).isin(valid_days)].copy()
    g_v['hour_of_day'] = (g_v['minute_of_day'] // 60).astype('int8')
    g_v['day_hour_idx'] = (g_v['PAXDAYM'].astype(int) * 24 +
                           g_v['hour_of_day'].astype(int))
    hourly = g_v.groupby('day_hour_idx')['PAXMTSM'].mean().sort_index()

    if len(hourly) >= 48:
        x  = hourly.values.astype(float)
        N  = len(x)
        gm = x.mean()
        means_by_hour = (g_v.groupby('hour_of_day')['PAXMTSM'].mean()
                            .reindex(range(24)))
        num_is = N * np.nansum((means_by_hour.values - gm) ** 2)
        den_is = 24 * np.sum((x - gm) ** 2)
        IS = float(num_is / den_is) if den_is > 0 else np.nan
        d = np.diff(x)
        num_iv = N * np.sum(d ** 2)
        den_iv = (N - 1) * np.sum((x - gm) ** 2)
        IV = float(num_iv / den_iv) if den_iv > 0 else np.nan
        s = pd.Series(x)
        M10 = float(s.rolling(10).mean().max())
        L5  = float(s.rolling(5).mean().min())
        RA  = (M10 - L5) / (M10 + L5) if (M10 + L5) > 0 else np.nan
        out.update({'IS': IS, 'IV': IV, 'RA': RA, 'M10': M10, 'L5': L5})
    else:
        out.update({k: np.nan for k in ['IS','IV','RA','M10','L5']})

    return out


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--paxmin', required=True)
    p.add_argument('--paxhd',  required=True)
    p.add_argument('--seqn',   required=True)
    p.add_argument('--out',    required=True)
    p.add_argument('--chunksize', type=int, default=500_000)
    args = p.parse_args()

    os.makedirs(args.out, exist_ok=True)
    fh = open(os.path.join(args.out, 'paxmin_features_log.txt'), 'w')
    log('=== Stage 2 PAXMIN feature extraction (patched sleep) ===', fh)
    log(f'  paxmin: {args.paxmin}', fh)
    log(f'  cut-point sets: {list(CUT_SETS.keys())} (primary=karas)', fh)
    log(f'  valid-day threshold: >={VALID_WEAR_MIN_PER_DAY}min wear; '
        f'>={MIN_VALID_DAYS} valid days; sleep merge gap {SLEEP_MERGE_GAP_MIN}min', fh)

    seqn_set = set(pd.read_csv(args.seqn)['SEQN'].astype('int64'))
    log(f'\nStage 1 analytic SEQNs: n={len(seqn_set)}', fh)

    paxhd = pd.read_sas(args.paxhd, format='xport')
    paxhd['SEQN'] = paxhd['SEQN'].astype('int64')
    paxhd = paxhd[paxhd['SEQN'].isin(seqn_set)][['SEQN','PAXSTS','PAXHAND']].copy()
    log(f'PAXHD rows for analytic SEQNs: {len(paxhd)}', fh)
    log(f'PAXSTS distribution: {paxhd["PAXSTS"].value_counts().to_dict()}', fh)

    log(f'\nStreaming PAXMIN_H in chunks of {args.chunksize:,}...', fh)
    keep = ['SEQN','PAXSSNMP','PAXMTSM','PAXPREDM','PAXDAYM','PAXDAYWM']
    chunks = []
    n_total = n_kept = 0
    t0 = time.time()
    for i, chunk in enumerate(pd.read_sas(args.paxmin, format='xport',
                                          chunksize=args.chunksize), 1):
        n_total += len(chunk)
        chunk['SEQN'] = chunk['SEQN'].astype('int64')
        chunk = chunk[chunk['SEQN'].isin(seqn_set)]
        if len(chunk):
            chunk = chunk[keep].copy()
            chunk['PAXSSNMP'] = chunk['PAXSSNMP'].astype('int64')
            chunk['PAXMTSM']  = chunk['PAXMTSM'].astype('float32')
            chunk['PAXPREDM'] = chunk['PAXPREDM'].astype('int8')
            chunk['PAXDAYM']  = chunk['PAXDAYM'].astype('int8')
            chunk['PAXDAYWM'] = chunk['PAXDAYWM'].astype('int8')
            chunks.append(chunk)
            n_kept += len(chunk)
        if i % 10 == 0:
            elapsed = time.time() - t0
            log(f'  chunk {i:>4}: total={n_total:>12,} kept={n_kept:>10,} '
                f'elapsed={elapsed:6.0f}s', fh)
    df = pd.concat(chunks, ignore_index=True); del chunks
    log(f'Streaming complete. Read {n_total:,}; kept {n_kept:,}.', fh)
    log(f'Filtered DF memory: {df.memory_usage(deep=True).sum()/1e9:.2f} GB', fh)

    df = derive_minute_index(df)

    log(f'\nComputing per-participant features for {df["SEQN"].nunique()} participants...', fh)
    rows = []
    t1 = time.time()
    for j, (seqn, g) in enumerate(df.groupby('SEQN', sort=False), 1):
        feats = compute_features_for_participant(g)
        if feats is not None:
            rows.append(feats)
        if j % 500 == 0:
            log(f'  processed {j} participants in {time.time()-t1:.0f}s', fh)
    feats_df = pd.DataFrame(rows).merge(paxhd, on='SEQN', how='left')

    out_csv = os.path.join(args.out, 'paxmin_features.csv')
    feats_df.to_csv(out_csv, index=False)
    n_incl = int(feats_df['meets_4day_inclusion'].sum())
    log(f'\nWrote {out_csv}: rows={len(feats_df)}, cols={feats_df.shape[1]}', fh)
    log(f'Meeting >={MIN_VALID_DAYS}-valid-day inclusion: {n_incl} '
        f'({100*n_incl/len(feats_df):.1f}%)', fh)

    log('\nDescriptive snapshot among included participants:', fh)
    incl = feats_df[feats_df['meets_4day_inclusion'] == 1]
    for c in ['mims_mean_per_day','sedentary_min_day_karas','mvpa_min_day_karas',
              'total_sleep_time_min','waso_min','sleep_efficiency','IS','IV','RA']:
        if c in incl.columns:
            s = incl[c]
            log(f'  {c:30s}  mean={s.mean():>8.3f}  sd={s.std():>8.3f}  '
                f'n_nonna={s.notna().sum():>5d}', fh)
    fh.close()


if __name__ == '__main__':
    main()
