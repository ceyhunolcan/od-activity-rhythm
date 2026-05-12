"""
Minute-level fragmentation features from PAXMIN_H.

For each participant, computes per-minute ASTP/SATP, bout-length distributions
for active and sedentary states, discrete-time hazard rates at 1/5/15/30/60 min,
and ASTP/SATP by time-of-day window (morning/afternoon/evening).

Karas et al. cut-points: sedentary < 10.558 MIMS/min, MVPA >= 37.5 MIMS/min.
PAXPREDM == 1 = wake (analysis state).

PAXSSNMP is samples since monitor start at 80 Hz (4800 samples/minute).
We use PAXFTIME from PAXHD_H to compute the participant's real clock-time
hour-of-day so the morning/afternoon/evening bins are correct.

Outputs fragmentation_features.csv with one row per SEQN.
"""

import argparse
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

warnings.filterwarnings('ignore')


SAMPLES_PER_MIN = 4800
SED_CUT  = 10.558
MVPA_CUT = 37.5
WAKE = 1
HAZ_DURATIONS = [1, 5, 15, 30, 60]


def to_int(series):
    """PAXPREDM can be bytes (b'1'), str ('1'), or float (1.0) -- normalize to int."""
    if series.dtype == 'O':
        def conv(x):
            if pd.isna(x):
                return np.nan
            if isinstance(x, bytes):
                try:
                    return int(x.decode('utf-8').strip())
                except (ValueError, AttributeError):
                    return np.nan
            try:
                return int(str(x).strip())
            except ValueError:
                return np.nan
        return series.apply(conv)
    return pd.to_numeric(series, errors='coerce')


def parse_ftime(x):
    """PAXFTIME from PAXHD_H: HH:MM:SS string, bytes, or numeric seconds-since-midnight.
    Returns minutes-since-midnight as float, or NaN."""
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


def find_bouts(state):
    if len(state) == 0:
        return pd.DataFrame(columns=['state', 'length'])
    state = np.asarray(state)
    change = np.where(np.diff(state) != 0)[0] + 1
    bounds = np.concatenate([[0], change, [len(state)]])
    bouts = []
    for i in range(len(bounds) - 1):
        s, e = bounds[i], bounds[i + 1]
        bouts.append({'state': int(state[s]), 'length': int(e - s)})
    return pd.DataFrame(bouts)


def hazards(bout_lengths, durations=HAZ_DURATIONS):
    bout_lengths = np.asarray(bout_lengths)
    if len(bout_lengths) == 0:
        return [np.nan] * len(durations)
    out = []
    for t in durations:
        at_risk = (bout_lengths >= t).sum()
        terminated = (bout_lengths == t).sum()
        out.append(terminated / at_risk if at_risk > 0 else np.nan)
    return out


def real_hour_of_day(paxssnmp, ftime_min):
    """Convert PAXSSNMP (80 Hz samples since monitor start) to real clock hour-of-day.

    samples / 4800 = minutes since monitor start
    (ftime_min + minutes_since_start) mod 1440 = minute-of-day
    minute-of-day // 60 = hour 0..23
    """
    min_since_start = paxssnmp // SAMPLES_PER_MIN
    minute_of_day = (ftime_min + min_since_start) % 1440
    return (minute_of_day // 60).astype(int)


def features_for_person(pdf, ftime_min):
    """Feature vector for one participant. Returns dict or None.

    Important: PAXPREDM==WAKE filtering drops sleep minutes, which would
    erroneously merge a Tuesday-night sedentary period with a Wednesday-morning
    sedentary period if we ran find_bouts on the full wake-only stream. To
    avoid this, we split the wake minutes into contiguous blocks (where the
    PAXSSNMP step between consecutive wake minutes equals SAMPLES_PER_MIN) and
    compute bouts per block, then concatenate the per-block bout lists.
    """
    out = {'SEQN': int(pdf['SEQN'].iloc[0])}
    wake = pdf[pdf['PAXPREDM'] == WAKE].copy()
    if len(wake) == 0:
        return None

    wake = wake.sort_values('PAXSSNMP').reset_index(drop=True)
    mims = wake['PAXMTSM'].values
    state3 = np.where(mims < SED_CUT, 0, np.where(mims >= MVPA_CUT, 2, 1))
    state  = np.where(state3 == 0, 0, 1)

    # identify contiguous wake blocks (gap > 1 minute indicates a sleep boundary)
    samp = wake['PAXSSNMP'].values.astype('int64')
    if len(samp) > 1:
        gaps = np.diff(samp)
        block_breaks = np.where(gaps > SAMPLES_PER_MIN)[0] + 1
    else:
        block_breaks = np.array([], dtype=int)
    block_bounds = np.concatenate([[0], block_breaks, [len(state)]])

    bout_dfs = []
    for i in range(len(block_bounds) - 1):
        lo, hi = block_bounds[i], block_bounds[i + 1]
        if hi - lo > 0:
            bout_dfs.append(find_bouts(state[lo:hi]))
    if bout_dfs:
        bouts = pd.concat(bout_dfs, ignore_index=True)
    else:
        bouts = pd.DataFrame(columns=['state', 'length'])

    if len(bouts) == 0:
        return None

    sed = bouts[bouts['state'] == 0]['length'].values
    act = bouts[bouts['state'] == 1]['length'].values
    n_min = len(state)

    # transitions: count within each contiguous block (don't span sleep)
    s_to_a = a_to_s = 0
    n_act_pred = n_sed_pred = 0
    for i in range(len(block_bounds) - 1):
        lo, hi = block_bounds[i], block_bounds[i + 1]
        if hi - lo <= 1:
            continue
        seg = state[lo:hi]
        s_to_a += int(((seg[:-1] == 0) & (seg[1:] == 1)).sum())
        a_to_s += int(((seg[:-1] == 1) & (seg[1:] == 0)).sum())
        n_act_pred += int((seg[:-1] == 1).sum())
        n_sed_pred += int((seg[:-1] == 0).sum())
    astp = a_to_s / n_act_pred if n_act_pred > 0 else np.nan
    satp = s_to_a / n_sed_pred if n_sed_pred > 0 else np.nan

    out['ASTP_minute']      = astp
    out['SATP_minute']      = satp
    out['n_transitions_total'] = int(a_to_s + s_to_a)
    out['n_act_to_sed_total']  = int(a_to_s)
    out['n_sed_to_act_total']  = int(s_to_a)

    wake_hrs = n_min / 60.0
    out['transitions_per_wake_hour'] = out['n_transitions_total'] / wake_hrs if wake_hrs > 0 else np.nan
    out['act_to_sed_per_wake_hour']  = a_to_s / wake_hrs if wake_hrs > 0 else np.nan
    out['sed_to_act_per_wake_hour']  = s_to_a / wake_hrs if wake_hrs > 0 else np.nan

    def desc(arr, prefix):
        return {
            f'{prefix}_n':      len(arr),
            f'{prefix}_mean':   float(np.mean(arr))   if len(arr) else np.nan,
            f'{prefix}_median': float(np.median(arr)) if len(arr) else np.nan,
            f'{prefix}_sd':     float(np.std(arr))    if len(arr) else np.nan,
            f'{prefix}_p90':    float(np.percentile(arr, 90)) if len(arr) else np.nan,
            f'{prefix}_max':    int(np.max(arr))      if len(arr) else np.nan,
        }
    out.update(desc(act, 'act_bout'))
    out.update(desc(sed, 'sed_bout'))

    if len(sed) > 0:
        out['frac_sed_in_long_bouts']      = float(sed[sed >= 30].sum()) / float(sed.sum())
        out['frac_sed_in_very_long_bouts'] = float(sed[sed >= 60].sum()) / float(sed.sum())
    else:
        out['frac_sed_in_long_bouts'] = np.nan
        out['frac_sed_in_very_long_bouts'] = np.nan

    if len(sed) >= 50:
        for d, h in zip(HAZ_DURATIONS, hazards(sed)):
            out[f'sed_hazard_at_{d}min'] = h
    else:
        for d in HAZ_DURATIONS:
            out[f'sed_hazard_at_{d}min'] = np.nan
    if len(act) >= 50:
        for d, h in zip(HAZ_DURATIONS, hazards(act)):
            out[f'act_hazard_at_{d}min'] = h
    else:
        for d in HAZ_DURATIONS:
            out[f'act_hazard_at_{d}min'] = np.nan

    # time-of-day decomposition using REAL clock hour (needs PAXFTIME).
    # Use the same per-block split as the global transition counting so we
    # don't introduce spurious transitions across sleep gaps when minutes
    # from different days are concatenated within a time-of-day window.
    if not np.isnan(ftime_min) and 'PAXSSNMP' in wake.columns:
        hr = real_hour_of_day(wake['PAXSSNMP'].values.astype('int64'), ftime_min)
        for label, lo, hi in [('morn', 6, 12), ('aft', 12, 18), ('eve', 18, 24)]:
            mask = (hr >= lo) & (hr < hi)
            if mask.sum() < 30:
                out[f'astp_{label}'] = np.nan
                out[f'satp_{label}'] = np.nan
                continue
            # accumulate transitions within each contiguous block restricted to window
            ats = sta = 0
            n_act = n_sed = 0
            for i in range(len(block_bounds) - 1):
                blo, bhi = block_bounds[i], block_bounds[i + 1]
                if bhi - blo <= 1:
                    continue
                block_mask = mask[blo:bhi]
                if block_mask.sum() < 2:
                    continue
                seg = state[blo:bhi][block_mask]
                ats += int(((seg[:-1] == 1) & (seg[1:] == 0)).sum())
                sta += int(((seg[:-1] == 0) & (seg[1:] == 1)).sum())
                n_act += int((seg[:-1] == 1).sum())
                n_sed += int((seg[:-1] == 0).sum())
            out[f'astp_{label}'] = ats / n_act if n_act > 0 else np.nan
            out[f'satp_{label}'] = sta / n_sed if n_sed > 0 else np.nan

    out['n_wake_min_total'] = n_min
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--paxmin', required=True)
    ap.add_argument('--paxhd',  required=True, help='for PAXFTIME (start clock time per SEQN)')
    ap.add_argument('--seqn',   required=True, help='analytic_seqn_list.csv')
    ap.add_argument('--out',    required=True)
    ap.add_argument('--chunksize', type=int, default=500_000)
    args = ap.parse_args()

    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / 'fragmentation_features.csv'

    targets = set(pd.read_csv(args.seqn)['SEQN'].astype(int))
    print(f'Target SEQNs: n = {len(targets)}')

    # Read PAXHD for per-SEQN start time
    print(f'Reading {args.paxhd}')
    ph = pd.read_sas(args.paxhd, format='xport')
    ph['SEQN'] = ph['SEQN'].astype('int64')
    ftime_col = 'PAXFTIME' if 'PAXFTIME' in ph.columns else 'PAXFTM'
    if ftime_col not in ph.columns:
        raise KeyError(f'No PAXFTIME/PAXFTM in PAXHD; columns: {list(ph.columns)}')
    ph['ftime_min'] = ph[ftime_col].apply(parse_ftime)
    seqn_to_ftime = dict(zip(ph['SEQN'], ph['ftime_min']))

    print(f'Streaming {args.paxmin}')
    needed = ['SEQN', 'PAXSSNMP', 'PAXMTSM', 'PAXPREDM']
    buffers = {}
    chunk_count = 0
    with pd.read_sas(args.paxmin, format='xport', chunksize=args.chunksize) as reader:
        for chunk in reader:
            chunk_count += 1
            chunk['PAXPREDM'] = to_int(chunk['PAXPREDM'])
            chunk['SEQN'] = pd.to_numeric(chunk['SEQN'], errors='coerce').astype('Int64')
            chunk = chunk[chunk['SEQN'].isin(targets)]
            if len(chunk) == 0:
                continue
            for seqn, person_chunk in chunk.groupby('SEQN'):
                buffers.setdefault(int(seqn), []).append(person_chunk[needed])
            if chunk_count % 10 == 0:
                print(f'  chunk {chunk_count}: {len(buffers)} SEQNs buffered')
    print(f'Read {chunk_count} chunks; {len(buffers)} SEQNs')

    results = []
    for i, (seqn, chunks) in enumerate(buffers.items(), 1):
        pdf = pd.concat(chunks, ignore_index=True).sort_values('PAXSSNMP').reset_index(drop=True)
        ftime_min = seqn_to_ftime.get(seqn, np.nan)
        f = features_for_person(pdf, ftime_min)
        if f is not None:
            results.append(f)
        if i % 500 == 0:
            print(f'  processed {i}/{len(buffers)} ({len(results)} successful)')

    print(f'\nDone: features computed for {len(results)} of {len(buffers)} participants')
    if not results:
        return

    df_out = pd.DataFrame(results)
    df_out.to_csv(out_path, index=False)
    print(f'Wrote {out_path} (shape: {df_out.shape})')


if __name__ == '__main__':
    main()
