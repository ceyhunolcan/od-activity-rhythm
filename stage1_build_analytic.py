"""
Build the analytic dataset for the OD x accelerometer paper.

Merges NHANES 2013-14 XPT files and applies inclusion criteria:
  - age >= 40 (PST eligibility)
  - PST completed (CSXEXSTS == 1)
  - not pregnant/breastfeeding at exam (CSQ241 != 1 in CSX_H)
  - no acute nasal symptoms at exam (CSQ260a/d/g/i/n in CSX_H)
  - >= 4 valid accelerometer days

Writes: analytic_full.csv, analytic_seqn_list.csv, attrition_log.csv

Run from the directory containing the NHANES XPT files:
    python stage1_build_analytic.py
"""

import pandas as pd
import numpy as np
from pathlib import Path
import sys


FILES = {
    'DEMO_H':   'DEMO_H.xpt',
    'BMX_H':    'BMX_H.xpt',
    'BPX_H':    'BPX_H.xpt',
    'CSX_H':    'CSX_H.xpt',     # chemosensory exam (PST + nasal symptoms + CSQ241 pregnancy)
    'CSQ_H':    'CSQ_H.xpt',     # household chemosensory Q (CSQ240 head injury, CSQ200 cold/flu)
    'SMQ_H':    'SMQ_H.xpt',
    'DIQ_H':    'DIQ_H.xpt',
    'GHB_H':    'GHB_H.xpt',
    'MCQ_H':    'MCQ_H.xpt',
    'BPQ_H':    'BPQ_H.xpt',
    'DPQ_H':    'DPQ_H.xpt',
    'HUQ_H':    'HUQ_H.xpt',
    'RXQ_RX_H': 'RXQ_RX_H.xpt',
}


# 8-item Pocket Smell Test: variable -> (label, correct response code)
# Verified against CSX_H codebook:
# https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2013/DataFiles/CSX_H.htm
PST_ITEMS = {
    'CSXCHOOD': ('chocolate',   2),   # 1=Lemon  2=Chocolate  3=Smoke   4=Black pepper
    'CSXSBOD':  ('strawberry',  1),   # 1=Strawberry 2=Garlic 3=Leather 4=Gasoline
    'CSXSMKOD': ('smoke',       3),   # 1=Garlic 2=Grass  3=Smoke 4=Peach
    'CSXLEAOD': ('leather',     3),   # 1=Mint   2=Flower 3=Leather 4=Apple
    'CSXSOAOD': ('soap',        1),   # 1=Soap   2=Black Pepper 3=Leather 4=Peanut
    'CSXGRAOD': ('grape',       2),   # 1=Gasoline 2=Grape 3=Rose 4=Peanut
    'CSXONOD':  ('onion',       3),   # 1=Chocolate 2=Strawberry 3=Onion 4=Fruit Punch
    'CSXNGSOD': ('natural_gas', 4),   # 1=Orange 2=Cinnamon 3=Cola 4=Natural Gas
}

FOOD_ITEMS      = ['CSXCHOOD', 'CSXSBOD', 'CSXGRAOD', 'CSXONOD']
WARNING_ITEMS   = ['CSXSMKOD', 'CSXNGSOD']
HOUSEHOLD_ITEMS = ['CSXLEAOD', 'CSXSOAOD']

# Nasal symptoms at chemosensory exam, in CSX_H (acute URI exclusion).
# NOTE: paper #1 excluded 237 for "current upper-respiratory illness". The exact
# subset of nasal-symptom variables that produces n=237 depends on the
# definition (any vs. specific subset). The default below is the union of the
# five non-"none" symptoms, which is likely *broader* than paper #1's criterion.
# Adjust NASAL_SYMPTOMS to match your inclusion-flow target.
NASAL_SYMPTOMS = ['CSQ260a', 'CSQ260d', 'CSQ260g', 'CSQ260i', 'CSQ260n']

# Comorbidity series in MCQ_H -- doctor ever told (1 = Yes)
COMORBID_VARS = ['MCQ160A', 'MCQ160B', 'MCQ160C', 'MCQ160D', 'MCQ160E', 'MCQ160F',
                 'MCQ160G', 'MCQ160K', 'MCQ160L', 'MCQ160M', 'MCQ160N', 'MCQ160P',
                 'MCQ220']

# PHQ-9 items in DPQ_H
PHQ_ITEMS = ['DPQ010', 'DPQ020', 'DPQ030', 'DPQ040', 'DPQ050',
             'DPQ060', 'DPQ070', 'DPQ080', 'DPQ090']


def read_xpt(name, path):
    if not Path(path).exists():
        sys.exit(f'missing {name}: {path}')
    return pd.read_sas(path, format='xport')


def score_pst(row):
    """Return # correctly identified odors (0-8), NaN if any item missing."""
    correct = 0
    for var, (_, expected) in PST_ITEMS.items():
        v = row.get(var)
        if pd.isna(v):
            return np.nan
        if v == expected:
            correct += 1
    return correct


def subtype_deficit(row, items, threshold):
    """1 if participant missed >= threshold of `items`; 0 otherwise; NaN if any missing."""
    missed = 0
    for var in items:
        v = row.get(var)
        if pd.isna(v):
            return np.nan
        expected = PST_ITEMS[var][1]
        if v != expected:
            missed += 1
    return int(missed >= threshold)


def main():
    print('Stage 1: build analytic dataset')

    df = {}
    for name, path in FILES.items():
        df[name] = read_xpt(name, path)
        print(f'  {name}: {len(df[name]):>6} rows')

    demo = df['DEMO_H'][['SEQN', 'RIAGENDR', 'RIDAGEYR', 'RIDRETH3', 'DMDEDUC2',
                         'INDFMPIR', 'WTMEC2YR', 'SDMVSTRA', 'SDMVPSU']].copy()
    demo['SEQN'] = demo['SEQN'].astype('int64')

    log = [('MEC examined', len(demo))]
    demo = demo[demo['RIDAGEYR'] >= 40].copy()
    log.append(('age >= 40', len(demo)))

    csx = df['CSX_H'].copy()
    csx['SEQN'] = csx['SEQN'].astype('int64')

    missing_pst = [v for v in PST_ITEMS if v not in csx.columns]
    if missing_pst:
        sys.exit(f'CSX_H is missing expected PST variables: {missing_pst}\n'
                 f'CSX columns present: '
                 f'{[c for c in csx.columns if c.startswith("CSX")][:20]}')

    csx['PST_correct']      = csx.apply(score_pst, axis=1)
    csx['food_deficit']     = csx.apply(lambda r: subtype_deficit(r, FOOD_ITEMS, 2), axis=1)
    csx['warning_deficit']  = csx.apply(lambda r: subtype_deficit(r, WARNING_ITEMS, 1), axis=1)
    csx['household_deficit']= csx.apply(lambda r: subtype_deficit(r, HOUSEHOLD_ITEMS, 1), axis=1)

    csx_done = csx[csx['CSXEXSTS'] == 1].copy()
    pst_cols = ['SEQN', 'PST_correct', 'food_deficit', 'warning_deficit', 'household_deficit']
    nasal_cols = [c for c in NASAL_SYMPTOMS if c in csx_done.columns]
    preg_col   = ['CSQ241'] if 'CSQ241' in csx_done.columns else []

    demo = demo.merge(csx_done[pst_cols + nasal_cols + preg_col], on='SEQN', how='inner')
    log.append(('completed PST', len(demo)))

    if 'CSQ241' in demo.columns:
        demo = demo[demo['CSQ241'] != 1].copy()
    log.append(('not pregnant/breastfeeding', len(demo)))

    if nasal_cols:
        has_symptom = pd.Series(False, index=demo.index)
        for c in nasal_cols:
            has_symptom |= (demo[c] == 1)
        demo = demo[~has_symptom].copy()
    log.append(('no acute nasal symptoms', len(demo)))

    # accelerometer valid-day flag (precomputed by an earlier minute-level stage)
    feats_path = Path('paxmin_features.csv')
    if feats_path.exists():
        feats = pd.read_csv(feats_path)
        feats['SEQN'] = feats['SEQN'].astype('int64')
        if 'meets_4day_inclusion' not in feats.columns:
            sys.exit('paxmin_features.csv lacks meets_4day_inclusion column')
        demo = demo.merge(feats[['SEQN', 'valid_days', 'meets_4day_inclusion']],
                          on='SEQN', how='left')
        demo = demo[demo['meets_4day_inclusion'] == 1].copy()
        log.append(('>=4 valid accel days', len(demo)))
    else:
        print('  WARNING: paxmin_features.csv not found; skipping >=4-day inclusion.')
        print('           Final n will not match the published n=2,327.')
        log.append(('>=4 valid accel days (SKIPPED)', len(demo)))

    keep = demo.copy()

    # BMI
    bmx = df['BMX_H'][['SEQN', 'BMXBMI']].copy()
    bmx['SEQN'] = bmx['SEQN'].astype('int64')
    keep = keep.merge(bmx.rename(columns={'BMXBMI': 'bmi'}), on='SEQN', how='left')

    # blood pressure -- mean of up to 3 readings
    bpx = df['BPX_H'][['SEQN', 'BPXSY1', 'BPXSY2', 'BPXSY3', 'BPXDI1', 'BPXDI2', 'BPXDI3']].copy()
    bpx['SEQN'] = bpx['SEQN'].astype('int64')
    bpx['sbp'] = bpx[['BPXSY1', 'BPXSY2', 'BPXSY3']].mean(axis=1)
    bpx['dbp'] = bpx[['BPXDI1', 'BPXDI2', 'BPXDI3']].mean(axis=1)
    keep = keep.merge(bpx[['SEQN', 'sbp', 'dbp']], on='SEQN', how='left')

    # smoking: SMQ020 ever-smoked-100-cigs (1=Yes, 2=No, 7/9 sentinel)
    #          SMQ040 currently smoke (1=Every day, 2=Some days, 3=Not at all, 7/9 sentinel)
    smq = df['SMQ_H'][['SEQN', 'SMQ020', 'SMQ040']].copy()
    smq['SEQN'] = smq['SEQN'].astype('int64')
    smq['smoker_status'] = np.nan
    smq.loc[smq['SMQ020'] == 2, 'smoker_status'] = 'never'
    smq.loc[(smq['SMQ020'] == 1) & (smq['SMQ040'] == 3),  'smoker_status'] = 'former'
    smq.loc[(smq['SMQ020'] == 1) & (smq['SMQ040'].isin([1, 2])), 'smoker_status'] = 'current'
    keep = keep.merge(smq[['SEQN', 'smoker_status']], on='SEQN', how='left')

    # diabetes: DIQ010 (1=Yes 2=No 7=Refused 9=DK -- treat 7/9/missing as NaN)
    diq = df['DIQ_H'][['SEQN', 'DIQ010']].copy()
    diq['SEQN'] = diq['SEQN'].astype('int64')
    diq['diabetes'] = np.where(diq['DIQ010'] == 1, 1,
                       np.where(diq['DIQ010'] == 2, 0, np.nan))
    keep = keep.merge(diq[['SEQN', 'diabetes']], on='SEQN', how='left')

    # HbA1c
    ghb = df['GHB_H'][['SEQN', 'LBXGH']].copy()
    ghb['SEQN'] = ghb['SEQN'].astype('int64')
    keep = keep.merge(ghb.rename(columns={'LBXGH': 'hba1c'}), on='SEQN', how='left')

    # comorbidity count
    mcq = df['MCQ_H'].copy()
    mcq['SEQN'] = mcq['SEQN'].astype('int64')
    present_comorbid = [c for c in COMORBID_VARS if c in mcq.columns]
    if not present_comorbid:
        sys.exit('No MCQ160* / MCQ220 variables found in MCQ_H')
    mcq['comorbidity_count'] = (mcq[present_comorbid] == 1).sum(axis=1).astype(int)
    keep = keep.merge(mcq[['SEQN', 'comorbidity_count']], on='SEQN', how='left')

    # hypertension: BPQ020 (1=Yes 2=No 7/9 sentinel)
    bpq = df['BPQ_H'][['SEQN', 'BPQ020']].copy()
    bpq['SEQN'] = bpq['SEQN'].astype('int64')
    bpq['hypertension'] = np.where(bpq['BPQ020'] == 1, 1,
                          np.where(bpq['BPQ020'] == 2, 0, np.nan))
    keep = keep.merge(bpq[['SEQN', 'hypertension']], on='SEQN', how='left')

    # PHQ-9: treat 7 and 9 as missing, require >= 7 of 9 items
    dpq = df['DPQ_H'].copy()
    dpq['SEQN'] = dpq['SEQN'].astype('int64')
    present_phq = [c for c in PHQ_ITEMS if c in dpq.columns]
    dpq_clean = dpq[present_phq].replace({7: np.nan, 9: np.nan})
    dpq['phq9'] = dpq_clean.sum(axis=1, min_count=7)
    keep = keep.merge(dpq[['SEQN', 'phq9']], on='SEQN', how='left')

    # head injury (CSQ240), chronic sinus (CSQ260), self-report smell (CSQ010) -- all in CSQ_H
    # All use 1=Yes 2=No 7/9 sentinel; preserve NaN.
    csq = df['CSQ_H'].copy()
    csq['SEQN'] = csq['SEQN'].astype('int64')
    csq_keep = ['SEQN']

    def _yes_no(series):
        return np.where(series == 1, 1, np.where(series == 2, 0, np.nan))

    if 'CSQ240' in csq.columns:
        csq['head_injury'] = _yes_no(csq['CSQ240'])
        csq_keep.append('head_injury')
    if 'CSQ260' in csq.columns:
        csq['sinus'] = _yes_no(csq['CSQ260'])
        csq_keep.append('sinus')
    if 'CSQ010' in csq.columns:
        csq['self_smell_problem'] = _yes_no(csq['CSQ010'])
        csq_keep.append('self_smell_problem')
    keep = keep.merge(csq[csq_keep], on='SEQN', how='left')

    # number of prescription medications
    rx = df['RXQ_RX_H'].copy()
    rx['SEQN'] = rx['SEQN'].astype('int64')
    rx_count = rx.groupby('SEQN').size().reset_index(name='nmedications')
    keep = keep.merge(rx_count, on='SEQN', how='left')
    keep['nmedications'] = keep['nmedications'].fillna(0).astype(int)

    # self-rated health (HUQ010) for frailty composite
    huq = df['HUQ_H'][['SEQN', 'HUQ010']].copy()
    huq['SEQN'] = huq['SEQN'].astype('int64')
    keep = keep.merge(huq, on='SEQN', how='left')

    # rename / derive analysis columns
    keep['female']    = (keep['RIAGENDR'] == 2).astype(int)
    keep['age']       = keep['RIDAGEYR']
    keep['race_eth']  = keep['RIDRETH3']
    keep['education'] = keep['DMDEDUC2']
    keep['pir']       = keep['INDFMPIR']
    # od_binary preserves NaN -- a participant with any missing PST item should not
    # be silently coded as non-OD. The PST_correct score itself returned NaN for
    # those cases (see score_pst).
    keep['od_binary'] = np.where(keep['PST_correct'].isna(), np.nan,
                                 (keep['PST_correct'] <= 5).astype(int))

    # smoker_status comes through as Python None / NaN mixed with strings; coerce
    # to a proper pandas Categorical so downstream R via reticulate/feather and
    # in-Python statsmodels handle it consistently.
    if 'smoker_status' in keep.columns:
        keep['smoker_status'] = pd.Categorical(
            keep['smoker_status'], categories=['never', 'former', 'current'])

    out_cols = ['SEQN', 'WTMEC2YR', 'SDMVSTRA', 'SDMVPSU',
                'age', 'female', 'race_eth', 'education', 'pir',
                'PST_correct', 'od_binary',
                'food_deficit', 'warning_deficit', 'household_deficit',
                'bmi', 'sbp', 'dbp', 'smoker_status', 'diabetes', 'hba1c',
                'comorbidity_count', 'hypertension', 'phq9',
                'head_injury', 'sinus', 'self_smell_problem',
                'nmedications', 'HUQ010']
    out = keep[[c for c in out_cols if c in keep.columns]].copy()

    out.to_csv('analytic_full.csv', index=False)
    pd.DataFrame({'SEQN': out['SEQN']}).to_csv('analytic_seqn_list.csv', index=False)
    pd.DataFrame(log, columns=['step', 'n']).to_csv('attrition_log.csv', index=False)

    print()
    print('Attrition:')
    for step, n in log:
        print(f'  {step:>30}: n = {n}')
    print()
    print(f'Final analytic n = {len(out)}')
    if (out['od_binary'] == 1).any():
        print(f'  OD (PST <= 5): {(out["od_binary"]==1).sum()} '
              f'({(out["od_binary"]==1).mean()*100:.1f}%)')
    print('Wrote analytic_full.csv, analytic_seqn_list.csv, attrition_log.csv')


if __name__ == '__main__':
    main()
