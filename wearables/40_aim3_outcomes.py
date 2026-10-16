"""Aim 3 step 1: clinical outcome table for linked CNV carriers (visit 1). Controls were not psychiatrically assessed.

Sources: G2MH REDCap export (Aug 2026; redcap_rarecnv.csv), CNB composites from A. Srinivasan's cnb_model_data.tsv
(cnb_acc_comp / cnb_spd_comp: equal weight per test, oriented higher = better), joined by rarecnv_id.
Codings (Feb-2026 dictionary): psy_asd 0/1; psy_id 0 none, 1-3 mild-severe, 4 GDD without IQ test, 5 borderline, 6 profound;
corepsych_*_dx 0 no / 1 yes / 2 subthreshold; psychosis_group_status Definitive / Putative Non-PS, Subthreshold, PS disorder.
Outcome definitions (binary: 1 vs 0, other levels excluded)
  asd        psy_asd 1 vs 0
  adhd       corepsych_adhd_dx 1 vs 0 (subthreshold excluded)
  id         psy_id in {1,2,3,4,6} vs 0 (borderline 5 excluded)
  psychosis  psychosis_group_status Subthreshold or PS disorder vs Definitive / Putative Non-PS
  anxiety    corepsych_anx_dx 1 vs 0;  depression corepsych_dep_dx 1 vs 0;  ocd corepsych_ocd_dx 1 vs 0
Continuous: fsiq (smry_iq_fsiq), cnb_acc, cnb_spd.
Usage: python -I 40_aim3_outcomes.py <analysis_root> <redcap_csv> <cnb_model_data.tsv> <out_csv>
"""
import sys, csv
import numpy as np, pandas as pd
ROOT, RC, CNB, OUT = sys.argv[1:5]
L = pd.read_csv(f'{ROOT}/phase0/linkage.csv')
L = L[(L.link_status == 'linked') & (L.visit == 1) & L.genotype.isin(['16pDel', '16pDup', '22qDel', '22qDup'])]
ids = set(L.rarecnv_id.dropna())
F = ['psy_asd', 'psy_id', 'corepsych_adhd_dx', 'corepsych_anx_dx', 'corepsych_dep_dx', 'corepsych_ocd_dx', 'psychosis_group_status', 'smry_iq_fsiq']
rows = []
with open(RC, newline='', encoding='utf-8', errors='replace') as fh:
    r = csv.reader(fh); h = next(r); ix = {c: h.index(c) for c in ['rarecnv_id'] + F}
    for x in r:
        if x[ix['rarecnv_id']] in ids: rows.append({c: x[i] for c, i in ix.items()})
R = pd.DataFrame(rows).replace({'': np.nan, 'NA': np.nan, 'NI': np.nan, 'UNK': np.nan, 'M': np.nan})
first = lambda s: s.dropna().iloc[0] if s.notna().any() else np.nan
R = R.groupby('rarecnv_id').agg(first)
num = lambda s: pd.to_numeric(s, errors='coerce')
O = pd.DataFrame(index=R.index)
bin10 = lambda s, yes, no: np.where(s.isin(yes), 1.0, np.where(s.isin(no), 0.0, np.nan))
O['asd'] = bin10(R.psy_asd, ['1'], ['0'])
O['adhd'] = bin10(R.corepsych_adhd_dx, ['1'], ['0'])
O['id'] = bin10(R.psy_id, ['1', '2', '3', '4', '6'], ['0'])
O['psychosis'] = bin10(R.psychosis_group_status, ['Subthreshold Psychotic', 'Psychosis Spectrum Disorder (PS)'], ['Definitive Non-PS', 'Putative Non-PS'])
O['anxiety'] = bin10(R.corepsych_anx_dx, ['1'], ['0'])
O['depression'] = bin10(R.corepsych_dep_dx, ['1'], ['0'])
O['ocd'] = bin10(R.corepsych_ocd_dx, ['1'], ['0'])
O['fsiq'] = num(R.smry_iq_fsiq)
C = pd.read_csv(CNB, sep='\t', usecols=['rarecnv_id', 'cnb_acc_comp', 'cnb_spd_comp']).drop_duplicates('rarecnv_id').set_index('rarecnv_id')
O = O.join(C.rename(columns={'cnb_acc_comp': 'cnb_acc', 'cnb_spd_comp': 'cnb_spd'}), how='left')
O = L[['empatica_id', 'rarecnv_id', 'genotype']].merge(O, left_on='rarecnv_id', right_index=True, how='left')
O.to_csv(OUT, index=False)
print('carriers', len(O))
for c in O.columns[3:]:
    s = O[c]
    print(f'{c:11s} n={s.notna().sum():3d}', (f'cases={int(s.sum())}' if s.dropna().isin([0, 1]).all() else f'mean={s.mean():.2f} sd={s.std():.2f}'),
          O.groupby('genotype')[c].apply(lambda x: f'{x.notna().sum()}' + (f'/{int(x.sum())}' if x.dropna().isin([0, 1]).all() else '')).to_dict())
