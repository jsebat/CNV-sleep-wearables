"""Master per-participant phenotype table (one row per participant with any device data) + data dictionary.

Columns
  identifiers / covariates : study_id, empatica_id (wrist), guid, rarecnv_id, family, genotype, locus, copy_number,
                             canonical_bp, link_status, visit, age_at_recording, sex, site, bmi_enrolment, medication flags
  hb_*      HEADBAND (Dreem EEG) core sleep phenotypes            (29_headband_core.py)        panel = core
  wr_*      WRIST BAND core phenotypes                            (20_phenotypes.py)           panel = core
  wr_p2_*   WRIST BAND second-panel phenotypes                    (25, 27, 28)                 panel = second (exploratory)
  mat_p2_*  MAT (Withings) second-panel phenotypes                (25)                         panel = second (exploratory)
  n_* columns give the nights/days behind each device's values.
R2 repeat visits are kept as separate rows (visit = 2); analyses use visit 1.
Usage: python -I 30_master_table.py <analysis_root> <out_dir>
"""
import os, sys, glob
import pandas as pd

ROOT, OUT = sys.argv[1:3]; os.makedirs(OUT, exist_ok=True)
P0, A4 = os.path.join(ROOT, 'phase0'), os.path.join(ROOT, 'aim4')
L = pd.read_csv(os.path.join(P0, 'linkage.csv')); L = L[~L.test_account]
C = pd.read_csv(os.path.join(P0, 'covariates.csv'))
ids = L[['empatica_id', 'base_id', 'visit', 'guid', 'rarecnv_id', 'rarecnv_family_id', 'genotype', 'locus', 'copy_number', 'canonical_bp',
         'link_status', 'sex_final', 'redcap_data_access_group']].rename(columns={'base_id': 'study_id', 'rarecnv_family_id': 'family',
                                                                                 'sex_final': 'sex', 'redcap_data_access_group': 'site'})
cv = C[['empatica_id', 'age_at_recording', 'bmi', 'meds_form_present', 'med_stimulant', 'med_adhd_nonstim', 'med_alpha2', 'med_antipsychotic',
        'med_antidepressant', 'med_lithium', 'med_anxiolytic_hypnotic', 'med_melatonin', 'med_antiepileptic']].rename(columns={'bmi': 'bmi_enrolment'})
M = ids.merge(cv, on='empatica_id', how='left')
# headband (keyed by study ID; headband-only participants added as extra rows)
HB = pd.read_csv(os.path.join(A4, 'headband_core_participant.csv'))
M = M.merge(HB, on='study_id', how='outer')
M.loc[M.empatica_id.isna(), 'visit'] = 1
# wrist core
core = pd.read_csv(os.path.join(A4, 'phenotypes_participant.csv'))
core = core.rename(columns={c: ('wr_' + c if c != 'empatica_id' else c) for c in core.columns})
M = M.merge(core, on='empatica_id', how='left')
# second panel
p2a = pd.read_csv(os.path.join(A4, 'phase2_hr_mat_participant.csv')).drop(columns=['study_id'], errors='ignore')
p2a = p2a.rename(columns={c: ('mat_p2_' + c[4:] if c.startswith('mat_') else ('mat_p2_' + c if c == 'n_mat_sessions' else 'wr_p2_' + c))
                          for c in p2a.columns if c != 'empatica_id'})
p2b = pd.read_csv(os.path.join(A4, 'phase2_blocks2to5_participant.csv'))
p2b = p2b.rename(columns={c: 'wr_p2_' + c for c in p2b.columns if c != 'empatica_id'})
M = M.merge(p2a, on='empatica_id', how='left').merge(p2b, on='empatica_id', how='left')
eda_f = glob.glob(os.path.join(A4, 'phase2_eda', '*.csv'))
if eda_f:
    E = pd.concat([pd.read_csv(f) for f in eda_f if os.path.getsize(f) > 0], ignore_index=True)
    E[E.columns[2:]] = E[E.columns[2:]].apply(pd.to_numeric, errors='coerce')
    Ep = E.groupby('empatica_id').mean(numeric_only=True).drop(columns=['night'], errors='ignore')
    Ep['n_eda_nights'] = E.groupby('empatica_id').size()
    Ep = Ep.rename(columns={c: 'wr_p2_' + c for c in Ep.columns}).reset_index()
    M = M.merge(Ep, on='empatica_id', how='left')
M = M.sort_values(['genotype', 'study_id', 'visit'], na_position='last')
M.to_csv(os.path.join(OUT, 'G2MH_wearable_phenotypes_per_participant.csv'), index=False)

def device(c):
    return ('headband (Dreem EEG)' if c.startswith('hb_') else 'mat (Withings)' if c.startswith('mat_p2_') else
            'wrist band (Empatica EmbracePlus)' if c.startswith('wr_') else 'identifier / covariate')
def panel(c):
    return 'second (exploratory)' if '_p2_' in c else ('core' if c[:3] in ('hb_', 'wr_') else '')
Dd = pd.DataFrame({'column': M.columns, 'device': [device(c) for c in M.columns], 'panel': [panel(c) for c in M.columns],
                   'n_nonmissing': [int(M[c].notna().sum()) for c in M.columns]})
Dd.to_csv(os.path.join(OUT, 'G2MH_wearable_phenotypes_dictionary_skeleton.csv'), index=False)
print('rows', len(M), '| columns', M.shape[1], '| headband rows', int(M.hb_TST.notna().sum()) if 'hb_TST' in M else 0,
      '| wrist core rows', int(M.wr_tst_h.notna().sum()) if 'wr_tst_h' in M else 0)
print(Dd.groupby(['device', 'panel']).size().to_string())
