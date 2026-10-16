"""Phase 0.5: per-participant covariates for the actigraphy analyses.

Source: G2MH REDCap export (Aug 2026, Expanse) via rarecnv_id from the linkage table.
  age at recording = (median date of the participant's wrist recording) - (1 July of year of birth)
                     [year of birth from dem_yob or dem_yob_fam; +/- 0.5 y precision]
                     fallback: age_final + (recording year - enroll_year) when no year of birth
  BMI, BMI Z (WHO/UK90 via whoref, as in G2MH_CNV_growth) from the enrolment measurements, with
                     the gap (years) between enrolment and recording - BMI is NOT measured at recording.
  medications       current medications (meds_current == 1) from the medications instrument. The export has
                     no ATC codes, so each entry is classified by drug name (meds_name, or the RxNorm name of
                     meds_code looked up in RxNav: inputs/rxcui_names.csv) against generic + brand keyword
                     lists below. Recorded at enrolment, not at recording - the gap is carried.
Usage: python -I 06_covariates.py <redcap_csv> <phase0_dir> <rxcui_names.csv>
Writes <phase0>/covariates.csv and prints coverage.
"""
import os, sys, csv
import numpy as np, pandas as pd

RC, P0, RX = sys.argv[1:4]
L = pd.read_csv(os.path.join(P0, 'linkage.csv'))
L = L[L.link_status == 'linked']
ids = set(L.rarecnv_id.dropna())
keep = ['rarecnv_id', 'redcap_event_name', 'redcap_repeat_instrument', 'dem_yob', 'dem_yob_fam', 'dem_enroll_date', 'dem_age',
        'dem_age_fam', 'med_height_m', 'med_height_cm', 'med_weight', 'med_height_m_fam', 'med_height_cm_fam', 'med_weight_fam',
        'meds_current', 'meds_name', 'meds_code']
rows = []
with open(RC, newline='', encoding='utf-8', errors='replace') as fh:
    r = csv.reader(fh); h = next(r); ix = {c: h.index(c) for c in keep if c in h}
    for x in r:
        if x[ix['rarecnv_id']] in ids:
            rows.append({c: x[i] for c, i in ix.items()})
R = pd.DataFrame(rows).replace({'': np.nan, 'NA': np.nan, 'UNK': np.nan, 'M': np.nan, 'NI': np.nan})
num = lambda s: pd.to_numeric(s, errors='coerce')
first = lambda s: s.dropna().iloc[0] if s.notna().any() else np.nan
g = R.groupby('rarecnv_id')
yob = g.dem_yob.agg(first).combine_first(g.dem_yob_fam.agg(first))
enroll = pd.to_datetime(g.dem_enroll_date.agg(first), errors='coerce')
def height(m, cm):
    m, cm = num(m), num(cm)
    h = m * 100 + cm.fillna(0)
    return h.where(m.notna(), cm.where(cm > 50))
H = (height(g.med_height_m.agg(first), g.med_height_cm.agg(first))
     .combine_first(height(g.med_height_m_fam.agg(first), g.med_height_cm_fam.agg(first))))
Wt = num(g.med_weight.agg(first)).combine_first(num(g.med_weight_fam.agg(first)))
ageE = num(g.dem_age.agg(first)).combine_first(num(g.dem_age_fam.agg(first)))
M = R[R.redcap_repeat_instrument == 'medications'].copy()
print('meds_current values:', M.meds_current.value_counts(dropna=False).to_dict())
M = M[M.meds_current.astype(str) == '1']
KW = {
 'med_stimulant': 'methylphenidate ritalin concerta focalin dexmethylphenidate biphentin foquest quillichew daytrana jornay metadate methylin amphetamine adderall vyvanse lisdexamfetamine dexedrine dextroamphetamine mydayis evekeo zenzedi elvanse',
 'med_adhd_nonstim': 'atomoxetine strattera viloxazine qelbree',
 'med_alpha2': 'guanfacine intuniv tenex clonidine kapvay catapres',
 'med_antipsychotic': 'risperidone risperdal aripiprazole abilify quetiapine seroquel olanzapine zyprexa clozapine clozaril haloperidol haldol ziprasidone geodon lurasidone latuda paliperidone invega pipamperone dipiperon brexpiprazole rexulti cariprazine vraylar asenapine chlorpromazine perphenazine loxapine fluphenazine periciazine levomepromazine methotrimeprazine pimozide amisulpride sulpiride',
 'med_lithium': 'lithium',
 'med_antidepressant': 'sertraline zoloft fluoxetine prozac citalopram celexa escitalopram lexapro cipralex paroxetine paxil fluvoxamine luvox venlafaxine effexor desvenlafaxine pristiq duloxetine cymbalta bupropion wellbutrin mirtazapine remeron trazodone desyrel amitriptyline nortriptyline imipramine clomipramine anafranil vortioxetine trintellix vilazodone',
 'med_anxiolytic_hypnotic': 'lorazepam ativan oxazepam clonazepam rivotril klonopin diazepam valium alprazolam xanax buspirone buspar hydroxyzine atarax vistaril zolpidem ambien zopiclone imovane eszopiclone lunesta temazepam',
 'med_melatonin': 'melatonin',
 'med_antiepileptic': 'valproate valproic depakote divalproex epival lamotrigine lamictal levetiracetam keppra topiramate topamax carbamazepine tegretol oxcarbazepine trileptal clobazam onfi frisium lacosamide vimpat ethosuximide phenytoin dilantin gabapentin neurontin pregabalin lyrica zonisamide perampanel brivaracetam epidiolex rufinamide'}
rx = {a: b for a, b in csv.reader(open(RX))}
M['drug'] = (M.meds_name.fillna('') + ' ' + M.meds_code.fillna('').map(lambda c: rx.get(str(c).strip(), ''))).str.lower()
med = pd.DataFrame(index=sorted(ids))
for k, words in KW.items():
    pat = r'\b(' + '|'.join(words.split()) + r')'
    hit = M[M.drug.str.contains(pat, regex=True)].rarecnv_id
    med[k] = med.index.isin(set(hit)).astype(int)
CLASSES = KW
M['classified'] = M.drug.str.contains(r'\b(' + '|'.join(' '.join(KW.values()).split()) + r')', regex=True)
print('current med entries', len(M), '| with any name', int((M.drug.str.strip() != '').sum()), '| in a psychotropic class', int(M.classified.sum()))
med['med_any_psychotropic'] = med[[c for c in CLASSES if c != 'med_melatonin']].max(axis=1)
med['n_current_meds'] = M.groupby('rarecnv_id').size().reindex(med.index).fillna(0).astype(int)
has_meds_form = set(R[R.redcap_repeat_instrument == 'medications'].rarecnv_id)
med['meds_form_present'] = med.index.isin(has_meds_form).astype(int)

# recording date per Empatica ID (median date with data)
days = pd.read_parquet(os.path.join(P0, 'qc', 'daily_qc.parquet'), columns=['empatica_id', 'date'])
rec = days.groupby('empatica_id').date.agg(lambda d: pd.to_datetime(pd.Series(d)).sort_values().iloc[len(d) // 2])
C = L[['empatica_id', 'guid', 'rarecnv_id', 'genotype', 'age_final', 'enroll_year', 'sex_final']].copy()
C['rec_date'] = C.empatica_id.map(rec)
C['yob'] = num(C.rarecnv_id.map(yob))
C['enroll_date'] = C.rarecnv_id.map(enroll)
C['age_at_recording'] = (C.rec_date - pd.to_datetime(C.yob.astype('Int64').astype(str) + '-07-01', errors='coerce')).dt.days / 365.25
fb = C.age_at_recording.isna() & C.age_final.notna()
C.loc[fb, 'age_at_recording'] = num(C.loc[fb, 'age_final']) + (C.loc[fb, 'rec_date'].dt.year - num(C.loc[fb, 'enroll_year']))
C['age_source'] = np.where(fb, 'age_final + years since enrolment', np.where(C.age_at_recording.notna(), 'year of birth', 'missing'))
C['years_enrol_to_recording'] = (C.rec_date - C.enroll_date).dt.days / 365.25
C['age_at_enrolment'] = C.rarecnv_id.map(ageE)
C['height_cm'] = C.rarecnv_id.map(H); C['weight_kg'] = C.rarecnv_id.map(Wt)
C['bmi'] = C.weight_kg / (C.height_cm / 100) ** 2
C.loc[~C.bmi.between(10, 80), 'bmi'] = np.nan
C = C.merge(med, left_on='rarecnv_id', right_index=True, how='left')
C.to_csv(os.path.join(P0, 'covariates.csv'), index=False)
pd.set_option('display.width', 200)
print('linked participants', len(C), '| age source:', C.age_source.value_counts().to_dict())
print('age at recording minus age_final (y): quantiles', (C.age_at_recording - num(C.age_final)).quantile([.05, .5, .95]).round(2).tolist())
print('years enrolment -> recording quantiles', C.years_enrol_to_recording.quantile([.1, .5, .9]).round(2).tolist())
print('BMI available', int(C.bmi.notna().sum()), '| meds form present', int(C.meds_form_present.sum()))
print(C.groupby('genotype')[['age_at_recording', 'med_stimulant', 'med_antipsychotic', 'med_antidepressant', 'med_melatonin', 'med_any_psychotropic']].agg(
    lambda x: round(x.mean(), 2)).to_string())
