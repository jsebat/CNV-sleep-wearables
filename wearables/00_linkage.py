"""Phase 0a: link every Empatica participant ID to GUID, genotype and covariates.

Empatica participant IDs are the G2MH wearable study IDs (site prefix + 4 digits, e.g. SITE0000),
the same IDs the Dreem headband uses; the manifest carries them in `dreem_id`.
`<ID>R2` folders are repeat recordings of <ID> (visit 2).

Genotype grouping follows the sleep-architecture analysis (01_ingest.R
normalise_genotype): 16pDel / 16pDup / 22qDel / 22qDup / Control by final_group_class,
any breakpoint type; double carriers and unknown classes excluded from the five groups.

Usage: python -I 00_linkage.py <empatica_root> <manifest.csv> <out_dir>
"""
import os, re, sys, glob, collections
import pandas as pd

ROOT, MANIFEST, OUT = sys.argv[1:4]
os.makedirs(OUT, exist_ok=True)

# ---- Empatica inventory: one row per participant-day folder ---------------------
rows = []
for day_dir in sorted(glob.glob(os.path.join(ROOT, '1', '1', 'participant_data', '*'))):
    date = os.path.basename(day_dir)
    for p in os.listdir(day_dir):
        pid, _, serial = p.partition('-')
        has_raw = os.path.isdir(os.path.join(day_dir, p, 'raw_data'))
        has_bio = os.path.isdir(os.path.join(day_dir, p, 'digital_biomarkers'))
        rows.append(dict(empatica_id=pid.strip(), date=date, device_sn=serial, has_raw=has_raw, has_bio=has_bio))
inv = pd.DataFrame(rows)
inv.to_csv(os.path.join(OUT, 'empatica_inventory_days.csv'), index=False)

per = inv.groupby('empatica_id').agg(n_days=('date', 'nunique'), first_date=('date', 'min'), last_date=('date', 'max'),
                                     n_raw_days=('has_raw', 'sum'), n_bio_days=('has_bio', 'sum'),
                                     devices=('device_sn', lambda s: ';'.join(sorted(set(s))))).reset_index()
per['test_account'] = per.empatica_id.str.contains('TEST|TROUBLE')
per['base_id'] = per.empatica_id.str.replace(r'R2$', '', regex=True)
per['visit'] = per.empatica_id.str.endswith('R2').map({True: 2, False: 1})
per['site_code'] = per.empatica_id.str.extract(r'^([A-Z]+?)(?:TEST|\d)')[0]

# ---- Manifest --------------------------------------------------------------------
m = pd.read_csv(MANIFEST, dtype=str)
m['dreem_id'] = m.dreem_id.str.strip()
mm = m.dropna(subset=['dreem_id'])
n_guid = mm.groupby('dreem_id').guid.nunique()


def geno5(cls):
    if not isinstance(cls, str):
        return None
    s = re.sub('[^A-Za-z0-9]', '', cls).upper()
    if ('16P' in s and '22Q' in s) or ('DEL' in s and 'DUP' in s) or '&' in cls:
        return None
    for k, lab in (('16PDEL', '16pDel'), ('16PDUP', '16pDup'), ('22QDEL', '22qDel'), ('22QDUP', '22qDup')):
        if k[:3] in s and k[3:] in s:
            return lab
    return 'Control' if s in ('NOCNV',) else None


cols = ['guid', 'rarecnv_id', 'rarecnv_family_id', 'final_group_class', 'final_group_16p_type', 'final_group_22q_type',
        'age_final', 'sex_final', 'redcap_data_access_group', 'enroll_year', 'proband', 'ancestry_FINAL']
out = []
for _, r in per.iterrows():
    rec = r.to_dict()
    hit = mm[mm.dreem_id == r.base_id]
    if r.test_account:
        rec['link_status'] = 'test account'
    elif hit.empty:
        rec['link_status'] = 'not in manifest'
    elif n_guid.get(r.base_id, 0) > 1:
        rec['link_status'] = 'AMBIGUOUS: dreem_id on %d GUIDs (%s)' % (n_guid[r.base_id], ', '.join(hit.guid))
    else:
        h = hit.iloc[0]
        rec.update({c: h.get(c) for c in cols})
        g = geno5(h.final_group_class)
        rec['genotype'] = g
        rec['link_status'] = 'linked' if g else ('linked, no genotype' if not isinstance(h.final_group_class, str)
                                                 else 'linked, excluded genotype (%s)' % h.final_group_class)
    out.append(rec)
L = pd.DataFrame(out)
L['locus'] = L.genotype.map({'16pDel': '16p', '16pDup': '16p', '22qDel': '22q', '22qDup': '22q'})
L['copy_number'] = L.genotype.map({'16pDel': 1, '22qDel': 1, '16pDup': 3, '22qDup': 3, 'Control': 2})
L['canonical_bp'] = ((L.genotype.isin(['16pDel', '16pDup']) & (L.final_group_16p_type == 'BP4-5')) |
                     (L.genotype.isin(['22qDel', '22qDup']) & (L.final_group_22q_type == 'A-D')) |
                     (L.genotype == 'Control'))
L.to_csv(os.path.join(OUT, 'linkage.csv'), index=False)

real = L[~L.test_account]
print('Empatica IDs:', len(L), '| real:', len(real), '| visits-2:', (real.visit == 2).sum())
print(real.link_status.str.replace(r'\(.*', '', regex=True).value_counts().to_string())
print('\nGenotype x site (linked, visit 1):')
v1 = real[(real.link_status == 'linked') & (real.visit == 1)]
print(pd.crosstab(v1.genotype, v1.site_code, margins=True).to_string())
print('\nUnlinked / problem IDs:')
print(real[real.link_status != 'linked'][['empatica_id', 'n_days', 'link_status']].to_string(index=False))
