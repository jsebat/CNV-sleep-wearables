"""Aim 3 prediction layer (WRIST BAND core phenotypes): does actigraphy improve prediction of clinical features beyond CNV group?
Carriers only. Baseline features = genotype (one-hot), age, age^2, sex, site; full = baseline + the 31 wrist core phenotypes.
Same nested CV as Aim 2: outer 5-fold x 20 repeats, inner 3-fold, family-grouped (stratified for binary outcomes);
XGBoost for both models (the inner loop picks depth 1-3 and 100/300 trees). Metric: AUC (binary) or R^2 (continuous).
Gain = full - baseline (mean over repeats). Significance: the wrist feature block is permuted across participants
(rows shuffled together, so the features' joint structure is kept) 500 times, each a full nested-CV run (1 repeat);
p = (1 + #null gain >= observed) / (1 + n_perm).
Env NO_MEDS=1: carriers on any psychotropic, or without a medication form, removed.
Usage: python -I 42_aim3_predict.py <analysis_root> <out_dir> <n_jobs> [n_perm]
"""
import os, sys, warnings
import numpy as np, pandas as pd, xgboost as xgb
from joblib import Parallel, delayed
from sklearn.model_selection import StratifiedGroupKFold, GroupKFold
from sklearn.metrics import roc_auc_score, r2_score
warnings.filterwarnings('ignore')
ROOT, OUT, NJ = sys.argv[1], sys.argv[2], int(sys.argv[3]); NPERM = int(sys.argv[4]) if len(sys.argv) > 4 else 500
os.makedirs(OUT, exist_ok=True)
NREP = int(os.environ.get('NREP', 20)); GRID = [dict(max_depth=d, n_estimators=n) for d in (1, 2, 3) for n in (100, 300)]
FIXED = dict(learning_rate=0.05, subsample=0.8, colsample_bytree=0.8, min_child_weight=1, reg_lambda=1.0, tree_method='hist', n_jobs=1, verbosity=0)
BIN = ['asd', 'adhd', 'id', 'psychosis', 'anxiety', 'depression']; CONT = ['fsiq', 'cnb_acc', 'cnb_spd']   # ocd: 10 cases, association only

M = pd.read_csv(f'{ROOT}/aim4/master/G2MH_wearable_phenotypes_per_participant.csv')
M = M[(M.visit == 1) & (M.link_status == 'linked') & M.genotype.isin(['16pDel', '16pDup', '22qDel', '22qDup']) & M.wr_tst_h.notna() & M.age_at_recording.notna()]
O = pd.read_csv(f'{ROOT}/aim3/outcomes.csv').drop(columns=['rarecnv_id', 'genotype'])
D = M.merge(O, on='empatica_id').reset_index(drop=True)
D['family'] = D.family.where(D.family.notna() & (D.family.astype(str) != ''), D.study_id).astype(str)
MEDS = ['med_stimulant', 'med_adhd_nonstim', 'med_alpha2', 'med_antipsychotic', 'med_antidepressant', 'med_lithium', 'med_anxiolytic_hypnotic', 'med_antiepileptic']
if os.environ.get('NO_MEDS') == '1':
    D = D[(D.meds_form_present == 1) & (D[MEDS].fillna(0).sum(1) == 0)].reset_index(drop=True); print('NO_MEDS', flush=True)
COUNTS = {'wr_n_valid_days', 'wr_n_valid_nights'}
WR = [c for c in D.columns if c.startswith('wr_') and '_p2_' not in c and c not in COUNTS]
a = D.age_at_recording - 20
BASE = pd.concat([pd.get_dummies(D.genotype, prefix='g', drop_first=True).astype(float), pd.DataFrame({'age': a, 'age2': a ** 2}),
                  (D.sex.astype(str) == sorted(D.sex.astype(str).unique())[0]).astype(float).rename('sex'),
                  pd.get_dummies(D.site.astype(str), prefix='s', drop_first=True).astype(float)], axis=1).to_numpy(float)
WRX = D[WR].to_numpy(float)
print('carriers', len(D), '| wrist features', len(WR), flush=True)

def splits(y, g, k, seed, binary):
    if binary: return StratifiedGroupKFold(k, shuffle=True, random_state=seed).split(y, y, g)
    fams = np.random.default_rng(seed).permutation(np.unique(g)); fold = {f: i % k for i, f in enumerate(fams)}
    fo = np.array([fold[x] for x in g])
    return [(np.where(fo != i)[0], np.where(fo == i)[0]) for i in range(k)]

def model(binary, p): return (xgb.XGBClassifier if binary else xgb.XGBRegressor)(**FIXED, **p)
def score(binary, y, p): return roc_auc_score(y, p) if binary else r2_score(y, p)
def predict(m, X, binary): return m.predict_proba(X)[:, 1] if binary else m.predict(X)

def nested(X, y, g, seed, binary):
    pr = np.full(len(y), np.nan)
    for tr, te in splits(y, g, 5, seed, binary):
        best, bs = 0, -np.inf
        inner = list(splits(y[tr], g[tr], 3, seed + 7, binary))
        for k, p in enumerate(GRID):
            s = [score(binary, y[tr][b], predict(model(binary, p).fit(X[tr][a_], y[tr][a_]), X[tr][b], binary))
                 for a_, b in inner if (not binary or len(np.unique(y[tr][b])) == 2)]
            if s and np.mean(s) > bs: bs, best = np.mean(s), k
        pr[te] = predict(model(binary, GRID[best]).fit(X[tr], y[tr]), X[te], binary)
    return score(binary, y, pr)

def data(outcome):
    ok = D[outcome].notna().to_numpy(); y = D.loc[ok, outcome].to_numpy(float)
    return BASE[ok], WRX[ok], y, D.family.to_numpy()[ok], outcome in BIN

def run(outcome, seed, perm=None):
    B, W, y, g, binary = data(outcome)
    if perm is not None: W = W[np.random.default_rng(10_000 + perm).permutation(len(W))]
    return outcome, seed, perm, nested(B, y, g, seed, binary), nested(np.hstack([B, W]), y, g, seed, binary)

OUTC = [o for o in BIN + CONT if D[o].notna().sum() >= 40 and (o not in BIN or min(D[o].sum(), (1 - D[o]).sum()) >= 8)]
obs = Parallel(n_jobs=NJ)(delayed(run)(o, r) for o in OUTC for r in range(NREP))
nul = Parallel(n_jobs=NJ)(delayed(run)(o, 900 + i, i) for o in OUTC for i in range(NPERM))
Ob = pd.DataFrame(obs, columns=['outcome', 'seed', 'perm', 'base', 'full']); Nu = pd.DataFrame(nul, columns=['outcome', 'seed', 'perm', 'base', 'full'])
Ob.to_csv(f'{OUT}/aim3_predict_repeats.csv', index=False); Nu.to_csv(f'{OUT}/aim3_predict_null.csv', index=False)
rows = []
for o in OUTC:
    b, n = Ob[Ob.outcome == o], Nu[Nu.outcome == o]; gain = (b.full - b.base).mean()
    B_, W_, y, g, binary = data(o)
    rows.append(dict(outcome=o, metric='AUC' if binary else 'R2', n=len(y), n_cases=int(y.sum()) if binary else np.nan,
                     base_mean=b.base.mean(), full_mean=b.full.mean(), full_lo=b.full.quantile(.025), full_hi=b.full.quantile(.975),
                     gain_mean=gain, gain_lo=(b.full - b.base).quantile(.025), gain_hi=(b.full - b.base).quantile(.975),
                     perm_p=(1 + ((n.full - n.base) >= gain).sum()) / (1 + len(n))))
R = pd.DataFrame(rows); R.to_csv(f'{OUT}/aim3_prediction.csv', index=False)
pd.set_option('display.width', 220); print(R.round(3).to_string(index=False))
