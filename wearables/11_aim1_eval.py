"""Aim 1, step 2: agreement of wrist methods with Dreem staging (reference = Dreem auto-staging).

Analysis set: Dreem recordings passing the sleep-architecture core QC with >= 80 % wrist
coverage. Epochs scored MT (movement) or U (unscored) by Dreem are excluded from epoch
metrics; 4-class labels: W / light (N1+N2) / deep (N3) / REM.

Methods
  A  Empatica sleep detection, binary: code 101 = sleep, else wake (102 = sleep as sensitivity).
     Codes are undocumented; mapping inferred from the Dreem cross-tab (10_aim1_pairs).
  C  van Hees 2015 sustained inactivity: sleep if no 5-s change in anglez > 5 deg for >= 5 min.
  D  XGBoost on wrist features + centred rolling context (1, 2.5, 5, 10, 30 min) of movement,
     angle change, HR, RMSSD, temperature, EDA; per-recording z-scores of HR/temp/EDA;
     elapsed fraction of the recording. Binary and 4-class models.
     Cross-validation: 5 folds grouped by FAMILY (no person or relative in both train and test).
  (B, Cole-Kripke / Sadeh, needs count rescaling for this device - next iteration.)
Metrics
  epoch: accuracy, sensitivity (sleep), specificity (wake), Cohen's kappa; 4-class accuracy,
  kappa, per-class F1.  night: TST/WASO/SOL bias and 95 % limits of agreement, ICC(A,1).
  By genotype: kappa and TST bias per group.
Usage: python -I 11_aim1_eval.py <phase0_dir> <aim1_dir> [n_threads]
"""
import os, sys, json
import numpy as np, pandas as pd
from sklearn.metrics import cohen_kappa_score, f1_score, confusion_matrix
from sklearn.model_selection import GroupKFold
import xgboost as xgb

P0, A1 = sys.argv[1:3]
NT = int(sys.argv[3]) if len(sys.argv) > 3 else 8
P = pd.read_parquet(os.path.join(A1, 'aim1_epoch_pairs.parquet'))
L = pd.read_csv(os.path.join(P0, 'linkage.csv'))
cov = P.groupby('recording_id').coverage.first()
keep = cov[cov >= 0.8].index
P = P[P.recording_id.isin(keep) & (P.qc_pass_core == True)].copy()
P = P.sort_values(['recording_id', 'epoch']).reset_index(drop=True)
P = P.merge(L[['empatica_id', 'base_id', 'rarecnv_family_id', 'genotype']], on='empatica_id', how='left')
P['family'] = P.rarecnv_family_id.fillna(P.base_id)
P['scored'] = ~P.stage.isin(['MT', 'U'])
P['sleep'] = P.stage.isin(['N1', 'N2', 'N3', 'R']).astype(int)
P['y4'] = P.stage.map({'W': 0, 'N1': 1, 'N2': 1, 'N3': 2, 'R': 3})

# ---- method A: vendor -----------------------------------------------------------------
P['A_sleep'] = (P.sleep_detection_stage == 101).astype(float).where(P.sleep_detection_stage.notna())
P['A2_sleep'] = P.sleep_detection_stage.isin([101, 102]).astype(float).where(P.sleep_detection_stage.notna())

# ---- method C: van Hees 2015 angle rule on 30-s epochs ------------------------------------
# an epoch is 'still' if max 5-s angle change in it <= 5 deg; sleep = inside a run of >= 10 still epochs (5 min)
def runs_ge(b, k):
    out = np.zeros(len(b), bool); i = 0
    while i < len(b):
        if b[i]:
            j = i
            while j < len(b) and b[j]:
                j += 1
            if j - i >= k:
                out[i:j] = True
            i = j
        else:
            i += 1
    return out
P['C_sleep'] = np.nan
for rid, idx in P.groupby('recording_id').groups.items():
    s = (P.loc[idx, 'dangle_max'].fillna(99) <= 5).values
    P.loc[idx, 'C_sleep'] = runs_ge(s, 10).astype(float)

# ---- method D features -----------------------------------------------------------------
base = ['enmo_mean', 'enmo_max', 'dangle_max', 'anglez_mean', 'hr_bpm', 'rmssd_ms', 'n_beats', 'temp_c', 'eda_us', 'vector_magnitude']
g = P.groupby('recording_id')
for c in ['hr_bpm', 'temp_c', 'eda_us', 'rmssd_ms']:
    P[c + '_z'] = (P[c] - g[c].transform('median')) / g[c].transform('std')
P['elapsed'] = g.epoch.transform(lambda e: e / max(e.max(), 1))
P['log_enmo'] = np.log10(P.enmo_mean.clip(lower=1e-4))
roll_src = ['log_enmo', 'dangle_max', 'hr_bpm_z', 'rmssd_ms_z', 'temp_c_z', 'eda_us_z']
feats = base + [c + '_z' for c in ['hr_bpm', 'temp_c', 'eda_us', 'rmssd_ms']] + ['elapsed', 'log_enmo']
for w in (2, 5, 10, 20, 60):                                  # 30-s epochs: 1, 2.5, 5, 10, 30 min
    for c in roll_src:
        P[f'{c}_m{w}'] = g[c].transform(lambda s: s.rolling(w, center=True, min_periods=1).mean())
        feats.append(f'{c}_m{w}')
    P[f'still_frac_m{w}'] = g.dangle_max.transform(lambda s: (s <= 5).astype(float).rolling(w, center=True, min_periods=1).mean())
    feats.append(f'still_frac_m{w}')

S = P[P.scored].reset_index(drop=True)
X = S[feats].astype('float32').values
gkf = GroupKFold(n_splits=5)
S['D_sleep'] = np.nan; S['D_p_sleep'] = np.nan; S['D4'] = np.nan
fold_info = []
for k, (tr, te) in enumerate(gkf.split(X, S.sleep, groups=S.family)):
    mb = xgb.XGBClassifier(n_estimators=400, max_depth=6, learning_rate=0.05, subsample=0.8, colsample_bytree=0.8,
                           tree_method='hist', n_jobs=NT, eval_metric='logloss')
    mb.fit(X[tr], S.sleep.values[tr])
    p = mb.predict_proba(X[te])[:, 1]
    S.loc[te, 'D_p_sleep'] = p; S.loc[te, 'D_sleep'] = (p >= 0.5).astype(float)
    m4 = xgb.XGBClassifier(n_estimators=400, max_depth=6, learning_rate=0.05, subsample=0.8, colsample_bytree=0.8,
                           tree_method='hist', n_jobs=NT, objective='multi:softprob', num_class=4, eval_metric='mlogloss')
    m4.fit(X[tr], S.y4.values[tr])
    S.loc[te, 'D4'] = m4.predict(X[te])
    fold_info.append(dict(fold=k, n_test_epochs=len(te), n_test_families=S.family.iloc[te].nunique()))
imp = pd.Series(mb.get_booster().get_score(importance_type='gain')).rename(lambda f: feats[int(f[1:])] if f.startswith('f') and f[1:].isdigit() else f)

# ---- epoch metrics -----------------------------------------------------------------------
def binm(y, yhat):
    ok = ~np.isnan(yhat); y, yhat = y[ok], yhat[ok].astype(int)
    tn, fp, fn, tp = confusion_matrix(y, yhat, labels=[0, 1]).ravel()
    return dict(n=int(ok.sum()), accuracy=(tp + tn) / len(y), sensitivity=tp / (tp + fn), specificity=tn / (tn + fp),
                kappa=cohen_kappa_score(y, yhat))
rows = []
for m in ['A_sleep', 'A2_sleep', 'C_sleep', 'D_sleep']:
    rows.append(dict(method=m, **binm(S.sleep.values, S[m].values.astype(float))))
E = pd.DataFrame(rows)
ok = S.D4.notna()
f1 = f1_score(S.y4[ok], S.D4[ok].astype(int), average=None, labels=[0, 1, 2, 3])
E4 = pd.DataFrame([dict(method='D_4class', n=int(ok.sum()), accuracy=(S.y4[ok] == S.D4[ok]).mean(),
                        kappa=cohen_kappa_score(S.y4[ok], S.D4[ok].astype(int)),
                        F1_wake=f1[0], F1_light=f1[1], F1_deep=f1[2], F1_REM=f1[3])])
cm4 = pd.crosstab(S.y4[ok].map({0: 'W', 1: 'light', 2: 'deep', 3: 'REM'}), S.D4[ok].map({0: 'W', 1: 'light', 2: 'deep', 3: 'REM'}),
                  rownames=['Dreem'], colnames=['model D'])

# ---- night metrics ----------------------------------------------------------------------
def night(df, col):
    s = df[col].fillna(0).values; first = np.argmax(s > 0) if s.any() else len(s)
    return pd.Series({'TST': s.sum() / 2, 'SOL': first / 2, 'WASO': (1 - s[first:]).sum() / 2 if first < len(s) else np.nan})
N = []
for rid, d in S.groupby('recording_id', sort=False):
    ref = night(d.assign(ref=d.sleep), 'ref')
    for m in ['A_sleep', 'C_sleep', 'D_sleep']:
        est = night(d, m)
        N.append(dict(recording_id=rid, empatica_id=d.empatica_id.iloc[0], genotype=d.genotype.iloc[0], method=m,
                      **{f'ref_{k}': v for k, v in ref.items()}, **{f'est_{k}': v for k, v in est.items()},
                      kappa=cohen_kappa_score(d.sleep, d[m].fillna(0).astype(int)) if d.sleep.nunique() > 1 else np.nan))
N = pd.DataFrame(N)

def icc_a1(a, b):                       # two-way random, absolute agreement, single measure (McGraw & Wong)
    Y = np.column_stack([a, b]); Y = Y[~np.isnan(Y).any(1)]; n, k = Y.shape
    gm = Y.mean(); msr = k * ((Y.mean(1) - gm) ** 2).sum() / (n - 1); msc = n * ((Y.mean(0) - gm) ** 2).sum() / (k - 1)
    mse = ((Y - Y.mean(1, keepdims=True) - Y.mean(0) + gm) ** 2).sum() / ((n - 1) * (k - 1))
    return (msr - mse) / (msr + (k - 1) * mse + k * (msc - mse) / n)
NA = []
for m, d in N.groupby('method'):
    for v in ['TST', 'SOL', 'WASO']:
        diff = d[f'est_{v}'] - d[f'ref_{v}']
        NA.append(dict(method=m, measure=v, n_nights=int(diff.notna().sum()), bias_min=diff.mean(), sd=diff.std(),
                       loa_low=diff.mean() - 1.96 * diff.std(), loa_high=diff.mean() + 1.96 * diff.std(),
                       icc=icc_a1(d[f'est_{v}'].values, d[f'ref_{v}'].values)))
NA = pd.DataFrame(NA)
G = N.assign(tst_diff=N.est_TST - N.ref_TST).groupby(['method', 'genotype']).agg(
    nights=('kappa', 'size'), kappa_median=('kappa', 'median'), tst_bias=('tst_diff', 'mean'), tst_sd=('tst_diff', 'std')).reset_index()

os.makedirs(os.path.join(A1, 'results'), exist_ok=True)
R = os.path.join(A1, 'results')
E.to_csv(f'{R}/epoch_binary.csv', index=False); E4.to_csv(f'{R}/epoch_4class.csv', index=False); cm4.to_csv(f'{R}/confusion_4class.csv')
N.to_csv(f'{R}/night_level.csv', index=False); NA.to_csv(f'{R}/night_agreement.csv', index=False); G.to_csv(f'{R}/by_genotype.csv', index=False)
imp.sort_values(ascending=False).head(30).to_csv(f'{R}/modelD_binary_importance_lastfold.csv')
S[['recording_id', 'empatica_id', 'genotype', 'family', 'epoch', 'ts_utc', 'stage', 'sleep', 'y4', 'A_sleep', 'C_sleep', 'D_sleep', 'D_p_sleep', 'D4']].to_parquet(f'{R}/epoch_predictions.parquet', index=False)
json.dump(dict(recordings=int(S.recording_id.nunique()), participants=int(S.empatica_id.nunique()), families=int(S.family.nunique()),
               epochs=int(len(S)), folds=fold_info), open(f'{R}/run_info.json', 'w'), indent=1)
pd.set_option('display.width', 200)
print(json.load(open(f'{R}/run_info.json')))
print(E.round(3).to_string(index=False)); print(E4.round(3).to_string(index=False)); print(cm4.to_string())
print(NA.round(2).to_string(index=False)); print(G.round(2).to_string(index=False))
print('top features (binary, last fold):', imp.sort_values(ascending=False).head(10).round(0).to_dict())
