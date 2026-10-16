"""Aim 2 feature weights (WRIST BAND, Empatica EmbracePlus) with bootstrap stability, joined to the univariate genotype effects.

For each comparison x feature set (core, extended_wrist), all ages:
  100 family-cluster bootstrap resamples; each: regress age, age^2, sex, site out of the features on the resample,
  fit XGBoost with the comparison's most-chosen nested-CV configuration (35_aim2_xgb.py, aim2/full and aim2/extended_wrist),
  TreeSHAP (xgboost pred_contribs) on ALL participants of the comparison (residualised with the resample's coefficients).
  Per feature: median over resamples of mean |SHAP|; share of resamples where the feature is in the top 5;
  direction = sign of Spearman(feature, SHAP) pooled over resamples.
Univariate: the matching Aim 4 / second-panel mixed-model estimate (y ~ group + age + age^2 + sex + site + random effects):
  carrier vs control -> 'group_vs_control' term for that group (5-group model);
  deletion vs duplication -> per-locus mirror (dosage) term, sign flipped so + = higher in deletions.
Usage: python -I 38_aim2_features.py <analysis_root> <out_dir> <n_jobs>
"""
import os, sys, warnings
import numpy as np, pandas as pd, xgboost as xgb
from joblib import Parallel, delayed
from scipy.stats import spearmanr
warnings.filterwarnings('ignore')
ROOT, OUT, NJ = sys.argv[1], sys.argv[2], int(sys.argv[3]); os.makedirs(OUT, exist_ok=True)
NB = 100
GRID = [dict(max_depth=d, n_estimators=n) for d in (1, 2, 3) for n in (100, 300)]
FIXED = dict(learning_rate=0.05, subsample=0.8, colsample_bytree=0.8, min_child_weight=1, reg_lambda=1.0, tree_method='hist', n_jobs=1, verbosity=0)
TASKS = [('16pDel', 'Control'), ('16pDup', 'Control'), ('22qDel', 'Control'), ('22qDup', 'Control'), ('16pDel', '16pDup'), ('22qDel', '22qDup')]

M = pd.read_csv(os.path.join(ROOT, 'aim4', 'master', 'G2MH_wearable_phenotypes_per_participant.csv'))
M = M[(M.visit == 1) & (M.link_status == 'linked') & M.genotype.isin(['16pDel', '16pDup', '22qDel', '22qDup', 'Control'])
      & M.wr_tst_h.notna() & M.age_at_recording.notna()].copy()
M['family'] = M.family.where(M.family.notna() & (M.family.astype(str) != ''), M.study_id).astype(str)
M['sex'] = M.sex.astype(str); M['site'] = M.site.astype(str)
COUNTS = {'wr_n_valid_days', 'wr_n_valid_nights', 'wr_p2_n_day_still_epochs', 'wr_p2_n_days_p2', 'wr_p2_n_eda_nights', 'mat_p2_n_mat_sessions'}  # data-volume counts, not phenotypes
CORE = [c for c in M.columns if c.startswith('wr_') and '_p2_' not in c and c not in COUNTS]
EXTW = CORE + [c for c in M.columns if c.startswith('wr_p2_') and c not in COUNTS]
SETS = {'core': CORE, 'extended_wrist': EXTW}

def conf(d, sexlev, sites):
    ac = d.age_at_recording.values - 20.0
    return np.column_stack([np.ones(len(d)), ac, ac ** 2, (d.sex.values == sexlev).astype(float)] + [(d.site.values == s).astype(float) for s in sites[1:]])

def best_cfg(run, case, ctrl, st):
    R = pd.read_csv(os.path.join(ROOT, 'aim2', run, 'aim2_repeats.csv'))
    R = R[(R.case == case) & (R.control == ctrl) & (R.set == st) & (R['sample'] == 'all')]
    return pd.Series(';'.join(R.chosen.astype(str)).split(';')).astype(int).mode()[0]

def one_boot(case, ctrl, st, b, cfg):
    d = M[M.genotype.isin([case, ctrl])].reset_index(drop=True); feats = SETS[st]
    sexlev, sites = sorted(d.sex.unique())[0], sorted(d.site.unique())
    X, C, y = d[feats].to_numpy(float), conf(d, sexlev, sites), (d.genotype == case).astype(int).to_numpy()
    rng = np.random.default_rng(b); fams = d.family.unique()
    ix = np.concatenate([np.where(d.family.values == f)[0] for f in rng.choice(fams, len(fams), replace=True)])
    R = np.full_like(X, np.nan)
    for j in range(X.shape[1]):
        ok = np.isfinite(X[ix, j])
        if ok.sum() < C.shape[1] + 3: continue
        coef = np.linalg.lstsq(C[ix][ok], X[ix, j][ok], rcond=None)[0]; R[:, j] = X[:, j] - C @ coef
    m = xgb.XGBClassifier(**FIXED, **GRID[cfg], random_state=b).fit(R[ix], y[ix])
    sh = m.get_booster().predict(xgb.DMatrix(R, missing=np.nan), pred_contribs=True)[:, :-1]
    return case, ctrl, st, b, np.nanmean(np.abs(sh), 0), R, sh

jobs = []
for case, ctrl in TASKS:
    for st, run in (('core', 'full'), ('extended_wrist', 'extended_wrist')):
        cfg = best_cfg(run, case, ctrl, st); jobs += [(case, ctrl, st, b, cfg) for b in range(NB)]
res = Parallel(n_jobs=NJ)(delayed(one_boot)(*j) for j in jobs)

# univariate genotype / dosage effects
A4 = os.path.join(ROOT, 'aim4')
U = pd.concat([pd.read_csv(os.path.join(A4, 'aim4_results.csv')).assign(prefix='wr_')] +
              [pd.read_csv(os.path.join(A4, 'phase2', f)).assign(prefix='wr_p2_') for f in ('block1_results.csv', 'blocks2to5_results.csv', 'eda_results.csv')])
U['feature'] = U.prefix + U['var']

rows = []
for (case, ctrl, st) in dict.fromkeys((r[0], r[1], r[2]) for r in res):
    rr = [r for r in res if r[:3] == (case, ctrl, st)]; feats = SETS[st]
    imp = np.vstack([r[4] for r in rr])                                   # boot x feature
    top5 = np.mean([np.isin(np.arange(len(feats)), np.argsort(-np.nan_to_num(v))[:5]) for v in imp], 0)
    for j, f in enumerate(feats):
        xs = np.concatenate([r[5][:, j] for r in rr]); ss = np.concatenate([r[6][:, j] for r in rr]); ok = np.isfinite(xs)
        rho = spearmanr(xs[ok], ss[ok]).correlation if ok.sum() > 10 and np.std(ss[ok]) > 0 else np.nan
        if ctrl == 'Control':
            u = U[(U.feature == f) & (U.model == 'group_vs_control') & (U.term == 'Genotype' + case)]; flip = 1
        else:
            u = U[(U.feature == f) & (U.model == 'mirror_single') & (U.cohort == case[:3])]; flip = -1
        u = u.iloc[0] if len(u) else None
        rows.append(dict(case=case, control=ctrl, set=st, feature=f, importance=np.nanmedian(imp[:, j]),
                         importance_lo=np.nanquantile(imp[:, j], .1), importance_hi=np.nanquantile(imp[:, j], .9), top5_freq=top5[j],
                         shap_direction=np.sign(rho) if np.isfinite(rho) else np.nan,
                         univ_std_effect=(flip * u.std_estimate) if u is not None else np.nan, univ_p=u.p_value if u is not None else np.nan,
                         univ_test=('group vs control' if ctrl == 'Control' else f'{case[:3]} dosage (sign: + = higher in del)') if u is not None else ''))
F = pd.DataFrame(rows)
F['rank'] = F.groupby(['case', 'control', 'set']).importance.rank(ascending=False, method='first').astype(int)
F['sign_agrees'] = np.where(F.univ_std_effect.notna() & F.shap_direction.notna(), np.sign(F.univ_std_effect) == F.shap_direction, np.nan)
F = F.sort_values(['case', 'control', 'set', 'rank']); F.to_csv(os.path.join(OUT, 'aim2_feature_weights.csv'), index=False)
pd.set_option('display.width', 230)
for k, t in F.groupby(['case', 'control', 'set'], sort=False):
    print('\n', k); print(t.head(8)[['feature', 'importance', 'top5_freq', 'shap_direction', 'univ_std_effect', 'univ_p', 'sign_agrees']].round(3).to_string(index=False))
# agreement summary: among each model's top 10, how many have a nominal univariate effect, and of those how many agree in sign
S = F[F['rank'] <= 10].groupby(['case', 'control', 'set']).agg(n_top=('feature', 'size'), n_univ_p05=('univ_p', lambda x: int((x < 0.05).sum())),
                                                               n_sign_agree=('sign_agrees', lambda x: int(np.nansum(x))), n_with_univ=('univ_p', lambda x: int(x.notna().sum()))).reset_index()
S.to_csv(os.path.join(OUT, 'aim2_feature_agreement.csv'), index=False); print('\n', S.to_string(index=False))
