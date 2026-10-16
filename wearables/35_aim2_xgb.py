"""Aim 2: XGBoost signature of each CNV group from WRIST BAND (Empatica EmbracePlus) participant-level phenotypes.

Unit = participant (visit 1, >=3 valid days). Tasks: each carrier group vs controls, and deletion vs duplication per locus.
Feature sets
  core     : the wr_ core phenotypes (20_phenotypes.py), confirmatory
  extended : core + wrist/mat second-panel phenotypes (wr_p2_*, mat_p2_*), exploratory; NaN left to XGBoost
Validation : nested CV (Varma & Simon 2006). Outer 5-fold, inner 3-fold, both StratifiedGroupKFold by family;
             outer loop repeated 20x. Inner loop picks hyperparameters by AUC.
Confounds  : age, age^2, sex, site (+ season as sin/cos of the participant's median recording day-of-year when SEASON=1) regressed out of every feature on the training fold, coefficients applied to the
             test fold (Snoek 2019). Baseline = logistic regression on the confounds alone, same folds.
Metrics    : out-of-fold AUC, balanced accuracy (threshold = training prevalence), Brier, calibration slope;
             mean and 2.5-97.5 % range over the 20 repeats.
Significance: permutation null, 1000 label shuffles within site, each a full nested-CV run (1 repeat);
             p = (1 + #null >= observed mean AUC) / (1 + n_perm). Conservative, since the null uses single repeats.
Interpretation: final model on all participants (confounds regressed on all), most-chosen hyperparameters;
             TreeSHAP via xgboost pred_contribs. Direction = sign of Spearman(feature, SHAP). Compared with the sign of the
             covariate-adjusted group difference (OLS: feature ~ group + age + age^2 + sex + site).
Samples    : all ages, and the pre-specified age 7-25 subsample.
Usage: python -I 35_aim2_xgb.py <master_csv> <out_dir> <n_jobs> [n_perm]
"""
import os, re, sys, itertools, warnings
import numpy as np, pandas as pd, xgboost as xgb
from joblib import Parallel, delayed
from scipy.stats import spearmanr
from sklearn.model_selection import StratifiedGroupKFold
from sklearn.metrics import roc_auc_score, balanced_accuracy_score, brier_score_loss
from sklearn.linear_model import LogisticRegression
warnings.filterwarnings('ignore')

MASTER, OUT, NJ = sys.argv[1], sys.argv[2], int(sys.argv[3]); NPERM = int(sys.argv[4]) if len(sys.argv) > 4 else 1000
os.makedirs(OUT, exist_ok=True)
NREP, OUTER, INNER = int(os.environ.get('NREP', 20)), 5, 3
GRID = [dict(max_depth=d, n_estimators=n) for d in (1, 2, 3) for n in (100, 300)]
FIXED = dict(learning_rate=0.05, subsample=0.8, colsample_bytree=0.8, min_child_weight=1, reg_lambda=1.0,
             tree_method='hist', n_jobs=1, eval_metric='logloss', verbosity=0)
TASKS = [('16pDel', 'Control'), ('16pDup', 'Control'), ('22qDel', 'Control'), ('22qDup', 'Control'),
         ('16pDel', '16pDup'), ('22qDel', '22qDup')]

M = pd.read_csv(MASTER)
M = M[(M.visit == 1) & (M.link_status == 'linked') & M.genotype.isin(['16pDel', '16pDup', '22qDel', '22qDup', 'Control'])
      & M.wr_tst_h.notna() & M.age_at_recording.notna()].copy()
if os.environ.get('NO_MEDS') == '1':
    # same rule as Aim 4 no_meds (23_hr_models.R): drop carriers on heart-rate-active psychotropics or without a meds form;
    # controls kept (family members were mostly not asked)
    HRMED = ['med_stimulant', 'med_adhd_nonstim', 'med_alpha2', 'med_antipsychotic', 'med_antidepressant', 'med_lithium']
    onmed = (M[HRMED].fillna(0).sum(1) > 0) | (M.meds_form_present != 1)
    M = M[~((M.genotype != 'Control') & onmed)]
    print('NO_MEDS: carriers on HR-active meds or without meds form removed', flush=True)
SEASON = os.environ.get('SEASON') == '1'
if SEASON:
    Dd = pd.read_csv(os.path.join(os.path.dirname(os.path.dirname(MASTER)), 'phenotypes_day.csv'), usecols=['empatica_id', 'date'])
    Dd['t'] = (pd.to_datetime(Dd.date) - pd.Timestamp('1970-01-01')).dt.days
    tm = Dd.groupby('empatica_id').t.median().round()
    med = pd.Series(pd.DatetimeIndex(pd.Timestamp('1970-01-01') + pd.to_timedelta(tm.values, unit='D')).dayofyear, index=tm.index)
    M['doy'] = M.empatica_id.map(med); M = M[M.doy.notna()]
    print('SEASON: sin/cos day-of-year added to confounds', flush=True)
M['family'] = M.family.where(M.family.notna() & (M.family.astype(str) != ''), M.study_id).astype(str)
COUNTS = {'wr_n_valid_days', 'wr_n_valid_nights', 'wr_p2_n_day_still_epochs', 'wr_p2_n_days_p2', 'wr_p2_n_eda_nights', 'mat_p2_n_mat_sessions'}  # data-volume counts, not phenotypes
CORE = [c for c in M.columns if c.startswith('wr_') and '_p2_' not in c and c not in COUNTS]
EXT = CORE + [c for c in M.columns if (c.startswith('wr_p2_') or c.startswith('mat_p2_')) and c not in COUNTS
              and c not in ('mat_p2_snore_frac',)]          # mat snoring is all zeros
SETS = {'core': CORE, 'extended': EXT,
        # wrist-only extension: drops the mat features, whose availability (and the undocumented breathing score) may differ by group
        'extended_wrist': [c for c in EXT if not c.startswith('mat_p2_')]}
# headband (Dreem EEG) per-participant sleep phenotypes (29_headband_core.py); TRT, stage minutes and MT fraction left out (redundant / QC)
HB = [c for c in M.columns if c.startswith('hb_') and not c.startswith(('hb_TRT', 'hb_min_', 'hb_frac_MT'))]
SETS['all_devices'] = EXT + HB
if os.environ.get('HB_SUBSET') == '1':
    # same-sample comparison: restrict every set to participants with headband data, so gains reflect features, not sample
    M = M[M.hb_TST.notna()]
    print('HB_SUBSET: participants with headband data only', flush=True)
SETS = {k: v for k, v in SETS.items() if k in re.split('[,:]', os.environ.get('FEATURE_SETS', 'core,extended'))}
print('participants', len(M), M.genotype.value_counts().to_dict(), {k: len(v) for k, v in SETS.items()}, flush=True)


def conf_matrix(d, levels):
    a = d.age_at_recording.values; ac = a - 20.0
    cols = [np.ones(len(d)), ac, ac ** 2, (d.sex.values == levels['sex'][0]).astype(float)]
    cols += [(d.site.values == s).astype(float) for s in levels['site'][1:]]
    if SEASON: cols += [np.sin(2 * np.pi * d.doy.values / 365.25), np.cos(2 * np.pi * d.doy.values / 365.25)]
    return np.column_stack(cols)


def residualise(Xtr, Ctr, Xte, Cte):
    Rtr, Rte = np.full_like(Xtr, np.nan), np.full_like(Xte, np.nan)
    for j in range(Xtr.shape[1]):
        ok = np.isfinite(Xtr[:, j])
        if ok.sum() < Ctr.shape[1] + 3: continue
        b = np.linalg.lstsq(Ctr[ok], Xtr[ok, j], rcond=None)[0]
        Rtr[:, j] = Xtr[:, j] - Ctr @ b; Rte[:, j] = Xte[:, j] - Cte @ b
    return Rtr, Rte


def fit(X, y, p):
    return xgb.XGBClassifier(**FIXED, **p).fit(X, y)


def nested_cv(X, C, y, g, seed):
    """One repeat of nested CV. Returns out-of-fold xgb prob, baseline prob, threshold per row, chosen configs."""
    pr, pb, th, chosen = np.full(len(y), np.nan), np.full(len(y), np.nan), np.full(len(y), np.nan), []
    for tr, te in StratifiedGroupKFold(OUTER, shuffle=True, random_state=seed).split(X, y, g):
        Xtr, Xte = residualise(X[tr], C[tr], X[te], C[te])
        best, bs = None, -1
        inner = list(StratifiedGroupKFold(INNER, shuffle=True, random_state=seed + 7).split(Xtr, y[tr], g[tr]))
        for k, p in enumerate(GRID):
            s = []
            for itr, ite in inner:
                if len(np.unique(y[tr][ite])) < 2: continue
                Ri, Re = residualise(X[tr][itr], C[tr][itr], X[tr][ite], C[tr][ite])
                s.append(roc_auc_score(y[tr][ite], fit(Ri, y[tr][itr], p).predict_proba(Re)[:, 1]))
            if np.mean(s) > bs: bs, best = np.mean(s), k
        chosen.append(best)
        pr[te] = fit(Xtr, y[tr], GRID[best]).predict_proba(Xte)[:, 1]
        sc = C[tr][:, 1:].std(0); sc[sc == 0] = 1
        lr = LogisticRegression(C=1.0, max_iter=2000).fit(C[tr][:, 1:] / sc, y[tr])
        pb[te] = lr.predict_proba(C[te][:, 1:] / sc)[:, 1]; th[te] = y[tr].mean()
    return pr, pb, th, chosen


def metrics(y, p, th):
    lp = np.log(np.clip(p, 1e-6, 1 - 1e-6) / (1 - np.clip(p, 1e-6, 1 - 1e-6)))
    slope = LogisticRegression(C=1e6, max_iter=2000).fit(lp[:, None], y).coef_[0, 0] if np.std(lp) > 0 else np.nan
    return dict(auc=roc_auc_score(y, p), bal_acc=balanced_accuracy_score(y, (p >= th).astype(int)),
                brier=brier_score_loss(y, p), cal_slope=slope)


def task_data(case, ctrl, feats, agerange):
    d = M[M.genotype.isin([case, ctrl])]
    if agerange: d = d[(d.age_at_recording >= agerange[0]) & (d.age_at_recording <= agerange[1])]
    d = d.reset_index(drop=True)
    lev = dict(sex=sorted(d.sex.astype(str).unique()), site=sorted(d.site.astype(str).unique()))
    d = d.assign(sex=d.sex.astype(str), site=d.site.astype(str))
    return d, d[feats].to_numpy(float), conf_matrix(d, lev), (d.genotype == case).astype(int).to_numpy(), d.family.to_numpy()


def one_perm(case, ctrl, setname, agerange, k):
    d, X, C, y, g = task_data(case, ctrl, SETS[setname], agerange)
    rng = np.random.default_rng(10_000 + k); yp = y.copy()
    for s in np.unique(d.site):
        ix = np.where(d.site.values == s)[0]; yp[ix] = rng.permutation(yp[ix])
    pr, _, _, _ = nested_cv(X, C, yp, g, seed=500 + k)
    return dict(case=case, control=ctrl, set=setname, sample=agerange_name(agerange), perm=k, auc=roc_auc_score(yp, pr))


def one_repeat(case, ctrl, setname, agerange, r):
    d, X, C, y, g = task_data(case, ctrl, SETS[setname], agerange)
    pr, pb, th, chosen = nested_cv(X, C, y, g, seed=r)
    row = dict(case=case, control=ctrl, set=setname, sample=agerange_name(agerange), rep=r, n_case=int(y.sum()), n_ctrl=int((1 - y).sum()),
               n_fam=len(np.unique(g)), chosen=';'.join(map(str, chosen)))
    row.update({f'xgb_{k}': v for k, v in metrics(y, pr, th).items()}); row.update({f'base_{k}': v for k, v in metrics(y, pb, th).items()})
    oof = pd.DataFrame(dict(study_id=d.study_id, case=case, control=ctrl, set=setname, sample=row['sample'], rep=r, y=y, p_xgb=pr, p_base=pb))
    return row, oof


def agerange_name(a): return 'age7to25' if a else 'all'


def explain(case, ctrl, setname, agerange, best):
    d, X, C, y, g = task_data(case, ctrl, SETS[setname], agerange)
    R, _ = residualise(X, C, X[:1], C[:1])
    m = fit(R, y, GRID[best])
    sh = m.get_booster().predict(xgb.DMatrix(R, missing=np.nan), pred_contribs=True)[:, :-1]
    rows = []
    for j, f in enumerate(SETS[setname]):
        ok = np.isfinite(X[:, j]); dd = d[ok]
        rho = spearmanr(R[ok, j], sh[ok, j]).correlation if ok.sum() > 5 and np.std(sh[ok, j]) > 0 else np.nan
        # covariate-adjusted group difference (OLS) for the direction check
        Z = np.column_stack([y[ok], C[ok]]); b = np.linalg.lstsq(Z, X[ok, j], rcond=None)[0][0]
        rows.append(dict(case=case, control=ctrl, set=setname, sample=agerange_name(agerange), feature=f, mean_abs_shap=np.abs(sh[ok, j]).mean(),
                         shap_direction=np.sign(rho) if np.isfinite(rho) else np.nan, spearman_feature_shap=rho,
                         adj_group_diff=b, sd_feature=np.nanstd(X[:, j]), n_nonmissing=int(ok.sum())))
    t = pd.DataFrame(rows).sort_values('mean_abs_shap', ascending=False); t['rank'] = np.arange(1, len(t) + 1)
    t['sign_agrees'] = np.where(t.mean_abs_shap > 0, np.sign(t.adj_group_diff) == t.shap_direction, np.nan)
    return t


combos = [(c, k, s, a) for (c, k) in TASKS for s in SETS for a in (None, (7, 25))]
res = Parallel(n_jobs=NJ, verbose=0)(delayed(one_repeat)(c, k, s, a, r) for (c, k, s, a) in combos for r in range(NREP))
P = pd.DataFrame([r for r, _ in res]); O = pd.concat([o for _, o in res]); O.to_csv(os.path.join(OUT, 'aim2_oof_predictions.csv'), index=False)
P.to_csv(os.path.join(OUT, 'aim2_repeats.csv'), index=False); print('repeats done', flush=True)
nulls = Parallel(n_jobs=NJ, verbose=0)(delayed(one_perm)(c, k, s, a, i) for (c, k, s, a) in combos for i in range(NPERM))
N = pd.DataFrame(nulls); N.to_csv(os.path.join(OUT, 'aim2_perm_null.csv'), index=False); print('permutations done', flush=True)

key = ['case', 'control', 'set', 'sample']
S = P.groupby(key).agg(n_case=('n_case', 'first'), n_ctrl=('n_ctrl', 'first'), n_fam=('n_fam', 'first'),
                       **{f'{m}_{s}': (m, f) for m in [c for c in P.columns if c.startswith(('xgb_', 'base_'))]
                          for s, f in [('mean', 'mean'), ('lo', lambda x: x.quantile(.025)), ('hi', lambda x: x.quantile(.975))]}).reset_index()
S['auc_gain_mean'] = S.xgb_auc_mean - S.base_auc_mean
S['perm_p'] = [(1 + (N[(N[key] == r[key].values).all(1)].auc >= r.xgb_auc_mean).sum()) / (1 + NPERM) for _, r in S.iterrows()]
S['null_auc_mean'] = [N[(N[key] == r[key].values).all(1)].auc.mean() for _, r in S.iterrows()]
S.to_csv(os.path.join(OUT, 'aim2_performance.csv'), index=False)

bestcfg = {tuple(k): pd.Series(';'.join(v.chosen).split(';')).astype(int).mode()[0] for k, v in P.groupby(key)}
E = pd.concat([explain(c, k, s, a, bestcfg[(c, k, s, agerange_name(a))]) for (c, k, s, a) in combos])
E.to_csv(os.path.join(OUT, 'aim2_shap.csv'), index=False)
pd.set_option('display.width', 250)
print(S[key + ['n_case', 'n_ctrl', 'xgb_auc_mean', 'xgb_auc_lo', 'xgb_auc_hi', 'base_auc_mean', 'auc_gain_mean', 'xgb_bal_acc_mean',
               'xgb_cal_slope_mean', 'perm_p']].round(3).to_string(index=False))
S0 = list(SETS)[-1]
for (c, k), t in E[(E.set == S0) & (E['sample'] == 'all')].groupby(['case', 'control'], sort=False):
    print(f'\n{c} vs {k}: top features ({S0}, all ages)')
    print(t.head(6)[['feature', 'mean_abs_shap', 'shap_direction', 'adj_group_diff', 'sign_agrees']].round(3).to_string(index=False))
