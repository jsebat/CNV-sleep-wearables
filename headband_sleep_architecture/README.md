# Sleep architecture × 16p11.2 / 22q11.2 copy number

Analysis pipeline for overnight Dreem EEG hypnograms from the SFARI / G2MH
cohort: dose-response effects of copy number on NREM depth, total sleep time
and sleep onset latency.

The model specification (mirror and U-shape copy-number terms, age, age², sex, site,
participant random effect) follows the original study design; ingest, QC, metrics,
covariate adjustment (partial residuals from one reference model) and figures are
implemented here.

## Running it

```bash
Rscript inspect_inputs.R    # sanity-check the two inputs -- ALWAYS run this first
Rscript run_all.R           # ingest -> QC -> adjust -> model -> figures -> sensitivity
Rscript make_report.R       # build the HTML report from results/
```

Everything you should need to change lives in **`config.R`**, and in practice
only one line: `DATA_ROOT`.

```
$DATA_ROOT/data/sfari_hypnograms/SFARI_Dreem_Prod_Hypno_<DEVICE>_<datetime>__1.txt
$DATA_ROOT/data/demographics_distributions/full_sample_manifest.csv
$DATA_ROOT/results/                                    <- created by run_all.R
```

`DATA_ROOT` points **outside** the Dropbox-synced tree
(`C:/Users/Jonathan Sebat/sleep_arch_data` on this machine). That is
deliberate: the extracted hypnograms are ~900 small files (~75 MB) that merely
re-extract `sfari_hypnograms.zip`, and `run_all.R` writes another ~25 MB of
output per run including a 9 MB `workspace.rds`. Syncing any of it is pure
churn. Only the 16 source files live in Dropbox; the curated deliverables go to
`G2MH/SleepArchitectureResults`.

Paths are used as-is if absolute and resolved against the project directory if
relative, so `DATA_ROOT <- "data"` also works for a self-contained checkout.
The `SLEEP_HYPNO_DIR` / `SLEEP_MANIFEST` / `SLEEP_OUT_DIR` environment
variables override `config.R` without editing it.

Requires R ≥ 4.1 with `lme4`, `lmerTest`, `ggplot2`; `writexl` optional.
`emmeans` is deliberately not used — model-implied means are computed directly
from the design matrix so the covariate reference profile is explicit and
provably matches the adjustment applied to the plotted points.

## Layout

| file | what it does |
|---|---|
| `config.R` | paths, QC thresholds, outcome definitions, model spec |
| `inspect_inputs.R` | prints what the pipeline sees in your inputs; writes nothing |
| `run_all.R` | the driver |
| `make_report.R` | builds a standalone HTML report from `results/` |
| `R/00_utils.R` | stage vocabulary, file sniffing, genotype label mapping |
| `R/01_ingest.R` | **the only format-dependent file** |
| `R/02_qc_metrics.R` | per-recording metrics and QC |
| `R/03_models.R` | mixed models + the covariate adjustment |
| `R/04_figures.R` | QC Venn, box/scatter panels, copy-number line plot |
| `R/05_sensitivity.R` | robustness battery |
| `tests/make_synthetic_data.R` | generates fake inputs of the same shape |
| `tests/check_adjustment.R` | asserts the adjustment's properties |

## Pipeline

**1. Ingest.** Dreem exports carry a metadata preamble followed by a
tab-delimited table introduced by `Sleep Stage⇥Time [hh:mm:ss]⇥Event⇥Duration[s]`.
The reader detects that header and reads from there; per-epoch `Duration[s]`
values are used for all time totals rather than assuming a fixed epoch length.
Stage labels `SLEEP-S0/S1/S2/S3/REM/MT` are normalised to `W/N1/N2/N3/R/MT`
(S0 is wake). A generic delimited/one-stage-per-line reader handles other
formats.

Recordings link to subjects by the **device ID** in the filename
(`SFARI_Dreem_Prod_Hypno_SITE0000_…`), looked up in the manifest's `dreem_id`
column. This covers 126/128 devices directly and needs no phenotype-file hop.

Manifest columns are mapped by name with exact-case priority — this manifest
has **both** `guid` and `GUID` holding different ID formats, so case must not
be collapsed. Genotype comes from `final_group_class`; `no_cnv` → Control.
Subjects carrying two CNVs (`22q11.2_Deletion & 16p11.2_Deletion`) or a
triplication are **excluded**, since they belong to no single dose series.

**2. QC.** Three exclusions, reproducing the original analysis:

| criterion | default |
|---|---|
| No sleep | `TST == 0` |
| Too short | `TRT < 120 min` |
| High movement | `MT > 50%` of epochs |

⚠️ The original Venn diagram used **TRT < 2 hours**, not 30 minutes.
`QC_MIN_TRT_MIN` defaults to 120 to reproduce it; set it to `30` for the
30-minute rule. `R/05_sensitivity.R` refits across a grid of both thresholds —
the depth-index result is unchanged at every setting tested.

**No outlier removal is applied.** Nothing is excluded on the basis of its
outcome value; these three criteria are the only exclusions, matching the
original analysis (782 recordings). `QC_MAX_PCT_N3` can apply a fourth filter
on implausible architecture (N3 above a given % of TST) but defaults to `Inf`,
i.e. off — retaining the extreme nights is the conservative choice, since
dropping them strengthens the 22q11.2 result. `R/05_sensitivity.R` reports what
that filter would do at every threshold from 40% to 70%.

Figures instead **clip the y-axis** via `FIG_CLIP_IQR` (default 3, meaning the
Tukey range Q1−3×IQR to Q3+3×IQR), applied with `coord_cartesian()` so points
outside are hidden but retained — boxplots, model fits and p-values all use
every observation, and each caption states how many points are off-axis.

**3. Outcomes** (per night).

| name | definition |
|---|---|
| `DI` | `N3 / (N1 + N2 + REM)` = `N3 / (TST − N3)`; wake excluded. A ratio, unbounded above. **Primary.** |
| `TST` | summed duration of all N1, N2, N3, REM epochs |
| `SOL` | recording start → first sleep epoch (no lights-out diary available) |
| `DI_over_NREM` | `N3 / (N1+N2+N3)` — sensitivity |
| `DI_over_TST` | `N3 / TST` — sensitivity |
| `SOL_log` | `log(1 + SOL)` — sensitivity; SOL is strongly right-skewed |
| `SOL_persist` | latency to 10 min of near-continuous sleep — sensitivity |

**4. Models.** Per outcome, per locus (copy number 1/2/3), controls shared:

| model | genotype term | tests |
|---|---|---|
| `mirror_single` | `copy_c` ∈ {−1, 0, +1} | linear dose / mirror effect |
| `ushape_single` | `copy_q` ∈ {+1, 0, +1} | both tails vs middle |
| `joint` | `copy_c + copy_q` | each adjusted for the other |
| `group_vs_control` | `Genotype` factor, Control = reference | Del and Dup separately |

All carry `+ Age_c + I(Age_c^2) + Sex + Site + (1 | GUID)`, where `Age_c` is
mean-centred age (numerically equivalent to raw age but avoids the Age/Age²
collinearity warning). A five-group model across the whole cohort is also
fitted. Constant covariates are dropped automatically.

Wald tests use REML with Satterthwaite d.f. (`lmerTest`); model comparisons use
LRTs on ML refits: *U-shape beyond mirror*, *mirror beyond U-shape*, and 2-df /
4-df omnibus tests. BH-FDR is applied within three families over the primary
outcomes.

**5. Covariate adjustment.** Plotted values are partial residuals with respect
to genotype:

```
y_adj[i] = y[i] − Σ_k (X[i,k] − mean(X[,k])) · β[k]     over covariate columns only
```

`β` comes from **one reference model per outcome**, fitted once on the full
five-group cohort with genotype as a saturated factor. Genotype is *in* the
model (so covariates aren't confounded by group) but *excluded from the
subtraction* (so group differences survive). Consequences:
`mean(y_adj) == mean(y)`; genotype, participant and residual structure are all
preserved; and a night plots at the same height in every panel. Model-implied
means are evaluated at the same reference — every covariate design column at
its sample mean, i.e. proportional weighting over sex and site.

`tests/check_adjustment.R` asserts these properties on every run. See
[`NOTES_on_original_code.md`](NOTES_on_original_code.md) for what this replaced.

**6. Sensitivity.** `R/05_sensitivity.R` refits the primary outcomes under
value caps (DI is an unbounded ratio), alternative outcome definitions, and a
grid of QC thresholds. Results go to `sensitivity_analyses.csv`.

## Outputs (`results/`)

| file | contents |
|---|---|
| `results_all_tests.csv` / `.xlsx` | **every statistical test**, one row each |
| `sensitivity_analyses.csv` | the robustness battery |
| `analysis_dataset_per_night.csv` | the cleaned per-night table, incl. `adj_*` columns |
| `per_recording_metrics_with_qc.csv` | all recordings with QC flags, pre-exclusion |
| `group_descriptives.csv` | N nights / subjects, age, sex by group |
| `qc_venn_counts.csv` | the seven Venn regions |
| `qc_summary.csv` | sequential QC filter counts |
| `figures/` | QC Venn, per-outcome panels, copy-number line plot |
| `report/report.html` | standalone report (figures base64-embedded) |
| `workspace.rds` | `d`, `tests`, `sens`, `fits`, `metrics` for interactive follow-up |

`singular = TRUE` on a test row means the random-intercept variance collapsed
to zero — that outcome shows no between-participant variance beyond residual
noise. Fixed effects remain interpretable; the flag is there to be checked.

## If ingest picks the wrong thing

`inspect_inputs.R` prints the raw first file, parsed stage counts, mapped
manifest columns and the linkage rate. Then:

| symptom | fix in `config.R` |
|---|---|
| `UNSCORED` dominates | set `HYPNO_STAGE_COL`, or switch `NUMERIC_STAGE_CODING` |
| REM scored as N3 | `NUMERIC_STAGE_CODING <- "AASM"` |
| no files found | set `HYPNO_GLOB` |
| low linkage rate | check `SUBJECT_ID_MODE` / `DEVICE_ID_REGEX` |
| wrong manifest column | add the name to the candidate list in `read_manifest()` |
| epochs aren't 30 s | `EPOCH_SEC` (ignored when the file carries durations) |

## Verification

Beyond the real analysis, the pipeline runs end to end on synthetic data of the
same shape with known injected effects, which it recovers:

```bash
Rscript tests/make_synthetic_data.R /tmp/synth
SLEEP_HYPNO_DIR=/tmp/synth/sfari_hypnograms \
SLEEP_MANIFEST=/tmp/synth/full_sample_manifest.csv \
SLEEP_OUT_DIR=results_synthetic_test Rscript run_all.R
Rscript tests/check_adjustment.R results_synthetic_test
```

The injected effect sizes exist only to give the pipeline a known signal to
recover; they are not estimates and carry no scientific claim.
