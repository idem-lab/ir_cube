# Cross-validation run: state, plan, and handover notes

Working document for the posterior predictive validation work on branch
`posterior-predictive-validation` (PR idem-lab/ir_cube#12, issue #10).
Written so that someone picking this up cold can continue without re-deriving
what has already been established or repeating what has already failed.

Last updated: 2026-09-04, at the launch of the corrected fold run.

---

## 1. What this is for

Cross-validation previously scored the model on its ability to predict the value
of an individual held-out bioassay, with the posterior collapsed to a mean
before scoring. A bioassay is a noisy, overdispersed measurement of the quantity
the model actually targets — the susceptibility fraction of the whole mosquito
population — so that comparison has an unreachable noise floor and discards
parameter uncertainty. The work replaces it with scoring of the full
out-of-sample posterior predictive distribution: calibration and coverage
alongside accuracy.

The plan is set out in full in the comment on issue #10. Nothing in the
statistical design has changed since; what follows is about executing it.

## 2. Environment: do not upgrade greta

Sampling this model fails on greta 0.6.0 — `mcmc()` raises a TensorFlow
`while_loop` shape error whenever the likelihood depends on
`iterate_dynamic_function()` output. Reported upstream as
greta-dev/greta.dynamics#45.

The working combination, verified for both master's model definition and
`fit_fold()`:

```r
remotes::install_version("tensorflow", version = "2.16.0")
remotes::install_github("njtierney/greta@4cc989f")            # 0.5.0.9000
remotes::install_github("greta-dev/greta.dynamics@db7df31")   # 0.2.2
```

Python side is a conda env built by `greta::install_greta_deps()`:
TensorFlow 2.15.1, TFP 0.23.0.

Two ordering constraints, both load-bearing, both encoded in `run_one_fold.R`:

- TensorFlow will not change its thread count after initialisation, so threads
  must be set before python comes up. greta exposes no interface for this;
  it must go through `reticulate::import("tensorflow")$config$threading$...`.
  Environment variables (`TF_NUM_INTRAOP_THREADS` etc.) are ignored.
- python must initialise **before** `terra` or `sf` are attached. Those load the
  system XML libraries, against which the conda environment's `pyexpat` is then
  resolved, and `tensorflow_probability` fails to import.

## 3. What is already done

| piece | file | state |
|---|---|---|
| Scoring functions | `R/validation_functions.R` | done, 24 checks pass |
| Checks on those | `R/check_validation_functions.R` | done |
| Overdispersion from replicates | `R/estimate_bioassay_rho.R` | done, run |
| Fold definitions | `R/validation_folds.R` | done, verified identical to master's |
| Null models | `R/null_models.R` | done, all 16 null folds saved |
| Covariates | `R/validation_covariates.R` | done |
| Model fit per fold | `R/fit_validation_fold.R` | done |
| Single fold runner | `R/run_one_fold.R` | done, smoke tested |
| Dispatcher | `R/run_validation_folds.R` | done, running |
| Metrics | `R/validation_metrics.R` | **needs adapting, see §6** |
| Figures | `R/fig_predictive_validation.R` | written, tested on synthetic draws |

External overdispersion estimates, from all 3,713 replicated
pixel-year-insecticide groups (`outputs/bioassay_rho.csv`): 0.155 overall;
0.118 organochlorines, 0.131 carbamates, 0.164 pyrethroids,
0.216 organophosphates. These set the noise floor and are independent of the
models being scored.

## 4. The run now in progress

Launched 2026-09-04 07:09 by `Rscript R/run_validation_folds.R`.

- 8 folds: 6 leave-one-country-out, plus spatial interpolation and temporal
  forecasting.
- **2 folds at a time, 4 chains and 4 threads each**, `Lmax = 30`,
  2,000 warmup, 500 initial samples extending in batches of 500 to a cap of
  5,000.
- One process per fold (`run_one_fold.R`), logging to
  `outputs/cv_logs/<experiment>__<fold>.log`.
- Draws saved to `outputs/cv_draws/dynamical__<experiment>__<fold>.rds`.
  Folds already on disk are skipped, so the run resumes after interruption.

Expected: roughly 40–50h per pair, so **150–200h in total**. The first fold
should land within ~50h.

Why this configuration, from measurements on one fold of this model:

- Chains are vectorised into one TensorFlow op, not run one per core, and that
  op scales poorly beyond ~4 threads: 8.39 s/iteration at 2 threads against
  6.96 s with the whole machine. So confining a fold to 4 threads costs little
  and leaves room for a second fold.
- Chain count costs more than linearly once cores are saturated: 6.96 s/iteration
  at 2 chains, 15.21 at 4, 70.61 at 8.
- But 2 chains adapt badly. greta pools information across chains during warmup,
  and the 2-chain run produced Rhat 7.6 on the Kenya fold with all 689
  parameters above 1.01, against 1.2–1.6 elsewhere. 4 chains is the compromise.

## 5. Monitoring, and what to check when a fold lands

Watch progress:

```bash
tail -f "outputs/cv_logs/spatial_extrapolation__Kenya.log"
grep -h "ESS" outputs/cv_logs/*.log | tail
```

Each fold logs, with timestamps: warmup and sampling progress bars, effective
sample size after every batch of extra samples, prediction ESS for `p` and
`rho`, and the worst Rhat.

**The first fold is the checkpoint that matters.** Within ~50h there should be
one `.rds` in `outputs/cv_draws/`. Check it before trusting the remaining 150h:

```r
d <- readRDS("outputs/cv_draws/dynamical__spatial_extrapolation__Côte d’Ivoire.rds")
names(d)                       # must include `draws`
class(d$draws)                 # greta_mcmc_list
dim(d$p_draws)                 # (n_chains * n_sampled) x n_test
median(d$ess_p); min(d$ess_p)  # ESS of the predicted fractions
median(d$ess_rho)
max(d$convergence[, 1])        # worst Rhat
```

What counts as acceptable: `ess_p` median in the hundreds is fine for the
metrics, which pool over ~1,000–1,900 held-out assays per fold. Rhat is the
harder question — see §7.

If the object lacks `draws`, or `p_draws` has the wrong number of rows, stop the
run; something in `fit_fold()` has regressed and the remaining folds will be
unusable too.

## 6. Work to do while the run proceeds

`R/validation_metrics.R` was written against the old saved format and needs
adapting. Specifically:

1. **It must not re-predict.** Use the `p_draws` and `rho_draws` in the saved
   object, which came from `calculate(values = draws)`. Do not write bespoke
   prediction code; if predictions at new points are ever needed, rebuild the
   model and use `calculate()` with the saved `draws`.
2. **Draw count has changed** from 1,000 to `n_chains * n_sampled`, up to
   20,000. `ppd_summary()` loops over observations building an
   `n_draws x (died + 1)` matrix each time, so this is ~20x slower than before
   and may need thinning. Thinning is cheap here: at ~50 draws per effective
   sample, taking every 10th draw loses almost nothing. Add a `thin` argument
   rather than silently subsampling.
3. **Use the per-class external `rho`** from `outputs/bioassay_rho.csv` for the
   noise floor, not the fitted values.
4. The null-model folds in `outputs/cv_draws/` have 1,000 draws and no `draws`
   object, by construction — they are analytic, not MCMC. The metrics code must
   handle both shapes.

Test it against the first completed fold as soon as it lands, together with the
16 null folds already on disk. That exercises the whole path end to end long
before the run finishes.

`R/fig_predictive_validation.R` reads only the tables the metrics script writes,
so it should need little change, but it has only ever been run on synthetic
draws.

## 7. Open questions for Nick

- **Convergence.** The production fit in `temporary/fitted_model.RData`
  (8 chains x 2,000) has worst Rhat 1.215 with 631 of 689 parameters above 1.01,
  and 45.7 draws per effective sample. The CV folds mix comparably. So
  non-convergence is a standing property of this model's geometry rather than
  something the validation introduced, and the options are to accept it with an
  explicit caveat in the supplement, or to treat reparameterisation as separate
  modelling work affecting the main fit too.
- **Fitted `rho` runs about double the external estimate** — 0.24–0.38 against
  0.118–0.216 — consistently across all four folds of the earlier run. If it
  holds, this is the diagnostic the plan anticipated: the model absorbing
  process misfit into the observation process, which should also show as
  over-coverage. Worth a supplementary panel.
- **ESS target.** 1,000 was set arbitrarily and is unreachable: the model needs
  ~50 draws per effective sample, so it would take ~50,000 draws per fold. The
  cap of 5,000 samples binds instead. Decide whether the achieved ESS on `p` and
  `rho` (the quantities that matter) is sufficient once the first fold reports.

## 8. Pitfalls already hit, worth not repeating

- `calculate(..., nsim = n)` returns an **independent resample** of the
  posterior. It preserves joint structure across quantities, so it is valid for
  prediction, but it destroys MCMC ordering, so effective sample size cannot be
  recovered from it, and it is not what greta's prediction interface expects to
  be handed. Use `calculate(values = draws)` without `nsim`.
- `coda::effectiveSize()` on 20 draws returns roughly 20. Any ESS computed from
  a short run is meaningless; several tuning conclusions were drawn from such
  numbers and had to be withdrawn.
- `future.callr` buffers worker stdout until the future resolves. A long run
  under it is invisible. Hence one process per fold.
- Functions called from a worker must have every dependency passed explicitly.
  `codetools::findGlobals(fit_fold, merge = FALSE)$variables` catches free
  variables; `n_unique_cells` was missing this way and cost a run.
- Verify that string edits to these scripts actually applied. A silently
  non-matching `str.replace` cost a second run.
- Always smoke test the full path with tiny settings before a long run:
  `Rscript R/run_one_fold.R spatial_extrapolation Kenya 2 4 5 5 5` takes a few
  minutes and exercises everything. Delete the resulting `.rds` afterwards, or
  the real fold will be skipped.

## 9. Superseded artefacts

`outputs/cv_draws_defunct/` holds four folds from the earlier 2-chain run. They
have no `draws` object, so they cannot be used with greta's prediction
interface, and their Kenya fold did not converge. Kept only for the diagnostics
quoted above; delete once the new run is complete.
