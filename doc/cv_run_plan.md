# Cross-validation: design, run recipe, and results

Written for whoever runs the next set of fits. Sections 1–3 are the recipe and
the constraints it exists to satisfy; sections 4–7 are the rebuild the review of
PR #12 asks for, and what it costs.

---

## 1. Environment: greta 0.6, no greta.dynamics

The selection recursion is computed in closed form by one greta op
(`closed_form_states()` in `R/dynamical_model.R`, #25), so the model no longer
uses `iterate_dynamic_function()`. That was what stopped greta 0.6.0 sampling
this model (a TensorFlow `while_loop` shape error, greta-dev/greta.dynamics#45),
and greta.dynamics is no longer needed.

greta must be at or after `43f9c52` (0.6.0.9000), which fills subassignments
into greta arrays column-major, as R does (greta-dev/greta#844, fixed in #847).
`check_greta_fill()` (`R/greta_setup.R`) checks this behaviour and stops
if it is wrong.

Tested with greta `282944f`, TensorFlow 2.21.0, TensorFlow Probability 0.25.0,
python 3.12, in a library and conda environment separate from the 0.5 setup:

```bash
Rscript -e 'remotes::install_github("greta-dev/greta@282944f", lib = "~/R/greta06-lib")'
~/.local/share/r-miniconda/bin/conda create -n greta06-env python=3.12
~/.local/share/r-miniconda/envs/greta06-env/bin/python -m pip install \
  "tensorflow==2.21.*" "tensorflow_probability[tf]==0.25.*"
```

and every script that uses greta is run with

```bash
export R_LIBS=~/R/greta06-lib
export RETICULATE_PYTHON=~/.local/share/r-miniconda/envs/greta06-env/bin/python
```

`R_LIBS` rather than `R_LIBS_USER`, which `~/.Renviron` sets. Set
`RETICULATE_PYTHON` explicitly: otherwise greta 0.6 finds the old
`greta-env-tf2` (TensorFlow 2.15), and its own uv-managed environment failed to
resolve here. The other R packages come from the default library, as before.

Under greta 0.6 the model gives the same log density as under 0.5.0.9000
(`4cc989f`, TensorFlow 2.15.1) at the same parameter values, full data and the
interpolation fold (difference 0, gradient to 6e-13), and samples.

## 2. Two ordering constraints, both load-bearing

Both are encoded at the top of `R/run_one_fold.R`, and the file will fail in
confusing ways if they are disturbed:

- TensorFlow will not change its thread count after initialisation, so threads
  must be set **before** python comes up. greta exposes no interface for this;
  it must go through `reticulate::import("tensorflow")$config$threading$...`.
  The environment variables (`TF_NUM_INTRAOP_THREADS` and friends) are ignored.
- python must initialise **before** `terra` or `sf` are attached. Those load the
  system XML libraries, against which the conda environment's `pyexpat` is then
  resolved, and `tensorflow_probability` fails to import.

## 3. Running the folds

```bash
Rscript R/fig_illustrate_bioassay_variability.R  # overdispersion per type; everything downstream needs it
Rscript R/run_validation_folds.R                 # nulls, then dispatch the model folds
Rscript R/validation_metrics.R                   # score everything on disk
Rscript R/validation_change.R                    # score predicted change, per forecast origin
Rscript R/validation_geometry.R                  # fold separation, and the leak check
Rscript R/variance_explained.R                   # variance explained and the noise ceiling
Rscript R/fig_variance_explained.R               # the bar figures
Rscript R/fig_predictive_validation.R            # figures and the table
```

The overdispersion fit comes first: `rho_lookup()` stops rather than falling
back, so nothing scores until `outputs/bioassay_rho_hierarchical.csv` and
`outputs/bioassay_rho_type_draws.rds` exist.

`run_validation_folds.R` fits the nulls in process (minutes) and then dispatches
each model fold as a separate `Rscript R/run_one_fold.R <experiment> <fold>`,
two at a time, logging to `outputs/cv_logs/<experiment>__<fold>.log`. A fold
whose `.rds` is already in `outputs/cv_draws/` is skipped, by `run_one_fold.R`
itself, so the run resumes after an interruption.

**Run from a frozen copy of the scripts.** R reads `--file=` incrementally, so
editing a script while a fold is running corrupts that fold — two folds were
lost this way after 62 h each, both crashing at the `saveRDS` block. Copy `R/`
to a scratch directory, symlink `data/`, `outputs/` and `temporary/` into it,
and launch from there.

Smoke test the whole path before committing to a long run:

```bash
Rscript R/run_one_fold.R spatial_blocks 1 2 4 5 5
```

takes a few minutes and exercises everything. Delete the resulting `.rds`
afterwards, or the real fold will be skipped.

### Sampling settings (#25)

`dynamical_mcmc_settings()` in `R/dynamical_model.R` holds them, and
`fit_fold()`, `run_one_fold.R`, `run_validation_folds.R` and `fit_model.R` all
take them from there: `windowed_hmc()` (`R/windowed_hmc.R`) with 30 to 60
leapfrog steps (`Lmin`, `Lmax`), redrawn every 10 iterations, target
acceptance 0.65, 4 chains, 2,000 warmup and 1,500 samples. The model samples the
countries' initial states centred and starts from `temporary/inits_refit.RDS`
(`dynamical_inits_file`).

**Centred selection hierarchy (#48).** Since #48 the model also samples the
class and type levels of the net use, IRS and population selection effects,
and the type level of the overdispersion, centred
(`centred_options_data_informed()`), which the data pin down; the crop
effects stay non-centred. The steps were 60 to 120 before. Full fits of the
default model (d_half 270, no floor) on the full data, the same data and
initial values, 4 chains, 2,000 + 3,000, 8 vCPU RunPod pods (cpu5c, except
cpu3c at 62 ms per gradient for 15-30); ESS per 1,000 gradients counts the
sampling iterations at the mean number of steps:

| hierarchy, steps | hours | worst rank Rhat | bulk ESS min / median | per 1,000 gradients min / median |
|---|---|---|---|---|
| non-centred, 60-120 | 6.18 | 1.01 | 400 / 1,999 | 1.48 / 7.41 |
| centred, 15-30 | 1.93 | 1.01 | 514 / 2,657 | 7.62 / 39.4 |
| **centred, 30-60** (these settings) | 3.00 | 1.00 | 1,208 / 7,217 | 8.95 / 53.5 |

The centred fits' posterior means of the 683 quantities all three share
(`R/centred_pilot.R`) differ from the non-centred fit's by a median 0.7
Monte Carlo standard errors (at most 3.8) and their sds by a median ratio
of 1.00. The worst-mixing quantity in all three is `sigma_class` for net use:
it is set by the type deviations the data barely inform (the
organophosphates, and bendiocarb and DDT alone in their classes), which
centred sit in a mild funnel with it.

At 30-60 steps the worst quantity had 0.40 effective samples per draw, so
1,500 samples (3,000 before #48) should give a minimum bulk ESS of about 600.
The warmup stays at 2,000: at 1,000 (+ 1,500, the same pod type), one of the
four chains ended warmup in that funnel's neck (`sigma_class` for net use at
0.13, against a posterior median of 0.27) and rejected every proposal while
sampling, at the step size the other chains adapted (rank Rhat up to 4.6).

At four threads the 2014 forecasting fold took 3.7 h at 3.0 s per iteration
(2,000 + 2,500, machine loaded to about 8 of 16 cores), so about 4.2 h at
3,000 samples. The interpolation fold ran 1.5 times slower per iteration than
the 2014 fold at 15-30 steps, so it should take about 6 h, and the full fit,
now at the same settings, about as long. It was roughly 62 h per fold with
the greta.dynamics loop.

**Result.** On the 2014 forecasting fold, for the default model, against
greta's `hmc()` with the non-centred hierarchy (on the interpolation fold, the
only fold that baseline was run on):

| | fold | draws | worst rank Rhat | params > 1.01 / > 1.05 | bulk ESS min / median |
|---|---|---|---|---|---|
| `hmc()`, non-centred, 15-30 steps | interpolation | 20,000 | 1.156 | 575 / 83 | 18 / 232 |
| windowed, centred, 15-30 steps | 2014 | 20,000 | 1.049 | 207 / 0 | 104 / 600 |
| **these settings** (2,500 samples) | 2014 | 10,000 | 1.013 | 1 / 0 | 361 / 2,723 |

At 3,000 samples the minimum should reach about 430. These settings were not
run on the interpolation fold, which mixed worse than the 2014 fold at every
setting tried (at 15-30 steps, minimum ESS 33 against 104). The worst-mixing
parameters are the countries' initial levels for type 7 (`init_country_level[,
7]`, ESS 361-630) and `init_country_sd[7]` (586); on the interpolation fold at
15-30 steps they were `mortality_floor`, the levels of countries with few data
for a type, and the levels of all countries of types 2 and 3 together.

**What limits mixing,** found from the posterior correlations, the principal
components of the draws and the per-chain means, not tried blind:

1. *greta's warmup.* `hmc()` sets the diagonal of the mass matrix between 10%
   and 40% of warmup from every warmup draw so far, transient included. On the
   interpolation fold its final `diag_sd` was 0.06 to 20 times the posterior
   sd; on 20 independent normals with sds 0.01-100 it was 0.04 to 9.3 times,
   with minimum ESS 14 of 4,000. `windowed_hmc()` estimates it in windows that
   start afresh, as Stan does (0.94-1.07 times the true sds, minimum ESS 933),
   pooled over chains.
2. *The non-centred hierarchy on the initial state.* Most countries have data
   for most types, which fixes each country's initial state, so the
   non-centred deviations of all countries in a region move together against
   their region's (correlations above 0.9, in the September fits too). The
   countries' levels are now sampled directly, around their region's level at
   the country's own mean initial-state covariates: an exact
   reparameterisation (the log densities match up to the Jacobian to 1e-8,
   and `R/check_dynamical_model.R` passes). Centring the deviations but not
   `logit_init_mean` moved the ridge to `logit_init_mean` (rank Rhat 1.41);
   centring the regions too put them in a funnel with `init_region_sd`, which
   5 regions barely identify (a chain sat still for 2,000 iterations at
   `init_region_sd` 0.06); levels at covariates of 0 traded off against the
   covariates' coefficients (-0.6), since the covariates are standardised over
   the whole mask and the data cells lie above its mean.
3. *Trajectory length.* The step size adapts to the acceptance target, but
   the number of leapfrog steps does not, and 15-30 steps were too few to
   move along the directions that remain slow: the mortality floor against
   the initial levels of many countries at once (correlations 0.4-0.5), the
   levels of one type against its selection and initial-state coefficients,
   and the levels of countries without data for a type, which in centred form
   sit in a mild funnel with `init_country_sd`. None is a linear ridge a
   reparameterisation removes exactly, and a dense mass matrix did not help.
   60-120 steps cost 2.2 times as much per iteration as 15-30 on the 2014
   fold and gave 7 times the minimum ESS per draw.

**Pilots on the full folds** (4 chains, 2,000 + 5,000, 4 threads; ESS per
core-hour is bulk ESS over threads times wall hours, warmup included; wall
time varied with the machine's load, 8 to 20 of 16 cores):

| run | fold | s/it | worst rank Rhat | > 1.05 | ESS min / median | per core-hour min / median |
|---|---|---|---|---|---|---|
| `hmc()`, non-centred | interp | 3.36 | 1.156 | 83 | 18 / 232 | 0.7 / 8.9 |
| windowed, regions and countries centred as deviations | interp | 2.68 | 1.407 | 267 | 9 / 134 | 0.4 / 6.5 |
| windowed, regions and countries centred as levels | 2014 | 1.20 | 1.233 | 480 | 14 / 164 | 1.5 / 17.6 |
| windowed, countries centred | interp | 2.00 | 1.168 | 53 | 22 / 254 | 1.4 / 16.3 |
| windowed, countries centred | 2014 | 1.19 | 1.059 | 7 | 56 / 539 | 6.0 / 58.2 |
| same, dense mass matrix | 2014 | 1.16 | 1.163 | 134 | 19 / 188 | 2.1 / 20.8 |
| + levels at the overall mean covariates | interp | 2.17 | 1.102 | 16 | 28 / 319 | 1.7 / 18.9 |
| + levels at the overall mean covariates | 2014 | 1.20 | 1.050 | 0 | 68 / 596 | 7.3 / 64.1 |
| + at each country's mean (**these settings**) | 2014 | 1.37 | 1.049 | 0 | 104 / 600 | 9.8 / 56.2 |
| + step size per chain | interp | 1.97 | 1.118 | 25 | 33 / 250 | 2.1 / 16.3 |
| + step size per chain | 2014 | 1.01 | 1.120 | 146 | 28 / 250 | 3.5 / 31.6 |
| shared step, failed proposals as rejections | 2014 | 1.93 | 1.066 | 19 | 44 / 374 | 2.9 / 24.9 |
| **60-120 steps** (2,000 + 2,500) | 2014 | 2.98 | 1.013 | 0 | 361 / 2,723 | 24.2 / 182.5 |

The per-chain step size freed a chain that moved once in 1,000 draws among 8
on the screening subset, but on the 2014 fold one chain then rejected 35% of
its proposals while sampling; both it and the change to failed proposals were
worse there, so the step size is adapted as `hmc()` does. Runs of one setting
vary by about as much as some of these differences.

**Chains, warmup and trajectory length,** screened on a smaller model: 12
countries from all 5 regions, 7,053 assays, with the model restricted to them.
(A subset of the data with all 46 countries in the model left whole regions
without data, a geometry the folds do not have.) 2 threads, 4 chains and
1,000 + 2,000 unless shown, 8,000 post-warmup draws in all:

| | s/it | worst rank Rhat | > 1.05 | ESS min / median | per core-hour min / median |
|---|---|---|---|---|---|
| 4 chains (two runs) | 0.77-0.84 | 1.051-1.061 | 1-7 | 55-82 / 341-464 | 43-59 / 265-331 |
| 4 chains, 2,500 warmup | 0.81 | 1.102 | 41 | 32 / 204 | 16 / 101 |
| 8 chains, 1,000 + 1,000 | 1.46 | 1.077 | 14 | 73 / 386 | 45 / 237 |
| 16 chains, 1,000 + 500 | 3.19 | 1.115 | 133 | 97 / 345 | 37 / 130 |
| L 5-15 | 0.77 | 1.272 | 174 | 12 / 95 | 9 / 74 |
| L 15-30 | 1.76 | 1.060 | 2 | 51 / 451 | 17 / 154 |
| L 30-60 | 3.10 | 1.050 | 1 | 79 / 961 | 15 / 186 |
| L 60-120 | 4.94 | 1.006 | 0 | 604 / 3,981 | 73 / 484 |
| L 120-240, 1,000 + 1,000 | 7.74 | 1.010 | 0 | 189 / 2,371 | 22 / 276 |

(The L runs ran together on a machine loaded to about 20 of 16 cores, so their
s/it are comparable with each other but not with the rows above.) Neither more
chains nor a longer warmup gave more effective samples per core-hour, and on
the full interpolation fold chains cost at least proportionally: 2.95, 5.46
and 15.14 s per iteration at 4, 8 and 16 chains (70 iterations, 4 threads).
Trajectory length did: 60-120 steps gave 4 times the minimum ESS per
core-hour of 15-30 on the subset, and 2.5 times on the 2014 fold; 120-240
gave less. Cost per iteration grows less than linearly with the steps
(overhead per iteration). Four chains rather than two because the metric is estimated from all
chains pooled: at two chains, under `hmc()`, the Kenya fold reached Rhat 7.6.
TensorFlow threads scale poorly beyond about four, which is why the budget
goes into concurrent folds rather than threads.

**Initial values.** `temporary/inits.RDS`, from the fits before the refit,
has no `logit_init_mean`, which the centred levels need (it started wherever
greta put it). `temporary/inits_refit.RDS` (not in git) holds the posterior
means of the default model on the interpolation fold from the `hmc()` run
above, before reversion was added (reversion starts at 0.01), and
`fit_model.R` rewrites it from the full fit. From the old inits the four
chains of that run agreed in their means: the stuck chains above came from
the sampler and the funnels, not from the starting values.

The September run reached 5,000 samples by taking 500 and topping up
with `extra_samples()` towards a 1,000 ESS target. That target was set on the
~690 raw hierarchical parameters, whose minimum ESS was 76–96 on every fold, so
it was never reachable and the cap always bound. The loop has been removed and
the samples are asked for directly; the two are statistically equivalent, since
`extra_samples()` continues the same chains without re-adapting.

### Production run, October 2026

The full fit and the five folds at the settings above, from frozen copies of
`R/` (commit `ecad4ba`), 4 threads each, the full fit and three folds at once
(16 threads, load 17-19 of 16 cores). Wall time includes prediction; the old
inits file is kept as `temporary/inits_refit_before_production.RDS`, and
`fit_model.R` rewrote `inits_refit.RDS` from this fit. Diagnostics over the
raw parameters (733), rank-normalised (`posterior`):

| fit | wall time | chains kept | worst rank Rhat | > 1.01 / > 1.05 | bulk ESS min / median | tail ESS min |
|---|---|---|---|---|---|---|
| full data | 20.3 h | 4 of 4 | 1.016 | 1 / 0 | 226 / 2,527 | 43 |
| spatial blocks 1 | 16.7 h | 3 of 4 | 1.007 | 0 / 0 | 740 / 2,369 | 416 |
| spatial blocks 2 | 17.7 h | 4 of 4 | 1.029 | 207 / 0 | 163 / 2,000 | 34 |
| interpolation | 19.0 h | 4 of 4 | 1.008 | 0 / 0 | 513 / 2,370 | 377 |
| forecasting 2014 | | | | | | |
| forecasting 2018 | | | | | | |

The full fit's worst parameter is `init_country_level[38, 6]`; blocks 2's are
`sigma_overall[11]` and the class-level coefficients of column 11, with one
chain's draws 20% repeated. The interpolation fold, which mixed worst at every
setting before, is now among the best (minimum ESS on `mortality_floor`).

One chain of blocks 1 repeated a single draw for all 3,000 samples (rank Rhat
1.53 on every parameter with it). `R/drop_stuck_chains.R` drops such a chain
from a fold's stored predictions and records it in `stuck_chains`, which
`paired_draw_index()` reads; the fold then keeps 1,500 of its 2,000 stored
draws. No other fit had a stuck chain (`stuck_chains()`: under half the draws
distinct).

Wall times are 3-4 times the pilots' estimate (4-6 h): with 16 threads on 8
physical cores the four fits ran at 12-15 s per iteration each. Stopping two
of the four for 7 minutes roughly doubled the other two's rate, so the
machine has about 8 cores' worth of TensorFlow throughput, not 16.

### Constrained refit, October 2026: two modes in the mortality floor

The refit with the initial-state coefficients constrained to be at most 0
(`12e51b4`) has two posterior modes: the mortality floor near 0.25-0.32 with
more susceptible initial country levels, or near 0.001-0.002. In the full fit
two chains found each; the low-floor chains had a log posterior about 59
higher, and the other two were dropped (`DROP_CHAINS=1,2 Rscript
R/drop_stuck_chains.R`; the original is kept as `draws_all_chains`, and in
`temporary/fitted_model_raw_constrained.RData`). The folds and the half-length
sensitivity fits (4 chains, 2,000 + 1,500), run on RunPod, checked the same
way: per-chain means of every parameter, rank Rhat (`posterior`), and where
chains disagree the log posterior per chain (60 draws each).

| fit | floor per chain | chains kept | worst rank Rhat before / after |
|---|---|---|---|
| full data | 0.269, 0.269, 0.0015, 0.0016 | 3, 4 | 1.74 / 1.03 |
| interpolation | 0.002 × 4 | all | 1.007 |
| spatial blocks 1 | 0.246 × 4 | all | 1.005 |
| spatial blocks 2 | 0.300 × 4 | all | 1.005 |
| forecasting 2014 | 0.0014 × 4 | all | 1.004 |
| forecasting 2018 | 0.0008 × 4 | all | 1.008 |
| d_half 5 | 0.002 × 4 | all | 1.037 |
| d_half 200 | 0.318 × 4 | all | 1.010 |
| crop trend 0.37 | 0.002 × 4 | all | 1.020 |
| crop and population trends 0.37 | 0.0012 × 4 | all | 1.023 |
| net weight 0.47 | 0.272 × 4 | all | 1.015 |
| legacy net use | 0.001 × 4 | all | 1.024 |
| no reversion | 0.268, 0.268, 0.0015, 0.269 | 3 | 1.54 / 1.06 (split, one chain) |

Without reversion the low-floor chain had a log posterior about 60 higher than
the other three (−94,284 against −94,342 to −94,350), as in the full fit, so it
alone is kept. Rhat after the drop is the split Rhat of that chain.

Every other fit has all four chains in one mode, so no chain can be dropped,
but five of them are in the high-floor mode: both block folds, d_half 200 and
net weight 0.47. Whether that is the minor mode in those fits is not known.
Evaluating each one's log posterior at the full fit's draws (which were not
fitted to that data) put both of the full fit's modes 65-1,160 below the
fit's own draws, and which of the two was higher varied: the low floor by 35
(blocks 1) and 100 (net weight 0.47), the high floor by 30 (blocks 2) and 420
(d_half 200). On the interpolation fold, a low-floor fit, the full fit's
low-floor draws were 45-95 below its own and the high-floor ones 140-150
below. So that check cannot settle it. The comparisons involving those fits — the block
folds against the previous rounds, and those two sensitivities against the
full fit — are partly a comparison between modes.

### Sensitivity fits

Each is the full fit with one change to `dynamical_model_options()`, at half
length (4 chains, 2,000 warmup + 1,500 samples), in both the constrained refit
and round 2 (floor Beta(1, 49)):

| fit | options |
|---|---|
| d_half 5 | `dynamical_model_options(selection_columns = selection_design(pop_d_half = 5))` |
| d_half 200 | `dynamical_model_options(selection_columns = selection_design(pop_d_half = 200))` |
| crop trend 0.37 | `dynamical_model_options(selection_columns = selection_design(trend_crops = linear_trend_matrix(0.37)))` |
| crop and population trends 0.37 | `dynamical_model_options(selection_columns = selection_design(trend_pop = linear_trend_matrix(0.37), trend_crops = linear_trend_matrix(0.37)))` |
| net weight 0.47 | `dynamical_model_options(selection_columns = selection_design(net_w = 0.47))` |
| legacy net use | `dynamical_model_options(selection_columns = selection_design(net_use_source = "legacy"))` |
| no reversion | `dynamical_model_options(reversion = FALSE)` |

`linear_trend_matrix(0.37)` (`R/model_covariates.R`) is g(t) = 0.37 in 1995
rising linearly to 1 in 2025 and on, the same in every region: FAOSTAT
insecticide use per ha of cropland in 1995 was about 1/2.7 of 2025's
(`R/faostat_trend_model.R`; #23), where the main fit's trends start at 0.

## 4. What is being changed, and why

The review of PR #12 found one blocking defect and one set of experiments that
measures something other than what it claims to. Both force refits; the rest of
the review is scoring and reporting, and has been done without refitting.

### 4.1 The temporal forecasting design: a leak, then a window with no power — rebuilt as a rolling origin

Two separate faults, found in that order.

**The leak.** The training set was "every record whose year is not 2020–2022".
The data run to 2024 while the covariates stop in 2022, so 888 records from 2023
and 2024 stayed in training: the model was fitted on both sides of the window it
was asked to forecast. 354 of the 1,461 held-out assays (24%) sat at pixels that
also carried post-horizon training data, and on those pixels the dynamical
model's MSE was 0.060 against 0.086 elsewhere. Fixed in `validation_folds.R`
(`year_start < cut_year`, with a `stopifnot`). The old fit is in
`outputs/cv_draws_leaky_forecast/` with a note; its scores are in the git
history of `outputs/cv_summary.csv` at `6b1ba1b`.

**The window.** The corrected 2020 three-year holdout turned out to be the worst
window in the record on both counts that matter. Sliding a three-year holdout
against the three years before it, the observed change at pixels assayed in both
is decisively negative at every origin from 2005 to 2017 and flattens only at
2018–2020: the 2020 origin reads +0.011 [−0.017, +0.040]. It is also the
thinnest, at 280 paired (pixel, insecticide) groups and 1,127 assays against
1,436 and 6,346 at a 2013 origin — enough to bound the mean change, nowhere near
enough to score any single site, which is why the direction test came out at
49%, a coin flip, and must not be reported as a finding.

The model is not misbehaving on the trend. It predicts −0.093 over that window,
and the slope fitted to the training years 2012–2019 implies −0.075 to −0.082.
It extrapolates the historical rate faithfully; the rate stopped.

**Five-year windows fix it.** At every origin they give 30–50% more paired
pixels and a signal roughly 5/3 larger, because the gap between window midpoints
is the window length, and a five-year window spans the 2018–2020 pause as well
as the decline either side of it. The 2018 origin's signal-to-resolvable ratio
goes from 0.2 at three years to 3.8 at five. The power analysis measures
this; the script and its figures are in the stub branch (§6).

**Two origins, fitted:**

| cut | training | % of data | holdout | paired pixels | observed change | ratio to resolvable |
|---|---|---|---|---|---|---|
| 2014 | 1995–2013, 14,285 | 52% | 2014–2018, 9,922 assays | 1,748 | −0.084 | 12.0 |
| 2018 | 1995–2017, 22,377 | 82% | 2018–2022, 4,096 assays | 881 | −0.042 | 3.8 |

Non-overlapping holdouts bar the shared endpoint, training at half and
four-fifths of the data, and a factor of two between their true rates of
decline. That contrast is the test: does the model track a slowing rate, or
carry a fixed slope forward? Earlier origins have a stronger signal still — a
2010 origin is the strongest in the record — but train on 17% of the data, so
they are not the model being deployed and their skill would not transfer.
Covariates end in 2022, so 2018 is the latest feasible five-year origin.

The 2020 three-year fold is deleted, not reported. Its training set was drawn
with `year_start <= cut`, so the cut year appeared in both training and holdout,
and its window landed on the one pause in twenty years of decline and was also
the thinnest. Its draws are parked in `outputs/cv_draws_defunct/`.

**The two origins are not independent.** The 2018 fold's before-window, 2013 to
2017, sits inside the 2014 fold's holdout, 2014 to 2018, and 2018 itself is in
both holdouts. Pooling them is still the right way to report a single forecasting
bar - the alternative is two bars whose difference is mostly window difficulty -
but the pooled interval carries roughly one and a bit folds' worth of
information rather than two, and should not be read as though the origins were
replicates. `outputs/cv_variance_explained_by_fold.csv` has them separately.

**How it is coded.** `validation_folds.R` exposes `forecasting_fold(cut_year,
window)` and a named list `temporal_forecasting_folds`, holding 2014 and 2018.
The single-fold alias `temporal_forecasting` is gone. Dispatch is `Rscript R/run_one_fold.R
temporal_forecasting 2014 4 4`. Each origin is stored under
`experiment = "temporal_forecasting_<cut>"` so that scoring never pools two
holdout windows whose true rates of decline differ by a factor of two; the file
name keeps the plain experiment name, so the origins sit together in
`outputs/cv_draws/`.

**Two changes the five-year window forced:**

- *The nearest neighbour null's lookback.* It searches the same year and
  `n_years_prior` earlier ones, intersected with training. With a five-year
  window and `n_years_prior = 3`, held-out years at lead 4 and 5 have no valid
  training year at all — and the old code then took `sort(...)[n]` of a vector
  of `Inf`, giving a threshold of `Inf`, and selected the *entire* training set:
  a global mean wearing a nearest-neighbour label, silently. `n_years_prior` is
  now the window length, and `predict_null_fixed_nn_counts()` stops rather than
  falling back. The three-year fold escaped this by one year, so its result is
  unaffected. The nearest-neighbour null is no longer tuned at all: it is
  reported at one neighbour, as a practice baseline, and separately at its best
  k on the held-out records, as an oracle bound.
- *Prediction volume.* The 2014 fold asks for 20,522 predictions against the
  three-year fold's 5,724. `fit_fold()` now thins the stored draws to 2,000 —
  which is what the scoring and the change score thin to anyway, and ESS is
  still measured on the unthinned ordered draws — and takes the before-window
  predictions in a second `calculate()` call, so peak memory tracks the larger
  window rather than both at once. `calculate(values = draws)` with no `nsim` is
  a deterministic function of the draws, so the pairing the change score needs
  survives the split.

### 4.2 Leave-one-country-out confounds spatial skill with an unidentified initial condition — removed

A held-out country's `init_country_raw` has no data, so its initial resistant
fraction reverts to the region prior, and that error is amplified through 15–29
years of deterministic selection before the comparison year. The review argued
this from the model structure; it was confirmed empirically before the folds
were dropped. Across the six folds the held-out bias correlated with that
country's fitted country effect at **r = −0.94**, and excess MSE against the
magnitude of that effect at **r = +0.91**. Côte d'Ivoire and Ethiopia have the
two largest negative country effects and were the two worst-predicted folds.

So these folds measure the difficulty of predicting an entirely unsampled
country, which is not a situation the deployed model faces — there are bioassays
in every country. The fold definitions and every code path that scored them have
been deleted, and the draws are parked in `outputs/cv_draws_defunct/`. The
sub-national blocks of §5 test spatial prediction without the confound.

### 4.3 A sub-national spatial block design replaces the national extrapolation concept — two new fits

Holding out blocks of cells within countries keeps every country intercept
identified, so the test isolates spatial prediction rather than compounding it
with the initial condition. See §5.

### 4.4 Scoring and reporting changes, all made without refitting

**One noise floor, one definition of variance explained.** The floor is
`noise_floor_mse()`: `y(1-y) k / (1-k)`, which since `E[y(1-y)] = p(1-p)(1-k)`
is exactly unbiased for a subset's mean sampling variance and assumes nothing
about how p is distributed. An empirical Bayes floor fitted per insecticide type
was tried and dropped — it was 1.4 to 9.6% high depending on rho, and is the
more fragile of the two where p is bimodal, as it is for DDT and
Alpha-cypermethrin. Variance explained is `1 - MSE/Var(y)`, floor-free, with the
noise share shown alongside as a band rather than divided out; the
intercept-null-referenced skill that also carried that name is gone, with the
second pixel-cluster bootstrap that supported it.

- **The PIT column was the mid-P value, not a randomised PIT.** `rowMeans()` over
  100 randomisations converges to `cdf_below + 0.5 * pmf_at`, which is not
  uniform under calibration for discrete data. Now `pit[, 1]`. The headline
  uniformity statistics used the full matrix and were never affected; the
  figures were. `check_validation_functions.R` now demonstrates the distortion
  on data with a large atom at 100% mortality: coverage 0.973 at nominal 0.95,
  CvM 10.58 against 0.03.
- **Skill is anchored on the intercept null**, not the nearest neighbour, which
  pinned an informative baseline at zero by construction. `excess` (MSE above
  the noise floor, in absolute mortality² units) and `rms_p` (its square root)
  are reported alongside every ratio.
- **Every model is scored at the external replicate-based overdispersion.**
  Letting each model fit its own made coverage a comparison of dispersion rather
  than of prediction: the intercept null reached 0.96 coverage by inflating rho
  to 0.45 against an external estimate of 0.12–0.22. Fitting a rho per null model
  is gone entirely — the nulls earn their place through mean squared error and
  variance explained, which need point predictions only, and the dynamical
  model's calibration is judged against held-out data directly. What the
  dynamical model's own residuals imply is still reported, in
  `cv_rho_comparison.csv` and its figure: a fitted 0.28 against an external 0.16
  is the model treating its own misfit as bioassay noise. The caveat to state in
  the paper: that posterior on p was fitted jointly with that rho, so scoring at
  the external value is not perfectly clean without a refit.
- **The metric surface is smaller.** Coverage, mean PIT and Cramér–von Mises are
  three functionals of one PIT distribution, which is enough; the
  Kolmogorov–Smirnov statistic and the PIT ECDF figure are gone. The CvM null
  band is gone from the reported tables, because it assumes independent PIT
  values and held-out records share a posterior, so it is too narrow — CvM is an
  ordering, not a test. The WHO-threshold block scored the sample quantity
  rather than the population quantity and averaged predictive quantiles across
  folds; gone. The pixel-year aggregation rung averaged 1.3 assays per group, so
  it was the unpooled comparison under another name; only the country-year rung
  (17–18 assays per group) is kept.
- **Per-fold and per-lead-year breakdowns are restored** (`cv_by_fold.csv`,
  `cv_by_year.csv`), which master reported and the first version of this
  pipeline dropped.
- **Skill is reported against distance and data volume**, in
  `R/validation_geometry.R`. For every held-out record: km to the nearest
  training record of the same insecticide, number of such records within 100 km,
  and years since the last observation at that pixel. This is what makes the
  arbitrary geometry of the folds matter less, and it is the answer to the
  question a user of the map actually has.

## 5. The sub-national spatial block design

**Which countries are split.** Only the six the national folds hold out: Côte
d'Ivoire, Ethiopia, Kenya, Nigeria, Senegal and Tanzania. `validation_folds.R`
selected those as the countries with more than ten bioassays of every one of the
nine insecticide types since 2010, and that criterion matters here for the same
reason. Keeping to them buys three things:

- the block experiment becomes directly comparable with the national one,
  differing in exactly the intended respect and no other — same countries, same
  insecticide coverage, and only whether each country's initial condition is
  identified, which is the whole mechanism under test;
- the held-out set keeps a controlled mix of insecticides, rather than becoming
  a weighted average over whatever the thinner countries happen to hold;
- each fold trains on 88% of the data rather than 67%, much closer to the
  production fit, so the result transfers to the deployed model more directly.

Splitting all 34 data-bearing countries would cost the same three fits and buy
precision that is not the binding constraint: the model differences are already
determined to a standard error of 0.003 to 0.006 in excess mean squared error.
The other 40 countries (17,588 records) stay wholly in training in every fold,
as do the other two blocks of each split country, so every country intercept is
identified throughout. That is the entire purpose of blocking rather than
holding out countries.

**Construction** (`R/validation_blocks.R`). Within each split country, its
cells are partitioned into two contiguous blocks that are **large and even in
area**, subject to a floor on the bioassays each one carries. Two families of
cut are tried — slabs (project the cells onto a direction and cut across it) and
sectors (cut on the bearing from the country's area centroid) — at 36
orientations each, and the winner maximises the smallest block's area. Slabs
win in all six countries.

Area means the number of land cells of the mask falling in the block, not the
convex hull of its data-bearing cells: the cut is applied to every land cell of
the country, so the three areas are directly comparable and sum to the country.
Because the three sum to a constant, maximising the smallest block's area is
exactly what "large and even" means — a cut that makes one block small
necessarily makes another large.

This is the third objective tried, and the reasoning behind it is worth keeping.

- *Balance the records, maximise the smallest block's area.* Does not work.
  Rotating a cut barely changes the areas, so the objective is nearly flat
  across orientations and chooses on density noise; the cuts it picked were thin
  slices, with a fifth of held-out records within 15 km of training data.
- *Balance the records, maximise separation directly.* Better, but bounded by
  the record constraint. Bioassay effort is wildly uneven in space, so cutting
  at the record terciles puts a boundary straight through the densest cluster,
  and the block holding a third of the records occupies a small area. Short
  separations then occur exactly where most of the data is.
- *Even areas with a record floor.* Lets a dense cluster sit whole inside one
  block, which keeps both blocks large and lengthens the separation for the
  sparser one. The cost is an uneven split of records between folds, which
  costs nothing overall: every record is still held out exactly once, so only the
  per-fold counts differ.

Floors, so that favouring area cannot leave a fold too thin to score: each block
must carry at least 150 bioassays or 12% of its country's records, whichever is
larger, and at least 10 data-bearing cells. The floor binds only in Kenya.

No buffer is applied. The blocks are large enough that a few cells near a
boundary cannot carry the result, and a buffer would remove training records
from exactly the countries whose intercepts this design exists to keep
identified.

**What the folds look like.** 5,382 / 3,312 held-out assays — 8,694 in total,
over the same six countries and the same "2010 and later" test window. All nine insecticide types
appear in every fold. Every blocked record from 2010 on is held out exactly
once; records at a held-out cell from before 2010 are dropped from the
experiment rather than returned to training, which is what the national folds
also do with pre-2010 records of a held-out country. No pixel is in both sets of
any fold.

Per country, with the smallest block's area, how even the three areas are, and
the worst departure of a block's record share from a third:

| country | family | smallest block | area evenness | worst record share |
|---|---|---|---|---|
| Senegal | slab | 68,900 km² | 1.00 | 0.19 |
| Côte d'Ivoire | slab | 109,300 km² | 1.00 | 0.12 |
| Kenya | slab | 130,600 km² | 0.43 | 0.41 |
| Nigeria | slab | 309,500 km² | 1.00 | 0.15 |
| Tanzania | slab | 319,100 km² | 1.00 | 0.19 |
| Ethiopia | slab | 383,700 km² | 1.00 | 0.16 |

Bioassay-weighted distance from a held-out pixel to the nearest training pixel,
pooled over the folds:

| | min | 10% | 25% | median | 75% | 90% | max |
|---|---|---|---|---|---|---|---|
| block folds | 5 | 19 | 46 | **79** | 132 | 178 | 375 |
| interpolation fold, for comparison | 14 | 21 | 25 | 36 | 48 | 66 | 466 |

7% of held-out assays sit within 15 km of training data and 14% within 25 km,
against 18% and 27% under the first objective. Separation is capped at a few
hundred km rather than the ~800 km the national folds reach, because a held-out
block is surrounded by training data in neighbouring countries. That is the
price of keeping the country intercepts identified, and it is the right price.

Per-country median separation: Nigeria 132 km, Ethiopia 98, Tanzania 93, Côte
d'Ivoire 83, Senegal 58, **Kenya 53**. Kenya is still the weakest, because it
holds the majority of its records in one cluster west of Lake Victoria; the area
objective at least keeps that cluster whole in one block, which lifted Kenya's
median separation from 26 km to 53 and halved the share within 15 km. The
distance-stratified reporting of §4.4 handles what remains: those records are
not discarded, they are read at the distance they actually represent.

One cosmetic note: a border cell is assigned to whichever country holds most of
its records, so fold 3's held-out set touches a seventh country through one such
cell.

**Sampling settings.** Unchanged from §3. Pinning each country's initial
condition in every fold removes the parameter that was previously unidentified,
so convergence should be better than the national folds, not worse. Check Rhat
on the first completed fold before launching the other two. If it has not
improved, a longer warmup is available for these folds, because they are a new
experiment and do not have to match the sampling settings of the retained
national folds.

## 6. Change-based scoring, and two analyses that are not here

`validation_change.R` scores predicted change between the before-window and the
holdout, per forecast origin, which is the quantity the forecasting experiment is
actually about. The before-window predictions it needs are saved with each fold
by `fit_validation_fold.R`; because sampling cannot be resumed across sessions
(§8), that has to be in place before a fit starts and cannot be added afterwards.

Two supporting analyses were written and are not in this branch, having served
their purpose:

- the power analysis over window length and origin, which chose the five-year
  window and showed 2014 to be the strongest feasible origin;
- skill against separation from the training data, which was a negative result:
  the joint test gave p = 0.058 for the dynamical model and 0.264 for the
  nearest survey. Distance has an intraclass correlation of 0.91 by pixel, so a
  pixel random effect absorbs the identifying contrast rather than controlling
  for it, and binned weighted least squares with a pixel-cluster bootstrap is
  what the design supports. The raw distance bins that preceded it were
  confounded in exactly that way and are also gone.

Both live in the stub branch if they are wanted again.

## 7. Cost

| fits | wall clock, two at a time |
|---|---|
| ~~temporal forecasting, corrected split~~ | done, 62 h |
| ~~two sub-national block folds~~ | done |
| forecast origins 2014 and 2018, five-year windows | ~62 h |

Per-fold cost is set by the dynamics graph, which is solved over all unique
cells × years × types regardless of how many rows enter the likelihood, so a
smaller K does not make each fold cheaper — it runs fewer of them. Four fits two
at a time is about 124 h against 250 h for the September run.

Running all four concurrently at two threads each is worth testing on one fold
first, given how weakly threads scale; memory is the constraint, since each
process holds the full dynamics graph.

The seven retained folds in `outputs/cv_draws/` are not refitted. They can be
slimmed offline — `rho_draws` is stored expanded to `n_draws × n_test` when only
`n_draws × n_classes` is distinct, about 1.8 GB across the folds — but the
scoring path already accepts either layout, and rewriting irreplaceable 62 h
files for disk space is not obviously worth the risk.

## 8. What survives a session, and what does not

Established by experiment, not assumption:

| operation on a reloaded `draws` object | result |
|---|---|
| `calculate(target, values = draws)` | **works** — targets are recoverable from `attr(draws, "model_info")` and the graph is re-traced |
| `extra_samples(draws, n_samples = ...)` | **fails** — `"object is from previous session and is now invalid"` |

The sampler state is bound to the session that created it, and redefining the
model produces new nodes the draws cannot attach to.

Consequences for planning:

- **A longer run must be requested up front**, through `warmup` and `n_samples`.
  Sampling cannot be topped up afterwards.
- **Folds to be compared must share their sampling settings.** Refitting one
  fold with longer warmup means refitting all of them.
- **Each saved fold therefore keeps `draws` and `prediction_arrays`.** Nothing in
  the committed pipeline reads them back, but they are what lets a finished fold
  produce a new prediction target — new cell-years, new aggregations, the
  before-window predictions of §6 — without re-running 62 h of MCMC. They are
  also most of each file's size. The four defunct folds in
  `outputs/cv_draws_defunct/` have no `draws` object and so cannot be used with
  greta's prediction interface at all, which is what that costs.

## 9. Pitfalls already hit, worth not repeating

- `calculate(..., nsim = n)` returns an **independent resample** of the
  posterior. It preserves joint structure across quantities, so it is valid for
  prediction, but it destroys MCMC ordering, so effective sample size cannot be
  recovered from it. Use `calculate(values = draws)` without `nsim`.
- `coda::effectiveSize()` on 20 draws returns roughly 20. Any ESS computed from
  a short run is meaningless; several tuning conclusions were drawn from such
  numbers and had to be withdrawn.
- `future.callr` buffers worker stdout until the future resolves, so a long run
  under it is invisible. Hence one process per fold.
- Functions called from a worker must have every dependency passed explicitly.
  `codetools::findGlobals(fit_fold, merge = FALSE)$variables` catches free
  variables; `n_unique_cells` was missing this way and cost a run.
- Verify that string edits to these scripts actually applied. A silently
  non-matching replacement cost a second run.
- Arm a monitor on the logs and check that the monitor itself is alive. A run
  died and went unnoticed for 14 h because the watcher had exited days earlier.
- `pgrep -f "validation_metrics.R"` matches the shell waiting on it as well as
  the R process. A completed run was reported as still running for six days on
  the strength of that.

## 10. Superseded artefacts

Parked in `outputs/cv_draws_defunct/`, not read by anything, kept as the record
of what was tried:

- the six leave-one-country-out folds (§4.2);
- the three-year 2020 forecasting fold, whose training set was drawn with
  `year_start <= cut`, so the cut year appeared in both training and holdout,
  and whose window landed on the one pause in twenty years of decline;
- `cv_draws_leaky_forecast/`, the same fold before the leak was fixed;
- four folds from the earlier two-chain run, with no `draws` object.

`validation_folds.R` defines `validation_experiments`, and both
`run_one_fold.R` and `validation_metrics.R` refuse anything outside it, so a
stray fold cannot be scored by accident.

Deleted on this branch: `predictive_validation.R` and
`dynamic_predictive_validation.R`, which held the plug-in deviance path and
three copies of the model definition; and `validation_metric_eval.R`, the
simulation study that supported the choice of metric.

## 11. What the validation found, and what follows

The dynamical model beats the nearest-recent-survey baseline on sub-national
spatial extrapolation (+8.8 [+3.3, +14.7] percentage points of variance
explained, paired within bootstrap replicates) but not on forecasting
(−3.6 [−8.5, +1.5]). Diagnosis, from the before-window predictions saved with
each fold:

- the model has the direction of local change right — 77–83% sign agreement at
  the 2014 origin among groups whose observed change exceeds its own noise — but
  predicts two to three times too much decline, and the signed error grows with
  forecast horizon;
- selection is linear in each covariate, and the net-use response saturates in
  the raw data (#23);
- resistance can only increase in the model: there is no fitness cost or decay
  term (#24).

Both fixes need a refit and are out of scope here.
