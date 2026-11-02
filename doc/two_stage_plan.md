# Two-stage correction: run sheet

Issue idem-lab/ir_cube#21. Stacked on #12 (`posterior-predictive-validation`). Code removed during development is kept at two tags:

- `two-stage-full-e31b2b5`: the diagnostics, the rejected variants (`omega_u`, `m_ref=loo`, the other meshes) and the separate simulation checks.
- `two-stage-options-6edcf0d`: the survey effect s, joint ρ and damped ξ.
- `archive/two-stage-is` (on goldingn/ir_cube): the importance-sampling check of the CCD integration over θ, reported on #58.

## Final model

Per insecticide type, on the logit scale, for assay i at pixel s and year t:

λ_i = m_i + ω(s) + ξ(s, t) + u(s, t) + p(s) + sy(study(i), t)

- m: the dynamical model's posterior mean logit at the training assays (`m_ref`).
- ω: a static Matérn field, correcting the initial conditions.
- ξ(s, t) = Σ η: the accumulation from ξ(·, t0 = 1995) = 0 of annual anomalies η, AR(1) in time (φ) with Matérn innovations. It is undamped, so a forecast holds the accumulated deviation and its mean plateaus at a rate set by φ.
- u (per pixel-year, SD τ), p (per pixel, SD σ_p) and sy (per study-year, SD σ_sy; #55, with the study label of #54): iid noise.
- Meshes (`build_correction_meshes()`):
  - In the data region, ω has ≤ 5000 nodes (15 km cutoff, inner edges ≤ 150 km) and ξ has ≤ 2500 nodes.
  - Around both, a coarse outer region (edges ≤ 800 km) extends 1500 km beyond the data and the prediction mask, so every map cell is inside the mesh. Callers pass the mask's coordinates (`prediction_mask_coords()`, from `data/clean/raster_mask.tif`).
- Priors (`correction_priors()`, #34): normal on log range for ω and η, with 95% of the prior in 50–1,500 km (the form of R-INLA's non-PC SPDE priors); normal on logit φ, with 99% of the prior for the persistence 1/(1 − φ) in 3–10 years; half-normal SDs, P(σ_ω > 1) = 0.05, P(σ_η > 0.05) = 0.05, and P(sd > 2) = 0.05 for the noise terms.
  - Persistence: without the study-year term, ξ fitted survey noise and φ ran to 0.01–0.07, so ξ was close to a random walk. With sy, the data identify φ: for Deltamethrin on the full data, a flat prior on logit φ gives φ = 0.46 (persistence 1.8 years; profile interval roughly 1.3–2.6 years), and the 3–10 year prior raises it to φ = 0.65 (2.9 years) at a cost of 2 nats. The prior is kept: resistance does not swing from year to year. The paper reports both.
- Fitting (`fit_correction()`, #34):
  - The hyperparameters θ at the penalised maximum of TMB's Laplace approximation to the marginal likelihood of the beta-binomial counts, at the external per-type ρ (`outputs/bioassay_rho_hierarchical.csv`), started from the priors' medians.
  - The posterior of θ by INLA's central composite design (81 points for the 8 hyperparameters); predictions mix the conditional latent posteriors N(mode(θ_j), H(θ_j)⁻¹) over the points.
- Cut posterior: the dynamical model is never updated. Each dynamical draw shifts the latent mode by H⁻¹A′W(m_ref − m_draw), linearised at the mode (W the beta-binomial weights), and is paired with its own latent draw.

**Target.** The prediction target is m + ω + ξ: that is what the maps, the published figures and the supplement figures show. u, p and sy are observation-level noise. They are in the fitted model, so that this noise stays out of ω and ξ. They enter predictions only of new assays, and only as fresh draws, N(0, τ²) per pixel-year, N(0, σ_p²) per pixel and N(0, σ_sy²) per study-year, never as their fitted values; fresh u and p also enter the maps' n\*_local. This reverses #21, which put u in the target, and an earlier version of this run sheet, which put p in it. Cross-validation scores the predictive distribution of held-out assays, so it includes fresh noise (`predict_correction()`).

## Headline results

These are the PQL fit of #21, without sy, on the bioassays before the deduplication of #63. The final model's results (Laplace + CCD, with sy) will replace them after the final run on #63's data (milestone 6).

Scored as in #12 (`R/validation_metrics.R` → `outputs/cv_scores.csv`, `cv_summary.csv`), where the final model is `two_stage`. Variance explained and its ceiling are from `R/validation_level.R` → `outputs/cv_level_skill.csv` and `_by_insecticide.csv`. The paired differences are from `R/two_stage_metrics.R` → `outputs/two_stage/cv_headline_two_stage.csv` and `cv_horizon_two_stage.csv`.

- Log score: mean beta-binomial log predictive density at the per-type ρ.
- CRPS: on the mortality scale.
- Variance explained: 100 (1 − MSE / Var(y)). The ceiling is the `noise_floor_mse()` bound.
- Coverage: the expected coverage of the central interval over the PIT randomisation (`expected_coverage()`), in every table. #12's tables used the mean over 100 PIT replicates (`cv_summary.csv`) or one replicate (`cv_by_fold.csv`); the switch changes #12's models' coverage by at most 0.002 in `cv_summary.csv` and 0.009 in `cv_by_fold.csv`.
- Intervals: 95% pixel bootstrap.
- NN: the nearest recent survey (k = 1). NN oracle: the nearest surveys at the k that minimises held-out error.

**Interpolation** (n = 1045; ceiling 78.0%)

| model | log score | CRPS | variance explained (%) | cover 50 / 95 |
|---|---|---|---|---|
| dynamical | −3.863 | 0.1353 | 37.0 [26.0, 46.3] | 0.298 / 0.787 |
| intercept | −4.233 | 0.1497 | 27.0 [22.7, 30.3] | 0.306 / 0.721 |
| NN | −4.010 | 0.1355 | 33.7 [20.5, 44.1] | 0.343 / 0.777 |
| NN oracle | −3.614 | 0.1107 | 52.0 [43.0, 59.7] | 0.390 / 0.833 |
| **two-stage** | **−3.372** | **0.1028** | **57.1 [49.2, 63.8]** | 0.448 / 0.921 |

**Spatial blocks 1+2** (n = 8694; ceiling 79.5%)

| model | log score | CRPS | variance explained (%) | cover 50 / 95 |
|---|---|---|---|---|
| dynamical | −4.327 | 0.1587 | 22.7 [17.9, 26.9] | 0.292 / 0.736 |
| intercept | −4.585 | 0.1640 | 19.3 [17.1, 21.3] | 0.304 / 0.709 |
| NN | −4.678 | 0.1685 | 13.8 [7.3, 19.8] | 0.293 / 0.691 |
| NN oracle | −4.189 | 0.1415 | 34.1 [29.3, 38.5] | 0.323 / 0.752 |
| **two-stage** | **−3.663** | **0.1327** | **36.9 [32.1, 41.3]** | 0.440 / 0.911 |

**Forecasting 2014+2018** (n = 14018; ceiling 82.8%)

| model | log score | CRPS | variance explained (%) | cover 50 / 95 |
|---|---|---|---|---|
| dynamical | −4.359 | 0.1619 | 26.2 [20.6, 31.4] | 0.311 / 0.732 |
| intercept | −4.703 | 0.1847 | 13.4 [10.8, 15.9] | 0.268 / 0.663 |
| NN | −4.322 | 0.1549 | 29.8 [25.7, 33.4] | 0.311 / 0.731 |
| NN oracle | −3.869 | **0.1281** | **48.3 [45.6, 50.9]** | 0.332 / 0.788 |
| **two-stage** | **−3.614** | 0.1315 | 42.7 [38.7, 46.6] | 0.457 / 0.900 |

**Two-stage minus dynamical**, paired:

| experiment | Δ log score | Δ CRPS (×10⁻³) | Δ variance explained (points) |
|---|---|---|---|
| interpolation | +0.491 [+0.366, +0.618] | −32.4 [−41.1, −23.4] | +20.1 [+13.5, +28.6] |
| blocks 1+2 | +0.664 [+0.597, +0.734] | −26.0 [−29.7, −22.3] | +14.2 [+11.6, +16.8] |
| forecasting 2014+2018 | +0.746 [+0.672, +0.825] | −30.4 [−34.6, −26.6] | +16.6 [+14.0, +19.2] |
| forecasting 2014 | +0.772 [+0.690, +0.860] | −29.7 [−34.3, −25.6] | +17.5 [+14.4, +20.8] |
| forecasting 2018 | +0.682 [+0.585, +0.784] | −32.3 [−38.4, −26.8] | +15.6 [+11.9, +19.0] |

**By forecast horizon** (both origins pooled):

| horizon (years) | n | variance explained | Δ log score vs dynamical | Δ explained vs dynamical | cover 95 | NN oracle explained |
|---|---|---|---|---|---|---|
| 1 | 4228 | 51.1 [46.0, 55.5] | +0.62 [+0.53, +0.71] | +18.5 [+14.8, +22.4] | 0.914 | 49.0 |
| 2 | 3176 | 51.9 [46.0, 57.1] | +0.61 [+0.51, +0.73] | +12.0 [+8.2, +16.0] | 0.914 | 55.5 |
| 3 | 2247 | 39.2 [27.5, 49.3] | +0.55 [+0.44, +0.66] | +10.3 [+6.2, +14.3] | 0.883 | 47.0 |
| 4 | 2222 | 30.4 [21.1, 38.6] | +1.05 [+0.86, +1.23] | +19.4 [+14.4, +24.8] | 0.882 | 45.0 |
| 5 | 2145 | 30.9 [21.1, 39.2] | +1.09 [+0.93, +1.26] | +23.5 [+19.4, +28.2] | 0.887 | 40.8 |

- The two-stage model has the best log score in every experiment.
- On variance explained it beats the NN oracle on interpolation (+5) and blocks (+3), and trails it on forecasting (−6), mostly at horizons of 3–5 years.
- Its 95% coverage is 0.89–0.92. Every other model's is 0.66–0.83.
- These are from the rerun of the five folds with #21's final code (meshes covering the prediction mask, `m_ref` unclamped, seeds per fold and type). Against the earlier draws, the posterior mean predictions correlate at 0.9998–0.9999 per fold (mean |difference| 0.2–0.4 points of mortality), and the headline scores move by at most 0.002 in log score, 0.2 points of variance explained and 0.002 in coverage.

## Decisions

Scores are from the cross-validation of the time. A decision's numbers can predate later choices; the tables above are #21's final model.

| option | result | adopted | code |
|---|---|---|---|
| ξ (`omega_xi_u` vs `omega_u`) | Forecasting +0.33 log score, +10 points explained; blocks +0.07, +6. Without ξ, forecasting gains nothing over the dynamical model in variance explained. | yes | full tag |
| Meshes omega5000_xi2500 (vs four coarser configurations) | Best everywhere. Blocks +0.055 log score and +2.1 points over the base mesh; ~11 GB per fit. | yes | full tag |
| `m_ref` = posterior mean (vs PSIS leave-out mean) | Leave-out shift 0.04–0.05 residual SD; double counting negligible. | posterior mean | full tag |
| Stage B, PQL on the beta-binomial (vs Gaussian empirical logit) | +0.06 to +0.19 log score; 95% coverage 0.86–0.88 → 0.90–0.92; removes the excess of 100% assays in the training residuals. | yes, then replaced by Laplace + CCD | — |
| Laplace + CCD over θ (vs PQL with θ plugged in), #34 | On blocks 1 and forecasting 2014 for Deltamethrin, Fenitrothion and Pirimiphos-methyl, every latent × θ method was within 0.06 nats per assay of the PQL fit, and Laplace with θ plugged in within 0.01 of Laplace with importance sampling over θ. CCD integrates over θ without Monte Carlo error. | yes | — |
| Pixel effect p | +0.003 to +0.005 log score; about 30% fewer ω hotspots. | yes, as observation noise | — |
| Survey effect s | +0.006 to +0.013 log score, all from a wider held-out predictive distribution; none in the mean or the map. | no; superseded by sy | options tag |
| Study-year noise sy, #55 | Variance explained against the current model's map, blocks 1: +3.8 points pooled, +6.0 for Deltamethrin. Deltamethrin σ_η 0.24 → 0.11–0.15, persistence 2.4 → 3.1–3.6 years; ξ's 95% band at a pixel without data 8.1–9.4 → 5.0–5.6 logit wide. | yes, as observation noise | — |
| Joint ρ | Estimated ρ is a median of 0.35 × the replicate estimate, with the noise moved into u (a PQL bias also seen in simulation); −0.04 to −0.07 log score at its own ρ. | no | options tag |
| Damped ξ (ψ) | Forecasting −0.038 log score, −0.12 and −0.16 at 4–5 years; spatial experiments unchanged. | no | options tag |
| u and p in the target | u and p are observation noise, as fresh draws for new assays only (see Final model). | target m + ω + ξ | — |
| Diagnostics: covariate and hotspot structure, training residuals at 0% and 100%, mesh comparison, stage-one effective parameters | Finished. They motivated stage B, p and the meshes. | removed | full tag |

## Open issues

- #28: the temporal panel of #12's variance-explained bar charts should score change, not level.
- Residual narrowness. Under the PQL fit, interior training residuals were narrower than their working variance implied (SD 0.70 vs 0.95), and a lower ρ did not explain it (joint ρ). Not rechecked under Laplace + CCD.

## How to run

Use OpenBLAS (`LD_PRELOAD=.../libopenblas.so.0 OPENBLAS_NUM_THREADS=3`), because reference BLAS is ~10× slower. Compile the template with TMBad, which `load_correction_template()` does.

1. **Simulation check:** `Rscript R/check_two_stage.R`.
2. **Folds** (`R/run_two_stage_folds.R`), for each of `spatial_interpolation all`, `spatial_blocks 1|2` and `temporal_forecasting 2014|2018`:
   - `Rscript R/run_two_stage_folds.R <experiment> <fold> <k>` for k = 1–9: one process per type. Each reads the paired dynamical logit draws that `R/run_one_fold.R` saved with the fold (`outputs/cv_draws/assays/`), not the fold. On 8-vCPU pods, two types at a time, a type took 15–45 min to fit and 5–35 min to predict, with peak memory up to 12 GB. A fit that fails or does not converge stops its type, and the fold is not assembled without it.
   - `... <experiment> <fold> assemble`: writes `outputs/cv_draws/two_stage__<experiment>__<fold>.rds` and `outputs/two_stage/fit_summary.csv`.
3. **Metrics:** `Rscript R/validation_metrics.R`, then `Rscript R/two_stage_metrics.R`; the rest of the scoring as in `doc/cv_run_plan.md`, section 3.
4. **Maps** (`R/two_stage_maps.R`):
   - `prepare`: loads `temporary/fitted_model.RData`, checks that its data and covariates are those the current scripts build, and saves per-type dynamical draws (2 min, 7 GB).
   - `fit <k>`: fits type k to all data and saves `outputs/two_stage/maps/<type>/fit.rds` and `hyperparameters.csv`. On 8-vCPU pods, two types at a time: 11–52 min and 6–11 GB per type (Deltamethrin 37 min, 11 GB).
   - `map <k>` for k = 1–9: writes the posterior mean and SD of mortality for every year 1995–2030 to `outputs/two_stage/ir_maps/<type>/ir_<year>_susceptibility(_sd).tif`, the layout of the dynamical model's `outputs/ir_maps` (`R/predict.R`), from 1000 paired draws in batches of 100 (Monte Carlo SE of a mean map: posterior SD / √1000, at most 1.3 percentage points where the SD reaches its maximum of 40). At the panel years it also writes the mean and SD with fresh u and p added to the same draws (`local_mean`, `local_sd`), which only n\*_local reads. Set `IR_CUBE_MAP_WORKERS` to map groups of cells in forked workers. On 8-vCPU, 64 GB pods a type took 41–110 min at 4 workers; each worker holds ~5 GB, and 4 workers ran out of memory at 32 GB.
   - `map llin_effective`, after the three pyrethroids: llin_effective's mean and SD as the sums of Alpha-cypermethrin's, Deltamethrin's and Permethrin's, weighted by `temporary/ingredient_weights.RDS`. Summing the SDs takes the pyrethroids' errors as perfectly correlated, the largest SD their marginal SDs allow (#40). A few minutes.
   - `figures`: writes `figures/two_stage/`. Per type and for llin_effective: the posterior uncertainty as n\*, the number of bioassays of 100 mosquitoes the prediction is worth (#40), with bioassay noise from ρ̂ only. u and p could be data error or real local variation the model cannot predict, and the data can't tell which, so n\* is mapped under both: n\*_pop treats them as data error (posterior of m + ω + ξ; an upper bound) and n\*_local as real (with fresh u and p; a lower bound). n\*_CV (`R/validation_level.R`) is the out-of-sample check: at held-out pixel-years it should be close to n\*_local, apart from survey error and any miscalibration. Per type also the correction ω + ξ and the difference from the dynamical model, and the all-types correction panel for 2020 and the hyperparameters.
5. **Supplement figures:** `Rscript R/fig_two_stage_components.R`. Reads the saved Deltamethrin fit and its draws, and writes `supp_components_realisations` and `supp_components_intervals` (captions: `doc/two_stage_supplement_captions.md`).
6. **Published figures.** The figure scripts that show predictions use the two-stage model: `R/fig_ir_maps.R`, `fig_ir_map_2025.R`, `fig_ir_map_2000_2030.R` and `fig_ento_epi_impact.R` read the rasters of step 4 (`ir_map_files()` in `R/functions.R`), and `R/fig_temporal_preds_data.R`, `fig_temporal_preds_net_use.R`, `summarise_model_fit.R` and `fig_internal_validation.R` draw from the saved fits at the data cells (`R/two_stage_predictions.R`). Each of the latter takes 4–11 min; `fig_ento_epi_impact.R` takes ~1 h.
