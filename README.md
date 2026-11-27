
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Spatiotemporal modelling of the spread of insecticide resistance in Africa

This repo contains models and code to model the evolution and spread of
insecticide resistance among malaria vectors in Africa. It uses a
semi-mechanistic model to predict future levels of resistance across the
continent, fitting closely to resistance data.

The published maps and time series are of the two-stage model: the
dynamical model, corrected by a second-stage geostatistical model of its
residuals (`R/two_stage_maps.R`, `doc/two_stage_plan.md`). The
cross-validation figures and the dynamical model’s own parameters
(e.g. the covariate effects) show the dynamical model separately.

Running order of scripts:

*To be tidied: this may not be the correct order*

1.  packages.R
2.  functions.R
3.  prep_admin.R
4.  prep_country_borders.R
5.  prep_bioassays.R
6.  prep_rasters.R
7.  prep_net_use_pyrethroid.R: also net use of PBO and dual-AI nets
8.  calculate_ingredient_fractions.R
9.  fit_model_glm.R
10. fit_model.R
11. drop_stuck_chains.R: leaves out the full fit's stuck chains
12. illustrate_validation.R
13. mtm_ir_explore.R
14. ploidy_demo.R
15. predict.R
16. two_stage_maps.R: `prepare`, `fit` 1-9, `map` 1-9, `map
    llin_effective`, then `figures`
    (`doc/two_stage_plan.md`). The figure scripts below that show
    predictions read its outputs
17. summarise_model_fit.R
18. visualise_colony_net_bioassay.R
19. visualise_data.R
20. fig_admin_maps.R
21. fig_baseline_susceptibility.R
22. fig_bioassay_maps.R
23. fig_covariate_effects.R
24. fig_covariate_maps.R
25. fig_data_distribution.R
26. fig_illustrate_bioassay_variability.R
27. fig_internal_validation.R
28. fig_ir_maps.R
29. fig_ir_map_2025.R
30. fig_ir_map_2000_2030.R
31. fig_ento_epi_impact.R: draws from the two-stage posterior (set
    `IR_CUBE_MAP_WORKERS` and preload OpenBLAS, as for
    `two_stage_maps.R`)
32. fig_temporal_preds_data.R
33. fig_temporal_preds_net_use.R
