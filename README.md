
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Spatiotemporal modelling of the spread of insecticide resistance in Africa

This repo contains models and code to model the evolution and spread of
insecticide resistance among malaria vectors in Africa. It uses a
semi-mechanistic model to predict future levels of resistance across the
continent, fitting closely to resistance data.

Running order of scripts:

*To be tidied: this may not be the correct order*

1.  packages.R
2.  functions.R
3.  prep_admin.R
4.  prep_country_borders.R
5.  prep_bioassays.R
6.  prep_rasters.R
7.  prep_net_use_pyrethroid.R
8.  calculate_ingredient_fractions.R
9.  fit_model_glm.R
10. fit_model.R
11. chain_mode_check.R, then drop_stuck_chains.R with DROP_CHAINS from
    its `.drop` file: leaves out chains in the minor mortality-floor mode
    (#37)
12. illustrate_validation.R
13. mtm_ir_explore.R
14. ploidy_demo.R
15. predict.R
16. summarise_model_fit.R
17. visualise_colony_net_bioassay.R
18. visualise_data.R
19. fig_admin_maps.R
20. fig_baseline_susceptibility.R
21. fig_bioassay_maps.R
22. fig_covariate_effects.R
23. fig_covariate_maps.R
24. fig_data_distribution.R
25. fig_illustrate_bioassay_variability.R
26. fig_internal_validation.R
27. fig_ir_maps.R
28. fig_temporal_preds_data.R
29. fig_temporal_preds_net_use.R
