# load packages

# # NOTE: patchwork is bugging out with recent ggplot, use 3.4.4:
# remotes::install_version("ggplot2",
#                          version = "3.4.4",
#                          repos = "http://cran.us.r-project.org")
# # and older tidyterra bc dependency
# remotes::install_version("tidyterra",
#                          version = "0.4.0",
#                          dependencies = FALSE,
#                          repos = "http://cran.us.r-project.org")
# # greta >= 43f9c52 (0.6.0.9000), no greta.dynamics; installed in a separate
# # library and conda environment as in doc/cv_run_plan.md, section 1:
# remotes::install_github("greta-dev/greta@282944f",
#                         lib = "~/R/greta06-lib")

library(tidyverse)
library(readxl)
# python is not started here; scripts that use greta call start_greta()
# (R/greta_setup.R) before sourcing this file
library(greta)

library(lme4)
library(terra)
library(tidyterra)
library(tidygeocoder)
library(future)
library(future.apply)
library(future.callr)
library(DHARMa)
library(Hmisc)
library(patchwork)
library(extraDistr)
library(ggtext)
library(geodata)
library(sf)
