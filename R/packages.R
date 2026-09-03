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
# # greta and greta.dynamics versions. The greta_2 branch this used to pin no
# # longer exists upstream; it was merged into main.
# #
# # Do not upgrade greta to 0.6.0: it cannot sample a model whose likelihood
# # depends on iterate_dynamic_function() output, which this model's does. The
# # failure is raised at mcmc(), not at model construction, and it is not
# # specific to this model. Holding greta.dynamics fixed and moving greta from
# # 0.5.0 to 0.6.0 flips a minimal example from sampling to failing.
# remotes::install_version("tensorflow", version = "2.16.0")
# remotes::install_github("njtierney/greta@4cc989f")            # 0.5.0.9000
# remotes::install_github("greta-dev/greta.dynamics@db7df31")   # 0.2.2

library(tidyverse)
library(readxl)
library(greta)
library(lme4)
library(terra)
library(tidyterra)
library(greta.dynamics)
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
