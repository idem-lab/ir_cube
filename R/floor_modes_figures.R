# Figures and tables for the two mortality-floor modes (#14), from
# outputs/floor_modes_<label>.rds (R/floor_modes.R).
#   Rscript R/floor_modes_figures.R round2
# Writes figures/floor_modes_{trajectories,residuals,loglik_contrib}.png,
# outputs/floor_modes_loglik_contrib.csv, outputs/floor_modes_places.csv,
# and prints the residual trends with cumulative net use.
args <- commandArgs(trailingOnly = TRUE)
label <- if (length(args)) args[1] else "round2"
suppressMessages({library(dplyr); library(tidyr); library(ggplot2)
  library(patchwork); library(mgcv)})
x <- readRDS(sprintf("outputs/floor_modes_%s.rds", label))
df <- x$df
mode_names <- c(low = sprintf("LOW (chains %s)", paste(x$chains$low, collapse = ",")),
                high = sprintf("HIGH (chain%s %s)",
                               ifelse(length(x$chains$high) > 1, "s", ""),
                               paste(x$chains$high, collapse = ",")))
mode_cols <- setNames(c("#0072B2", "#E69F00"), mode_names)
floor_high <- mean(x$modes$high$floor)
floor_low <- mean(x$modes$low$floor)

a <- bind_rows(lapply(x$modes, `[[`, "assay")) |>
  select(mode, row, pred, loglik, pearson) |>
  pivot_wider(names_from = mode, values_from = c(pred, loglik, pearson))
d <- bind_cols(df, a[order(a$row), ]) |>
  mutate(dll = loglik_low - loglik_high,
         country = x$countries[country_id], type = x$types[type_id],
         obs = died / mosquito_number,
         period = cut(year_start, c(0, 2009, 2016, 3000),
                      c("pre-2010", "2010-2016", "2017+")))

# 1. informative places
places <- d |> group_by(country_id, country) |>
  summarise(assays = n(), pyrethroid_assays = sum(group == "Pyrethroids"),
            ddt_assays = sum(group == "DDT"),
            net_use_2010_on = mean(nets[year_start >= 2010]),
            max_cum_nets = max(cum_nets),
            dll_low_minus_high = sum(dll), .groups = "drop") |>
  mutate(high_mode_floor_correlate = country_id %in% c(22, 7, 27, 10, 30, 2),
         rank_assays = rank(-assays), rank_net_use = rank(-net_use_2010_on)) |>
  arrange(-abs(dll_low_minus_high))
write.csv(places, "outputs/floor_modes_places.csv", row.names = FALSE)

# 2. trajectories
show <- c("Burkina Faso", "Benin", "Nigeria", "Ethiopia",
          "Madagascar", "Rwanda", "Tanzania", "Kenya")
traj <- bind_rows(lapply(x$modes, `[[`, "traj")) |>
  mutate(country = x$countries[country_id], year = x$years[year_id],
         mode = mode_names[mode]) |>
  filter(country %in% show, year <= 2025)
obs_cy <- d |> filter(country %in% show, group != "other") |>
  group_by(country, group, year = year_start) |>
  summarise(mortality = sum(died) / sum(mosquito_number),
            mosquitoes = sum(mosquito_number), .groups = "drop")
panel <- function(dd) factor(paste(dd$country, dd$group, sep = ": "),
                             levels = paste(rep(show, each = 2),
                                            c("Pyrethroids", "DDT"), sep = ": "))
traj$panel <- panel(traj); obs_cy$panel <- panel(obs_cy)
p_traj <- ggplot(traj, aes(year)) +
  geom_hline(yintercept = floor_high, colour = mode_cols[2],
             linetype = "dashed", linewidth = 0.4) +
  geom_ribbon(aes(ymin = lo, ymax = hi, fill = mode), alpha = 0.25) +
  geom_line(aes(y = mean, colour = mode), linewidth = 0.7) +
  geom_point(data = obs_cy, aes(y = mortality, size = mosquitoes),
             shape = 21, fill = "grey30", colour = "white", stroke = 0.3,
             alpha = 0.8) +
  facet_wrap(~ panel, ncol = 4, drop = FALSE) +
  scale_colour_manual(values = mode_cols, name = NULL) +
  scale_fill_manual(values = mode_cols, name = NULL) +
  scale_size_area(max_size = 5, name = "mosquitoes tested") +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = NULL, y = "mortality",
       title = sprintf("Fitted mortality under the two floor modes, %s main fit", label),
       caption = paste0(
         "Lines: posterior mean of predicted mortality averaged over the pixels with assays of that group in the country, ",
         "weighted by their assay counts;\nbands: 90% posterior interval of that average. ",
         "Points: observed country-year mortality (pooled over assays), area by mosquitoes tested. ",
         sprintf("Dashed line: HIGH-mode floor (mean %.2f); LOW-mode floor %.4f.", floor_high, floor_low),
         "\nBurkina Faso, Benin, Nigeria, Ethiopia: West and Horn; Madagascar, Rwanda, Tanzania, Kenya: East.")) +
  theme_bw(base_size = 9) +
  theme(legend.position = "top", strip.background = element_blank(),
        plot.caption = element_text(hjust = 0))
ggsave("figures/floor_modes_trajectories.png", p_traj, width = 11, height = 9,
       dpi = 150)

# 3. residuals, per mode, pyrethroids pooled and DDT
east <- c("Eastern Africa", "Southern Africa")
r <- d |> filter(group != "other") |>
  select(group, region, year_start, nets, cum_nets, pred_low, pred_high,
         pearson_low, pearson_high) |>
  mutate(west_east = ifelse(region %in% east, "East and Southern", "West, Central, North")) |>
  pivot_longer(c(pred_low, pred_high, pearson_low, pearson_high),
               names_to = c(".value", "mode"), names_sep = "_") |>
  mutate(mode = factor(mode_names[mode], mode_names))
rl <- r |> pivot_longer(c(year_start, pred, nets, cum_nets),
                        names_to = "x", values_to = "value") |>
  mutate(x = factor(x, c("year_start", "pred", "nets", "cum_nets"),
                    c("year", "predicted mortality (that mode)",
                      "net use that year", "cumulative net use since 1995")))
set.seed(1)
p_res <- ggplot(rl, aes(value, pearson, colour = mode)) +
  geom_hline(yintercept = 0, colour = "grey50") +
  geom_point(data = slice_sample(rl, prop = 0.25), alpha = 0.05, size = 0.4) +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 8), linewidth = 0.8) +
  facet_grid(group ~ x, scales = "free_x") +
  coord_cartesian(ylim = c(-1.5, 1.5)) +
  scale_colour_manual(values = mode_cols, name = NULL) +
  labs(x = NULL, y = "Pearson residual (observed above predicted > 0)",
       title = sprintf("Residuals under the two floor modes, %s main fit", label),
       caption = paste0(
         "Pearson residual: (died - n p) / sqrt(n p (1 - p) (1 + (n - 1) rho)), p the mode's posterior mean prediction, ",
         "rho its posterior mean overdispersion for the type.\n",
         "Lines and bands: GAM smooth of the residuals with its 95% interval. Points: a random quarter of assays; ",
         "y axis cut at +-1.5.")) +
  theme_bw(base_size = 9) +
  theme(legend.position = "top", strip.background = element_blank(),
        plot.caption = element_text(hjust = 0))
ggsave("figures/floor_modes_residuals.png", p_res, width = 11, height = 6,
       dpi = 150)

# 4. saturation: residual trend with cumulative net use, by mode, group and
# region; GAM of the residual on s(cum_nets), weighted equally per assay
cat("\nmean Pearson residual by cumulative net use bin:\n")
sat <- r |> mutate(cum_bin = cut(cum_nets, c(-0.01, 1, 2, 3, 5, 11))) |>
  group_by(group, west_east, mode, cum_bin) |>
  summarise(n = n(), resid = round(mean(pearson), 2), .groups = "drop") |>
  pivot_wider(names_from = cum_bin, values_from = c(resid, n))
print(as.data.frame(sat), width = 200)
cat("\nGAM pearson ~ s(cum_nets), fitted at cum_nets 1 and 6, and deviance explained:\n")
for (g in unique(r$group)) for (we in unique(r$west_east)) for (m in levels(r$mode)) {
  rr <- filter(r, group == g, west_east == we, mode == m)
  fit <- gam(pearson ~ s(cum_nets, k = 6), data = rr)
  pr <- predict(fit, data.frame(cum_nets = c(1, 6)))
  cat(sprintf("%-12s %-22s %-20s n=%5d  at 1: %+.2f  at 6: %+.2f  dev expl %.1f%%  p=%.2g\n",
              g, we, m, nrow(rr), pr[1], pr[2], 100 * summary(fit)$dev.expl,
              summary(fit)$s.table[1, 4]))
}

# log likelihood contributions
tab <- bind_rows(
  d |> group_by(dimension = "country", level = country) |>
    summarise(n = n(), dll = sum(dll), .groups = "drop"),
  d |> group_by(dimension = "type", level = type) |>
    summarise(n = n(), dll = sum(dll), .groups = "drop"),
  d |> group_by(dimension = "group x period",
                level = paste(group, period)) |>
    summarise(n = n(), dll = sum(dll), .groups = "drop"),
  d |> group_by(dimension = "predicted mortality bin (LOW)",
                level = as.character(cut(pred_low, c(0, .1, .3, .5, .7, .9, .98, 1)))) |>
    summarise(n = n(), dll = sum(dll), .groups = "drop"),
  d |> group_by(dimension = "observed mortality bin",
                level = as.character(cut(obs, c(-.01, .02, .1, .3, .5, .7, .9, .98, 1)))) |>
    summarise(n = n(), dll = sum(dll), .groups = "drop"),
  d |> group_by(dimension = "country x type x period",
                level = paste(country, type, period, sep = " | ")) |>
    summarise(n = n(), dll = sum(dll), .groups = "drop")) |>
  mutate(fit = label, dll = round(dll, 2)) |>
  rename(loglik_low_minus_high = dll)
write.csv(tab, "outputs/floor_modes_loglik_contrib.csv", row.names = FALSE)
cat(sprintf("\ntotal LOW - HIGH expected log likelihood: %.1f\n", sum(d$dll)))
print(as.data.frame(tab |> filter(dimension == "country x type x period") |>
                      arrange(-abs(loglik_low_minus_high)) |> head(15)))

bar <- function(dd, title) {
  ggplot(dd, aes(loglik_low_minus_high, reorder(level, loglik_low_minus_high),
                 fill = loglik_low_minus_high > 0)) +
    geom_col(width = 0.7) + geom_vline(xintercept = 0) +
    scale_fill_manual(values = c(`TRUE` = mode_cols[[1]], `FALSE` = mode_cols[[2]]),
                      labels = c(`TRUE` = "LOW fits better", `FALSE` = "HIGH fits better"),
                      name = NULL) +
    labs(x = "LOW - HIGH log likelihood", y = NULL, title = title) +
    theme_bw(base_size = 9) + theme(legend.position = "none")
}
p1 <- bar(tab |> filter(dimension == "country") |>
            slice_max(abs(loglik_low_minus_high), n = 20), "by country (top 20)")
p2 <- bar(tab |> filter(dimension == "group x period"), "by insecticide group and period")
p3 <- bar(tab |> filter(dimension == "observed mortality bin") |>
            mutate(level = factor(level, level)), "by observed mortality") +
  aes(y = level)
p4 <- bar(tab |> filter(dimension == "country x type x period") |>
            slice_max(abs(loglik_low_minus_high), n = 15), "by country x type x period (top 15)")
p_ll <- (p1 | (p2 / p3)) / p4 +
  plot_layout(heights = c(1.4, 1)) +
  plot_annotation(
    title = sprintf("Where the log likelihood difference between the floor modes comes from, %s main fit (total LOW - HIGH %.1f)",
                    label, sum(d$dll)),
    caption = paste0("Each assay's log likelihood is its betabinomial log density averaged over the mode's posterior draws ",
                     "(150 per chain); bars sum the LOW - HIGH difference.\nBlue: LOW fits better; orange: HIGH fits better."),
    theme = theme(plot.caption = element_text(hjust = 0)))
ggsave("figures/floor_modes_loglik_contrib.png", p_ll, width = 11, height = 10,
       dpi = 150)
cat("done\n")
