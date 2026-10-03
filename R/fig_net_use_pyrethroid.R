# Summaries and a figure of the pyrethroid-only net share (#26;
# R/prep_net_use_pyrethroid.R): the net-crop-weighted share by year for the
# continent and each region, the change in the net covariate at the bioassay
# pixel-years, and figures/net_use_pyrethroid_share.png (the share by year per
# region, and a 2024 map of the share at each cell).

source("R/packages.R")
source("R/functions.R")

net_type_w <- c(0.25, 0.47)
years <- 2000:2024
summary_years <- c(2000, 2005, 2010, 2015, 2020, 2024)

netcrop <- read_csv("data/raw/itn/net_use_20260929/netcrop_multitype_timeseries.csv",
                    show_col_types = FALSE) %>%
  filter(year %in% years) %>%
  mutate(country_name = printable_country_name(
    if_else(country == "Côte d'Ivoire", "Côte d’Ivoire", country))) %>%
  left_join(country_region_lookup(), by = "country_name")
stopifnot(!anyNA(netcrop$region))

# net-crop-weighted share, the sum of w cITN + LLIN over all nets
crop_share <- function(data, ...) {
  data %>%
    group_by(..., year) %>%
    summarise(citn = sum(cITN), llin = sum(LLIN), total = sum(`Total Nets`),
              .groups = "drop") %>%
    cross_join(tibble(w = net_type_w)) %>%
    mutate(share = (w * citn + llin) / total)
}
shares <- bind_rows(
  crop_share(netcrop) %>% mutate(region = "Africa"),
  crop_share(netcrop, region)
)

cat("pyrethroid-only share of the net crop (w cITN + LLIN over all nets):\n")
shares %>%
  filter(year %in% summary_years) %>%
  select(w, region, year, share) %>%
  pivot_wider(names_from = year, values_from = share) %>%
  arrange(w, region) %>%
  as.data.frame() %>%
  print(digits = 2)


# the net covariate at bioassay pixel-years -----------------------------------

mask <- rast("data/clean/raster_mask.tif")
bioassays <- readRDS("data/clean/all_gambiae_complex_data.RDS") %>%
  mutate(cell = terra::cellFromXY(mask, cbind(longitude, latitude)),
         year = pmin(pmax(year_start, min(years)), max(years))) %>%
  filter(!is.na(cell), year_start >= 1995) %>%
  distinct(cell, year)
cube_at <- function(cube) {
  x <- as.matrix(terra::extract(cube, bioassays$cell))
  x[cbind(seq_len(nrow(x)), bioassays$year - min(years) + 1)]
}
old <- cube_at(rast("data/clean/net_use_cube.tif"))
ratios <- bind_rows(lapply(net_type_w, function(w) {
  new <- cube_at(rast(sprintf("data/clean/net_use_pyrethroid_cube_w%.2f.tif",
                              w)))
  tibble(w = w, year = bioassays$year, old = old, new = new)
}))
cat("\nbioassay pixel-years:", nrow(bioassays), "(years before 2000 as 2000);",
    sum(is.na(old)), "with no net use\n")
cat("pyrethroid-only over all net use at bioassay pixel-years",
    "(mean of the ratio where use > 0.01; ratio of means):\n")
ratios %>%
  filter(!is.na(old)) %>%
  group_by(w, year) %>%
  summarise(n = n(),
            mean_ratio = mean((new / old)[old > 0.01]),
            ratio_of_means = mean(new) / mean(old),
            mean_old = mean(old),
            mean_new = mean(new),
            .groups = "drop") %>%
  filter(year %in% summary_years) %>%
  as.data.frame() %>%
  print(digits = 3)


# figure ----------------------------------------------------------------------

region_colours <- c(Africa = "#222222",
                    "Eastern Africa" = "#2a78d6",
                    "Middle Africa" = "#eb6834",
                    "Northern Africa" = "#e87ba4",
                    "Southern Africa" = "#1baf7a",
                    "Western Africa" = "#eda100")
shares_main <- filter(shares, w == net_type_w[1])
p_share <- shares_main %>%
  ggplot(aes(year, share, colour = region)) +
  geom_line(linewidth = 0.7) +
  ggrepel::geom_text_repel(data = filter(shares_main, year == max(years)),
                           aes(label = region), hjust = 0, nudge_x = 0.6,
                           direction = "y", size = 3, segment.colour = NA,
                           show.legend = FALSE) +
  scale_colour_manual(values = region_colours) +
  scale_x_continuous(expand = expansion(mult = c(0.02, 0.35))) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = NULL, y = "Pyrethroid-only share of nets",
       colour = NULL,
       title = sprintf("Net-crop share, w = %.2f", net_type_w[1]),
       caption = "Northern Africa: Sudan only") +
  theme_minimal() +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())

share_2024 <- rast(sprintf("data/clean/net_pyrethroid_share_cube_w%.2f.tif",
                           net_type_w[1]))[["share_2024"]]
p_map <- ggplot() +
  geom_spatraster(data = share_2024) +
  scale_fill_distiller(palette = "Blues", direction = 1, limits = c(0, 1),
                       na.value = "transparent", name = "Share") +
  labs(title = "2024") +
  theme_void()

dir.create("figures", showWarnings = FALSE)
ggsave("figures/net_use_pyrethroid_share.png", p_share + p_map,
       width = 11, height = 5, bg = "white")
