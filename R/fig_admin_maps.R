# plot maps of admin units used in model

# load packages and functions
source("R/packages.R")
source("R/functions.R")

# make the colours and order match the main reigon figure
country_borders <- readRDS("data/clean/country_borders.RDS") %>%
  mutate(
    region = factor(region,
                    levels = c(
                      "Southern Africa",
                      "Western Africa",
                      "Middle Africa", 
                      "Eastern Africa",
                      "Northern Africa"
                    ))
  )

region_map <- ggplot() +
  geom_sf(
    aes(
      fill = region
    ),
    data = country_borders,
    col = "black"
  ) +
  scale_fill_brewer(
    name = "",
    palette = "Accent"
  ) +
  theme_ir_maps()

region_map

ggsave(
  "figures/region_map.png",
  bg = "white",
  width = 4.5,
  height = 3
)
