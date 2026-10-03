# FAOSTAT pesticide use (domain RP) and pesticide trade (domain RT) for the
# African countries and FAO's African aggregates, as candidate data for the
# time trend g(t) of non-vector-control insecticide use.
#
# Bulk downloads (normalized CSVs), fetched 2026-09-29; the CSVs inside are
# dated 2026-07-29:
#   https://fenixservices.fao.org/faostat/static/bulkdownloads/Inputs_Pesticides_Use_E_All_Data_(Normalized).zip
#   https://fenixservices.fao.org/faostat/static/bulkdownloads/Inputs_Pesticides_Trade_E_All_Data_(Normalized).zip
# unzipped into data/raw/faostat/.
#
# Flags (from the *_Flags.csv files): A official; E estimated; I imputed by a
# receiving agency; X from an external organisation; L missing (data exist but
# were not collected); B time series break (trade only). The use data also
# carry a Note giving the estimation method.

source("R/packages.R")
source("R/functions.R")

faostat_dir <- "data/raw/faostat"

read_faostat <- function(file) {
  read_csv(file.path(faostat_dir, file),
           col_types = cols(.default = "c"),
           locale = locale(encoding = "UTF-8")) %>%
    mutate(
      m49 = str_remove(`Area Code (M49)`, "'"),
      year = as.integer(Year),
      value = as.numeric(Value)
    )
}

# African countries, with their model region and printable name, by M49 code
africa_m49 <- read_delim("data/raw/UNSD — Methodology.csv",
                         delim = ";",
                         col_types = cols(.default = "c")) %>%
  filter(`Region Name` == "Africa") %>%
  select(country = `Country or Area`, m49 = `M49 Code`) %>%
  mutate(country_name = printable_country_name(country)) %>%
  left_join(country_region_lookup(), by = join_by(country_name))

# FAOSTAT's pre-2011 Sudan (736) is not a UNSD country; give it Sudan's region
africa_m49 <- bind_rows(
  africa_m49,
  tibble(country = "Sudan (former)", m49 = "736",
         country_name = "Sudan (former)", region = "Northern Africa")
)

# the modelled countries: those in the bioassay data or the prediction mask
model_countries <- union(
  readRDS("data/clean/all_gambiae_complex_data.RDS")$country_name,
  readRDS("data/clean/country_borders.RDS")$country_name
)

# FAO's aggregates for Africa and its five subregions
africa_aggregates <- c("Africa", "Eastern Africa", "Middle Africa",
                       "Northern Africa", "Southern Africa", "Western Africa")

tidy_faostat <- function(raw) {
  countries <- raw %>%
    inner_join(africa_m49, by = join_by(m49)) %>%
    mutate(aggregate = FALSE,
           modelled = country_name %in% model_countries)
  aggregates <- raw %>%
    filter(Area %in% africa_aggregates) %>%
    mutate(country_name = Area,
           region = Area,
           aggregate = TRUE,
           modelled = FALSE)
  bind_rows(countries, aggregates) %>%
    select(country_name, region, aggregate, modelled,
           item_code = `Item Code`, item = Item,
           element = Element, unit = Unit,
           year, value, flag = Flag,
           any_of(c(note = "Note")))
}

pesticide_use <- read_faostat(
  "Inputs_Pesticides_Use_E_All_Data_(Normalized).csv"
) %>%
  tidy_faostat()

pesticide_trade <- read_faostat(
  "Inputs_Pesticides_Trade_E_All_Data_(Normalized).csv"
) %>%
  tidy_faostat()

saveRDS(pesticide_use, "data/clean/faostat_pesticide_use.RDS")
saveRDS(pesticide_trade, "data/clean/faostat_pesticide_trade.RDS")

missing <- setdiff(model_countries, pesticide_use$country_name)
cat("modelled countries with no FAOSTAT use rows:",
    paste(missing, collapse = ", "), "\n")
