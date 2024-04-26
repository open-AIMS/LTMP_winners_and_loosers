## ---- load_data
load("../data/primary/fish.benthic.2024.RData")
load("../data/primary/240417orig.ltmp.w.RData")
lookup <- readxl::read_xlsx("../data/primary/Species by decade glmm.xlsx",
  sheet = "Species table"
)
## ----end

## ---- extract_benthic
benthic <- fish.benthic.2024 |>
  dplyr::select(P_CODE.y:AB) |>
  distinct() |>
  group_by(P_CODE.y, AIMS_REEF_NAME, cREPORT_YEAR, Decade, REEFSITE, REEFSITETRANSECT) |>
  summarise(
    HC = sum(HC, na.rm = TRUE),
    n.points = sum(n.points),
    total.points = unique(total.points)
  ) |>
  ungroup() |>
  mutate(REPORT_YEAR = as.numeric(as.character(cREPORT_YEAR))) |>
  dplyr::select(AIMS_REEF_NAME, REEFSITE, REEFSITETRANSECT, REPORT_YEAR, Decade, HC) |>
  mutate(
    ## SITE_NO = str_replace(REEFSITE, paste0(toupper(AIMS_REEF_NAME), ".*([0-9]$)"), "\\1"),
    ## TRANSECT_NO = str_replace(REEFSITETRANSECT, paste0(toupper(AIMS_REEF_NAME), ".*([0-9]$)"), "\\1")
    SITE_NO = str_replace(REEFSITE, ".*([0-9]$)", "\\1"),
    TRANSECT_NO = str_replace(REEFSITETRANSECT, ".*([0-9]$)", "\\1")
  ) |>
  filter(!is.na(AIMS_REEF_NAME)) |>
  droplevels()

rm("fish.benthic.2024")
gc()
saveRDS(benthic, file = "../data/processed/benthic.RData")
## ----end

## get a list of species
## ---- get_species
species <- orig.ltmp.w |>
  colnames() |>
  str_subset("[A-Z][a-z]+\\_[a-z]+")
saveRDS(species, file = "../data/processed/species.RData")
## ----end


## ---- join_benthic_and_fish_data
data <- orig.ltmp.w |>
  filter(REPORT_YEAR < 2024) |>
  droplevels() |> 
  full_join(benthic, by = c("AIMS_REEF_NAME", "SITE_NO", "TRANSECT_NO", "REPORT_YEAR"))
## ----end

## ---- add_nested_effects
data <- data |>
  mutate(
    Site = paste(REEF_NAME, SITE_NO),
    Transect = paste(Site, TRANSECT_NO),
    Decade = factor(Decade)
  ) |>
  filter(!is.na(SITE_NO), !is.na(HC)) |>
  droplevels() 
## ----end

## ---- nest_data
data <- data |>
  pivot_longer(cols = species, names_to = "Species") |>
  nest(.by = Species)
  ## nest_by(Species) # seems this might be getting depricated??
saveRDS(data, file = "../data/processed/data.RData")
## ----end

