## ---- load_libraries_30
library(INLA)
library(readxl)
library(glmmTMB)
library(brms)
library(emmeans)
library(tidyverse)
library(tidybayes)
## ----end

source("../R/helper_functions.R")

## ---- load_data_40
data_pred2 <- readRDS("../data/modelled/data_pred2.RData")
benthic_pred2 <- readRDS("../data/modelled/benthic_pred2.RData")
## ----end

data_pred_old <- data_pred2

data_pred2 <- data_pred_old

## ---- load_model_lookup
model_lookup <- read_csv("../data/primary/All species averages.csv")
model_lookup <- model_lookup |>
  mutate(`Best model` = ifelse(Species == "Lutjanus_adetii", NA, `Best model`)) |> 
  mutate(`Best model` = ifelse(Species == "Chromis_agilis", "Zero-inflated poisson 1", `Best model`)) |> 
  mutate(`Best model` = ifelse(Species == "Chromis_viridis", "Zero-inflated poisson 1", `Best model`)) |> 
  mutate(`Best model` = ifelse(Species == "Chromis_amboinensis", "Zero-inflated poisson 1", `Best model`)) |> 
  mutate(`Best model` = ifelse(Species == "Chromis_amboinensis", "Zero-inflated poisson 1", `Best model`)) |> 
  mutate(`Best model` = ifelse(Species == "Lutjanus_kasmira", "Zero-inflated poisson 1", `Best model`)) |> 
  mutate(`Best model` = ifelse(Species == "Chaetodon_speculum", "Zero-inflated poisson 1", `Best model`)) |> 
  mutate(`Best model` = ifelse(Species == "Chaetodon_bennetti", "Zero-inflated poisson 1", `Best model`)) |> 
  mutate(`Best model` = ifelse(Species == "Plectropomus_areolatus", "Zero-inflated poisson 1", `Best model`)) 
## ----end



## ---- apply_model_lookup
data_pred2 <- data_pred2 |>
  left_join(model_lookup |> dplyr::select(Trophic, Family, Species, Best_model = `Best model`, Exclude),
    by = c("Species")
  ) |>
  mutate(Best_model = ifelse(is.na(Best_model), "Zero-inflated Negative Binomial 0", Best_model)) |>
    filter(Exclude != "TRUE") |>
  droplevels()
## ----end


## Arranged by Family level
## ---- arrange_species_by_family_level
data_pred2_fam <- data_pred2 |>
  dplyr::select(Family, Species, contr, Best_model) |>
  arrange(Family, Species) |>
  mutate(
    Family = factor(Family, levels = unique(Family)),
    Species = factor(Species, levels = unique(Species))
  )
## ----end

data_pred2_fam
data_pred2_fam[1, "contr"][[1]][[1]]
data_pred2_fam[1, "contr"][[1]][[1]][1, "pred"][[1]][[1]]

## ---- compile_family_arranged_plot_data
data_pred3_fam <-
  data_pred2_fam |>
  mutate(plot_data = pmap(
    .l = list(contr, Best_model),
    .f = ~ {
      Best_model <- ..2
      which_model <- case_when(
        Best_model == "Zero-inflated poisson 1" ~ "zeroinflatedpoisson1",
        Best_model == "Zero-inflated poisson 0" ~ "zeroinflatedpoisson0",
        Best_model == "Zero-inflated negative binomial 1" ~ "zeroinflatednbinomial1",
        Best_model == "Negative binomial" ~ "nbinomial",
        Best_model == "Negative binomial 1" ~ "nbinomial2",
        Best_model == "Poisson" ~ "poisson",
        .default = "zeroinflatednbinomial0"
      )
      df <- ..1 |>
        filter(name == which_model, mod == "inla") |>
        unnest(pred) |>
        dplyr::select(-mod, -name, -fam, -type, -zi, -nm, -data)
      df
    }
  )) |>
  unnest(plot_data)
## ----end

## ---- compile_hcc_plot_data
benthic_pred3 <- benthic_pred2[1, "contr"][[1]][[1]] |>
  filter(mod == "inla", name == "betabinomial") |>
  pull(pred) |>
  _[[1]] |>
  filter(variable == "Ratio")
## ----end



## ---- generate_family_arranged_plots
p1 <- 
  data_pred3_fam |>
  filter(variable == "Ratio") |>
  filter(Species != "Plectropomus_areolatus") |>
  droplevels() |> 
  mutate(Effect = case_when(
    lower > 1 ~ "Positive",
    upper < 1 ~ "Negative",
    .default = "None"
  )) |> 
  group_by(Family, Decade) |>
  arrange(mean) |>
  mutate(Species = factor(Species, levels = unique(Species))) |>
  ## filter(Family == "Pomacentridae") |>
  ## filter(Family == "Serranidae") |>
  droplevels() |>
  ggplot(aes(y = Species, x = mean)) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  geom_vline(data = benthic_pred3, aes(xintercept = mean), colour = "blue") +
  geom_ribbon(data = benthic_pred3,
    aes(x = mean, xmin = lower, xmax = upper, ymin = -Inf, ymax = Inf), colour = "blue", fill = "blue",
    inherit.aes = FALSE) +
  geom_pointrange(aes(xmin = lower, xmax = upper, colour = Effect)) +
  scale_y_discrete("", expand = c(0, 1.5), labels = function(x) str_replace(x, "_", " ")) +
  scale_x_continuous("Percentage change in fish abundance (100 x focal decade /1990's)",
    trans = scales::log2_trans(),
    label = function(x) 100 * (x -1)) +
  scale_colour_manual(
    breaks = c("Negative", "None", "Positive"),
    values = c("red", "grey", "darkgreen")
  ) +
  facet_grid(Family ~ Decade,
    space = "free",
    scales = "free") +
  theme_bw() +
  theme(axis.text.y = element_text(size = rel(0.75)), strip.text.y.right = element_text(angle = 0))

ggsave(
  file = "../data/modelled/family_level_plot.png",
  plot = p1,
  width = 15, height = 15
)
## ----end



## Arranged by Trophic level
## ---- arrange_species_by_trophic_level
data_pred2_trop <- data_pred2 |>
  dplyr::select(Trophic, Species, contr, Best_model) |>
  arrange(Trophic, Species) |>
  mutate(
    Tropic = factor(Trophic, levels = unique(Trophic)),
    Species = factor(Species, levels = unique(Species))
  )
## ----end

data_pred2_trop
data_pred2_trop[1, "contr"][[1]][[1]]
data_pred2_trop[1, "contr"][[1]][[1]][1, "pred"][[1]][[1]]

## ---- compile_trophic_arranged_plot_data
data_pred3_trop <-
  data_pred2_trop |>
  mutate(plot_data = pmap(
    .l = list(contr, Best_model),
    .f = ~ {
      Best_model <- ..2
      which_model <- case_when(
        Best_model == "Zero-inflated poisson 1" ~ "zeroinflatedpoisson1",
        Best_model == "Zero-inflated poisson 0" ~ "zeroinflatedpoisson0",
        Best_model == "Zero-inflated negative binomial 1" ~ "zeroinflatednbinomial1",
        Best_model == "Negative binomial" ~ "nbinomial",
        Best_model == "Negative binomial 1" ~ "nbinomial2",
        Best_model == "Poisson" ~ "poisson",
        .default = "zeroinflatednbinomial0"
      )
      df <- ..1 |>
        filter(name == which_model, mod == "inla") |>
        unnest(pred) |>
        dplyr::select(-mod, -name, -fam, -type, -zi, -nm, -data)
      df
    }
  )) |>
  unnest(plot_data)
## ----end

## ---- generate_trophic_arranged_plots
p2 <-
  data_pred3_trop |>
  filter(variable == "Ratio") |>
  filter(Species != "Plectropomus_areolatus") |>
  droplevels() |> 
  mutate(Effect = case_when(
    lower > 1 ~ "Positive",
    upper < 1 ~ "Negative",
    .default = "None"
  )) |> 
  group_by(Trophic, Decade) |>
    arrange(mean) |>
    mutate(Species = factor(Species, levels = unique(Species))) |>
    ## filter(Family == "Pomacentridae") |>
    ## filter(Family == "Serranidae") |>
    droplevels() |>
    ggplot(aes(y = Species, x = mean)) +
    geom_vline(xintercept = 1, linetype = "dashed") +
    geom_pointrange(aes(xmin = lower, xmax = upper, colour = Effect)) +
    scale_y_discrete("", expand = c(0, 1.5), labels = function(x) str_replace(x, "_", " ")) +
  scale_x_continuous("Percentage change in fish abundance (100 x focal decade /1990's)",
    trans = scales::log2_trans(),
    label = function(x) 100 * (x -1)) +
  scale_colour_manual(
    breaks = c("Negative", "None", "Positive"),
    values = c("red", "grey", "darkgreen")
  ) +
  facet_grid(Trophic ~ Decade,
    space = "free",
    scales = "free") +
    theme_bw() +
    theme(axis.text.y = element_text(size = rel(0.75)), strip.text.y.right = element_text(angle = 0))

ggsave(
  file = "../data/modelled/trophic_level_plot.png",
  plot = p2,
  width = 15, height = 15
)
## ----end
