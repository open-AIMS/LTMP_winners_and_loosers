## ---- load_libraries_20
library(INLA)
library(readxl)
library(glmmTMB)
library(brms)
library(emmeans)
library(tidyverse)
## ----end

source("../R/helper_functions.R")

## ---- read_processed_benthic_data
load("../data/primary/fish.benthic.2024.RData")
## ----end

## ---- process_benthic_data
coral_fams <-  fish.benthic.2024 |> filter(HC>0) |> pull(FAMILY_2021) |> unique()
benthic <- fish.benthic.2024 |>
  filter(FAMILY_2021 %in% coral_fams) |> 
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
  dplyr::select(AIMS_REEF_NAME, REEFSITE, REEFSITETRANSECT, REPORT_YEAR, Decade, HC, n.points, total.points) |>
  mutate(
    SITE_NO = str_replace(REEFSITE, ".*([0-9]$)", "\\1"),
    TRANSECT_NO = str_replace(REEFSITETRANSECT, ".*([0-9]$)", "\\1")
  ) |>
  filter(!is.na(AIMS_REEF_NAME)) |>
  droplevels()

rm("fish.benthic.2024")
gc()

## ----end



## ---- compile_inla_input_data
benthic <- benthic |>
  nest() |> 
  mutate(data = map(
    .x = data,
    .f = ~ .x |>
      mutate(
        Site = REEFSITE,
        Transect = REEFSITETRANSECT
      ) |> 
      mutate(Obs = 1:n()) |> 
      droplevels()
  )) |> 
  mutate(newdata = list({
    newdata <- data.frame(Decade = c("1990s", "2000s", "2010s", "2020s"), total.points = 1)
    newdata
  })) |> 
  mutate(data.1 = map2(
    .x = data, .y = newdata,
    .f = ~ .x |> bind_rows(.y)
  )) |>
  mutate(inewdata = map2(
    .x = data,
    .y = data.1,
    .f = ~ (nrow(.x) + 1):nrow(.y)
  ))
## ----end


## ---- fit_models_hcc

##for binomial inla models, need to pass the formula into the wrapper
## functions as a string, otherwise the environment of the formula
## messes it all up!

benthic <- benthic |>
  mutate(model = map(
    .x = data.1,
    .f = ~ {
      nm <- fit_models_inla_hcc(data.1 = .x,
        form = 'n.points ~ Decade +
          f(AIMS_REEF_NAME, model = "iid") +
          f(Site, model = "iid") +
          f(Transect, model = "iid")'
      )
      nm1 <- fit_models_glmmTMB_hcc(data.1 = .x,
        form = 'cbind(n.points, total.points - n.points) ~ Decade  +
          (1 | AIMS_REEF_NAME) +
          (1 | Site) +
          (1 | Transect)'
        )
      c(nm, nm1)
    }
  ))
saveRDS(data, file = "../data/modelled/benthic.RData")
## ----end

## ---- model_lookup_hcc
benthic_model_lookup <- tribble(
  ~mod,  ~name, ~fam, ~type, ~zi,
  "inla",  "binomial", "Binomial", "0", FALSE,
  "inla",  "betabinomial", "Beta-Binomial", "0", FALSE,
  "glmmTMB",    "binomial", "Binomial", "0", FALSE,
  "glmmTMB",    "betabinomial", "Beta-Binomial", "0", FALSE,
  )
## ----end

## ---- posteriors_hcc
benthic_pred <- benthic |>
  mutate(pred = pmap(
    .l = list(newdata = newdata, inewdata = inewdata, model = model),
    .f = ~ {
      nms <- ..3
      newdata <- ..1
      inewdata <- ..2
      models <- benthic_model_lookup |>
        mutate(nm = nms) |> 
        ## Had to replace the following line with the other two, because
        ## it was not supported by the R version on the HPC
        ## nest(.by = c(mod, name, fam, type, zi, nm)) |>
        group_by(mod, name, fam, type, zi, nm) %>%
        summarise(data = list(cur_data_all()), .groups = "drop") %>% 
        mutate(pred = map(
          .x = nm,
         .f =  ~ {
           print(.x)
           mod <- readRDS(.x)
           cellmeans_nm <- calc_cellmeans(mod, newdata, inewdata, .x)
           contrasts_nm <- calc_contrasts(mod, newdata, .x)
           dharma_nm <- calc_dharma(mod, .x, inewdata)
           data.frame(Type = c("Cellmeans", "Contrasts", "DHARMa"),
                      File = c("cellmeans_nm" = cellmeans_nm,
                               "contrasts_nm" = contrasts_nm,
                               "dharma_nm" = dharma_nm))
         }
        ))
      models
    }
  ))

saveRDS(benthic_pred, file = "../data/modelled/benthic_pred.RData")
## ----end

## ---- collate_posteriors_hcc
benthic_pred1 <- benthic_pred |>
  ## slice(1) |>
  mutate(cellmeans = map(
    .x = pred,
    .f = ~ {
      cellmeans_nm <- str_replace(.x$nm, "mod_", "pred_")
      .x |>
        mutate(nm = str_replace(nm, "mod_", "pred_")) |> 
        mutate(pred = map(
          .x = nm,
          .f = ~ {
            nm <- .x
            cellmeans <- readRDS(.x)
            tibble(cellmeans)
          }
        ))
    }
  )) |> 
  mutate(contr = map(
    .x = pred,
    .f = ~ {
      contr_nm <- str_replace(.x$nm, "mod_", "pred_")
      .x |>
        mutate(nm = str_replace(nm, "mod_", "contr_")) |> 
        mutate(pred = map(
          .x = nm,
          .f = ~ {
            nm <- .x
            contr <- readRDS(.x)
            tibble(contr)
          }
        ))
    }
  )) 

saveRDS(benthic_pred1, file = "../data/modelled/benthic_pred1.RData")
## ----end

## ---- generate_plots
benthic_pred2 <- benthic_pred1 |>
  mutate(cellmeans_plot = map(
    .x = cellmeans,
    .f = ~ {
      nm <- paste0("../data/modelled/plot_cellmeans_hcc.png")
      p <- .x |>
        unnest(c(pred)) |>
        ggplot(aes(y = mean, x = Decade, colour = mod)) +
        geom_pointrange(aes(ymin = lower, ymax = upper), position = position_dodge(width = 0.3)) +
        facet_grid(fam ~ zi) +
        scale_x_discrete("Decade") +
        scale_y_continuous("Hard coral cover") +
        theme_bw()
      print(nm)
      ggsave(filename = nm, plot = p, width = 6, height = 6)
      nm
    }
  )) |> 
  mutate(contr_plot = map(
    .x = contr,
    .f = ~ {
      nm <- paste0("../data/modelled/plot_contr_hcc.png")
      p1 <- .x |>
        unnest(c(pred)) |>
        filter(variable == "Ratio") |> 
        ggplot(aes(x = mean, y = Decade, colour = mod)) +
        geom_vline(xintercept = 1, linetype = "dashed") +
        geom_pointrange(aes(xmin = lower, xmax = upper), position = position_dodge(width = 0.3)) +
        facet_grid(fam ~ zi, scales = "free") +
        scale_x_continuous(
"Percentage change in percentage coral cover (100 x focal decade / 1990's)",
          trans = scales::log2_trans(),
labels = function(x) sprintf("%0.1f%%", 100 * (x - 1)),
breaks = c(1, 0.9, 0.8, 0.7, 0.6, 0.5)
        ) +
          ## breaks = c(0.5, 1/1.5, 1, 1.5, 2),
          ## labels = function(x) sprintf("%0.1f%%", 100 * (x - 1))
        ##) +
        scale_y_discrete("Decade") +
      theme_bw()
      p2 <- .x |>
        unnest(c(pred)) |>
        filter(variable == "Values") |>
        ggplot(aes(x = mean, y = Decade, colour = mod)) +
        geom_vline(xintercept = 0, linetype = "dashed") +
        geom_pointrange(aes(xmin = lower, xmax = upper), position = position_dodge(width = 0.3)) +
        facet_grid(fam ~ zi, scales = "free") +
        scale_y_discrete("Decade") +
        scale_x_continuous("Absolute change in percentage coral cover (focal decade - 1990's)",
          labels = function(x) x*100) +
      theme_bw()
      print(nm)
      p <- p1/p2
      ggsave(filename = nm, plot = p, width = 6, height = 10)
      nm
    }
  ))
saveRDS(benthic_pred2, file = "../data/modelled/benthic_pred2.RData")
## ----end

## ---- load_posteriors_2_hcc
benthic_pred2 <- readRDS("../data/modelled/benthic_pred2.RData")
## ----end

## ---- loop_posteriors_hcc
## cat("::: {.panel-tabset}\n\n")

pwalk(
  .l = with(benthic_pred2,
    list(cellmeans_plot, contr_plot)),
  .f = ~ {
    cat(":::: {.panel-tabset}\n\n")
    cat("### Cellmeans\n\n")
    nm1 <- ..1
    #print(nm1)
    capt <- "Estimated percentage hard coral cover per decades for each of the candidate models.  Top and bottom facets represent beta-binomial and binomial models respectively.  Red and blue colours represent frequentist (glmmTMB) and Bayesian (INLA) models."
    cat(paste0("![", capt, "](", nm1, "){#fig-", str_replace_all(nm1, "..", "_"), "}\n\n"))
    cat("### Contrasts\n\n")
    nm1 <- ..2
    #print(nm1)
    capt <- "Effect sizes (in percent change units in the top figure and absolute change in the bottom figure) between each of the decades and 1990's for each of the candidate models.  Top and bottom facets represent beta-binomial and binomial models respectively.  Red and blue colours represent frequentist (glmmTMB) and Bayesian (INLA) models."
    cat(paste0("![",capt,"](", nm1, "){#fig-", str_replace_all(nm1, "..", "_"), "}\n\n"))
    cat("### DHARMa residuals\n\n")
    cat("::::\n\n")
  }
)
## cat(":::\n\n")
## ----end
