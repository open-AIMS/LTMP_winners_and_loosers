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

## ---- load_data_30
data <- readRDS(file = "../data/modelled/data.RData")
## ----end

## ---- model_lookup
model_lookup <- tribble(
  ~mod,  ~name, ~fam, ~type, ~zi,
  "inla",  "zeroinflatednbinomial0", "Negative Binomial", "0", TRUE,
  "inla",   "zeroinflatednbinomial1", "Negative Binomial", "1", TRUE,
  "inla",   "zeroinflatednbinomial2", "Negative Binomial", "2", TRUE,
  "inla",   "zeroinflatedpoisson0", "Poisson", "0", TRUE,
  "inla",   "zeroinflatedpoisson1", "Poisson", "1",TRUE,
  "inla",   "nbinomial","Negative Binomial", "", FALSE,
  "inla",   "nbinomial2","Negative Binomial", "2",FALSE,
  "inla",   "poisson", "Poisson", "",FALSE,
  "glmmTMB",    "nbinom1", "Negative Binomial", "1", TRUE,
  "glmmTMB",    "nbinom2", "Negative Binomial", "2",TRUE,
  "glmmTMB",    "poisson","Poisson", "", TRUE,
  "glmmTMB",    "nbinom1", "Negative Binomial", "1", FALSE,
  "glmmTMB",    "nbinom2", "Negative Binomial",  "2", FALSE,
  "glmmTMB",    "poisson","Poisson", "", FALSE,
  )
## ----end



Species <- data[1,'Species'][[1]][[1]]
newdata <- data[1,'newdata'][[1]][[1]]
inewdata <- data[1,'inewdata'][[1]][[1]]
nms <- model <- data[1,'model'][[1]][[1]]

.x <- nm <- nms[3]
## ---- posteriors 
data_pred <- data |>
  slice(1) |> 
  mutate(pred = pmap(
    .l = list(Species = Species, newdata = newdata, inewdata = inewdata, model = model),
    .f = ~ {
      nms <- ..4
      newdata <- ..2
      inewdata <- ..3
      models <- model_lookup |>
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

saveRDS(data_pred, file = "../data/modelled/data_pred.RData")
## ----end


cellmeans <- readRDS(cellmeans_nm)
contr <- readRDS(contrasts_nm)


## newdata <- data[1, "newdata"][[1]][[1]]
## inewdata <- data[1, "inewdata"][[1]][[1]]
## pred <- readRDS("../data/modelled/pred_inla_Acanthochromis_polyacanthus__zeroinflatednbinomial0.RData")


data_pred <- readRDS(file = "../data/modelled/data_pred.RData")
## data_pred <- readRDS(file = "../data/modelled/data_pred_test.RData")

.x <- data_pred[1, "pred"][[1]][[1]]
data_pred[1, "pred"][[1]][[1]][7, "pred"][[1]][[1]]
data_pred[1, "pred"][[1]][[1]][1, "nm"]

## ---- collate_posteriors
data_pred1 <- data_pred |>
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

saveRDS(data_pred1, file = "../data/modelled/data_pred1.RData")
## ----end



data_pred1[1, "cellmeans"][[1]][[1]]
data_pred1[1, "cellmeans"][[1]][[1]][1, "pred"][[1]][[1]]
data_pred1[1, "contr"][[1]][[1]][1, "pred"][[1]][[1]]


## ---- load_posteriors
data_pred1 <- readRDS("../data/modelled/data_pred1.RData")
## ----end

## ---- generate_plots
data_pred2 <- data_pred1 |>
  mutate(cellmeans_plot = map2(
    .x = cellmeans,
    .y = Species,
    .f = ~ {
      nm <- paste0("../data/modelled/plot_cellmeans_", .y, ".png")
      p <- .x |>
        unnest(c(pred)) |>
        ggplot(aes(y = mean, x = Decade, colour = mod, shape = type)) +
        geom_pointrange(aes(ymin = lower, ymax = upper), position = position_dodge(width = 0.3)) +
        facet_grid(fam ~ zi) +
        scale_x_discrete("Decade") +
        scale_y_continuous("Fish abundance") +
        theme_bw()
      print(nm)
      ggsave(filename = nm, plot = p, width = 10, height = 6)
      nm
    }
  )) |> 
  mutate(contr_plot = map2(
    .x = contr,
    .y = Species,
    .f = ~ {
      nm <- paste0("../data/modelled/plot_contr_", .y, ".png")
      p1 <- .x |>
        unnest(c(pred)) |>
        filter(variable == "Ratio") |> 
        ggplot(aes(x = mean, y = Decade, colour = mod, shape = type)) +
        geom_vline(xintercept = 1, linetype = "dashed") +
        geom_pointrange(aes(xmin = lower, xmax = upper), position = position_dodge(width = 0.3)) +
        facet_grid(fam ~ zi, scales = "free") +
        scale_x_continuous(
          trans = scales::log2_trans(),
          labels = function(x) sprintf("%0.1f%%", 100 * (x - 1))
        ) +
          ## breaks = c(0.5, 1/1.5, 1, 1.5, 2),
          ## labels = function(x) sprintf("%0.1f%%", 100 * (x - 1))
        ##) +
        scale_y_discrete("Decade") +
        scale_x_continuous("Percentage change in fish abundance (100 x focal decade / 1990's)") +
      theme_bw()
      p2 <- .x |>
        unnest(c(pred)) |>
        filter(variable == "Values") |>
        ggplot(aes(x = mean, y = Decade, colour = mod, shape = type)) +
        geom_vline(xintercept = 1, linetype = "dashed") +
        geom_pointrange(aes(xmin = lower, xmax = upper), position = position_dodge(width = 0.3)) +
        facet_grid(fam ~ zi, scales = "free") +
        scale_y_discrete("Decade") +
        scale_x_continuous("Absolute change in fish abundance (focal decade - 1990's)") +
      theme_bw()
      print(nm)
      p <- p1/p2
      ggsave(filename = nm, plot = p, width = 10, height = 10)
      nm
    }
  ))
saveRDS(data_pred2, file = "../data/modelled/data_pred2.RData")
## ----end

## ---- load_posteriors_2
data_pred2 <- readRDS("../data/modelled/data_pred2.RData")
## ----end

## ---- loop_posteriors
## cat("::: {.panel-tabset}\n\n")

pwalk(
  ## .l = with(slice(data_pred2, 1),
  .l = with(data_pred2,
    list(Species, cellmeans_plot, contr_plot)),
  .f = ~ {
    cat(paste0("## ", ..1, "\n"))
    cat(":::: {.panel-tabset}\n\n")
    cat("### Cellmeans\n\n")
    nm1 <- ..2
    #print(nm1)
    capt <- "Estimated figh abundances per decades for each of the candidate models.  Left and right facets represent without and with zero-inflation and top and bottom facets represent negative binomial and poisson models respectively.  Red and blue colours represent frequentist (glmmTMB) and Bayesian (INLA) models.  Point shapes represent alternative likelihood formulation types."
    cat(paste0("![", capt, "](", nm1, "){#fig-", str_replace_all(nm1, "..", "_"), "}\n\n"))
    cat("### Contrasts\n\n")
    nm1 <- ..3
    #print(nm1)
    capt <- "Effect sizes (in percent change units in the top figure and absolute change in the bottom figure) between each of the decades and 1990's for each of the candidate models.  Left and right facets represent without and with zero-inflation and top and bottom facets represent negative binomial and poisson models respectively.  Red and blue colours represent frequentist (glmmTMB) and Bayesian (INLA) models.  Point shapes represent alternative likelihood formulation types."
    cat(paste0("![",capt,"](", nm1, "){#fig-", str_replace_all(nm1, "..", "_"), "}\n\n"))
    cat("### DHARMa residuals\n\n")
    cat("::::\n\n")
  }
)
## cat(":::\n\n")
## ----end


## pwalk(
##   .l = with(slice(data_pred1, 1),
##     list(Species, cellmeans, cellmeans_plot, contr_plot)),
##   .f = ~ {
##     cat(paste0("## ", ..1, "\n"))
##     cat(":::: {.panel-tabset}\n\n")
##     cellmeans <- ..2
##     pwalk(
##       .l = with(cellmeans, list(mod, fam, type, zi, pred, nm)),
##       .f = ~ {
##         mod <- ..1
##         fam <- ..2
##         type <- ..3
##         zi <- ..4
##         pred <- ..5
##         nm <- ..6
##         cat(paste0(
##           "### ", ifelse(zi, "Zero-inflated ", ""),
##           mod, " ", fam, " ", type, " ", "\n\n"
##         ))
##         cat("::::: {.panel-tabset}\n\n")
##         cat("#### DHARMa residuals\n\n")
##         nm <- str_replace(nm, "pred_", "dharma_")
##         nm <- str_replace(nm, "RData", "pdf")
##         nm1 <- str_replace(nm, "pdf", "png")
##         print(nm1)
##         system(paste0("convert ", nm, " ", nm1))
##         cat(paste0("![](", nm1, ")\n\n"))
##         cat("#### Cellmeans\n\n")
##         cat("#### Contrasts\n\n")
##         cat(":::::\n\n")
##       }
##     )
##     cat("::::\n\n")
##   }
## )



## ## Cellmeans
## dat_sum <- data_pred1 |>
##   dplyr::select(Species, cellmeans) |>
##   unnest(cellmeans) |>
##   dplyr::select(-data) |>
##   unnest(pred)

## p1 <-
##   dat_sum |>
##   filter(Species == "Acanthochromis_polyacanthus") |> 
##   ggplot(aes(y = mean, x = Decade, colour = mod, shape = type)) +
##   geom_pointrange(aes(ymin = lower, ymax = upper), position = position_dodge(width = 0.3)) +
##   facet_grid(fam ~ zi)
## p1


## ## Contrasts
## dat_sum <- data_pred1 |>
##   dplyr::select(Species, contr) |>
##   unnest(contr) |>
##   dplyr::select(-data) |>
##   unnest(pred)

## p2 <-
##   dat_sum |>
##   filter(Species == "Acanthochromis_polyacanthus", variable == "Values") |> 
##   ggplot(aes(x = mean, y = Decade, colour = mod, shape = type)) +
##   geom_vline(xintercept = 0, linetype = "dashed")+
##   geom_pointrange(aes(xmin = lower, xmax = upper), position = position_dodge(width = 0.3)) +
##   facet_grid(fam ~ zi)
## p2

## p3 <-
##   dat_sum |>
##   filter(Species == "Acanthochromis_polyacanthus", variable == "Ratio") |>
##   ggplot(aes(x = mean, y = Decade, colour = mod, shape = type)) +
##   geom_vline(xintercept = 1, linetype = "dashed") +
##   geom_pointrange(aes(xmin = lower, xmax = upper), position = position_dodge(width = 0.3)) +
##   facet_grid(fam ~ zi) +
##   scale_x_continuous(
##     trans = scales::log2_trans(),
##     breaks = c(0.5, 1/1.5, 1, 1.5, 2),
##     labels = function(x) sprintf("%0.1f%%", 100 * (x - 1))
##   )
## p3
