## ---- load_libraries_20
library(INLA)
library(readxl)
library(glmmTMB)
library(brms)
library(emmeans)
library(tidyverse)
## ----end

source("../R/helper_functions.R")

## ---- read_processed_data
data <- readRDS(data, file = "../data/processed/data.RData")
## ----end

## ---- define_newdata
data <- data |>
  mutate(newdata = list({
    newdata <- data.frame(Decade = c("1990s", "2000s", "2010s", "2020s"))
    newdata
  }))
## ----end

## ---- compile_inla_input_data
data <- data |>
  mutate(data = map(
    .x = data,
    .f = ~ .x |>
      mutate(Obs = 1:n()) |> 
      filter(!is.na(value)) |>
      droplevels()
  )) |> 
  mutate(data = map(.x = data,
    .f =  ~ {
      s <- .x |>
        group_by(Site) |>
        summarise(Sum = sum(value)) |>
        filter(Sum == 0) |>
        pull(Site) |>
        unique()
      .x |>
        filter(!Site %in% s) |>
        droplevels()
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


inla.setOption(inla.timeout = 30)

## ---- fit_models
data <- data |>
  mutate(model = map2(
    .x = data.1,
    .y = Species,
    .f = ~ {
      cat(paste0(.y, "================\n\n"))
      nm <- fit_models_inla(.x,
        form = value ~ Decade * scale(log(HC + 0.01)) +
          f(AIMS_REEF_NAME, model = "iid") +
          f(Site, model = "iid") +
          f(Transect, model = "iid"),
        Sp = .y
      )
      nm1 <- fit_models_glmmTMB(.x,
        form = value ~ Decade * scale(log(HC + 0.01)) +
          (1 | AIMS_REEF_NAME) +
          (1 | Site) +
          (1 | Transect),
        Sp = .y
      )
      c(nm, nm1)
    }
  ))
saveRDS(data, file = "../data/modelled/data.RData")
## ----end





## data_1 <- data |>
##   slice(-1:-11) |> 
##   mutate(model = map2(
##     .x = data.1,
##     .y = Species,
##     .f = ~ {
##       cat(paste0(.y, "================\n\n"))
##       nm1 <- fit_models_glmmTMB(.x,
##         form = value ~ Decade * scale(log(HC + 0.01)) +
##           (1 | AIMS_REEF_NAME) +
##           (1 | Site) +
##           (1 | Transect),
##         Sp = .y
##       )
##       nm1
##     }
##   ))


## data1 <- data |>
##   mutate(model = map(
##     .x = Species,
##     .f = ~ {
##       fams <- c(
##         "zeroinflatednbinomial0",
##         "zeroinflatednbinomial1",
##         "zeroinflatednbinomial2",
##         "zeroinflatedpoisson0",
##         "zeroinflatedpoisson1",
##         "nbinomial",
##         "nbinomial2",
##         "poisson"
##       )
##       nm <- paste0("../data/modelled/mod_inla_", .x, "__", fams, ".RData")
##       fams <- c(
##         "nbinom1",
##         "nbinom2",
##         "poisson",
##         "nbinom1",
##         "nbinom2",
##         "poisson"
##       )
##       zi_nm <- c("zeroinflated", "zeroinflated", "zeroinflated", "", "", "")
##       nm1 <- paste0("../data/modelled/mod_glmmTMB_", .x, "__", zi_nm, fams, ".RData")
##       c(nm, nm1)
##     }
##   ))

## data <- data1


























## data_m <- data_m |>
##   mutate(pred_inla = pmap(.l =  list(mod = model, nd = newdata, i = inewdata),
##     .f = ~ {
##       mod_inla <- readRDS(gsub("glmmTMB", "inla", ..1))
##       predict_inla(mod_inla, ..2, ..3)
##     }
##   )) |> 
##   mutate(pred_glmmTMB = pmap(.l =  list(mod = model),
##     .f = ~ {
##       mod_glmmTMB <- readRDS(gsub("inla", "glmmTMB", ..1))
##       predict_glmmTMB(mod_glmmTMB)
##     }
##   ))


## mod_inla <- readRDS(data_m[1, "model"][[1]][[1]])
## nd <- data[1, "newdata"][[1]][[1]]
## i <- data[1, "inewdata"][[1]][[1]]
## predict_inla(mod_inla, nd, i)

## mod_inla <- fit_model_inla(data.1,
##   form = value ~ Decade * scale(log(HC + 0.01)) +
##     f(AIMS_REEF_NAME, model = "iid") +
##     f(Site, model = "iid") +
##     f(Transect, model = "iid"),
##   family = "zeroinflatednbinomial0"
## )
## mod_inla_p <- predict_inla(mod_inla, newdata, inewdata)
## mod_inla1 <- mod_inla
## -sum(log(mod_inla1$cpo$cpo[!mod_inla1$cpo$fail>0 & !is.na(mod_inla1$cpo$fail)]))
## -sum(log(mod_inla$cpo$cpo[!mod_inla$cpo$fail>0 & !is.na(mod_inla$cpo$fail)]))

## mod_glmmTMB <- fit_model_glmmTMB(data.2,
##   form = value ~ Decade * scale(log(HC + 0.01)) +
##       (1 | AIMS_REEF_NAME) +
##       (1 | Site) +
##     (1 | Transect),
##   family = "nbinom2",
##   ziformula = ~ Decade +
##       (1 | AIMS_REEF_NAME) +
##       (1 | Site) +
##       (1 | Transect)
## )
## mod_glmmTMB_p <- predict_glmmTMB(mod_glmmTMB)



## data.1 <- data[1, "data.1"][[1]][[1]]
## newdata <- data[1, "newdata"][[1]][[1]]
## inewdata <- data[1, "inewdata"][[1]][[1]]
## Species <- data[1, "Species"][[1]][[1]]
## data.2 <- data[1, "data"][[1]][[1]]



## ############################################
## data.1 <- data[12, "data.1"][[1]][[1]]
## newdata <- data[12, "newdata"][[1]][[1]]
## inewdata <- data[12, "inewdata"][[1]][[1]]
## Species <- data[12, "Species"][[1]][[1]]
## data.2 <- data[12, "data"][[1]][[1]]

## data.2 |>
##   filter(!is.na(value)) |> 
##   group_by(Decade) |>
##   summarise(Sum = sum(value), nZeros = sum(value == 0), N = n()) |>
##   mutate(zero_fraction = nZeros / N)


## mod_inla <- fit_model_inla(data.1, form = value~Decade, family = "poisson")
## mod_inla_p <- predict_inla(mod_inla, newdata, inewdata)

## mod_glmmTMB <- fit_model_glmmTMB(data.2, form = value~Decade, family = "poisson", ziformula = NULL)
## mod_inla_p <- predict_glmmTMB(mod_glmmTMB)

## ############
## mod_inla <- fit_model_inla(data.1,
##   form = value ~ Decade,
##   family = "nbinomial2"
## )
## mod_inla_p <- predict_inla(mod_inla, newdata, inewdata)

## mod_glmmTMB <- fit_model_glmmTMB(data.2,
##   form = value ~ Decade,
##   family = "nbinom2",
##   ziformula = NULL
## )
## mod_glmmTMB_p <- predict_glmmTMB(mod_glmmTMB)

## ############
## mod_inla <- fit_model_inla(data.1,
##   form = value ~ Decade + 
##     f(AIMS_REEF_NAME, model = "iid") +
##     f(Site, model = "iid") +
##     f(Transect, model = "iid"),
##   family = "poisson"
## )
## mod_inla_p <- predict_inla(mod_inla, newdata, inewdata)

## mod_glmmTMB <- fit_model_glmmTMB(data.2,
##   form = value ~ Decade + 
##       (1 | AIMS_REEF_NAME) +
##       (1 | Site) +
##       (1 | Transect),
##   family = "poisson",
##   ziformula = NULL
## )
## mod_glmmTMB_p <- predict_glmmTMB(mod_glmmTMB)


## ggplot() +
##   geom_pointrange(data = mod_glmmTMB_p, aes(y = rate, x = Decade, ymin = asymp.LCL, ymax = asymp.UCL, colour = "glmmTMB")) +
##   geom_pointrange(
##     data = mod_inla_p, aes(y = mean, x = Decade, ymin = `0.025quant`, ymax = `0.975quant`, color = "INLA"),
##     position = position_nudge(x = -0.4)
##   )

## ############
## mod_inla <- fit_model_inla(data.1,
##   form = value ~ Decade + 
##     f(AIMS_REEF_NAME, model = "iid") +
##     f(Site, model = "iid") +
##     f(Transect, model = "iid"),
##   family = "nbinomial"
## )
## mod_inla_p <- predict_inla(mod_inla, newdata, inewdata)

## mod_glmmTMB <- fit_model_glmmTMB(data.2,
##   form = value ~ Decade + 
##       (1 | AIMS_REEF_NAME) +
##       (1 | Site) +
##       (1 | Transect),
##   family = "nbinom2",
##   ziformula = NULL
## )
## mod_glmmTMB_p <- predict_glmmTMB(mod_glmmTMB)


## ggplot() +
##   geom_pointrange(data = mod_glmmTMB_p, aes(y = response, x = Decade, ymin = asymp.LCL, ymax = asymp.UCL, colour = "glmmTMB")) +
##   geom_pointrange(
##     data = mod_inla_p, aes(y = mean, x = Decade, ymin = `0.025quant`, ymax = `0.975quant`, color = "INLA"),
##     position = position_nudge(x = -0.4)
##   )

## ############
## mod_inla <- fit_model_inla(data.1,
##   form = value ~ Decade + 
##     f(AIMS_REEF_NAME, model = "iid") +
##     f(Site, model = "iid") +
##     f(Transect, model = "iid"),
##   family = "zeroinflatedpoisson2"
## )
## mod_inla_p <- predict_inla(mod_inla, newdata, inewdata)
## mod_inla1 <- mod_inla
## -sum(log(mod_inla1$cpo$cpo[!mod_inla1$cpo$fail>0 & !is.na(mod_inla1$cpo$fail)]))
## -sum(log(mod_inla$cpo$cpo[!mod_inla$cpo$fail>0 & !is.na(mod_inla$cpo$fail)]))

## mod_glmmTMB <- fit_model_glmmTMB(data.2,
##   form = value ~ Decade + 
##       (1 | AIMS_REEF_NAME) +
##       (1 | Site) +
##       (1 | Transect),
##   family = "poisson",
##   ziformula = ~ Decade +
##       (1 | AIMS_REEF_NAME) +
##       (1 | Site) +
##       (1 | Transect)
## )
## mod_glmmTMB_p <- predict_glmmTMB(mod_glmmTMB)


## ggplot() +
##   geom_pointrange(data = mod_glmmTMB_p, aes(y = rate, x = Decade, ymin = asymp.LCL, ymax = asymp.UCL, colour = "glmmTMB")) +
##   geom_pointrange(
##     data = mod_inla_p, aes(y = mean, x = Decade, ymin = `0.025quant`, ymax = `0.975quant`, color = "INLA"),
##     position = position_nudge(x = -0.4)
##   )


## ############
## mod_inla <- fit_model_inla(data.1,
##   form = value ~ Decade * scale(log(HC + 0.01)) +
##     f(AIMS_REEF_NAME, model = "iid") +
##     f(Site, model = "iid") +
##     f(Transect, model = "iid"),
##   family = "zeroinflatedpoisson2"
## )
## mod_inla_p <- predict_inla(mod_inla, newdata, inewdata)
## mod_inla1 <- mod_inla
## -sum(log(mod_inla1$cpo$cpo[!mod_inla1$cpo$fail>0 & !is.na(mod_inla1$cpo$fail)]))
## -sum(log(mod_inla$cpo$cpo[!mod_inla$cpo$fail>0 & !is.na(mod_inla$cpo$fail)]))

## mod_glmmTMB <- fit_model_glmmTMB(data.2,
##   form = value ~ Decade * scale(log(HC + 0.01)) +
##       (1 | AIMS_REEF_NAME) +
##       (1 | Site) +
##       (1 | Transect),
##   family = "poisson",
##   ziformula = ~ Decade +
##       (1 | AIMS_REEF_NAME) +
##       (1 | Site) +
##       (1 | Transect)
## )
## mod_glmmTMB_p <- predict_glmmTMB(mod_glmmTMB)


## ggplot() +
##   geom_pointrange(data = mod_glmmTMB_p, aes(y = rate, x = Decade, ymin = asymp.LCL, ymax = asymp.UCL, colour = "glmmTMB")) +
##   geom_pointrange(
##     data = mod_inla_p, aes(y = mean, x = Decade, ymin = `0.025quant`, ymax = `0.975quant`, color = "INLA"),
##     position = position_nudge(x = -0.4)
##   )

## ############
## mod_inla <- fit_model_inla(data.1,
##   form = value ~ Decade * scale(log(HC + 0.01)) +
##     f(AIMS_REEF_NAME, model = "iid") +
##     f(Site, model = "iid") +
##     f(Transect, model = "iid"),
##   family = "zeroinflatednbinomial0"
## )
## mod_inla_p <- predict_inla(mod_inla, newdata, inewdata)
## mod_inla1 <- mod_inla
## -sum(log(mod_inla1$cpo$cpo[!mod_inla1$cpo$fail>0 & !is.na(mod_inla1$cpo$fail)]))
## -sum(log(mod_inla$cpo$cpo[!mod_inla$cpo$fail>0 & !is.na(mod_inla$cpo$fail)]))

## mod_glmmTMB <- fit_model_glmmTMB(data.2,
##   form = value ~ Decade * scale(log(HC + 0.01)) +
##       (1 | AIMS_REEF_NAME) +
##       (1 | Site) +
##     (1 | Transect),
##   family = "nbinom2",
##   ziformula = ~ Decade +
##       (1 | AIMS_REEF_NAME) +
##       (1 | Site) +
##       (1 | Transect)
## )
## mod_glmmTMB_p <- predict_glmmTMB(mod_glmmTMB)


## ggplot() +
##   geom_pointrange(data = mod_glmmTMB_p, aes(y = response, x = Decade, ymin = asymp.LCL, ymax = asymp.UCL, colour = "glmmTMB")) +
##   geom_pointrange(
##     data = mod_inla_p, aes(y = mean, x = Decade, ymin = `0.025quant`, ymax = `0.975quant`, color = "INLA"),
##     position = position_nudge(x = -0.4)
##   )




## fit_glmmTMB <- function(data.2) {
##   mod <- glmmTMB(
##     value ~ 1 + Decade * scale(log(HC+0.01)) + 
##       (1 | AIMS_REEF_NAME) +
##       (1 | Site) +
##       (1 | Transect),
##     ziformula = ~1,
##     data = data.2,
##     family = "nbinom2"
##   )
## }

## summary(mod)


## mod <- inla(
##   value ~ Decade * scale(log(HC + 0.01)) +
##     f(AIMS_REEF_NAME, model = "iid") +
##     f(Site, model = "iid") +
##     f(Transect, model = "iid"),
##   data = data.2,
##   family = "zeroinflatednbinomial2",
##   ## family = "zeroinflatednbinomial1",
##   control.predictor = list(link = 1, compute = TRUE),
##   control.compute = list(
##     dic = TRUE, cpo = TRUE, waic = TRUE,
##     config = TRUE)
## )
## mod <- inla(
##   value ~ Decade * scale(log(HC + 0.01)) +
##     f(AIMS_REEF_NAME, model = "iid") +
##     f(Site, model = "iid") +
##     f(Transect, model = "iid"),
##   data = data.2,
##   family = "nbinomial2",
##   ## family = "zeroinflatednbinomial1",
##   control.predictor = list(link = 1, compute = TRUE),
##   control.compute = list(
##     dic = TRUE, cpo = TRUE, waic = TRUE,
##     config = TRUE)
## )
## summary(mod)

##   mod_g <- glmmTMB(
##     value ~ 1 + Decade * scale(log(HC+0.01)) + 
##       (1 | AIMS_REEF_NAME) +
##       (1 | Site) +
##       (1 | Transect),
##     data = data.2,
##     family = "nbinom2"
##   )

## summary(mod_g)


## data <-
##   data |>
##   slice(12) |> 
##   mutate(model = map2(
##     .x = data.1,
##     .y = Species,
##     .f =  ~ {
##       cat(.y, "\n")
##       .x <- .x |>
##         filter(!is.na(value)) 
        
##       fit_model(.x, .y)
##     }
##   ))
## #6,7,11,12
## #8





## mod <- inla(
##   value ~ 1,
##   data = data.2,
##   family = "nbinomial2",
##   ## family = "zeroinflatednbinomial1",
##   control.predictor = list(link = 1, compute = TRUE),
##   control.compute = list(
##     dic = TRUE, cpo = TRUE, waic = TRUE,
##     config = TRUE)
## )
## summary(mod)

##   mod_g <- glmmTMB(
##     value ~ 1,
##     data = data.2,
##     family = "nbinom2"
##   )

## summary(mod_g)



## library(brms)


## #################################

## data.2 |>
##   filter(!is.na(value)) |>
##   group_by(Decade) |>
##   summarise(
##     Mean = log(mean(value)),
##     Median = log(median(value)),
##     SD = log(sd(value))
##   )

## priors <- prior(normal(-2.3, 2), class = "Intercept") +
##   prior(normal(0, 3), class = "b") 
## form <- bf(value ~ Decade, family = "negbinomial")
## mod_b <- brm(form,
##   data = data.2,
##   prior = priors,
##   sample_prior = "yes",
##   iter = 1000,
##   warmup = 500,
##   thin = 2,
##   chains = 3,
##   cores = 3
## )
## summary(mod_b)


## mod <- inla(
##   value ~ Decade,
##   data = data.2,
##   ## family = "nbinomial2",
##   family = "nbinomial",
##   ## family = "zeroinflatednbinomial1",
##   control.predictor = list(link = 1, compute = TRUE),
##   control.compute = list(
##     dic = TRUE, cpo = TRUE, waic = TRUE,
##     config = TRUE)
## )
## summary(mod)

## mod_g <- glmmTMB(
##   value ~ Decade,
##   data = data.2,
##   family = "nbinom2"
## )

## summary(mod_g)

## #################################

## data.2 |>
##   filter(!is.na(value)) |>
##   group_by(Decade) |>
##   summarise(
##     Mean = log(mean(value)),
##     Median = log(median(value)),
##     SD = log(sd(value))
##   )

## priors <- prior(normal(-2.3, 2), class = "Intercept") +
##   prior(normal(0, 3), class = "b") 
## form <- bf(value ~ Decade * scale(log(HC + 0.01)), family = "negbinomial")
## mod_b <- brm(form,
##   data = data.2,
##   prior = priors,
##   sample_prior = "yes",
##   iter = 1000,
##   warmup = 500,
##   thin = 2,
##   chains = 3,
##   cores = 3
## )
## summary(mod_b)


## mod <- inla(
##   value ~ Decade,
##   data = data.2,
##   ## family = "nbinomial2",
##   family = "nbinomial",
##   ## family = "zeroinflatednbinomial1",
##   control.predictor = list(link = 1, compute = TRUE),
##   control.compute = list(
##     dic = TRUE, cpo = TRUE, waic = TRUE,
##     config = TRUE)
## )
## summary(mod)

## mod_g <- glmmTMB(
##   value ~ Decade,
##   data = data.2,
##   family = "nbinom2"
## )

## summary(mod_g)

## #################################

## data.2 |>
##   filter(!is.na(value)) |>
##   group_by(Decade) |>
##   summarise(
##     Mean = log(mean(value)),
##     Median = log(median(value)),
##     SD = log(sd(value))
##   )

## priors <- prior(normal(-2.3, 2), class = "Intercept") +
##   prior(normal(0, 3), class = "b") 
## form <- bf(value ~ Decade * scale(log(HC + 0.01)) +
##              (1 | AIMS_REEF_NAME),
##   family = "negbinomial")
## mod_b <- brm(form,
##   data = data.2,
##   prior = priors,
##   sample_prior = "yes",
##   iter = 1000,
##   warmup = 500,
##   thin = 2,
##   chains = 3,
##   cores = 3
## )
## summary(mod_b)


## mod <- inla(
##   value ~ Decade* scale(log(HC + 0.01)) +
##              f(AIMS_REEF_NAME, model = "iid"),
##   data = data.2,
##   ## family = "nbinomial2",
##   family = "nbinomial",
##   ## family = "zeroinflatednbinomial1",
##   control.predictor = list(link = 1, compute = TRUE),
##   control.compute = list(
##     dic = TRUE, cpo = TRUE, waic = TRUE,
##     config = TRUE)
## )
## summary(mod)

## mod_g <- glmmTMB(
##   value ~ Decade * scale(log(HC + 0.01)) +
##              (1 | AIMS_REEF_NAME),
##   data = data.2,
##   family = "nbinom2"
## )

## summary(mod_g)

## #################################

## data.2 |>
##   filter(!is.na(value)) |>
##   group_by(Decade) |>
##   summarise(
##     Mean = log(mean(value)),
##     Median = log(median(value)),
##     SD = log(sd(value))
##   )

## priors <- prior(normal(-2.3, 2), class = "Intercept") +
##   prior(normal(0, 3), class = "b") 
## form <- bf(value ~ Decade * scale(log(HC + 0.01)) +
##              (1 | AIMS_REEF_NAME) + (1|Site) + (1|Transect),
##   family = "negbinomial")
## mod_b <- brm(form,
##   data = data.2,
##   prior = priors,
##   sample_prior = "yes",
##   iter = 1000,
##   warmup = 500,
##   thin = 2,
##   chains = 3,
##   cores = 3
## )
## summary(mod_b)


## mod <- inla(
##   value ~ Decade* scale(log(HC + 0.01)) +
##     f(AIMS_REEF_NAME, model = "iid") +
##     f(Site, model = "iid") +
##     f(Transect, model = "iid"),
##   data = data.2,
##   ## family = "nbinomial2",
##   family = "nbinomial",
##   ## family = "zeroinflatednbinomial1",
##   control.predictor = list(link = 1, compute = TRUE),
##   control.compute = list(
##     dic = TRUE, cpo = TRUE, waic = TRUE,
##     config = TRUE)
## )
## summary(mod)

## mod_g <- glmmTMB(
##   value ~ Decade * scale(log(HC + 0.01)) +
##              (1 | AIMS_REEF_NAME) + (1|Site) + (1|Transect),
##   data = data.2,
##   family = "nbinom2"
## )

## summary(mod_g)

## #################################

## data.2 |>
##   filter(!is.na(value)) |>
##   group_by(Decade) |>
##   summarise(
##     Mean = log(mean(value)),
##     Median = log(median(value)),
##     SD = log(sd(value))
##   )

## priors <- prior(normal(-2.3, 2), class = "Intercept") +
##   prior(normal(0, 3), class = "b") 
## form <- bf(value ~ Decade * scale(log(HC + 0.01)) +
##              (1 | AIMS_REEF_NAME) + (1|Site) + (1|Transect),
##   zi ~ 1,
##   family = "zero_inflated_negbinomial")
## mod_b <- brm(form,
##   data = data.2,
##   prior = priors,
##   sample_prior = "yes",
##   iter = 1000,
##   warmup = 500,
##   thin = 2,
##   chains = 3,
##   cores = 3
## )
## summary(mod_b)
## saveRDS(mod_b, file = "../data/modelled/test_brms.RData")

## mod_b_p <- emmeans(mod_b, ~Decade, type = "response") |>
##   as.data.frame()



## mod <- inla(
##   value ~ Decade* scale(log(HC + 0.01)) +
##     f(AIMS_REEF_NAME, model = "iid") +
##     f(Site, model = "iid") +
##     f(Transect, model = "iid"),
##   data = data.1,
##   ## family = "nbinomial2",
##   ## family = "nbinomial",
##   ## family = "zeroinflatednbinomial0",
##   ## family = "zeroinflatedpoisson1",
##   family = "zeroinflatedpoisson0",
##   control.predictor = list(link = 1, compute = TRUE),
##   control.compute = list(
##     dic = TRUE, cpo = TRUE, waic = TRUE,
##     config = TRUE)
## )
## summary(mod)

## mod_p <- mod$summary.fitted.values[(nrow(data.2) + 1):nrow(data.1), ] |>
##   bind_cols(newdata)

## mod_g <- glmmTMB(
##   value ~ Decade * scale(log(HC + 0.01)) +
##              (1 | AIMS_REEF_NAME) + (1|Site) + (1|Transect),
##   ziformula = ~1,
##   data = data.2,
##   family = "nbinom2"
## )

## summary(mod_g)
## mod_g_p <- emmeans(mod_g, ~Decade, type = "response") |>
##   summary(infer = TRUE) |>
##   as.data.frame()

## raw_p <- data.2 |>
##   group_by(Decade) |>
##   reframe(mean_cl_normal(value))

##   reframe(t(Hmisc::smean.cl.normal(value)))


## ggplot() +
##   geom_pointrange(data = mod_g_p, aes(y = response, x = Decade, ymin = asymp.LCL, ymax = asymp.UCL, colour = "glmmTMB")) +
##   geom_pointrange(
##     data = mod_b_p, aes(y = prob, x = Decade, ymin = lower.HPD, ymax = upper.HPD, colour = "brms"),
##     position = position_nudge(x = 0.2)
##   ) +
##   ## geom_pointrange(
##   ##   data = raw_p, aes(y = y, x = Decade, ymin = ymin, ymax = ymax, color = "Raw"),
##   ##   position = position_nudge(x = -0.2)
##   ## ) +
##   geom_pointrange(
##     data = mod_p, aes(y = mean, x = Decade, ymin = `0.025quant`, ymax = `0.975quant`, color = "INLA"),
##     position = position_nudge(x = -0.4)
##   )
