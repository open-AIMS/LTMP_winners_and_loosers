
## Model types
## 1. zero inflated negative binomial 2
## 2. zero inflated negative binomial 1
## 3. zero inflated negative binomial 0
## 4. zero inflated poisson
## 5. zero inflated poisson 1
## 6. nbinomial
## 7. nbinomial 2
## 8. poisson

fit_models_inla <- function(data.1, form = value ~ 1, Sp) {
  fams <- c(
    "zeroinflatednbinomial0",
    "zeroinflatednbinomial1",
    "zeroinflatednbinomial2",
    "zeroinflatedpoisson0",
    "zeroinflatedpoisson1",
    "nbinomial",
    "nbinomial2",
    "poisson"
  )
  nms <- NULL
  for (f in fams) {
    mod_inla <- fit_model_inla(data.1, form = form, family = f)
    nm <- paste0("../data/modelled/mod_inla_", Sp, "__", f, ".RData")
    saveRDS(mod_inla, file = nm)
    nms <- c(nms, nm)
  }
  nms
}

fit_models_glmmTMB <- function(data.1, form = value ~ 1, Sp) {
  fams <- c(
    "nbinom1",
    "nbinom2",
    "poisson",
    "nbinom1",
    "nbinom2",
    "poisson"
  )
  zi <- c(
    ~ Decade + (1 | AIMS_REEF_NAME) + (1 | Site) + (1 | Transect),
    ~ Decade + (1 | AIMS_REEF_NAME) + (1 | Site) + (1 | Transect),
    ~ Decade + (1 | AIMS_REEF_NAME) + (1 | Site) + (1 | Transect),
    ~ 0,
    ~0,
    ~0
  )
  zi_nm <- c("zeroinflated", "zeroinflated", "zeroinflated", "", "", "")
  nms <- NULL
  for (i in 1:length(fams)) {
    mod_glmmTMB <- fit_model_glmmTMB(data.1, form = form, family = fams[[i]], ziformula = zi[[i]])
    nm <- paste0("../data/modelled/mod_glmmTMB_", Sp, "__", zi_nm[[i]], fams[[i]], ".RData")
    saveRDS(mod_glmmTMB, file = nm)
    nms <- c(nms, nm)
  }
  nms
}




fit_model_inla <- function(data.1, form = value ~ 1, family = "poisson") {
  e <- tryCatch({
    cat(paste0("trying ", family, " ======\n"))
    mod <- inla(
      form,
      data = data.1,
      family = family,
      control.predictor = list(link = 1, compute = TRUE),
      control.compute = list(
        dic = TRUE, cpo = TRUE, waic = TRUE,
        config = TRUE)
    )
  }, error =  function(e) e)

  if (!("inla" %in% class(e))) {
    cat("Model failed ===\n\n")
    return(NULL)
  } else {
    return(mod)
  }
}



fit_model_glmmTMB <- function(data.1, form = value ~ 1, family = "poisson", ziformula = ~ 0) {
  e <- tryCatch({
    cat(paste0("trying ", family, " ======\n"))
    mod <- glmmTMB(
      form,
      data = data.1,
      family = family,
      ziformula = ziformula,
      REML = TRUE
    )
  }, error =  function(e) e)

  if (!("glmmTMB" %in% class(e))) {
    cat("Model failed ===\n\n")
    return(NULL)
  } else {
    return(mod)
  }
}


## Cellmeans
calc_cellmeans <- function(mod, newdata, inewdata, .x) {
  pnm <- str_replace_all(.x, "mod_", "pred_")
  if (is.null(mod)) {
    pred <- NULL
  } else {
    if (str_detect(.x, "inla")) {
      pred <- predict_inla(mod, newdata, inewdata)
    } else {
      pred <- predict_glmmTMB(mod)
    }
  }
  saveRDS(pred, file = pnm)
  pnm
}

predict_inla <- function(mod, newdata, inewdata) {
  mod$summary.fitted.values[inewdata, ] |>
    bind_cols(newdata) |>
    dplyr::rename(mean = mean, lower = `0.025quant`, upper = `0.975quant`) |>
    dplyr::select(-`0.5quant`, -sd, -mode)
}

predict_glmmTMB <- function(mod) {
  e <- tryCatch(
  {
    lookup <- c("mean" = "rate", "mean" = "response", "mean" = "prob", "lower" = "asymp.LCL", "upper" = "asymp.UCL",
                "lower" = "lower.CL", "upper" = "upper.CL")
    pred <- emmeans(mod, ~Decade, type = "response") |>
      summary(infer = TRUE) |>
      as.data.frame() |>
      dplyr::rename(any_of(lookup)) |>
      dplyr::select(Decade, mean, lower, upper)
  },
    error = function(e) e
  )
  if ("simpleError" %in% class(e)) {
    return(NULL)
  } else {
    return(pred)
  }
}

## Contrasts (compare each to the first
calc_contrasts <- function(mod, newdata, .x) {
  cnm <- str_replace_all(.x, "mod_", "contr_")
  if (is.null(mod)) {
    contr <- NULL
  } else {
    if (str_detect(.x, "inla")) {
      contr <- contr_inla(mod, newdata)
    } else {
      contr <- contr_glmmTMB(mod, newdata)
    }
  }
  saveRDS(contr, file = cnm)
  cnm
}

## ---- my_ilink function
my_ilink <- function(x, link) {
  switch(link, identity = x, log = exp(x), logm1 = expp1(x), 
         log1p = expm1(x), inverse = 1/x, sqrt = x^2, `1/mu^2` = 1/sqrt(x), 
         tan_half = 2 * atan(x), logit = standist::inv_logit(x), probit = pnorm(x), 
         cauchit = pcauchy(x), cloglog = inv_cloglog(x), probit_approx = pnorm(x), 
        softplus = log1p_exp(x), stop2("Link '", link, "' not supported."))
  }
## ----end

posterior_fitted.inla <- function(object, newdata = NULL, ndraws = 1000, form = NULL) {
  draws <- inla.posterior.sample(n=ndraws, result = object)
  contents <- object$misc$configs$contents

  if (is.null(form)) form <- object$.args$formula
  gf <- INLA:::inla.interpret.formula(form)
  Xmat <- model.matrix(update(gf$fixf, NULL ~ .), newdata)
  nms <- colnames(Xmat)
  i_lp <- contents$start[contents$tag %in% nms]
  lp <- t(sapply(draws, function(x) x$latent[i_lp]))
  b <- tcrossprod(Xmat, lp)
  link <- object$misc$linkfunctions$names
  out <- my_ilink(b, link) 
  cellmeans.full <- newdata |> bind_cols(as.data.frame(out)) |>
    pivot_longer(cols=matches("^V[0-9]*$"),
      names_to='.draws',
      values_to='Values') |> 
    group_by(across(all_of(colnames(newdata))))
  cellmeans <-
          cellmeans.full |>
          posterior::as_draws() |>
          dplyr::select(-.draw) |>
          dplyr::mutate(.draw = as.integer(str_replace(.draws, "V", ""))) |> 
          dplyr::select(-.draws) 
  
}

contr_inla <- function(mod, newdata) {
  ## newdata <- newdata |>
  ##   filter(Decade != first(Decade)) |>
  ##   droplevels()
  contr <- posterior_fitted.inla(mod, newdata, form = value ~ Decade) |>
    group_by(.draw) |>
    ## reframe(
    summarise(
      Decade = Decade[-1],
      Ratio = exp(log(Values[-1]) - log(Values[1])),
      Values = Values[-1] - Values[1]
    ) |>
    group_by(Decade) |>
    tidybayes::summarise_draws(mean, HDInterval::hdi,
      pl1 = ~ mean(.x < 1),
      pg1 = ~ mean(.x > 1),
      pl0 = ~ mean(.x < 0),
      pg0 = ~ mean(.x > 0)
      )
  contr
}

contr_glmmTMB <- function(mod, newdata) {
  e <- tryCatch(
  {
    lookup <- c("mean" = "estimate", "mean" = "ratio", "mean" = "odds.ratio", "lower" = "asymp.LCL", "upper" = "asymp.UCL",
                "lower" = "lower.CL", "upper" = "upper.CL")
    cmat <- cbind(c(-1, 1, 0, 0), c(-1, 0, 1, 0), c(-1, 0, 0, 1))
    contr1 <- emmeans(mod, ~Decade, type = "response") |>
      contrast(method = list(Decade = cmat)) |>
      summary(infer = TRUE) |>
      as.data.frame() |>
      dplyr::rename(any_of(lookup)) |>
      mutate(Decade = newdata$Decade[-1]) |>
      dplyr::select(Decade, mean, lower, upper) |>
      mutate(variable = "Ratio")
    
    contr2 <- emmeans(mod, ~Decade) |>
      regrid() |> 
      contrast(method = list(Decade = cmat)) |>
      summary(infer = TRUE) |>
      as.data.frame() |>
      dplyr::rename(any_of(lookup)) |>
      mutate(Decade = newdata$Decade[-1]) |>
      dplyr::select(Decade, mean, lower, upper) |>
      mutate(variable = "Values")
    contr <- rbind(contr1, contr2)
  },
    error = function(e) e
  )
  if ("simpleError" %in% class(e)) {
    return(NULL)
  } else {
    return(contr)
  }
}


## DHARMa diagnostics
calc_dharma <- function(mod, .x, inewdata) {
  dnm <- str_replace_all(.x, "mod_", "dharma_")
  dnm <- str_replace_all(dnm, ".RData", ".pdf")
  if (is.null(mod)) {
    dharma <- NULL
  } else {
    if (str_detect(.x, "inla")) {
      dharma <- dharma_inla(mod, inewdata)
    } else {
      dharma <- dharma_glmmTMB(mod)
    }
  }
  ggsave(filename = dnm, dharma, width = 10, height = 4)
  dnm
}


dharma_inla <- function(mod, inewdata) {
  e <- tryCatch(
  {
  draws <- inla.posterior.sample(n=1000, result = mod)
  contents <- mod$misc$configs$contents

  cnt <- contents$start[1]:(inewdata[1]-1)
  if (mod$.args$family %in% c('binomial', 'betabinomial')) {
    lp <- plogis(sapply(draws, function(x) x$latent[cnt]))
    response <- mod$.args$data$HC[cnt]/100
    fitted_means <- apply(lp, 1, mean)
    
  } else {   #fishes
    lp <- exp(sapply(draws, function(x) x$latent[cnt]))
    response <- mod$.args$data$value[cnt]
    fitted_means <- apply(lp, 1, mean)
  }
  mod_resids <- DHARMa::createDHARMa(
    simulatedResponse = lp,
    observedResponse = response,
    fittedPredictedResponse = fitted_means,
    integerResponse = TRUE
  )
  
  p <- patchwork::wrap_elements(~DHARMa::testUniformity(mod_resids)) +
    patchwork::wrap_elements(~DHARMa::plotResiduals(mod_resids)) +
    patchwork::wrap_elements(~DHARMa::testDispersion(mod_resids)) 

  },
  error = function(e) e
  )
  if ("simpleError" %in% class(e)) {
    return(NULL)
  } else {
    return(p)
  }
}


dharma_glmmTMB <- function(mod) {
  e <- tryCatch(
  {
    mod_resids <- DHARMa::simulateResiduals(mod, plot = F)
    p <- patchwork::wrap_elements(~DHARMa::testUniformity(mod_resids)) +
        patchwork::wrap_elements(~DHARMa::plotResiduals(mod_resids)) +
        patchwork::wrap_elements(~DHARMa::testDispersion(mod_resids)) 

  },
    error = function(e) e
  )
  if ("simpleError" %in% class(e)) {
    return(NULL)
  } else {
    return(p)
  }
}








## zero inflated negative binomial 2
fit_model_type1_inla <- function(data.1, form, family = "poisson") {
  e <- tryCatch({
    cat("trying zeroinflatednbinomial2 ======\n")
    if (is.null(form)) {
      form <- value ~ Decade * scale(log(HC + 0.01)) +
        f(AIMS_REEF_NAME, model = "iid") +
        f(Site, model = "iid") +
        f(Transect, model = "iid")
    }
    mod <- inla(
      form,
      data = data.1,
      family = "zeroinflatednbinomial2",
      control.predictor = list(link = 1, compute = TRUE),
      control.compute = list(
        dic = TRUE, cpo = TRUE, waic = TRUE,
        config = TRUE)
    )
  }, error =  function(e) e)

  if (!("inla" %in% class(e))) {
    cat("Model failed ===\n\n")
    return(NULL)
  } else {
    return(mod)
  }
}

## zero inflated negative binomial 1
fit_model_type2_inla <- function(data.1) {
  mod <- inla(
    value ~ Decade * scale(log(HC + 0.01)) +
      f(AIMS_REEF_NAME, model = "iid") +
      f(Site, model = "iid") +
      f(Transect, model = "iid"),
    data = data.1,
    family = "zeroinflatednbinomial1",
    control.predictor = list(link = 1, compute = TRUE),
    control.compute = list(
      dic = TRUE, cpo = TRUE, waic = TRUE,
      config = TRUE)
  )
  mod
}

## zero inflated negative binomial 0
fit_model_type3_inla <- function(data.1) {
  mod <- inla(
    value ~ Decade * scale(log(HC + 0.01)) +
      f(AIMS_REEF_NAME, model = "iid") +
      f(Site, model = "iid") +
      f(Transect, model = "iid"),
    data = data.1,
    family = "zeroinflatednbinomial0",
    control.predictor = list(link = 1, compute = TRUE),
    control.compute = list(
      dic = TRUE, cpo = TRUE, waic = TRUE,
      config = TRUE)
  )
  mod
}

## zero inflated poisson 0
fit_model_type4_inla <- function(data.1) {
  mod <- inla(
    value ~ Decade * scale(log(HC + 0.01)) +
      f(AIMS_REEF_NAME, model = "iid") +
      f(Site, model = "iid") +
      f(Transect, model = "iid"),
    data = data.1,
    family = "zeroinflatedpoisson0",
    control.predictor = list(link = 1, compute = TRUE),
    control.compute = list(
      dic = TRUE, cpo = TRUE, waic = TRUE,
      config = TRUE)
  )
  mod
}

## zero inflated poisson 1
fit_model_type5_inla <- function(data.1) {
  mod <- inla(
    value ~ Decade * scale(log(HC + 0.01)) +
      f(AIMS_REEF_NAME, model = "iid") +
      f(Site, model = "iid") +
      f(Transect, model = "iid"),
    data = data.1,
    family = "zeroinflatedpoisson1",
    control.predictor = list(link = 1, compute = TRUE),
    control.compute = list(
      dic = TRUE, cpo = TRUE, waic = TRUE,
      config = TRUE)
  )
  mod
}

## negative binomial
fit_model_type6_inla <- function(data.1) {
  mod <- inla(
    value ~ Decade * scale(log(HC + 0.01)) +
      f(AIMS_REEF_NAME, model = "iid") +
      f(Site, model = "iid") +
      f(Transect, model = "iid"),
    data = data.1,
    family = "nbinomial",
    control.predictor = list(link = 1, compute = TRUE),
    control.compute = list(
      dic = TRUE, cpo = TRUE, waic = TRUE,
      config = TRUE)
  )
  mod
}

## negative binomial 2
fit_model_type7_inla <- function(data.1) {
  mod <- inla(
    value ~ Decade * scale(log(HC + 0.01)) +
      f(AIMS_REEF_NAME, model = "iid") +
      f(Site, model = "iid") +
      f(Transect, model = "iid"),
    data = data.1,
    family = "nbinomial2",
    control.predictor = list(link = 1, compute = TRUE),
    control.compute = list(
      dic = TRUE, cpo = TRUE, waic = TRUE,
      config = TRUE)
  )
  mod
}

## poisson
fit_model_type8_inla <- function(data.1) {
  mod <- inla(
    value ~ Decade * scale(log(HC + 0.01)) +
      f(AIMS_REEF_NAME, model = "iid") +
      f(Site, model = "iid") +
      f(Transect, model = "iid"),
    data = data.1,
    family = "poisson",
    control.predictor = list(link = 1, compute = TRUE),
    control.compute = list(
      dic = TRUE, cpo = TRUE, waic = TRUE,
      config = TRUE)
  )
  mod
}


## fit_model_inla <- function(data.1) {

##   e <- tryCatch({
##     cat("trying zeroinflatednbinomial2 ======\n")
##     mod <- fit_model_type1_inla(data.1)
##   }, error =  function(e) e)
##   print(class(e))
  
## }

fit_model_old <- function(data.1, Species) {
  e <- tryCatch({
    cat("trying zeroinflatednbinomial2 ======\n")
    mod <- inla(
      value ~ Decade + scale(log(HC)) +
        f(AIMS_REEF_NAME, model = "iid") +
        f(Site, model = "iid") +
        f(Transect, model = "iid"),
      data = data.1,
      family = "zeroinflatednbinomial2",
      ## family = "zeroinflatednbinomial1",
      control.predictor = list(link = 1, compute = TRUE),
      control.compute = list(
        dic = TRUE, cpo = TRUE, waic = TRUE,
        config = TRUE)
      )
  }, error =  function(e) e)
  print(class(e))
  if (!("inla" %in% class(e))) {
    e <- tryCatch({
    cat("trying zeroinflatednbinomial1 ======\n")
      mod <- inla(
        value ~ Decade + scale(log(HC)) +
          f(AIMS_REEF_NAME, model = "iid") +
          f(Site, model = "iid") +
          f(Transect, model = "iid"),
        data = data.1,
        family = "zeroinflatednbinomial1",
        control.predictor = list(link = 1, compute = TRUE),
        control.compute = list(
          dic = TRUE, cpo = TRUE, waic = TRUE,
          config = TRUE)
      )
    }, error =  function(e) e)
  }
  print(class(e))
  if (!("inla" %in% class(e))) {
    e <- tryCatch({
      cat("trying zeroinflatedpoisson1 ======\n")
      mod <- inla(
        value ~ Decade + scale(log(HC)) +
          f(AIMS_REEF_NAME, model = "iid") +
          f(Site, model = "iid") +
          f(Transect, model = "iid"),
        data = data.1,
        family = "zeroinflatedpoisson1",
        control.predictor = list(link = 1, compute = TRUE),
        control.compute = list(
          dic = TRUE, cpo = TRUE, waic = TRUE,
          config = TRUE)
      )
    }, error =  function(e) e)
  }
  print(class(e))
  if (!("inla" %in% class(e))) {
    e <- tryCatch({
      cat("trying nbinomial2 ======\n")
      mod <- inla(
        value ~ Decade + scale(log(HC)) +
          f(AIMS_REEF_NAME, model = "iid") +
          f(Site, model = "iid") +
          f(Transect, model = "iid"),
        data = data.1,
        family = "nbinomial2",
        control.predictor = list(link = 1, compute = TRUE),
        control.compute = list(
          dic = TRUE, cpo = TRUE, waic = TRUE,
          config = TRUE)
      )
    }, error =  function(e) e)
  }
  print(class(e))
  summary(mod)
  nm <- paste0("../data/modelled/mod_", Species, ".RData")
  save(mod, file = nm)
  nm
}


fit_models_inla_hcc <- function(data.1 = data.1, form = n.points ~ 1) {
  fams <- c(
    "binomial",
    "betabinomial"
  )
  nms <- NULL
  for (f in fams) {
    mod_inla <- fit_model_inla_hcc(data.1, form = form, family = f)
    nm <- paste0("../data/modelled/mod_inla_hcc_", f, ".RData")
    saveRDS(mod_inla, file = nm)
    nms <- c(nms, nm)
  }
  nms
}

fit_model_inla_hcc <- function(data.1 = data.1, form = n.points ~ 1, family = "binomial") {
  e <- tryCatch({
    cat(paste0("trying ", family, " ======\n"))

    form <- as.formula(form)
    mod <- inla(
      form,
      data = data.1,
      family = family,
      Ntrials = data.1$total.points,
      control.predictor = list(link = 1, compute = TRUE),
      control.compute = list(
        dic = TRUE, cpo = TRUE, waic = TRUE,
        config = TRUE)
    )
  }, error =  function(e) e)

  if (!("inla" %in% class(e))) {
    cat("Model failed ===\n\n")
    return(NULL)
  } else {
    return(mod)
  }
}

fit_models_glmmTMB_hcc <- function(data.1 = data.1, form = value ~ 1) {
  fams <- c(
    "binomial",
    "betabinomial"
  )
  nms <- NULL
  for (i in 1:length(fams)) {
    mod_glmmTMB <- fit_model_glmmTMB_hcc(data.1 = data.1, form = form, family = fams[[i]])
    nm <- paste0("../data/modelled/mod_glmmTMB_hcc_", fams[[i]], ".RData")
    saveRDS(mod_glmmTMB, file = nm)
    nms <- c(nms, nm)
  }
  nms
}


fit_model_glmmTMB_hcc <- function(data.1, form = n.points ~ 1, family = "binomial") {
  e <- tryCatch({
    cat(paste0("trying ", family, " ======\n"))
    form <- as.formula(form)
    mod <- glmmTMB(
      form,
      data = data.1,
      family = family,
      REML = TRUE
    )
  }, error =  function(e) e)

  if (!("glmmTMB" %in% class(e))) {
    cat("Model failed ===\n\n")
    return(NULL)
  } else {
    return(mod)
  }
}
