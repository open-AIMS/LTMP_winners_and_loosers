fit_model <- function(data.1, Species) {
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

