library(tidyverse)
library(INLA)
library(readxl)

## use annual_report.sif

## Hi Murray,

## This is a question relating to the same paper Mike was recently
## running PCAs and adonis results for me – the last thing I’d like to
## do is identify species-level ‘winners and losers’ across decades
## (1990s, 2000s, 2010s, 2020s). I’ve started by running frequentist
## glmms with the orig.ltmp.w dataset (see results in the ‘Species
## table’ tab for result) because brms was taking so long I basically
## never saw a result. I thought as well that I should probably add
## coral cover so I’ve gone as far as creating the second attached
## data file. What I’d really like is a caterpillar plot with a
## Bayesian hierarchical model using Decade and HC as fixed factors
## and Reef, Site, Transect as random factors, displayed as a
## caterpillar plot with species split by trophic group (the orig
## column) – or a multi-panel figure, each panel a different trophic
## group. Would you have time to do something like that for me?

## Cheers

## Dani



load("../data/primary/fish.benthic.2024.RData")
load("../data/primary/240417orig.ltmp.w.RData")

lookup <- readxl::read_xlsx("../data/primary/Species by decade glmm.xlsx",
  sheet = "Species table"
)

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


## get a list of species
species <- orig.ltmp.w |>
  colnames() |>
  str_subset("[A-Z][a-z]+\\_[a-z]+")


data <- orig.ltmp.w |>
  filter(REPORT_YEAR < 2024) |>
  droplevels() |> 
  full_join(benthic, by = c("AIMS_REEF_NAME", "SITE_NO", "TRANSECT_NO", "REPORT_YEAR"))
  

data <- data |>
  mutate(
    Site = paste(REEF_NAME, SITE_NO),
    Transect = paste(Site, TRANSECT_NO)
  ) |>
  filter(!is.na(SITE_NO), !is.na(HC)) |>
  droplevels() 

data.1 <- data |>
  pivot_longer(cols = species, names_to = "Species") |>
  nest_by(Species)
data.1

data.1 <- data.1 |>
  mutate(newdata = list({
    newdata <- data.frame(Decade = c("1990s", "2020s"))
    newdata
  }))

data.1 <- data.1 |>
  mutate(data.1 = list({
    data |> bind_rows(newdata)
  }))

data.1 <- data.1 |>
  mutate(lincomb = list({
    nd <- data.frame(Decade = c(1, 0, 0, 1, 0))
    inla.make.lincombs(as.data.frame(nd))
  }))
  
fit_model <- function(data.1, Species) {
  e <- tryCatch({
    cat("trying zeroinflatednbinomial1 ======\n")
    mod <- inla(
      value ~ Decade + scale(log(HC)) +
        f(AIMS_REEF_NAME, model = "iid") +
        f(Site, model = "iid") +
        f(Transect, model = "iid"),
      data = data.1,
      ## family = "zeroinflatednbinomial2",
      family = "zeroinflatednbinomial1",
      control.predictor = list(link = 1, compute = TRUE),
      control.compute = list(
        dic = TRUE, cpo = TRUE, waic = TRUE,
        config = TRUE)
      )
  }, error =  function(e) e)
  print(class(e))
  if (!("inla" %in% class(e))) {
    e <- tryCatch({
    cat("trying zeroinflatednbinomial2 ======\n")
      mod <- inla(
        value ~ Decade + scale(log(HC)) +
          f(AIMS_REEF_NAME, model = "iid") +
          f(Site, model = "iid") +
          f(Transect, model = "iid"),
        data = data.1,
        family = "zeroinflatednbinomial2",
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

data.1 <- data.1 |>
  mutate(model = list({
    cat(Species, "\n")
    fit_model(data.1, Species)
  }))
data.2 <- data.1[-1:-17,] |>
  mutate(model = list({
    cat(Species, "\n")
    fit_model(data.1, Species)
  }))
data.2 <- data.1[-1:-22,] |>
  mutate(model = list({
    cat(Species, "\n")
    fit_model(data.1, Species)
  }))
data.2 <- data.1[-1:-65,] |>
  mutate(model = list({
    cat(Species, "\n")
    fit_model(data.1, Species)
  }))

d <- data.1[24,'data.1'][[1]][[1]]
S <- data.1[24,'Species'][[1]][[1]]

fit_model(d, S)

d <- data.1[1,'data.1'][[1]][[1]]
S <- data.1[1,'Species'][[1]][[1]]

fit_model(d, S)

## dd <- data.frame(y =  rnorm(10), x = rnorm(10))
## inla(y ~ x, data =  dd)

## mod <- inla(
##   value ~ Decade + scale(log(HC)) +
##       f(AIMS_REEF_NAME, model = "iid") +
##       f(Site, model = "iid") +
##       f(Transect, model = "iid"),
##   data = d,
##   family = "zeroinflatednbinomial1",
##   ## lincomb =  lincomb,
##   ## family = "nbinomial2",
##   control.predictor = list(link = 1, compute = TRUE),
##   control.compute = list(
##     dic = TRUE, cpo = TRUE, waic = TRUE,
##     config = TRUE)
## )

## summary(mod)


##   mod <- inla(
##     value ~ Decade + scale(log(HC)) +
##       f(AIMS_REEF_NAME, model = "iid") +
##       f(Site, model = "iid") +
##       f(Transect, model = "iid"),
##     data = d,
##     ## family = "zeroinflatednbinomial2",
##     ## lincomb =  lincomb,
##     family = "zeroinflatednbinomial1",
##     control.predictor = list(link = 1, compute = TRUE),
##     control.compute = list(
##       dic = TRUE, cpo = TRUE, waic = TRUE,
##       config = TRUE)
##   )
