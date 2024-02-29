## Retrieve the data
load(file = "../data/processed/data_q2.R")


## Strip out all unnecessary variables

full_data <- data
data <- full_data

## Switch to zero-inflated model (intercept only)
{
  ## Data preparation
  {
    data <- full_data
  }
  ## Raw means
  {
    cellmeans_summ_raw <- data |>
      reframe(
        type = c("mean", "median"),
        Mean = c(mean(n.points / total.points), median(n.points / total.points)),
        SD = c(sd(n.points / total.points), MAD = mad(n.points / total.points)),
        N = c(n(), n())
      ) |>
      mutate(
        lower = Mean - 2 * (SD / sqrt(N)),
        upper = Mean + 2 * (SD / sqrt(N))
      )
    save(cellmeans_summ_raw, file = "../data/modelled/cellmeans_summ_raw_1.2.RData")
  }
  ## glmmTMB
  {
    ## Fit model
    {
      if (rerun_models) {
        mod_glmmTMB <- glmmTMB(
          cbind(n.points, total.points - n.points) ~ 1 +
            (1 | AIMS_REEF_NAME) +
            (1 | Site) +
            (1 | Transect),
          ziformula = ~1,
          data = data,
          family = "binomial",
          REML = TRUE
        )

        summary(mod_glmmTMB)
        save(mod_glmmTMB, file = "../data/modelled/mod_glmmTMB_1.2.RData")
      }
      ## Partial effects
      {
        load(file = "../data/modelled/mod_glmmTMB_1.2.RData")
        cellmeans_summ_glmmTMB <- mod_glmmTMB |>
          emmeans(~1, type = "response") |>
          as.data.frame()
        save(cellmeans_summ_glmmTMB, file = "../data/modelled/cellmeans_summ_glmmTMB_1.2.RData")
      }
    }
  }
  ## brms
  {
    ## Fit model
    {
      if (rerun_models) {
        form <- bf(
          n.points | trials(total.points) ~ 1 +
            (1 | AIMS_REEF_NAME) +
            (1 | Site) +
            (1 | Transect),
          zi = ~1,
          family = "zero_inflated_binomial"
        )
        priors <- prior(normal(0, 1), class = "Intercept") +
          prior(student_t(3, 0, 1), class = "sd") +
          prior(logistic(0, 1), class = "Intercept", dpar = "zi")
        mod_brm <- brm(form,
          data = data,
          prior = priors,
          iter = 5000, warmup = 1000,
          chains = 3, cores = 3,
          thin = 4,
          backend = "cmdstanr",
          control = list(adapt_delta = 0.99),
          silent = 0 # ,
          ## refresh = 100
        )
        summary(mod_brm)
        save(mod_brm, file = "../data/modelled/mod_brm_1.2.RData")
        load(file = "../data/modelled/mod_brm_1.2.RData")
      }
    }
    ## Partial effects
    {
      load(file = "../data/modelled/mod_brm_1.2.RData")
      cellmeans_summ_brm <- mod_brm |>
        emmeans(~1, type = "response") |>
        as.data.frame()
      save(cellmeans_summ_brm, file = "../data/modelled/cellmeans_summ_brm_1.2.RData")
    }
  }
  ## INLA
  {
    ## Prepare data
    {
      data.pred <- data |>
        dplyr::select(
          n.points, total.points, AIMS_REEF_NAME,
          Site, Transect
        )
      newdata <- data.frame(AIMS_REEF_NAME = NA, Site = NA, Transect = NA)
      data_pred <- data |> bind_rows(newdata)
      i_newdata <- (nrow(data) + 1):nrow(data_pred)
    }
    ## Fit model
    {
      if (rerun_models) {
        mod <- inla(
          n.points ~ 1 +
            f(model = "iid", AIMS_REEF_NAME) +
            f(model = "iid", Site) +
            f(model = "iid", Transect),
          data = data.pred,
          Ntrials = data.pred$total.points,
          family = "zeroinflatedbinomial2",
          control.predictor = list(link = 1, compute = TRUE),
          control.compute = list(
            config = TRUE,
            dic = TRUE, waic = TRUE, cpo = TRUE
          )
        )
        summary(mod)
        ## autoplot(mod)
        save(mod, file = "../data/modelled/mod_1.2.RData")
      }
    }
    ## Partial effects
    {
      load(file = "../data/modelled/mod_1.2.RData")
      draws <- inla.posterior.sample(n = 1000, result = mod)
      contents <- mod$misc$configs$contents
      i_newdata <- contents$start[contents$tag == "(Intercept)"]
      cellmeans <- newdata |> cbind(t(sapply(draws, function(x) x$latent[i_newdata])))
      cellmeans_summ_inla <- cellmeans |>
        pivot_longer(cols = matches("^[0-9]*$"), names_to = ".draws") |>
        posterior::as_draws() |>
        mutate(.draw = .draws) |>
        dplyr::select(-.draws, -AIMS_REEF_NAME, -Site, -Transect) |>
        mutate(value = plogis(value)) |>
        summarise_draws(median, HDInterval::hdi)

      save(cellmeans_summ_inla, file = "../data/modelled/cellmeans_summ_inla_1.2.RData")
    }
  }
  ## Comparisons
  {
    
    load(file = "../data/modelled/cellmeans_summ_raw_1.2.RData")
    load(file = "../data/modelled/cellmeans_summ_glmmTMB_1.2.RData")
    load(file = "../data/modelled/cellmeans_summ_brm_1.2.RData")
    load(file = "../data/modelled/cellmeans_summ_inla_1.2.RData")

    cellmeans_summ <- bind_rows(
      cellmeans_summ_raw |>
        mutate(method = "raw") |>
        dplyr::rename(median = Mean, lower = lower, upper = upper),
      cellmeans_summ_glmmTMB |>
        mutate(method = "glmmTMB", type = "median") |>
        dplyr::rename(median = prob, lower = asymp.LCL, upper = asymp.UCL),
      cellmeans_summ_brm |>
        mutate(method = "brm", type = "median") |>
        dplyr::rename(median = prob, lower = lower.HPD, upper = upper.HPD),
      cellmeans_summ_inla |>
        mutate(method = "inla", type = "median") |> 
        dplyr::rename(median = median, lower = lower, upper = upper)
    )
    cellmeans_summ

    cellmeans_summ |> ggplot(aes(y = median, x = method, shape = type)) +
            geom_pointrange(aes(ymin = lower, ymax = upper))

  }
}
