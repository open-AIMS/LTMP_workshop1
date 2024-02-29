## Retrieve the data
load(file = "../data/processed/data_q2.R")


## Strip out all unnecessary variables

full_data <- data
data <- full_data


## Start with an intercept only model
{
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
    save(cellmeans_summ_raw, file = "../data/modelled/cellmeans_summ_raw_1.1.RData")
  }
  ## glmmTMB
  {
    ## fit model
    {
      if (rerun_models) {
        mod_glmmTMB <- glmmTMB(cbind(n.points, total.points-n.points) ~ 1 +
                                 (1|AIMS_REEF_NAME) +
                                 (1|Site) +
                                 (1|Transect),
          data = data,
          family = "binomial", 
          REML = TRUE
        )

        summary(mod_glmmTMB)
        resids <- simulateResiduals(mod_glmmTMB, plot = TRUE)
        testDispersion(resids)
        testZeroInflation(resids)
        save(mod_glmmTMB, file = "../data/modelled/mod_glmmTMB_1.1.RData")
      }
    }
    ## Partial effects
    {
      load(file = "../data/modelled/mod_glmmTMB_1.1.RData")
      cellmeans_summ_glmmTMB <- mod_glmmTMB |>
        emmeans(~1, type = "response") |>
        as.data.frame()
      save(cellmeans_summ_glmmTMB, file = "../data/modelled/cellmeans_summ_glmmTMB_1.1.RData")
    }
  }
  ## brms
  {
    ## fit model
    {
      if (rerun_models) {
        mod_brm <- brm(
          bf(
            n.points | trials(total.points) ~ 1 +
              (1 | AIMS_REEF_NAME) +
              (1 | Site) +
              (1 | Transect),
            family = "binomial"
          ),
          data = data,
          iter = 5000, warmup = 1000,
          chains = 3, cores = 3,
          thin = 4,
          backend = "cmdstanr",
          control = list(adapt_delta = 0.99)
        )
        summary(mod_brm)
        save(mod_brm, file = "../data/modelled/mod_brm_1.1.RData")
      }
    }
    ## Partial effects
    {
      load(file = "../data/modelled/mod_brm_1.1.RData")
      cellmeans_summ_brm <- mod_brm |>
        emmeans(~1, type = "response") |>
        as.data.frame()
      save(cellmeans_summ_brm, file = "../data/modelled/cellmeans_summ_brm_1.1.RData")
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
        mod <- inla(n.points ~ 1 +
                      f(model = "iid", AIMS_REEF_NAME) +
                      f(model = "iid", Site) +
                      f(model = "iid", Transect),
          data = data.pred,
          Ntrials = data.pred$total.points,
          family = "binomial", 
          control.predictor = list(link = 1, compute = TRUE),
          control.compute = list(config = TRUE,
            dic = TRUE, waic = TRUE, cpo = TRUE
          )
        )
        summary(mod)
        ## autoplot(mod)
        save(mod, file = "../data/modelled/mod_1.1.RData")
      }
    }
    ## Diagnostics
    {
      load(file = "../data/modelled/mod_1.1.RData")
      wrap_plots(
        pit_qq_plot(mod, i.mod = 1:nrow(data), logit_scale = TRUE),
        pit_plot(mod, i.mod = 1:nrow(data)),
        pit_resid_plot(mod, i.mod = 1:nrow(data)),
        pit_hist_plot(mod, i.mod = 1:nrow(data))
      ) +
      plot_layout(ncol = 3)

      mod_cpo <- data.frame(cpo = mod$cpo$cpo,
        pit = mod$cpo$pit,
        failure = mod$cpo$failure) |>
        filter(failure == 0) |>
        mutate(row = 1:n())

      mod$cpo$failure

      g1 <- mod_cpo |>
        ggplot(aes(y = cpo, x = row)) +
        geom_point()
      g2 <- mod_cpo |>
        ggplot(aes(x = cpo)) +
        geom_histogram()
      g3 <- mod_cpo |>
        ggplot(aes(y = pit, x = row)) +
        geom_point()
      g4 <- mod_cpo |>
        ggplot(aes(x = pit)) +
        geom_histogram()

      (g1 + g2)/(g3 + g4)

   ## Ideally, we want to see that no CPO or PIT is very different in
   ## value from any others (often not informative in binomial models)
   ## and that the histograms are relatively flat (not really the case
   ## here).

   ## ggplot_inla_residuals(mod, observed = data$n.points)
   ## ggplot_inla_residuals2(mod, observed = data$n.points)
    }
    ## DHARMa 
    {
      load(file = "../data/modelled/mod_1.1.RData")
      ## preds <- posterior_predict.inla(mod, newdata = data)
      ## mod_resids <- createDHARMa(simulatedResponse = t(preds),
      ##   observedResponse = data$n.points,
      ##   fittedPredictedResponse = apply(preds, 2, mean),
      ##   integerResponse = TRUE)
      ## wrap_elements(~testUniformity(mod_resids)) +
      ##   wrap_elements(~plotResiduals(mod_resids)) +
      ##   wrap_elements(~testDispersion(mod_resids)) 

      ## testDispersion(mod_resids)
      ## testZeroInflation(mod_resids)
    }
    ## Partial effects
    {
      load(file = "../data/modelled/mod_1.1.RData")
      draws <- inla.posterior.sample(n=1000, result = mod)
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

      save(cellmeans_summ_inla, file = "../data/modelled/cellmeans_summ_inla_1.1.RData")
    }
  }
  ## Comparisons
  {
    load(file = "../data/modelled/cellmeans_summ_raw_1.1.RData")
    load(file = "../data/modelled/cellmeans_summ_glmmTMB_1.1.RData")
    load(file = "../data/modelled/cellmeans_summ_brm_1.1.RData")
    load(file = "../data/modelled/cellmeans_summ_inla_1.1.RData")

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
