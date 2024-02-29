## Retrieve the data
load(file = "../data/processed/data_q2.R")


## Strip out all unnecessary variables

full_data <- data
data <- full_data

## Add Before/After
{
  ## Data preparation
  {
    data <- full_data
  }
  ## Raw means
  {
    ## Cellmeans
    {
      cellmeans_summ_raw <- data |>
        group_by(Dist.time) |> 
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
      save(cellmeans_summ_raw, file = "../data/modelled/cellmeans_summ_raw_2.1.RData")

      eff_summ_raw <- cellmeans_summ_raw |>
        group_by(type) |>
        summarise(Mean = Mean[Dist.time == "After"] / Mean[Dist.time == "Before"])
      save(eff_summ_raw, file = "../data/modelled/eff_summ_raw_2.1.RData")
    }
    ## Effects
    {
      eff_summ_raw <-
        data |>
        ## filter(Transect == "Feather Reef 1 5", Dist.number == 2) |> 
        mutate(cover = n.points / total.points) |> 
        dplyr::select(Dist.time, AIMS_REEF_NAME, Site, Transect, Dist.number,
          s, c, d, b, u, cover) |> 
        mutate(Dist = case_when(s == 1 ~ "s",
          c == 1 ~ "c",
          d == 1 ~ "d",
          b == 1 ~ "b",
          u == 1 ~ "u")) |>
        dplyr::select(-any_of(c("s", "c", "d", "b", "u"))) |>
        group_by(AIMS_REEF_NAME, Site, Transect, Dist.number, Dist) |>
        reframe(Value = mean(cover[Dist.time == "After"]) - mean(cover[Dist.time == "Before"])) |> 
        ungroup() |> 
        reframe(
          type = c("mean", "median"),
          Mean = c(mean(Value, na.rm = TRUE), median(Value, na.rm = TRUE)),
          SD = c(sd(Value, na.rm = TRUE), MAD = mad(Value, na.rm = TRUE)),
          N = c(n(), n())
        ) |>
        mutate(
          lower = Mean - 2 * (SD / sqrt(N)),
          upper = Mean + 2 * (SD / sqrt(N))
        )
      save(eff_summ_raw, file = "../data/modelled/eff_summ_raw_2.1a.RData")
    }
  }
  ## glmmTMB
  {
    ## Fit model
    {
      if (rerun_models) {
        mod_glmmTMB <- glmmTMB(
          cbind(n.points, total.points - n.points) ~ Dist.time +
            (1 | AIMS_REEF_NAME) +
            (1 | Site) +
            (1 | Transect),
          ziformula = ~ 1 + (1 | AIMS_REEF_NAME) + (1 | Site) + (1 | Transect),
          data = data,
          family = "binomial",
          REML = TRUE
        )

        summary(mod_glmmTMB)
        save(mod_glmmTMB, file = "../data/modelled/mod_glmmTMB_2.1.RData")
      }
    }
    ## DHARMa
    {
      resids <- simulateResiduals(mod_glmmTMB, plot = TRUE)
      wrap_elements(~testUniformity(resids)) +
        wrap_elements(~plotResiduals(resids)) +
        wrap_elements(~testDispersion(resids)) 
      testDispersion(resids)
      testZeroInflation(resids)
    }
    ## Partial plots
    {
      load(file = "../data/modelled/mod_glmmTMB_2.1.RData")
      cellmeans_summ_glmmTMB <- mod_glmmTMB |>
        emmeans(~Dist.time, type = "response") |>
        as.data.frame()
      save(cellmeans_summ_glmmTMB, file = "../data/modelled/cellmeans_summ_glmmTMB_2.1.RData")

      cellmeans_summ_glmmTMB |>
        ggplot(aes(y = prob, x = Dist.time)) +
        geom_pointrange(aes(ymin = asymp.LCL, ymax = asymp.UCL))
    }
    ## Contrasts
    {
      load(file = "../data/modelled/mod_glmmTMB_2.1.RData")
      eff_summ_glmmTMB <-
        mod_glmmTMB |>
        emmeans(~Dist.time, type = "response") |>
        contrast(method = list(Dist.time = c(-1, 1))) |>
        summary(infer = TRUE) |>
        as.data.frame()
      save(eff_summ_glmmTMB, file = "../data/modelled/eff_summ_glmmTMB_2.1.RData")

      eff_summ_glmmTMB |>
        ggplot(aes(x = odds.ratio, y = contrast)) +
        geom_vline(xintercept = 1, linetype = "dashed") +
        geom_pointrange(aes(xmin = asymp.LCL, xmax = asymp.UCL)) +
        scale_x_continuous("Effect (Before - After) on a fold scale",
          trans = "log2", breaks = scales::extended_breaks(n = 8)
        )
    }
  }
  ## brms
  {
    ## Fit model
    {
      if (rerun_models) {
        form <- bf(
          n.points | trials(total.points) ~ Dist.time +
            (1 | AIMS_REEF_NAME) +
            (1 | Site) +
            (1 | Transect),
          zi = ~ 1 + (1 | AIMS_REEF_NAME) + (1 | Site) + (1 | Transect),
          family = "zero_inflated_binomial"
        )
        priors <- prior(normal(0, 1), class = "Intercept") +
          prior(normal(0, 1), class = "b") +
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
        save(mod_brm, file = "../data/modelled/mod_brm_2.1.RData")
        load(file = "../data/modelled/mod_brm_2.1.RData")
      }
    }
    ## Partial plots
    {
      load(file = "../data/modelled/mod_brm_2.1.RData")
      cellmeans_summ_brm <- mod_brm |>
        emmeans(~Dist.time, type = "response") |>
        as.data.frame()
      save(cellmeans_summ_brm, file = "../data/modelled/cellmeans_summ_brm_2.1.RData")

      cellmeans_summ_brm |>
        ggplot(aes(y = prob, x = Dist.time)) +
        geom_pointrange(aes(ymin = lower.HPD, ymax = upper.HPD))
    }
    ## Contrasts
    {
      load(file = "../data/modelled/mod_brm_2.1.RData")
      eff_summ_brm <-
        mod_brm |>
        emmeans(~Dist.time, type = "response") |>
        contrast(method = list(Dist.time = c(-1, 1))) |>
        summary(infer = TRUE) |>
        as.data.frame()
      save(eff_summ_brm, file = "../data/modelled/eff_summ_brm_2.1.RData")

      eff_summ_brm |>
        ggplot(aes(x = odds.ratio, y = contrast)) +
        geom_vline(xintercept = 1, linetype = "dashed") +
        geom_pointrange(aes(xmin = lower.HPD, xmax = upper.HPD)) +
        scale_x_continuous("Effect (Before - After) on a fold scale",
          trans = "log2", breaks = scales::extended_breaks(n = 8)
        )
    }
  }
  ## INLA
  {
    ## Prepare data
    {
      ## Cellmeans
      newdata <- data.frame(Dist.time = factor(c("Before", "After"),
        levels = c("Before", "After")))

      data_pred <- data |>
        dplyr::select(
          n.points, total.points, AIMS_REEF_NAME,
          Site, Transect, Dist.time
        ) |> 
        bind_rows(newdata)
      i_newdata <- (nrow(data) + 1):nrow(data_pred)
      ## Contrasts
      Xmat <- model.matrix(~Dist.time, data = newdata)
      nd <- newdata |>
        bind_cols(Xmat) |>
        summarise(across(where(is.numeric), \(x) x[2] - x[1])) |>
        as.matrix()
      lincomb <- inla.make.lincombs(as.data.frame(nd))
    }
    ## Fit model
    {
      if (rerun_models) {
        mod <- inla(
          n.points ~ Dist.time +
            f(model = "iid", AIMS_REEF_NAME) +
            f(model = "iid", Site) +
            f(model = "iid", Transect),
          data = data_pred,
          Ntrials = data_pred$total.points,
          lincomb = lincomb,
          family = "zeroinflatedbinomial1",
          control.predictor = list(link = 1, compute = TRUE),
          control.compute = list(
            config = TRUE,
            dic = TRUE, waic = TRUE, cpo = TRUE,
            return.marginals.predictor = TRUE
          )
        )
        ## summary(mod)
        ## autoplot(mod)
        save(mod, file = "../data/modelled/mod_2.1.RData")
      }
    }
    ## Diagnostics
    {
      load(file = "../data/modelled/mod_2.1.RData")
      ggplot_inla_residuals(mod,
        observed = data$n.points / data$total.points, CI = TRUE
      )
      ggplot_inla_residuals2(mod,
        observed = data$n.points / data$total.points, CI = TRUE
      )
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
      load(file = "../data/modelled/mod_2.1.RData")
      preds <- posterior_predict.inla(mod, newdata = data)
      mod_resids <- createDHARMa(simulatedResponse = t(preds),
        observedResponse = data$n.points,
        fittedPredictedResponse = apply(preds, 2, mean),
        integerResponse = TRUE)
      mod_resids |> plot()
      mod_resids |> testDispersion()
      mod_resids |> testZeroInflation()
    }
    ## Partial plots - version 1
    {
      load(file = "../data/modelled/mod_2.1.RData")
      cellmeans_summ_inla <- newdata |>
        bind_cols(mod$summary.fitted.values[i_newdata, ])
      cellmeans_summ_inla
      save(cellmeans_summ_inla, file = "../data/modelled/cellmeans_summ_inla_2.1.RData")

      cellmeans_summ_inla |> ggplot(aes(y = `0.5quant`, x = Dist.time)) +
        geom_pointrange(aes(ymin = `0.025quant`, ymax = `0.975quant`))
    }
    ## Partial plots - version 2
    {
      newdata_fitted <- mod |>
        posterior_fitted.inla(newdata)

      newdata_fitted <- newdata_fitted |>
        group_by(Dist.time) |>
        posterior::summarise_draws(median,
          HDInterval::hdi)

      newdata_fitted |>
        ggplot(aes(y = median, x = Dist.time)) +
        geom_pointrange(aes(ymin = lower, ymax = upper))
    }
    ## Partial effects - version 3
    {
      load(file = "../data/modelled/mod_2.1.RData")
      draws <- inla.posterior.sample(n = 1000, result = mod)
      contents <- mod$misc$configs$contents
      cellmeans <- newdata |> cbind(sapply(draws, function(x) x$latent[i_newdata]))
      cellmeans_inla <- cellmeans |>
              pivot_longer(cols = matches("^[0-9]*$"), names_to = ".draws") |>
              posterior::as_draws() |>
              mutate(.draw = .draws) |>
              dplyr::select(-.draws) |>
        mutate(value = plogis(value))
      cellmeans_summ_inla <- cellmeans_inla |> 
        group_by(Dist.time) |> 
        summarise_draws(median, HDInterval::hdi)

      cellmeans_summ_inla |> ggplot(aes(y = median, x = Dist.time)) +
        geom_pointrange(aes(ymin = lower, ymax = upper), position = position_dodge(width = 0.5))

      save(cellmeans_summ_inla, file = "../data/modelled/cellmeans_summ_inla_2.1.RData")
    }

    ## Contrasts - version 1
    {
      nd_pred <- 
        bind_cols(mod$summary.lincomb.derived) |>
        mutate(across(where(is.numeric), exp))
      nd_pred |> ggplot(aes(x = `0.5quant`, y = ID)) +
        geom_vline(xintercept = 1, linetype = "dashed") +
        geom_pointrange(aes(xmin = `0.025quant`, xmax = `0.975quant`)) +
        scale_x_continuous("Effect (Before - After) on a fold scale",
          trans = "log2", breaks = scales::extended_breaks(n = 8)
        )
    }
    ## Contrasts - version 2
    {
      newdata_fitted <- mod |>
        posterior_fitted.inla(newdata)

      newdata_fitted <- newdata_fitted |>
              ungroup() |>
              group_by(.draw) |>
              summarise(Values = exp(diff(log(Values)))) |>
        ungroup() |>
        posterior::summarise_draws(median,
          HDInterval::hdi)

      newdata_fitted |> ggplot(aes(x = median, y = variable)) +
        geom_vline(xintercept = 1, linetype = "dashed") +
        geom_pointrange(aes(xmin = lower, xmax = upper)) +
        scale_x_continuous("Effect (Before - After) on a fold scale",
          trans = "log2", breaks = scales::extended_breaks(n = 8)
        )
    }
    ## Contrasts - version 3
    {
      
      eff_inla <- cellmeans_inla |>
        group_by(.draw) |>
        summarise(value = exp(log(value[Dist.time == "After"]) - log(value[Dist.time == "Before"])))
      eff_summ_inla <- eff_inla |> 
        ungroup() |> 
        summarise_draws(median, HDInterval::hdi)

      eff_summ_inla |> ggplot(aes(x = median, y = 1)) +
        geom_vline(xintercept = 1, linetype = "dashed") +
        geom_pointrange(aes(xmin = lower, xmax = upper), position = position_dodge(width = 0.5))

      save(eff_summ_inla, file = "../data/modelled/eff_summ_inla_2.1.RData")
    }
  }
  ## Comparisons
  {
    
    load(file = "../data/modelled/cellmeans_summ_raw_2.1.RData")
    load(file = "../data/modelled/cellmeans_summ_glmmTMB_2.1.RData")
    load(file = "../data/modelled/cellmeans_summ_brm_2.1.RData")
    load(file = "../data/modelled/cellmeans_summ_inla_2.1.RData")

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

    cellmeans_summ |> ggplot(aes(y = median, x = method, shape = type, colour = Dist.time)) +
            geom_pointrange(aes(ymin = lower, ymax = upper), position = position_dodge(width = 0.5))

    load(file = "../data/modelled/eff_summ_raw_2.1.RData")
    load(file = "../data/modelled/eff_summ_glmmTMB_2.1.RData")
    load(file = "../data/modelled/eff_summ_brm_2.1.RData")
    load(file = "../data/modelled/eff_summ_inla_2.1.RData")

    eff_summ <- bind_rows(
      eff_summ_raw |>
        mutate(method = "raw") |>
        dplyr::rename(median = Mean),
      eff_summ_glmmTMB |>
        mutate(method = "glmmTMB") |>
        dplyr::rename(median = odds.ratio, lower = asymp.LCL, upper = asymp.UCL),
      eff_summ_brm |>
        mutate(method = "brm") |>
        dplyr::rename(median = odds.ratio, lower = lower.HPD, upper = upper.HPD),
      eff_summ_inla |>
        mutate(method = "inla") |> 
        dplyr::rename(median = median, lower = lower, upper = upper)
    )
    eff_summ

    eff_summ |>
      ggplot(aes(x = median, y = method)) +
        geom_vline(xintercept = 1, linetype = "dashed") +
        geom_pointrange(aes(xmin = lower, xmax = upper)) +
        scale_x_continuous("Effect (Before - After) on a fold scale",
          trans = "log2", breaks = scales::extended_breaks(n = 8)
        )
  }
}

