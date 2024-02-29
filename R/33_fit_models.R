rerun_models <- FALSE

## Model 1.1  β₀ + γ₀ᵢ + γ₁ⱼ + γ₂ₖ  (Binomial)
source("33_fit_models_1.1.R")

## Model 1.2  β₀ + γ₀ᵢ + γ₁ⱼ + γ₂ₖ  (Zero-inflated Binomial)
source("33_fit_models_1.2.R")

## data <- data |>
##   dplyr::select(n.points, total.points, A_SECTOR, SHELF, SecShelf,
##     AIMS_REEF_NAME, Site, Transect, Dist.time)


##     ## Prepare data
##     {
##       newdata <- data.frame(Dist.time = factor(c("Before", "After"),
##         levels = c("Before", "After"))) |>
##         crossing(A_SECTOR = data$A_SECTOR, SHELF = data$SHELF) |>
##         mutate(SecShelf = factor(paste(A_SECTOR, SHELF)))

##       data.pred <- data |> bind_rows(newdata)
##       i.newdata <- (nrow(data) + 1):nrow(data.pred)
      
##     }



## Switch to zero-inflated model (RE)
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
    save(cellmeans_summ_raw, file = "../data/modelled/cellmeans_summ_raw_1.3.RData")
  }
  ## glmmTMB
  {
    ## Fit model
    {
    mod_glmmTMB <- glmmTMB(cbind(n.points, total.points-n.points) ~ 1 +
                             (1|AIMS_REEF_NAME) +
                             (1|Site) +
                             (1|Transect),
      ziformula =  ~ 1 + (1|AIMS_REEF_NAME) + (1|Site) + (1|Transect),
      data = data,
      family = "binomial", 
      REML = TRUE
    )

    summary(mod_glmmTMB)
    save(mod_glmmTMB, file = "../data/modelled/mod_glmmTMB_1.3.RData")
    }
    ## Partial effects
    {
      load(file = "../data/modelled/mod_glmmTMB_1.3.RData")
      cellmeans_summ_glmmTMB <- mod_glmmTMB |>
        emmeans(~1, type = "response") |>
        as.data.frame()
      save(cellmeans_summ_glmmTMB, file = "../data/modelled/cellmeans_summ_glmmTMB_1.3.RData")
    }
  }
  ## brms
  {
    ## Fit model
    {
      form <- bf(
        n.points | trials(total.points) ~ 1 +
          (1 | AIMS_REEF_NAME) +
          (1 | Site) +
          (1 | Transect),
        zi = ~ 1 + (1 | AIMS_REEF_NAME) + (1 | Site) + (1 | Transect),
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
        silent =  0#,
        ## refresh = 100
      )
      summary(mod_brm)
      save(mod_brm, file = "../data/modelled/mod_brm_1.3.RData")
    }
    ## Partial effects
    {
      load(file = "../data/modelled/mod_brm_1.3.RData")
      cellmeans_summ_brm <- mod_brm |>
        emmeans(~1, type = "response") |>
        as.data.frame()
      save(cellmeans_summ_brm, file = "../data/modelled/cellmeans_summ_brm_1.3.RData")
    }
  }
  ## INLA
  {
    ## Fit model
    {
      data.pred <- data |>
        dplyr::select(
          n.points, total.points, AIMS_REEF_NAME,
          Site, Transect
        )
      newdata <- data.frame(AIMS_REEF_NAME = NA, Site = NA, Transect = NA)
      data_pred <- data |> bind_rows(newdata)
      i_newdata <- (nrow(data) + 1):nrow(data_pred)
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
      save(mod, file = "../data/modelled/mod_1.3.RData")
    }
    ## Partial effects
    {
      load(file = "../data/modelled/mod_1.3.RData")
      draws <- inla.posterior.sample(n = 1000, result = mod)
      contents <- mod$misc$configs$contents
      i_newdata <- contents$start[contents$tag == "(Intercept)"]
      cellmeans <- newdata |> bind_cols(sapply(draws, function(x) x$latent[i_newdata]))
      cellmeans_summ_inla <- cellmeans |>
        pivot_longer(cols = matches("^[\\.]{3}[0-9]*"), names_to = ".draws") |>
        posterior::as_draws() |>
        mutate(.draw = .draws) |>
        dplyr::select(-.draws, -AIMS_REEF_NAME, -Site, -Transect) |>
        mutate(value = plogis(value)) |>
        summarise_draws(median, HDInterval::hdi)

      cellmeans_summ_inla |> ggplot(aes(y = median, x = 1)) +
        geom_pointrange(aes(ymin = lower, ymax = upper), position = position_dodge(width = 0.5))

      save(cellmeans_summ_inla, file = "../data/modelled/cellmeans_summ_inla_1.3.RData")
    }
  }
  ## Comparisons
  {
    
    load(file = "../data/modelled/cellmeans_summ_raw_1.3.RData")
    load(file = "../data/modelled/cellmeans_summ_glmmTMB_1.3.RData")
    load(file = "../data/modelled/cellmeans_summ_brm_1.3.RData")
    load(file = "../data/modelled/cellmeans_summ_inla_1.3.RData")

    cellmeans_summ <- bind_rows(
      cellmeans_summ_raw |>
        mutate(method = "raw") |>
        dplyr::rename(median = Mean, lower = lower, upper = upper),
      cellmeans_summ_glmmTMB |>
        mutate(method = "glmmTMB") |>
        dplyr::rename(median = prob, lower = asymp.LCL, upper = asymp.UCL),
      cellmeans_summ_brm |> mutate(method = "brm") |>
        dplyr::rename(median = prob, lower = lower.HPD, upper = upper.HPD),
      cellmeans_summ_inla |>
        mutate(method = "inla") |> 
        dplyr::rename(median = median, lower = lower, upper = upper)
    )
    cellmeans_summ

    cellmeans_summ |> ggplot(aes(y = median, x = method)) +
            geom_pointrange(aes(ymin = lower, ymax = upper))

  }
}




## Add Before/After
{
  ## Data preparation
  {
    data <- full_data
  }
  ## Raw means
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
  ## glmmTMB
  {
    ## Fit model
    {
      mod_glmmTMB <- glmmTMB(cbind(n.points, total.points-n.points) ~ Dist.time +
                               (1|AIMS_REEF_NAME) +
                               (1|Site) +
                               (1|Transect),
        ziformula =  ~ 1 + (1|AIMS_REEF_NAME) + (1|Site) + (1|Transect),
        data = data,
        family = "binomial", 
        REML = TRUE
      )

      summary(mod_glmmTMB)
      save(mod_glmmTMB, file = "../data/modelled/mod_glmmTMB_2.1.RData")
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
        silent =  0#,
        ## refresh = 100
      )
      summary(mod_brm)
      save(mod_brm, file = "../data/modelled/mod_brm_2.1.RData")
      load(file = "../data/modelled/mod_brm_2.1.RData")
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
      mod <- inla(n.points ~ Dist.time +
                    f(model = "iid", AIMS_REEF_NAME) +
                    f(model = "iid", Site) +
                    f(model = "iid", Transect),
        data = data_pred,
        Ntrials = data_pred$total.points,
        lincomb =  lincomb,
        family = "zeroinflatedbinomial1", 
        control.predictor = list(link = 1, compute = TRUE),
        control.compute = list(config = TRUE,
          dic = TRUE, waic = TRUE, cpo = TRUE,
          return.marginals.predictor = TRUE
        )
      )
      ## summary(mod)
      ## autoplot(mod)
      save(mod, file = "../data/modelled/mod_2.1.RData")
    }

    ## Diagnostics
    {
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
      cellmeans_summ_inla <- newdata |>
        bind_cols(mod$summary.fitted.values[i_newdata, ])
      newdata_pred 
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
      cellmeans <- newdata |> bind_cols(sapply(draws, function(x) x$latent[i_newdata]))
      cellmeans_inla <- cellmeans |>
              pivot_longer(cols = matches("^[\\.]{3}[0-9]*"), names_to = ".draws") |>
              posterior::as_draws() |>
              mutate(.draw = .draws) |>
              dplyr::select(-.draws) |>
        mutate(value = plogis(value))
      cellmeans_summ_inla <- cellmeans_inla |> 
        group_by(Dist.time) |> 
        summarise_draws(median, HDInterval::hdi)

      cellmeans_summ_inla |> ggplot(aes(y = median, x = 1)) +
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
        mutate(method = "glmmTMB") |>
        dplyr::rename(median = prob, lower = asymp.LCL, upper = asymp.UCL),
      cellmeans_summ_brm |> mutate(method = "brm") |>
        dplyr::rename(median = prob, lower = lower.HPD, upper = upper.HPD),
      cellmeans_summ_inla |>
        mutate(method = "inla") |> 
        dplyr::rename(median = median, lower = lower, upper = upper)
    )
    cellmeans_summ

    cellmeans_summ |> ggplot(aes(y = median, x = method, colour = Dist.time)) +
            geom_pointrange(aes(ymin = lower, ymax = upper))

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

## Before/After and Sector/Shelf
{
  ## Prepare data
  {
      data <- full_data

      data <- data |>
        dplyr::select(n.points, total.points, A_SECTOR, SHELF, SecShelf,
          AIMS_REEF_NAME, Site, Transect, Dist.time)
      ## what combinations are missing
      data |>
              group_by(Dist.time, SecShelf) |>
              summarise(Points = sum(n.points)) |>
              filter(Points == 0)
      ## we need to exclude these
      data <- data |>
              filter(!SecShelf %in% c("CG M", "PC M", "PC O")) |>
              droplevels()
  }
  ## Raw means
  {
    cellmeans_summ_raw <- data |>
      group_by(Dist.time, A_SECTOR, SHELF) |> 
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
    save(cellmeans_summ_raw, file = "../data/modelled/cellmeans_summ_raw_3.1.RData")

    eff_summ_raw <- cellmeans_summ_raw |>
            group_by(type, A_SECTOR, SHELF) |>
            summarise(Mean = Mean[Dist.time == "After"] / Mean[Dist.time == "Before"])
    save(eff_summ_raw, file = "../data/modelled/eff_summ_raw_3.1.RData")
  }
  ## glmmTMB
  {
    ## Fit model
    {
      mod_glmmTMB <- glmmTMB(cbind(n.points, total.points-n.points) ~ Dist.time*SecShelf +
                               (1|AIMS_REEF_NAME) +
                               (1|Site) +
                               (1|Transect),
        ziformula =  ~ 1 + (1|AIMS_REEF_NAME) + (1|Site) + (1|Transect),
        data = data,
        family = "binomial", 
        REML = TRUE
      )

      summary(mod_glmmTMB)
      save(mod_glmmTMB, file = "../data/modelled/mod_glmmTMB_3.1.RData")
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
      load(file = "../data/modelled/mod_glmmTMB_3.1.RData")
      cellmeans_summ_glmmTMB <- mod_glmmTMB |>
        emmeans(~ Dist.time | SecShelf, type = "response") |>
        as.data.frame() |>
        separate(SecShelf, into = c("A_SECTOR", "SHELF"), remove = FALSE) 
      save(cellmeans_summ_glmmTMB, file = "../data/modelled/cellmeans_summ_glmmTMB_3.1.RData")

      ## cellmeans_summ_glmmTMB |>
      ##   ggplot(aes(y = prob, x = Dist.time)) +
      ##   geom_pointrange(aes(ymin = asymp.LCL, ymax = asymp.UCL))
      ## mod_glmmTMB |>
      ##   emmeans(~ Dist.time * SecShelf, type = "response") |>
      ##   as.data.frame() |>
      ##   separate(SecShelf, into = c("A_SECTOR", "SHELF"), remove = FALSE) |> 
      ##   ggplot(aes(y = prob, x = A_SECTOR, colour = Dist.time)) +
      ##   geom_pointrange(aes(ymin = asymp.LCL, ymax = asymp.UCL), position = position_dodge(width = 0.5)) +
      ##   facet_grid(~SHELF)
    }
    ## Contrasts
    {
      load(file = "../data/modelled/mod_glmmTMB_3.1.RData")
      eff_summ_glmmTMB <-
        mod_glmmTMB |>
        emmeans(~Dist.time|SecShelf, type = "response") |>
        contrast(method = list(Dist.time = c(-1, 1))) |>
        summary(infer = TRUE) |>
        as.data.frame() |> 
        separate(SecShelf, into = c("A_SECTOR", "SHELF"), remove = FALSE) 
      save(eff_summ_glmmTMB, file = "../data/modelled/eff_summ_glmmTMB_3.1.RData")

      ## eff_summ_glmmTMB |>
      ##   ggplot(aes(x = odds.ratio, y = contrast)) +
      ##   geom_vline(xintercept = 1, linetype = "dashed") +
      ##   geom_pointrange(aes(xmin = asymp.LCL, xmax = asymp.UCL)) +
      ##   scale_x_continuous("Effect (Before - After) on a fold scale",
      ##     trans = "log2", breaks = scales::extended_breaks(n = 8)
      ##   )
      ## mod_glmmTMB |>
      ##         emmeans(~ Dist.time | SecShelf, type = "response") |>
      ##         contrast(method = list(Dist.time = c(-1, 1))) |>
      ##         summary(infer = TRUE) |>
      ##         as.data.frame() |>
      ##         separate(SecShelf, into = c("A_SECTOR", "SHELF"), remove = FALSE) |>
      ##         ggplot(aes(x = odds.ratio, y = A_SECTOR)) +
      ##         geom_vline(xintercept = 1, linetype = "dashed") +
      ##         geom_pointrange(aes(xmin = asymp.LCL, xmax = asymp.UCL)) +
      ##   scale_x_continuous("Effect (Before - After) on a fold scale", trans = "log2", breaks = scales::extended_breaks(n = 8)) +
      ##   facet_grid(~SHELF)
    }
  }
  ## brms
  {
    ## Fit model
    {
      form <- bf(
        n.points | trials(total.points) ~ Dist.time*SecShelf +
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
        silent =  0#,
        ## refresh = 100
      )
      summary(mod_brm)
      save(mod_brm, file = "../data/modelled/mod_brm_3.1.RData")
    }
    ## Partial plots
    {
      load(file = "../data/modelled/mod_brm_3.1.RData")
      cellmeans_summ_brm <- mod_brm |>
        emmeans(~ Dist.time | SecShelf, type = "response") |>
        as.data.frame() |>
        separate(SecShelf, into = c("A_SECTOR", "SHELF"), remove = FALSE) 
      save(cellmeans_summ_brm, file = "../data/modelled/cellmeans_summ_brm_3.1.RData")
    }
    ## Contrasts
    {
      load(file = "../data/modelled/mod_brm_3.1.RData")
      eff_summ_brm <-
        mod_brm |>
        emmeans(~Dist.time|SecShelf, type = "response") |>
        contrast(method = list(Dist.time = c(-1, 1))) |>
        summary(infer = TRUE) |>
        as.data.frame() |> 
        separate(SecShelf, into = c("A_SECTOR", "SHELF"), remove = FALSE) 
      save(eff_summ_brm, file = "../data/modelled/eff_summ_brm_3.1.RData")
    }
  }
  ## INLA
  {
    ## Prepare data
    {
      data <- full_data

      data <- data |>
        dplyr::select(n.points, total.points, A_SECTOR, SHELF, SecShelf,
          AIMS_REEF_NAME, Site, Transect, Dist.time)
      ## what combinations are missing
      data |>
              group_by(Dist.time, SecShelf) |>
              summarise(Points = sum(n.points)) |>
              filter(Points == 0)
      ## we need to exclude these
      data <- data |>
              filter(!SecShelf %in% c("CG M", "PC M", "PC O")) |>
              droplevels()
      ## Cellmeans
      newdata <- data.frame(Dist.time = factor(c("Before", "After"),
        levels = c("Before", "After"))) |> 
        crossing(SecShelf = data$SecShelf) |>
        separate(SecShelf, into = c("A_SECTOR", "SHELF"), sep = " ", remove = FALSE)
        ## crossing(A_SECTOR = data$A_SECTOR, SHELF = data$SHELF) |>
        ## mutate(SecShelf = factor(paste(A_SECTOR, SHELF)))

      data_pred <- data |> bind_rows(newdata)
      i_newdata <- (nrow(data) + 1):nrow(data_pred)
      ## Contrasts
      Xmat <- model.matrix(~Dist.time * SecShelf, data = newdata)
      nd <- newdata |>
        bind_cols(Xmat) |>
        group_by(A_SECTOR, SHELF, SecShelf) |>
        summarise(across(where(is.numeric), \(x) x[2] - x[1])) |>
        ungroup() |>
        dplyr::select(-A_SECTOR, -SHELF, -SecShelf) |>
        as.matrix()
      lincomb <- inla.make.lincombs(as.data.frame(nd))
      nd1 <- newdata |>
        bind_cols(Xmat) |>
        group_by(A_SECTOR, SHELF, SecShelf) |>
        summarise(across(where(is.numeric), \(x) x[2] - x[1])) |> 
        ungroup() |>
        dplyr::select(A_SECTOR, SHELF, SecShelf)
    }
    ## Fit model
    {
      mod <- inla(n.points ~ Dist.time * SecShelf +
                    f(model = "iid", AIMS_REEF_NAME) +
                    f(model = "iid", Site) +
                    f(model = "iid", Transect),
        data = data_pred,
        Ntrials = data_pred$total.points,
        lincomb =  lincomb,
        family = "zeroinflatedbinomial2", #"binomial", 
        control.predictor = list(link = 1, compute = TRUE),
        control.compute = list(config = TRUE,
          dic = TRUE, waic = TRUE, cpo = TRUE,
          return.marginals.predictor = TRUE
        )
      )
      ## summary(mod)
      ## autoplot(mod)
      save(mod, file = "../data/modelled/mod_3.1.RData")
    }

    ## diagnostics
    {
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
    ## DHARMa version 1
    {
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
      newdata_pred <- newdata |>
        bind_cols(mod$summary.fitted.values[i_newdata, ])
      newdata_pred 

      newdata_pred |>
        ggplot(aes(y = `0.5quant`, x = A_SECTOR, colour = Dist.time)) +
        geom_pointrange(aes(ymin = `0.025quant`, ymax = `0.975quant`),
          position = position_dodge(width = 0.5)) +
        facet_grid(~SHELF)
    }
    ## Partial plots - version 2
    {
      newdata_fitted <- mod |>
        posterior_fitted.inla(newdata)

      newdata_fitted <- newdata_fitted |>
              group_by(Dist.time, A_SECTOR, SHELF) |>
              dplyr::select(-SecShelf) |> 
        posterior::summarise_draws(median,
          HDInterval::hdi)

      newdata_fitted |>
        ggplot(aes(y = median, x = A_SECTOR, colour = Dist.time)) +
        geom_pointrange(aes(ymin = lower, ymax = upper),
          position = position_dodge(width = 0.5)) +
        facet_grid(~SHELF)
    }
    ## Partial effects - version 3
    {
      load(file = "../data/modelled/mod_3.1.RData")
      draws <- inla.posterior.sample(n = 1000, result = mod)
      contents <- mod$misc$configs$contents
      cellmeans <- newdata |> bind_cols(sapply(draws, function(x) x$latent[i_newdata]))
      cellmeans_inla <- cellmeans |>
              pivot_longer(cols = matches("^[\\.]{3}[0-9]*"), names_to = ".draws") |>
              posterior::as_draws() |>
              mutate(.draw = .draws) |>
              dplyr::select(-.draws, -SecShelf) |>
        mutate(value = plogis(value))
      cellmeans_summ_inla <- cellmeans_inla |> 
        group_by(Dist.time, A_SECTOR, SHELF) |> 
        summarise_draws(median, HDInterval::hdi)

      save(cellmeans_summ_inla, file = "../data/modelled/cellmeans_summ_inla_3.1.RData")
    }
    ## Contrasts
    {
      nd.pred <- nd1 |>
        bind_cols(mod$summary.lincomb.derived) |>
        mutate(across(where(is.numeric), exp))
      nd.pred |> ggplot(aes(x = `0.5quant`, y = A_SECTOR)) +
              geom_vline(xintercept = 1, linetype = "dashed") +
              geom_pointrange(aes(xmin = `0.025quant`, xmax = `0.975quant`)) +
        scale_x_continuous("Effect (Before - After) on a fold scale",
          trans = "log2", breaks = scales::extended_breaks(n = 8)) +
        facet_grid(~SHELF)
    }
    ## Contrasts - version 2
    {
      newdata_fitted <- mod |>
        posterior_fitted.inla(newdata)

      newdata_fitted <- newdata_fitted |>
        ungroup() |>
        dplyr::select(-Dist.time, -SecShelf) |>
        group_by(.draw, A_SECTOR, SHELF) |>
        summarise(Values = exp(diff(log(Values)))) |>
        ungroup() |>
        group_by(A_SECTOR, SHELF) |>
        posterior::summarise_draws(
          median,
          HDInterval::hdi,
          Pl = ~ mean(.x < 1),
          Pg = ~ mean(.x > 1)
        )

      newdata_fitted |> ggplot(aes(x = median, y = A_SECTOR)) +
        geom_vline(xintercept = 1, linetype = "dashed") +
        geom_pointrange(aes(xmin = lower, xmax = upper)) +
        scale_x_continuous("Effect (Before - After) on a fold scale",
          trans = "log2", breaks = scales::extended_breaks(n = 8)
        ) +
        facet_grid(~SHELF)

      newdata_fitted <- mod |>
        posterior_fitted.inla(newdata)

      newdata_fitted <- newdata_fitted |>
              ungroup() |>
              dplyr::select(-Dist.time, -SecShelf) |>
              group_by(.draw, A_SECTOR, SHELF) |>
        summarise(Values = exp(diff(log(Values))))

      newdata_fitted |>
        ungroup() |>
        dplyr::select(-SHELF) |>
        group_by(A_SECTOR) |>
        posterior::summarise_draws(
          median,
          HDInterval::hdi,
          Pl = ~ mean(.x < 1),
          Pg = ~ mean(.x > 1)
        )

      newdata_fitted |>
        ungroup() |>
        dplyr::select(-SHELF, -A_SECTOR) |>
        posterior::summarise_draws(
          median,
          HDInterval::hdi,
          Pl = ~ mean(.x < 1),
          Pg = ~ mean(.x > 1)
        )
    }
    ## Contrasts - version 3
    {
      eff_inla <- cellmeans_inla |>
        group_by(.draw, A_SECTOR, SHELF) |>
        summarise(value = exp(log(value[Dist.time == "After"]) - log(value[Dist.time == "Before"])))
      eff_summ_inla <- eff_inla |> 
        ungroup() |> 
        group_by(A_SECTOR, SHELF) |>
        summarise_draws(median, HDInterval::hdi)

      save(eff_summ_inla, file = "../data/modelled/eff_summ_inla_3.1.RData")
    }
  }
  ## Comparisons
  {
    
    load(file = "../data/modelled/cellmeans_summ_raw_3.1.RData")
    load(file = "../data/modelled/cellmeans_summ_glmmTMB_3.1.RData")
    load(file = "../data/modelled/cellmeans_summ_brm_3.1.RData")
    load(file = "../data/modelled/cellmeans_summ_inla_3.1.RData")

    cellmeans_summ <- bind_rows(
      cellmeans_summ_raw |>
        mutate(method = "raw") |>
        dplyr::rename(median = Mean, lower = lower, upper = upper),
      cellmeans_summ_glmmTMB |>
        mutate(method = "glmmTMB", type = "median") |>
        dplyr::rename(median = prob, lower = asymp.LCL, upper = asymp.UCL),
      cellmeans_summ_brm |> mutate(method = "brm", type = "median") |>
        dplyr::rename(median = prob, lower = lower.HPD, upper = upper.HPD),
      cellmeans_summ_inla |>
        mutate(method = "inla", type = "median") |> 
        dplyr::rename(median = median, lower = lower, upper = upper)
    )
    cellmeans_summ

    cellmeans_summ |> ggplot(aes(y = median, x = method, shape = type, colour = Dist.time)) +
      geom_pointrange(aes(ymin = lower, ymax = upper),
        position = position_dodge(width = 0.5)) +
      facet_wrap(~A_SECTOR + SHELF, scales = "free")
      ## facet_grid(A_SECTOR ~ SHELF, scales = "free")

    load(file = "../data/modelled/eff_summ_raw_3.1.RData")
    load(file = "../data/modelled/eff_summ_glmmTMB_3.1.RData")
    load(file = "../data/modelled/eff_summ_brm_3.1.RData")
    load(file = "../data/modelled/eff_summ_inla_3.1.RData")

    eff_summ <- bind_rows(
      eff_summ_raw |>
        mutate(method = "raw") |>
        dplyr::rename(median = Mean),
      eff_summ_glmmTMB |>
        mutate(method = "glmmTMB", type = "median") |>
        dplyr::rename(median = odds.ratio, lower = asymp.LCL, upper = asymp.UCL),
      eff_summ_brm |>
        mutate(method = "brm", type = "median") |>
        dplyr::rename(median = odds.ratio, lower = lower.HPD, upper = upper.HPD),
      eff_summ_inla |>
        mutate(method = "inla", type = "median") |> 
        dplyr::rename(median = median, lower = lower, upper = upper)
    )
    eff_summ

    eff_summ |>
      ggplot(aes(x = median, y = method, shape = type)) +
        geom_vline(xintercept = 1, linetype = "dashed") +
        geom_pointrange(aes(xmin = lower, xmax = upper)) +
        scale_x_continuous("Effect (Before - After) on a fold scale",
          trans = "log2", breaks = scales::extended_breaks(n = 8)
        ) +
      facet_grid(A_SECTOR ~ SHELF, scales = "free")
  }
}

## Disturbances (zi intercept only)
{
  ## Prepare data
  {
    data <- full_data

    ## Focus on only the necessary variables
    data <- data |>
      dplyr::select(
        n.points, total.points, Dist.time, s, c, d, b, u,
        AIMS_REEF_NAME, Site, Transect, Dist.number
      ) |>
      mutate(across(c(s, c, d, b, u), \(x) ifelse(Dist.time == "Before", 0, x)))
  }
  ## Raw means
  {
    cellmeans_summ_raw <- data |>
        filter((Dist.time == "Before" & s == 0 & d == 0 & c == 0 & b == 0 & u == 0) |
                 (Dist.time == "After" & s == 1 & d == 0 & c == 0 & b == 0 & u == 0) |
                 (Dist.time == "After" & s == 0 & d == 1 & c == 0 & b == 0 & u == 0) |
                 (Dist.time == "After" & s == 0 & d == 0 & c == 1 & b == 0 & u == 0) |
                 (Dist.time == "After" & s == 0 & d == 0 & c == 0 & b == 1 & u == 0) |
                 (Dist.time == "After" & s == 0 & d == 0 & c == 0 & b == 0 & u == 1) 
        ) |> 
        mutate(Dist = case_when(s == 1 ~ "s",
          c == 1 ~ "c",
          d == 1 ~ "d",
          b == 1 ~ "b",
          u == 1 ~ "u",
          .default = "Before")) |>
        dplyr::select(-any_of(c("s", "c", "d", "b", "u"))) |> 
      group_by(Dist.time, Dist) |> 
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
    save(cellmeans_summ_raw, file = "../data/modelled/cellmeans_summ_raw_4.1.RData")

    eff_summ_raw <- cellmeans_summ_raw |>
            mutate(Values = Mean) |>
            nest(.by = type) |>
            mutate(eff = map(
                    .x = data,
                    .f = ~ before_vs_afters(.x)
            )) |>
      dplyr::select(-data) |> 
      unnest(c(eff))
    save(eff_summ_raw, file = "../data/modelled/eff_summ_raw_4.1.RData")
  }
  ## glmmTMB
  {
    ## Fit model
    {
      mod_glmmTMB <- glmmTMB(
        cbind(n.points, total.points - n.points) ~ Dist.time + (s + c + d + b + u) +
          (1 | AIMS_REEF_NAME) +
          (1 | Site) +
          (1 | Transect),
        ziformula = ~ 1,
        data = data,
        family = "binomial",
        REML = TRUE
      )

      summary(mod_glmmTMB)
      save(mod_glmmTMB, file = "../data/modelled/mod_glmmTMB_4.1.RData")
    }
    ## DHARMa
    {
      resids <- simulateResiduals(mod_glmmTMB, plot = TRUE)
      wrap_elements(~ testUniformity(resids)) +
        wrap_elements(~ plotResiduals(resids)) +
        wrap_elements(~ testDispersion(resids))
      testDispersion(resids)
      testZeroInflation(resids)
    }
    ## Partial plots
    {
      load(file = "../data/modelled/mod_glmmTMB_4.1.RData")

      ## newdata <- crossing(Dist.time = data$Dist.time) |>
      ##   crossing(u = c(0, 1), b = c(0, 1), d = c(0, 1), c = c(0, 1), s = c(0,1)) |>
      ##   mutate(AIMS_REEF_NAME = NA, Site = NA, Transect = NA, Dist.number =  NA)
      ## newdata_1 <- newdata |> mutate(D = paste(s, c, d, b, u, Dist.time))
      ## newdata_1 <- newdata |>
      ##   mutate(
      ##   ## D = paste(Dist.time, u, b, d, c, s),
      ##   D = paste(s, c, d, b, u, Dist.time),
      ##   D = factor(D, levels = D)
      ## )
      ## Xmat <- model.matrix(~ -1 + D, data = newdata_1)
      ## Xmat[1, ]
      ## newdata_1
      ## colnames(Xmat)
      ## Before <- str_which(colnames(Xmat), "0 0 0 0 0 Before")
      ## S <- str_which(colnames(Xmat), "1 . . . . After")
      ## Xmat_s <- colMeans(Xmat[S, ])
      ## S <- str_which(colnames(Xmat), "1 0 0 0 0 After")
      ## Xmat_s <- Xmat[S, ]
      ## C <- str_which(colnames(Xmat), ". 1 . . . After")
      ## Xmat_c <- colMeans(Xmat[C, ])
      ## C <- str_which(colnames(Xmat), "0 1 0 0 0 After")
      ## Xmat_c <- Xmat[C, ]
      ## D <- str_which(colnames(Xmat), ". . 1 . . After")
      ## Xmat_d <- colMeans(Xmat[D, ])
      ## D <- str_which(colnames(Xmat), "0 0 1 0 0 After")
      ## Xmat_d <- Xmat[D, ]
      ## B <- str_which(colnames(Xmat), ". . . 1 . After")
      ## Xmat_b <- colMeans(Xmat[B, ])
      ## B <- str_which(colnames(Xmat), "0 0 0 1 0 After")
      ## Xmat_b <- Xmat[B, ]
      ## U <- str_which(colnames(Xmat), ". . . . 1 After")
      ## Xmat_u <- colMeans(Xmat[U, ])
      ## U <- str_which(colnames(Xmat), "0 0 0 0 1 After")
      ## Xmat_u <- Xmat[U, ]
      ## Xmat_all <- rbind(
      ##         Xmat[Before, ],
      ##         Xmat_s, Xmat_c, Xmat_d, Xmat_b, Xmat_u
      ## )
      
      ## mod_glmmTMB |>
      ##         emmeans(~ s + c + d + b + u + Dist.time, type = "response") 
      ## mod_glmmTMB |>
      ##         emmeans(~ s + c + d + b + u + Dist.time, type = "response") |>
      ##         contrast(method = list(t(Xmat_all))) |> 
      ##         ## contrast(method = list(Xmat[1, ]))
      ##         ## contrast(method = list(c(1, rep(0, 63)))) |>
      ##         summary(infer = TRUE) |>
      ##         as.data.frame() |>
      ##         mutate(across(c(estimate, asymp.LCL, asymp.UCL), plogis))
      
      newdata <- crossing(Dist.time = data$Dist.time) |>
        crossing(s = c(0,1), d = c(0, 1), c = c(0, 1), b = c(0, 1), u = c(0, 1)) |>
        mutate(AIMS_REEF_NAME = NA, Site = NA, Transect = NA, Dist.number =  NA)
      newdata <- newdata |>
        filter((Dist.time == "Before" & s == 0 & d == 0 & c == 0 & b == 0 & u == 0) |
                 (Dist.time == "After" & s == 1 & d == 0 & c == 0 & b == 0 & u == 0) |
                 (Dist.time == "After" & s == 0 & d == 1 & c == 0 & b == 0 & u == 0) |
                 (Dist.time == "After" & s == 0 & d == 0 & c == 1 & b == 0 & u == 0) |
                 (Dist.time == "After" & s == 0 & d == 0 & c == 0 & b == 1 & u == 0) |
                 (Dist.time == "After" & s == 0 & d == 0 & c == 0 & b == 0 & u == 1) 
        )

      p <- predict(mod_glmmTMB, newdata = newdata, se.fit = TRUE)
      newdata <- newdata |>
        bind_cols(fit = p$fit, se = p$se.fit) |>
        mutate(Pred = plogis(fit), lower = plogis(fit - 2 * se), upper = plogis(fit + 2 * se))
      newdata <- newdata |>
        mutate(Dist = case_when(
          s == 1 ~ "s",
          c == 1 ~ "c",
          d == 1 ~ "d",
          b == 1 ~ "b",
          u == 1 ~ "u",
          .default = "Before"
        )) 

      cellmeans_summ_glmmTMB <- newdata |> dplyr::select(Dist.time, Dist, Pred, lower, upper)
      save(cellmeans_summ_glmmTMB, file = "../data/modelled/cellmeans_summ_glmmTMB_4.1.RData")
      
      ## newdata |>
      ##   ggplot(aes(y = Pred, x = Dist, colour = factor(Dist.time))) +
      ##   geom_pointrange(aes(ymin = lower, ymax = upper),
      ##     position = position_dodge(width = 0.5)) 
    }
    ## Contrasts
    {
      eff_summ_glmmTMB <- 
        cellmeans_summ_glmmTMB |>
        mutate(Values = Pred) |> 
        before_vs_afters()
      save(eff_summ_glmmTMB, file = "../data/modelled/eff_summ_glmmTMB_4.1.RData")
    }
  }
  ## brms
  {
    ## Fit model
    {
      form <- bf(
        n.points | trials(total.points) ~ Dist.time + (s + c + d + b + u) + 
          (1 | AIMS_REEF_NAME) +
          (1 | Site) +
          (1 | Transect),
        zi = ~ 1 ,
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
        silent =  0#,
        ## refresh = 100
      )
      summary(mod_brm)
      save(mod_brm, file = "../data/modelled/mod_brm_4.1.RData")
    }
    ## Partial plot
    {
      load(file = "../data/modelled/mod_brm_4.1.RData")
      newdata <- crossing(Dist.time = data$Dist.time) |>
        crossing(s = c(0,1), d = c(0, 1), c = c(0, 1), b = c(0, 1), u = c(0, 1)) |>
        mutate(AIMS_REEF_NAME = NA, Site = NA, Transect = NA, Dist.number =  NA)
      newdata <- newdata |>
        filter((Dist.time == "Before" & s == 0 & d == 0 & c == 0 & b == 0 & u == 0) |
                 (Dist.time == "After" & s == 1 & d == 0 & c == 0 & b == 0 & u == 0) |
                 (Dist.time == "After" & s == 0 & d == 1 & c == 0 & b == 0 & u == 0) |
                 (Dist.time == "After" & s == 0 & d == 0 & c == 1 & b == 0 & u == 0) |
                 (Dist.time == "After" & s == 0 & d == 0 & c == 0 & b == 1 & u == 0) |
                 (Dist.time == "After" & s == 0 & d == 0 & c == 0 & b == 0 & u == 1)) |>
        mutate(total.points = 1)

      p <- t(brms::posterior_epred(mod_brm, newdata = newdata, se.fit = TRUE))
      cellmeans_brm <-
        newdata |>
        cbind(p) |>
        pivot_longer(cols = matches("^[0-9]*$"), names_to = ".draws") |>
        mutate(Dist = case_when(
          s == 1 ~ "s",
          c == 1 ~ "c",
          d == 1 ~ "d",
          b == 1 ~ "b",
          u == 1 ~ "u",
          .default = "Before"
        )) |>
        dplyr::select(-any_of(c(
          "s", "d", "c", "b", "u", "AIMS_REEF_NAME",
          "Dist.time", "Site", "Transect", "Dist.number", "total.points"
        )))

      cellmeans_summ_brm <- cellmeans_brm |>
              dplyr::select(-.draws) |> 
        group_by(Dist) |>
        summarise_draws(median, HDInterval::hdi) |> 
        dplyr::select(Dist.time, Dist, Pred, lower, upper)
      save(cellmeans_summ_brm, file = "../data/modelled/cellmeans_summ_brm_4.1.RData")
      
    }
    ## Contrasts
    {
      eff_brm <- cellmeans_brm |>
        mutate(Values = value) |>
        nest(.by = .draws) |>
        mutate(eff = map(
          .x = data,
          .f = ~ before_vs_afters(.x)
        )) |>
        dplyr::select(-data) |>
        unnest(c(eff))

      eff_summ_brm <- eff_brm |> 
        dplyr::select(-.draws) |>
        group_by(Dist) |> 
        summarise_draws(
          median,
          HDInterval::hdi
        )
      
      save(eff_summ_brm, file = "../data/modelled/eff_summ_brm_4.1.RData")
    }
  }
  ## INLA
  {
    ## Prepare data
    {
      data <- full_data

      data <- data |>
        dplyr::select(n.points, total.points, s, c, d, b, u,
          AIMS_REEF_NAME, Site, Transect, Dist.time)
      ## data <- data |>
      ##         filter(!SecShelf %in% c("CG M", "PC M", "PC O")) |>
      ##         droplevels()
      ## Cellmeans
      newdata <- crossing(Dist.time = data$Dist.time) |>
        crossing(s = c(0,1), d = c(0, 1), c = c(0, 1), b = c(0, 1), u = c(0, 1)) |>
        mutate(AIMS_REEF_NAME = NA, Site = NA, Transect = NA, Dist.number =  NA)
      newdata <- newdata |>
        filter((Dist.time == "Before" & s == 0 & d == 0 & c == 0 & b == 0 & u == 0) |
                 (Dist.time == "After" & s == 1 & d == 0 & c == 0 & b == 0 & u == 0) |
                 (Dist.time == "After" & s == 0 & d == 1 & c == 0 & b == 0 & u == 0) |
                 (Dist.time == "After" & s == 0 & d == 0 & c == 1 & b == 0 & u == 0) |
                 (Dist.time == "After" & s == 0 & d == 0 & c == 0 & b == 1 & u == 0) |
                 (Dist.time == "After" & s == 0 & d == 0 & c == 0 & b == 0 & u == 1)) |>
        mutate(total.points = 1)

      ## newdata <- data.frame(Dist.time = factor(c("Before", "After"),
      ##   levels = c("Before", "After"))) |> 
      ##   crossing(s = c(0,1), d = c(0, 1), c = c(0, 1), b = c(0, 1), u = c(0, 1)) |>
      ##   mutate(AIMS_REEF_NAME = NA, Site = NA, Transect = NA, Dist.number =  NA)

      data_pred <- data |> bind_rows(newdata)
      i_newdata <- (nrow(data) + 1):nrow(data_pred)
    }
    ## Fit model
    {
      mod <- inla(n.points ~ Dist.time + (s + c + d + b + u) +
                    f(model = "iid", AIMS_REEF_NAME) +
                    f(model = "iid", Site) +
                    f(model = "iid", Transect),
        data = data_pred,
        Ntrials = data_pred$total.points,
        family = "zeroinflatedbinomial2", #"binomial", 
        control.predictor = list(link = 1, compute = TRUE),
        control.compute = list(config = TRUE,
          dic = TRUE, waic = TRUE, cpo = TRUE,
          return.marginals.predictor = TRUE
        )
      )
      summary(mod)
      ## autoplot(mod)
      save(mod, file = "../data/modelled/mod_4.1.RData")
    }

    ## Partial plots - version 1
    {
      newdata_pred <- newdata |>
        bind_cols(mod$summary.fitted.values[i_newdata, ])
      newdata_pred |> 
        mutate(Dist = case_when(
          s == 1 ~ "s",
          c == 1 ~ "c",
          d == 1 ~ "d",
          b == 1 ~ "b",
          u == 1 ~ "u",
          .default = "Before"
        )) |>
        dplyr::select(-any_of(c(
          "s", "d", "c", "b", "u", "AIMS_REEF_NAME",
          "Dist.time", "Site", "Transect", "Dist.number", "total.points"
        )))  
    }
    ## Partial plots - version 2
    {
      ## newdata_fitted <- mod |>
      ##   posterior_fitted.inla(newdata)

      ## newdata_fitted <- newdata_fitted |>
      ##         group_by(Dist.time, A_SECTOR, SHELF) |>
      ##         dplyr::select(-SecShelf) |> 
      ##   posterior::summarise_draws(median,
      ##     HDInterval::hdi)
    }
    ## Partial effects - version 3
    {
      load(file = "../data/modelled/mod_4.1.RData")
      draws <- inla.posterior.sample(n = 1000, result = mod)
      contents <- mod$misc$configs$contents
      ## cellmeans <- newdata |> bind_cols(sapply(draws, function(x) x$latent[i_newdata]))
      cellmeans <- newdata |> cbind(sapply(draws, function(x) x$latent[i_newdata]))
      cellmeans_inla <- cellmeans |>
              ## pivot_longer(cols = matches("^[\\.]{3}[0-9]*"), names_to = ".draws") |>
              pivot_longer(cols = matches("^[0-9]*$"), names_to = ".draws") |>
              posterior::as_draws() |>
              mutate(.draw = .draws) |>
        mutate(Dist = case_when(
          s == 1 ~ "s",
          c == 1 ~ "c",
          d == 1 ~ "d",
          b == 1 ~ "b",
          u == 1 ~ "u",
          .default = "Before"
        )) |>
        dplyr::select(-any_of(c(
          "s", "d", "c", "b", "u", "AIMS_REEF_NAME",
          "Dist.time", "Site", "Transect", "Dist.number", "total.points"
        ))) |> 
        mutate(value = plogis(value))

      cellmeans_summ_inla <- cellmeans_inla |> 
      dplyr::select(-.draws) |>
        group_by(Dist) |> 
        summarise_draws(median, HDInterval::hdi)

      save(cellmeans_summ_inla, file = "../data/modelled/cellmeans_summ_inla_4.1.RData")
    }
    ## Contrasts
    {
      nd.pred <- nd1 |>
        bind_cols(mod$summary.lincomb.derived) |>
        mutate(across(where(is.numeric), exp))
      nd.pred |> ggplot(aes(x = `0.5quant`, y = A_SECTOR)) +
              geom_vline(xintercept = 1, linetype = "dashed") +
              geom_pointrange(aes(xmin = `0.025quant`, xmax = `0.975quant`)) +
        scale_x_continuous("Effect (Before - After) on a fold scale",
          trans = "log2", breaks = scales::extended_breaks(n = 8)) +
        facet_grid(~SHELF)
    }
    ## Contrasts - version 2
    {
      newdata_fitted <- mod |>
        posterior_fitted.inla(newdata)

      newdata_fitted <- newdata_fitted |>
        ungroup() |>
        dplyr::select(-Dist.time, -SecShelf) |>
        group_by(.draw, A_SECTOR, SHELF) |>
        summarise(Values = exp(diff(log(Values)))) |>
        ungroup() |>
        group_by(A_SECTOR, SHELF) |>
        posterior::summarise_draws(
          median,
          HDInterval::hdi,
          Pl = ~ mean(.x < 1),
          Pg = ~ mean(.x > 1)
        )

      newdata_fitted |> ggplot(aes(x = median, y = A_SECTOR)) +
        geom_vline(xintercept = 1, linetype = "dashed") +
        geom_pointrange(aes(xmin = lower, xmax = upper)) +
        scale_x_continuous("Effect (Before - After) on a fold scale",
          trans = "log2", breaks = scales::extended_breaks(n = 8)
        ) +
        facet_grid(~SHELF)

      newdata_fitted <- mod |>
        posterior_fitted.inla(newdata)

      newdata_fitted <- newdata_fitted |>
              ungroup() |>
              dplyr::select(-Dist.time, -SecShelf) |>
              group_by(.draw, A_SECTOR, SHELF) |>
        summarise(Values = exp(diff(log(Values))))

      newdata_fitted |>
        ungroup() |>
        dplyr::select(-SHELF) |>
        group_by(A_SECTOR) |>
        posterior::summarise_draws(
          median,
          HDInterval::hdi,
          Pl = ~ mean(.x < 1),
          Pg = ~ mean(.x > 1)
        )

      newdata_fitted |>
        ungroup() |>
        dplyr::select(-SHELF, -A_SECTOR) |>
        posterior::summarise_draws(
          median,
          HDInterval::hdi,
          Pl = ~ mean(.x < 1),
          Pg = ~ mean(.x > 1)
        )
    }
    ## Contrasts - version 3
    {

      eff_inla <- cellmeans_inla |>
        mutate(Values = value) |>
        nest(.by = .draws) |>
        mutate(eff = map(
          .x = data,
          .f = ~ before_vs_afters(.x)
        )) |>
        dplyr::select(-data) |>
        unnest(c(eff))

      eff_summ_inla <- eff_inla |> 
        dplyr::select(-.draws) |>
        group_by(Dist) |> 
        summarise_draws(
          median,
          HDInterval::hdi
        )

      ## eff_inla <- cellmeans_inla |>
      ##   group_by(.draw, A_SECTOR, SHELF) |>
      ##   summarise(value = exp(log(value[Dist.time == "After"]) - log(value[Dist.time == "Before"])))
      ## eff_summ_inla <- eff_inla |> 
      ##   ungroup() |> 
      ##   group_by(A_SECTOR, SHELF) |>
      ##   summarise_draws(median, HDInterval::hdi)

      save(eff_summ_inla, file = "../data/modelled/eff_summ_inla_4.1.RData")
    }
  }
  ## Comparisons
  {
    
    load(file = "../data/modelled/cellmeans_summ_raw_4.1.RData")
    load(file = "../data/modelled/cellmeans_summ_glmmTMB_4.1.RData")
    load(file = "../data/modelled/cellmeans_summ_brm_4.1.RData")
    load(file = "../data/modelled/cellmeans_summ_inla_4.1.RData")

    cellmeans_summ <- bind_rows(
      cellmeans_summ_raw |>
        mutate(method = "raw") |>
        dplyr::rename(median = Mean, lower = lower, upper = upper),
      cellmeans_summ_glmmTMB |>
        mutate(method = "glmmTMB", type = "median") |>
        dplyr::rename(median = Pred),
      cellmeans_summ_brm |> mutate(method = "brm", type = "median") |>
        dplyr::rename(median = median, lower = lower, upper = upper),
      cellmeans_summ_inla |>
        mutate(method = "inla", type = "median") |> 
        dplyr::rename(median = median, lower = lower, upper = upper)
    )
    cellmeans_summ

    cellmeans_summ |>
      ggplot(aes(y = median, x = method, shape = type, colour = Dist)) +
      geom_pointrange(aes(ymin = lower, ymax = upper),
        position = position_dodge(width = 0.5)) 

    load(file = "../data/modelled/eff_summ_raw_4.1.RData")
    load(file = "../data/modelled/eff_summ_glmmTMB_4.1.RData")
    load(file = "../data/modelled/eff_summ_brm_4.1.RData")
    load(file = "../data/modelled/eff_summ_inla_4.1.RData")

    eff_summ <- bind_rows(
      eff_summ_raw |>
        mutate(method = "raw") |>
        dplyr::rename(median = Values),
      eff_summ_glmmTMB |>
        mutate(method = "glmmTMB", type = "median") |>
        dplyr::rename(median = Values),
      eff_summ_brm |>
        mutate(method = "brm", type = "median") |>
        dplyr::rename(median = median, lower = lower, upper = upper),
      eff_summ_inla |>
        mutate(method = "inla", type = "median") |> 
        dplyr::rename(median = median, lower = lower, upper = upper)
    )
    eff_summ

    eff_summ |>
            ggplot(aes(x = median, y = Dist, colour = method, shape = type)) +
            geom_vline(xintercept = 1, linetype = "dashed") +
      geom_pointrange(aes(xmin = lower, xmax = upper),
        position = position_dodge(width = 0.5)) +
            scale_x_continuous("Effect (Before - After) on a fold scale",
                    trans = "log2", breaks = scales::extended_breaks(n = 8)
            )
  }
}


## Disturbances (full zi RE)
{
  ## Prepare data
  {
    data <- full_data

    ## Focus on only the necessary variables
    data <- data |>
      dplyr::select(
        n.points, total.points, Dist.time, s, c, d, b, u,
        AIMS_REEF_NAME, Site, Transect, Dist.number
      ) |>
      mutate(across(c(s, c, d, b, u), \(x) ifelse(Dist.time == "Before", 0, x))) 
  }
  ## Raw means
  {
    cellmeans_summ_raw <- data |>
        filter((Dist.time == "Before" & s == 0 & d == 0 & c == 0 & b == 0 & u == 0) |
                 (Dist.time == "After" & s == 1 & d == 0 & c == 0 & b == 0 & u == 0) |
                 (Dist.time == "After" & s == 0 & d == 1 & c == 0 & b == 0 & u == 0) |
                 (Dist.time == "After" & s == 0 & d == 0 & c == 1 & b == 0 & u == 0) |
                 (Dist.time == "After" & s == 0 & d == 0 & c == 0 & b == 1 & u == 0) |
                 (Dist.time == "After" & s == 0 & d == 0 & c == 0 & b == 0 & u == 1) 
        ) |> 
        mutate(Dist = case_when(s == 1 ~ "s",
          c == 1 ~ "c",
          d == 1 ~ "d",
          b == 1 ~ "b",
          u == 1 ~ "u",
          .default = "Before")) |>
        dplyr::select(-any_of(c("s", "c", "d", "b", "u"))) |> 
      group_by(Dist.time, Dist) |> 
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
    save(cellmeans_summ_raw, file = "../data/modelled/cellmeans_summ_raw_5.1.RData")

    eff_summ_raw <- cellmeans_summ_raw |>
            mutate(Values = Mean) |>
            nest(.by = type) |>
            mutate(eff = map(
                    .x = data,
                    .f = ~ before_vs_afters(.x)
            )) |>
      dplyr::select(-data) |> 
      unnest(c(eff))
    save(eff_summ_raw, file = "../data/modelled/eff_summ_raw_5.1.RData")
  }
  ## glmmTMB
  {
    ## Fit model
    {
      library(glmmTMB)
      mod_glmmTMB <- glmmTMB(cbind(n.points, total.points-n.points) ~ Dist.time*(s+c+d+b+u) +
                               (1|AIMS_REEF_NAME) +
                               (1|Site) +
                               (1|Transect),
        ziformula = ~ 1 + (1|AIMS_REEF_NAME) + (1|Site) + (1|Transect),
        data = data,
        family = "binomial", 
        REML = TRUE
      )

      summary(mod_glmmTMB)
      save(mod_glmmTMB, file = "../data/modelled/mod_glmmTMB_5.1.RData")
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
      load(file = "../data/modelled/mod_glmmTMB_5.1.RData")
      newdata <- crossing(Dist.time = data$Dist.time) |>
        crossing(s = c(0,1), d = c(0, 1), c = c(0, 1), b = c(0, 1), u = c(0, 1)) |>
        mutate(AIMS_REEF_NAME = NA, Site = NA, Transect = NA, Dist.number =  NA)
      newdata <- newdata |>
        filter((Dist.time == "Before" & s == 0 & d == 0 & c == 0 & b == 0 & u == 0) |
                 (Dist.time == "After" & s == 1 & d == 0 & c == 0 & b == 0 & u == 0) |
                 (Dist.time == "After" & s == 0 & d == 1 & c == 0 & b == 0 & u == 0) |
                 (Dist.time == "After" & s == 0 & d == 0 & c == 1 & b == 0 & u == 0) |
                 (Dist.time == "After" & s == 0 & d == 0 & c == 0 & b == 1 & u == 0) |
                 (Dist.time == "After" & s == 0 & d == 0 & c == 0 & b == 0 & u == 1) 
        )

      p <- predict(mod_glmmTMB, newdata = newdata, se.fit = TRUE)
      newdata <- newdata |>
        bind_cols(fit = p$fit, se = p$se.fit) |>
        mutate(Pred = plogis(fit), lower = plogis(fit - 2 * se), upper = plogis(fit + 2 * se))
      newdata <- newdata |>
        mutate(Dist = case_when(
          s == 1 ~ "s",
          c == 1 ~ "c",
          d == 1 ~ "d",
          b == 1 ~ "b",
          u == 1 ~ "u",
          .default = "Before"
        )) 

      cellmeans_summ_glmmTMB <- newdata |> dplyr::select(Dist.time, Dist, Pred, lower, upper)
      save(cellmeans_summ_glmmTMB, file = "../data/modelled/cellmeans_summ_glmmTMB_5.1.RData")
    }
    ## Contrasts
    {
      eff_summ_glmmTMB <- 
        cellmeans_summ_glmmTMB |>
        mutate(Values = Pred) |> 
        before_vs_afters()
      save(eff_summ_glmmTMB, file = "../data/modelled/eff_summ_glmmTMB_5.1.RData")
    }
  }
  ## brms
  {
    ## Fit model
    {
      form <- bf(
        n.points | trials(total.points) ~ Dist.time*(s+c+d+b+u) +
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
        silent =  0#,
        ## refresh = 100
      )
      summary(mod_brm)
      save(mod_brm, file = "../data/modelled/mod_brm_5.1.RData")
      load(file = "../data/modelled/mod_brm_5.1.RData")
    }
    ## Partial plot
    {
      load(file = "../data/modelled/mod_brm_5.1.RData")
      newdata <- crossing(Dist.time = data$Dist.time) |>
        crossing(s = c(0,1), d = c(0, 1), c = c(0, 1), b = c(0, 1), u = c(0, 1)) |>
        mutate(AIMS_REEF_NAME = NA, Site = NA, Transect = NA, Dist.number =  NA)
      newdata <- newdata |>
        filter((Dist.time == "Before" & s == 0 & d == 0 & c == 0 & b == 0 & u == 0) |
                 (Dist.time == "After" & s == 1 & d == 0 & c == 0 & b == 0 & u == 0) |
                 (Dist.time == "After" & s == 0 & d == 1 & c == 0 & b == 0 & u == 0) |
                 (Dist.time == "After" & s == 0 & d == 0 & c == 1 & b == 0 & u == 0) |
                 (Dist.time == "After" & s == 0 & d == 0 & c == 0 & b == 1 & u == 0) |
                 (Dist.time == "After" & s == 0 & d == 0 & c == 0 & b == 0 & u == 1)) |>
        mutate(total.points = 1)

      p <- t(brms::posterior_epred(mod_brm, newdata = newdata, se.fit = TRUE))
      cellmeans_brm <-
        newdata |>
        cbind(p) |>
        pivot_longer(cols = matches("^[0-9]*$"), names_to = ".draws") |>
        mutate(Dist = case_when(
          s == 1 ~ "s",
          c == 1 ~ "c",
          d == 1 ~ "d",
          b == 1 ~ "b",
          u == 1 ~ "u",
          .default = "Before"
        )) |>
        dplyr::select(-any_of(c(
          "s", "d", "c", "b", "u", "AIMS_REEF_NAME",
          "Dist.time", "Site", "Transect", "Dist.number", "total.points"
        )))

      cellmeans_summ_brm <- cellmeans_brm |>
              dplyr::select(-.draws) |> 
        group_by(Dist) |>
        summarise_draws(median, HDInterval::hdi) |> 
        dplyr::select(Dist, median, lower, upper)
      save(cellmeans_summ_brm, file = "../data/modelled/cellmeans_summ_brm_5.1.RData")
      
    }
    ## Contrasts
    {
      eff_brm <- cellmeans_brm |>
        mutate(Values = value) |>
        nest(.by = .draws) |>
        mutate(eff = map(
          .x = data,
          .f = ~ before_vs_afters(.x)
        )) |>
        dplyr::select(-data) |>
        unnest(c(eff))

      eff_summ_brm <- eff_brm |> 
        dplyr::select(-.draws) |>
        group_by(Dist) |> 
        summarise_draws(
          median,
          HDInterval::hdi
        )
      
      save(eff_summ_brm, file = "../data/modelled/eff_summ_brm_5.1.RData")
    }
  }
  ## INLA
  {
    ## Prepare data
    {
      data <- full_data

      ## Focus on only the necessary variables
      data <- data |>
        dplyr::select(
          n.points, total.points, Dist.time, s, c, d, b, u,
          AIMS_REEF_NAME, Site, Transect, Dist.number
        ) |>
        mutate(across(c(s, c, d, b, u), \(x) ifelse(Dist.time == "Before", 0, x))) 
      ## Cellmeans
      newdata <- data.frame(Dist.time = unique(data$Dist.time)) |>
              crossing(s = data$s, c = data$c, d = data$d, b = data$b, u = data$u)

      ## Restrict this to only where the sum of the rows is one so
      ## that our newdata is one row per disturbance (for Before and
      ## After)
      newdata <- newdata |>
              rowwise() |>
              filter(sum(c_across(where(is.numeric))) == 1)
      ## Further compress this to just a single Before (no disturbances)
      ## and single disturbance types for each After
      newdata <- newdata |>
        filter(Dist.time == "After") |>
        bind_rows(data.frame(
          Dist.time = "Before",
          s = 0, c = 0, d = 0, b = 0, u = 0
        ))

      data_pred <- data |>
              bind_rows(newdata) 
      i_newdata <- (nrow(data) + 1):nrow(data_pred)

      ## Pre-defined contrasts
      ## Compare Afters (for each disturbance) to Before (no disturbance)
      ## for each sector/shelf
      Xmat <- model.matrix(~ Dist.time + (s + b + c + d + u), data = newdata)

      nd <- newdata |>
        mutate(Dist = case_when(s == 1 ~ "s",
          c == 1 ~ "c",
          d == 1 ~ "d",
          b == 1 ~ "b",
          u == 1 ~ "u")) |>
        dplyr::select(-any_of(c("s", "c", "d", "b", "u"))) |>
        bind_cols(Xmat) |>
        filter(!is.na(Dist)) |>
        mutate(`(Intercept)` = 0,
          Dist.timeBefore = -1) |>
        dplyr::select(-Dist.time, -Dist) |>
        as.matrix()

      lincomb <- inla.make.lincombs(as.data.frame(nd))

      nd1 <- newdata |>
        mutate(Dist = case_when(s == 1 ~ "s",
          c == 1 ~ "c",
          d == 1 ~ "d",
          b == 1 ~ "b",
          u == 1 ~ "u")) |>
        dplyr::select(-any_of(c("s", "c", "d", "b", "u"))) |> 
        filter(!is.na(Dist)) 
    }
    ## Fit model
    {
      mod <- inla(n.points ~ Dist.time + (s + b + c + d + u) +
                    f(model = "iid", Dist.number) +
                    f(model = "iid", AIMS_REEF_NAME) +
                    f(model = "iid", Site) +
                    f(model = "iid", Transect),
        data = data_pred,
        Ntrials = data_pred$total.points,
        lincomb =  lincomb,
        family = "zeroinflatedbinomial1", #"binomial" 
        control.predictor = list(link = 1, compute = TRUE),
        control.compute = list(config = TRUE,
          dic = TRUE, waic = TRUE, cpo = TRUE,
          return.marginals.predictor = TRUE
        )
      )
      ## summary(mod)
      ## autoplot(mod)
      save(mod, file = "../data/modelled/mod_inla_5.1.RData")
    }
    ## Diagnostics
    {
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

      mod.cpo <- data.frame(cpo = mod$cpo$cpo,
        pit = mod$cpo$pit,
        failure = mod$cpo$failure) |>
        filter(failure == 0) |>
        mutate(row = 1:n())

      mod$cpo$failure

      g1 <- mod.cpo |>
        ggplot(aes(y = cpo, x = row)) +
        geom_point()
      g2 <- mod.cpo |>
        ggplot(aes(x = cpo)) +
        geom_histogram()
      g3 <- mod.cpo |>
        ggplot(aes(y = pit, x = row)) +
        geom_point()
      g4 <- mod.cpo |>
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
    ## DHARMa version 1
    {
      preds <- posterior_predict.inla(mod, newdata = data)
      mod_resids <- createDHARMa(simulatedResponse = t(preds),
        observedResponse = data$n.points,
        fittedPredictedResponse = apply(preds, 2, mean),
        integerResponse = TRUE)
      mod_resids |> plot()
      mod_resids |> testDispersion()
      mod_resids |> testZeroInflation()
    }
    ## Partial plots
    {
      newdata_pred <- newdata |>
        bind_cols(mod$summary.fitted.values[i_newdata, ])
      newdata_pred 

      newdata_pred <- newdata_pred |>
        mutate(Dist = case_when(
          s == 1 ~ "s",
          c == 1 ~ "c",
          d == 1 ~ "d",
          b == 1 ~ "b",
          u == 1 ~ "u"
        ))
      ## newdata_pred |>
      ##   ggplot(aes(y = `0.5quant`, x = factor(Dist), color = Dist.time)) +
      ##   geom_pointrange(
      ##     aes(
      ##       ymin = `0.025quant`,
      ##       ymax = `0.975quant`
      ##     ),
      ##     position = position_dodge(width = 0.5)
      ##   )
    }
    ## Partial plot - version 3
    {

      draws <- inla.posterior.sample(n = 1000, result = mod)
      contents <- mod$misc$configs$contents
      ## cellmeans <- newdata |> bind_cols(sapply(draws, function(x) x$latent[i_newdata]))
      cellmeans <- newdata |> cbind(sapply(draws, function(x) x$latent[i_newdata]))
      cellmeans_inla <- cellmeans |>
              ## pivot_longer(cols = matches("^[\\.]{3}[0-9]*"), names_to = ".draws") |>
              pivot_longer(cols = matches("^[0-9]*$"), names_to = ".draws") |>
              posterior::as_draws() |>
              mutate(.draw = .draws) |>
        mutate(Dist = case_when(
          s == 1 ~ "s",
          c == 1 ~ "c",
          d == 1 ~ "d",
          b == 1 ~ "b",
          u == 1 ~ "u",
          .default = "Before"
        )) |>
        dplyr::select(-any_of(c(
          "s", "d", "c", "b", "u", "AIMS_REEF_NAME",
          "Dist.time", "Site", "Transect", "Dist.number", "total.points"
        ))) |> 
        mutate(value = plogis(value))

      cellmeans_summ_inla <- cellmeans_inla |> 
      dplyr::select(-.draws) |>
        group_by(Dist) |> 
        summarise_draws(median, HDInterval::hdi)

      save(cellmeans_summ_inla, file = "../data/modelled/cellmeans_summ_inla_5.1.RData")
      
    }
    ## Contrasts
    {
      nd.pred <- nd1 |>
        bind_cols(mod$summary.lincomb.derived) |>
        mutate(across(where(is.numeric), exp))
      nd.pred |> ggplot(aes(x = `0.5quant`, y = Dist)) +
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
        mutate(Dist = case_when(s == 1 ~ "s",
          c == 1 ~ "c",
          d == 1 ~ "d",
          b == 1 ~ "b",
          u == 1 ~ "u")) |>
        dplyr::select(-any_of(c("s", "c", "d", "b", "u"))) |>
                ungroup() |>
                dplyr::select(-Dist.time) |>
                group_by(.draw) |>
                nest() |> 
      ## newdata_fitted[104, ] |> 
        mutate(data1 = map(
          .x = data,
          .f = ~ {
            .x <- .x |> slice(6:1)
            xmat <- cbind(-1, 1 * contr.treatment(6, base = 1, contrast = TRUE))
            xmat <- xmat[-1, ]
            x <- log(as.vector(as.vector(.x$Values)))
            data.frame(
                    Dist = .x$Dist[-1],
                    Values = exp(as.vector(x %*% t(xmat)))
            )
           
          }
        )) |>
        dplyr::select(-data) |>
                unnest(c(data1)) |>
                ungroup() |>
                group_by(Dist) |> 
        posterior::summarise_draws(
          median,
          HDInterval::hdi,
          Pl = ~ mean(.x < 1),
          Pg = ~ mean(.x > 1)
        )

      newdata_fitted |> ggplot(aes(x = median, y = Dist)) +
              geom_vline(xintercept = 1, linetype = "dashed") +
              geom_pointrange(aes(xmin = lower, xmax = upper)) +
              scale_x_continuous("Effect (Before - After) on a fold scale",
                      trans = "log2", breaks = scales::extended_breaks(n = 8)
              )

    }
    ## Contrasts - version 3
    {

      eff_inla <- cellmeans_inla |>
        mutate(Values = value) |>
        nest(.by = .draws) |>
        mutate(eff = map(
          .x = data,
          .f = ~ before_vs_afters(.x)
        )) |>
        dplyr::select(-data) |>
        unnest(c(eff))

      eff_summ_inla <- eff_inla |> 
        dplyr::select(-.draws) |>
        group_by(Dist) |> 
        summarise_draws(
          median,
          HDInterval::hdi
        )

      save(eff_summ_inla, file = "../data/modelled/eff_summ_inla_5.1.RData")
    }
  }
  ## Comparisons
  {
    
    load(file = "../data/modelled/cellmeans_summ_raw_5.1.RData")
    load(file = "../data/modelled/cellmeans_summ_glmmTMB_5.1.RData")
    load(file = "../data/modelled/cellmeans_summ_brm_5.1.RData")
    load(file = "../data/modelled/cellmeans_summ_inla_5.1.RData")

    cellmeans_summ <- bind_rows(
      cellmeans_summ_raw |>
        mutate(method = "raw") |>
        dplyr::rename(median = Mean, lower = lower, upper = upper),
      cellmeans_summ_glmmTMB |>
        mutate(method = "glmmTMB", type = "median") |>
        dplyr::rename(median = Pred),
      cellmeans_summ_brm |> mutate(method = "brm", type = "median") |>
        dplyr::rename(median = median, lower = lower, upper = upper),
      cellmeans_summ_inla |>
        mutate(method = "inla", type = "median") |> 
        dplyr::rename(median = median, lower = lower, upper = upper)
    )
    cellmeans_summ

    cellmeans_summ |>
      ggplot(aes(y = median, x = method, shape = type, colour = Dist)) +
      geom_pointrange(aes(ymin = lower, ymax = upper),
        position = position_dodge(width = 0.5)) 

    load(file = "../data/modelled/eff_summ_raw_5.1.RData")
    load(file = "../data/modelled/eff_summ_glmmTMB_5.1.RData")
    load(file = "../data/modelled/eff_summ_brm_5.1.RData")
    load(file = "../data/modelled/eff_summ_inla_5.1.RData")

    eff_summ <- bind_rows(
      eff_summ_raw |>
        mutate(method = "raw") |>
        dplyr::rename(median = Values),
      eff_summ_glmmTMB |>
        mutate(method = "glmmTMB", type = "median") |>
        dplyr::rename(median = Values),
      eff_summ_brm |>
        mutate(method = "brm", type = "median") |>
        dplyr::rename(median = median, lower = lower, upper = upper),
      eff_summ_inla |>
        mutate(method = "inla", type = "median") |> 
        dplyr::rename(median = median, lower = lower, upper = upper)
    )
    eff_summ

    eff_summ |>
            ggplot(aes(x = median, y = Dist, colour = method, shape = type)) +
            geom_vline(xintercept = 1, linetype = "dashed") +
      geom_pointrange(aes(xmin = lower, xmax = upper),
        position = position_dodge(width = 0.5)) +
            scale_x_continuous("Effect (Before - After) on a fold scale",
                    trans = "log2", breaks = scales::extended_breaks(n = 8)
            )
  }
}



## SecShelf, Dist.time and Disturbances (zi intercept only) 
{
  ## Prepare data
  {
    data <- full_data

    ## Focus on only the necessary variables
    data <- data |>
      dplyr::select(
        n.points, total.points, Dist.time, s, c, d, b, u,
        AIMS_REEF_NAME, Site, Transect, SecShelf, Dist.number,
        SecShelf
      ) |>
      filter(!SecShelf %in% c("CG M", "PC M", "PC O")) |>
              mutate(across(c(s, c, d, b, u), \(x) ifelse(Dist.time == "Before", 0, x))) 
  }
  ## Raw means
  {
    data_summ_s <- data |>
      dplyr::select(n.points, total.points, Dist.time, SecShelf, s, c, d, b, u) |>
      group_by(Dist.time, SecShelf, s) |>
      summarise(
        Median = median(n.points / total.points),
        Mean = mean(n.points / total.points),
        SD = sd(n.points / total.points),
        N = n()
      ) |>
      mutate(
        lower = Mean - 2 * (SD / N),
        upper = Mean + 2 * (SD / N)
      ) |>
      separate(SecShelf, into = c("A_SECTOR", "SHELF"), sep = " ") |>
              filter((Dist.time == "Before" & s == 0) | Dist.time == "After" & s == 1)

    data_summ_c <- data |>
      dplyr::select(n.points, total.points, Dist.time, SecShelf, s, c, d, b, u) |>
      group_by(Dist.time, SecShelf, c) |>
      summarise(
        Median = median(n.points / total.points),
        Mean = mean(n.points / total.points),
        SD = sd(n.points / total.points),
        N = n()
      ) |>
      mutate(
        lower = Mean - 2 * (SD / N),
        upper = Mean + 2 * (SD / N)
      ) |>
      separate(SecShelf, into = c("A_SECTOR", "SHELF"), sep = " ") |>
              filter((Dist.time == "Before" & c == 0) | Dist.time == "After" & c == 1)

    data_summ_d <- data |>
      dplyr::select(n.points, total.points, Dist.time, SecShelf, s, c, d, b, u) |>
      group_by(Dist.time, SecShelf, d) |>
      summarise(
        Median = median(n.points / total.points),
        Mean = mean(n.points / total.points),
        SD = sd(n.points / total.points),
        N = n()
      ) |>
      mutate(
        lower = Mean - 2 * (SD / N),
        upper = Mean + 2 * (SD / N)
      ) |>
      separate(SecShelf, into = c("A_SECTOR", "SHELF"), sep = " ") |>
              filter((Dist.time == "Before" & d == 0) | Dist.time == "After" & d == 1)

    data_summ_b <- data |>
      dplyr::select(n.points, total.points, Dist.time, SecShelf, s, c, d, b, u) |>
      group_by(Dist.time, SecShelf, b) |>
      summarise(
        Median = median(n.points / total.points),
        Mean = mean(n.points / total.points),
        SD = sd(n.points / total.points),
        N = n()
      ) |>
      mutate(
        lower = Mean - 2 * (SD / N),
        upper = Mean + 2 * (SD / N)
      ) |>
      separate(SecShelf, into = c("A_SECTOR", "SHELF"), sep = " ") |>
              filter((Dist.time == "Before" & b == 0) | Dist.time == "After" & b == 1)
    
    data_summ_u <- data |>
      dplyr::select(n.points, total.points, Dist.time, SecShelf, s, c, d, b, u) |>
      group_by(Dist.time, SecShelf, u) |>
      summarise(
        Median = median(n.points / total.points),
        Mean = mean(n.points / total.points),
        SD = sd(n.points / total.points),
        N = n()
      ) |>
      mutate(
        lower = Mean - 2 * (SD / N),
        upper = Mean + 2 * (SD / N)
      ) |>
      separate(SecShelf, into = c("A_SECTOR", "SHELF"), sep = " ") |>
              filter((Dist.time == "Before" & u == 0) | Dist.time == "After" & u == 1)

    data_summ <- bind_rows(data_summ_s, data_summ_c, data_summ_d,
      data_summ_b, data_summ_u) |> 
      mutate(Dist = case_when(
        s == 1 ~ "s",
        c == 1 ~ "c",
        d == 1 ~ "d",
        b == 1 ~ "b",
        u == 1 ~ "u",
        .default = "Before"
      ))
      
    ## data_summ |>
    ##   ggplot(aes(y = Mean, x = Dist.time)) +
    ##   geom_pointrange(aes(ymin = lower, ymax = upper)) +
    ##   facet_grid(A_SECTOR ~ SHELF, scales = "free")

    data_summ |> ggplot(aes(y = Mean, x = Dist, colour = factor(Dist.time))) +
      geom_pointrange(aes(ymin = lower, ymax = upper), position = position_dodge(width = 0.5)) +
      facet_grid(A_SECTOR ~ SHELF, scales = "free")
  
  }
  ## glmmTMB
  {
    ## Fit model
    {
      mod_glmmTMB <- glmmTMB(cbind(n.points, total.points-n.points) ~ 1+
                               (s + c + d + b + u) +
                               (1 | Dist.time:SecShelf:(s + c + d + b + u)) +
                               (1|AIMS_REEF_NAME) +
                               (1|Site) +
                               (1|Transect),
        ziformula = ~1 + (1|AIMS_REEF_NAME) +
                               (1|Site) +
                               (1|Transect),
        data = data,
        family = "binomial", 
        REML = TRUE
      )

      summary(mod_glmmTMB)

    }
    ## DHARMa
    {
      resids <- simulateResiduals(mod_glmmTMB, plot = TRUE, integerResponse = TRUE)
      wrap_elements(~testUniformity(resids)) +
        wrap_elements(~plotResiduals(resids)) +
        wrap_elements(~testDispersion(resids)) 
      testDispersion(resids)
      testZeroInflation(resids)
    }
    ## Partial plots
    {
      ## Since this model operates mainly through random effects, we cannot
      ## use emmeans to produce the cell means or contrasts
      ## Instead we just have to perform predictions.
      
      newdata <- crossing(Dist.time = data$Dist.time) |>
        crossing(s = c(0,1), d = c(0, 1), c = c(0, 1), b = c(0, 1), u = c(0, 1),
          SecShelf = as.character(unique(data$SecShelf))) |>
        mutate(AIMS_REEF_NAME = NA, Site = NA, Transect = NA, Dist.number =  NA)
      newdata <- newdata |>
        filter((Dist.time == "Before" & s == 0 & d == 0 & c == 0 & b == 0 & u == 0) |
                 (Dist.time == "After" & s == 1 & d == 0 & c == 0 & b == 0 & u == 0) |
                 (Dist.time == "After" & s == 0 & d == 1 & c == 0 & b == 0 & u == 0) |
                 (Dist.time == "After" & s == 0 & d == 0 & c == 1 & b == 0 & u == 0) |
                 (Dist.time == "After" & s == 0 & d == 0 & c == 0 & b == 1 & u == 0) |
                 (Dist.time == "After" & s == 0 & d == 0 & c == 0 & b == 0 & u == 1) 
        )

      p <- predict(mod_glmmTMB, newdata = newdata, se.fit = TRUE)
      newdata <- newdata |>
        bind_cols(fit = p$fit, se = p$se.fit) |>
        mutate(Pred = plogis(fit), lower = plogis(fit - 2 * se), upper = plogis(fit + 2 * se))
      newdata <- newdata |>
        mutate(Dist = case_when(
          s == 1 ~ "s",
          c == 1 ~ "c",
          d == 1 ~ "d",
          b == 1 ~ "b",
          u == 1 ~ "u",
          .default = "Before"
        )) |>
        separate(SecShelf, into = c("A_SECTOR", "SHELF"), sep = " ")

      newdata |> dplyr::select(Dist.time, Dist, A_SECTOR, SHELF, Pred, lower, upper)
      
      newdata |>
        ggplot(aes(y = Pred, x = Dist, colour = factor(Dist.time))) +
        geom_pointrange(aes(ymin = lower, ymax = upper),
          position = position_dodge(width = 0.5)) +
        facet_grid(A_SECTOR ~ SHELF, scales = "free")
    }
    ## Contrasts
    {
      ## These are not really possible - at least not with confidence intervals
    }
    
  }
  ## brms
  {
    form <- bf(
      n.points | trials(total.points) ~ 1+ (s + c + d + b + u) + 
        (1 | Dist.time:SecShelf:(s + c + d + b + u)) +
        (1 | AIMS_REEF_NAME) +
        (1 | Site) +
        (1 | Transect),
      zi = ~ 1 +,
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
      silent =  0#,
      ## refresh = 100
    )
    summary(mod_brm)
    save(mod_brm, file = "../data/modelled/mod_brm_6.1.RData")
    load(file = "../data/modelled/mod_brm_6.1.RData")
  }
  ## INLA
  {
    ## Prepare data OLD
    {
      if (1 != 2) {
        data <- full_data

        ## Focus on only the necessary variables
        data <- data |>
          dplyr::select(
            n.points, total.points, Dist.time, s, c, d, b, u,
            AIMS_REEF_NAME, Site, Transect, SecShelf, Dist.number
          ) |>
          mutate(across(c(s, c, d, b, u), \(x) ifelse(Dist.time == "Before", 0, x)))
        ## Cellmeans
        newdata <- data.frame(Dist.time = unique(data$Dist.time)) |>
          crossing(
            SecShelf = data$SecShelf,
            s = data$s, c = data$c, d = data$d, b = data$b, u = data$u
          )

        ## Restrict this to only where the sum of the rows is one so
        ## that our newdata is one row per disturbance (for Before and
        ## After)
        newdata <- newdata |>
          rowwise() |>
          filter(sum(c_across(where(is.numeric))) == 1)
        ## Further compress this to just a single Before (no disturbances)
        ## and single disturbance types for each After
        newdata <- newdata |>
          filter(Dist.time == "After") |>
          bind_rows(data.frame(
            Dist.time = "Before",
            s = 0, c = 0, d = 0, b = 0, u = 0
          ) |> crossing(SecShelf = newdata$SecShelf))

        data_pred <- data |>
          bind_rows(newdata) |>
          mutate(SecShelf = forcats::fct_relevel(SecShelf, "CA M"))
        i_newdata <- (nrow(data) + 1):nrow(data_pred)

        ## Pre-defined contrasts
        ## Compare Afters (for each disturbance) to Before (no disturbance)
        ## for each sector/shelf
        Xmat <- model.matrix(~ SecShelf*Dist.time * (s + b + c + d + u), data = newdata)

        nd <- newdata |>
          mutate(Dist = case_when(s == 1 ~ "s",
            c == 1 ~ "c",
            d == 1 ~ "d",
            b == 1 ~ "b",
            u == 1 ~ "u")) |>
          dplyr::select(-any_of(c("s", "c", "d", "b", "u"))) |>
          bind_cols(Xmat) |>
          filter(!is.na(Dist)) |>
          mutate(`(Intercept)` = 0) |>
          dplyr::select(-Dist.time, -Dist, -SecShelf) |>
          as.matrix()

        lincomb <- inla.make.lincombs(as.data.frame(nd))

        nd1 <- newdata |>
          mutate(Dist = case_when(s == 1 ~ "s",
            c == 1 ~ "c",
            d == 1 ~ "d",
            b == 1 ~ "b",
            u == 1 ~ "u")) |>
          dplyr::select(-any_of(c("s", "c", "d", "b", "u"))) |> 
          filter(!is.na(Dist)) 
      }
    }
    ## Prepare data
    {
      newdata <- crossing(Dist.time = data$Dist.time) |>
        crossing(s = c(0,1), d = c(0, 1), c = c(0, 1), b = c(0, 1), u = c(0, 1),
          SecShelf = as.character(unique(data$SecShelf))) |>
        separate(SecShelf, into = c("A_SECTOR", "SHELF"), sep = " ", remove = FALSE) |> 
        mutate(AIMS_REEF_NAME = NA, Site = NA, Transect = NA, Dist.number =  NA)
      newdata <- newdata |>
        filter((Dist.time == "Before" & s == 0 & d == 0 & c == 0 & b == 0 & u == 0) |
                 (Dist.time == "After" & s == 1 & d == 0 & c == 0 & b == 0 & u == 0) |
                 (Dist.time == "After" & s == 0 & d == 1 & c == 0 & b == 0 & u == 0) |
                 (Dist.time == "After" & s == 0 & d == 0 & c == 1 & b == 0 & u == 0) |
                 (Dist.time == "After" & s == 0 & d == 0 & c == 0 & b == 1 & u == 0) |
                 (Dist.time == "After" & s == 0 & d == 0 & c == 0 & b == 0 & u == 1) 
        )

      data_pred <- data |>
        bind_rows(newdata |>
                    dplyr::select(Dist.time, SecShelf, s, b, c, d, u)) |> 
        mutate(SS = paste(Dist.time, SecShelf, s, b, c, d, u))
      i_newdata <- (nrow(data) + 1):nrow(data_pred)
    }
    ## Fit model
    {
      mod <- inla(n.points ~ Dist.time + (s + b + c + d + u) +
                    f(model = "iid", SS) +
                    f(model = "iid", AIMS_REEF_NAME) +
                    f(model = "iid", Site) +
                    f(model = "iid", Transect),
        data = data_pred,
        Ntrials = data_pred$total.points,
        family = "zeroinflatedbinomial1", #"binomial" 
        control.predictor = list(link = 1, compute = TRUE),
        control.compute = list(config = TRUE,
          dic = TRUE, waic = TRUE, cpo = TRUE
        )
      )
      summary(mod)
      ## autoplot(mod)
    }
    ## Fit model
    {
      mod <- inla(n.points ~ SecShelf * Dist.time * (s + b + c + d + u) +
                    f(model = "iid", Dist.number) +
                    f(model = "iid", AIMS_REEF_NAME) +
                    f(model = "iid", Site) +
                    f(model = "iid", Transect),
        data = data_pred,
        Ntrials = data_pred$total.points,
        ## lincomb =  lincomb,
        family = "zeroinflatedbinomial1", #"binomial" 
        control.predictor = list(link = 1, compute = TRUE),
        control.compute = list(config = TRUE,
          dic = TRUE, waic = TRUE, cpo = TRUE
        )
      )
      summary(mod)
      ## autoplot(mod)
    }
    ## Diagnostics
    {
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

      mod.cpo <- data.frame(cpo = mod$cpo$cpo,
        pit = mod$cpo$pit,
        failure = mod$cpo$failure) |>
        filter(failure == 0) |>
        mutate(row = 1:n())

      mod$cpo$failure

      g1 <- mod.cpo |>
        ggplot(aes(y = cpo, x = row)) +
        geom_point()
      g2 <- mod.cpo |>
        ggplot(aes(x = cpo)) +
        geom_histogram()
      g3 <- mod.cpo |>
        ggplot(aes(y = pit, x = row)) +
        geom_point()
      g4 <- mod.cpo |>
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
    ## DHARMa version 1
    {
      preds <- posterior_predict.inla(mod, newdata = data)
      mod_resids <- createDHARMa(simulatedResponse = t(preds),
        observedResponse = data$n.points,
        fittedPredictedResponse = apply(preds, 2, mean),
        integerResponse = TRUE)
      mod_resids |> plot()
      mod_resids |> testDispersion()
      mod_resids |> testZeroInflation()
    }
    ## Partial plots
    {
      newdata_pred <- newdata |>
        bind_cols(mod$summary.fitted.values[i_newdata, ]) |> 
        mutate(Dist = case_when(
          s == 1 ~ "s",
          c == 1 ~ "c",
          d == 1 ~ "d",
          b == 1 ~ "b",
          u == 1 ~ "u",
          .default = "Before"
        )) |>
        separate(SecShelf, into = c("A_SECTOR", "SHELF"), sep = " ")
      newdata_pred 
      
      newdata_pred |>
        ggplot(aes(y = `0.5quant`, x = Dist, colour = factor(Dist.time))) +
        geom_pointrange(aes(ymin = `0.025quant`, ymax = `0.975quant`),
          position = position_dodge(width = 0.5)) +
        facet_grid(A_SECTOR ~ SHELF, scales = "free")


      newdata_pred |> ggplot(aes(y = `0.5quant`, x = Dist, colour = factor(Dist.time))) +
        geom_pointrange(aes(ymin = `0.025quant`, ymax = `0.975quant`, shape = "INLA"), position = position_dodge(width = 0.5)) +
        geom_pointrange(data = data_summ, aes(y = Mean, ymin = lower, ymax = upper, shape = "Raw"), position = position_nudge(x = 0.2)) +
        geom_pointrange(data =  newdata, aes(y =  Pred, ymin = lower, ymax = upper, shape = "glmmTMB"), position = position_nudge(x = -0.2)) +
        facet_wrap(A_SECTOR ~ SHELF, scales = "free")

      
      ## newdata_pred <- newdata_pred |>
      ##   mutate(Dist = case_when(
      ##     s == 1 ~ "s",
      ##     c == 1 ~ "c",
      ##     d == 1 ~ "d",
      ##     b == 1 ~ "b",
      ##     u == 1 ~ "u"
      ##   ))
      ## newdata_pred <- newdata_pred |>
      ##   separate(SecShelf, into = c("A_SECTOR", "SHELF"), remove = FALSE)
      ## ## newdata_before <- data.frame(Dist = unique(newdata_pred$Dist)) |>
      ## ##   crossing(
      ## ##     A_SECTOR = newdata_pred$A_SECTOR,
      ## ##     SHELF = newdata_pred$SHELF
      ## ##   ) |>
      ## ##   bind_cols(newdata_pred |> filter(Dist.time == "Before") |> dplyr::select(-SecShelf, -A_SECTOR, -SHELF, -s, -c, -d, -b, -u, -Dist))
      
      ##   newdata_pred |>
      ##           ## bind_rows(newdata_before) |> 
      ##         separate(SecShelf, into = c("A_SECTOR", "SHELF"), remove = FALSE) |> 
      ##         ggplot(aes(y = `0.5quant`, x = factor(Dist), color = Dist.time)) +
      ##         geom_pointrange(
      ##                 aes(
      ##                         ymin = `0.025quant`,
      ##                         ymax = `0.975quant`
      ##                 ),
      ##                 position = position_dodge(width = 0.5)
      ##         ) +
      ## facet_grid(A_SECTOR ~ SHELF)
    }
    ## Partial plots - version 2
    {
      draws <- inla.posterior.sample(n=1000, result = mod)
      cellmeans <- newdata |>
        mutate(Dist = case_when(
          s == 1 ~ "s",
          c == 1 ~ "c",
          d == 1 ~ "d",
          b == 1 ~ "b",
          u == 1 ~ "u",
          .default = "Before"
        )) |>
        separate(SecShelf, into = c("A_SECTOR", "SHELF"), sep = " ") |> 
        bind_cols(sapply(draws, function(x) x$latent[i_newdata]))

      cellmeans <- cellmeans |>
              pivot_longer(cols = matches("^[\\.]{3}[0-9]*"), names_to = ".draws") |>
              dplyr::select(Dist.time, A_SECTOR, SHELF, Dist, .draws, value) |>
              posterior::as_draws() |>
              mutate(.draw = .draws) |>
        mutate(value = plogis(value)) |> 
        dplyr::select(-.draws)
      cellmeans_summ <- cellmeans |> 
        group_by(Dist.time, A_SECTOR, SHELF, Dist) |>
        summarise_draws(median, HDInterval::hdi)

      cellmeans_summ |>
        ggplot(aes(y = median, x = Dist, colour = factor(Dist.time))) +
        geom_pointrange(aes(ymin = lower, ymax = upper),
          position = position_dodge(width = 0.5)) +
        facet_grid(A_SECTOR ~ SHELF, scales = "free")
    }
    ## Contrasts - version 2
    {
      eff <-
              cellmeans |>
        ## filter(.draw == "...14", A_SECTOR == "CA", SHELF == "I") |>
              nest(.by = c(A_SECTOR, SHELF, .draw)) |>
              mutate(data1 = map(
                      .x = data,
                      .f = ~ {
                              .x <- .x |>
                                      mutate(Dist = factor(Dist,
                                              levels = c("Before", "s", "c", "d", "b", "u")
                                      )) |>
                                      arrange(Dist)
                              xmat <- cbind(-1, 1 * contr.treatment(6, base = 1, contrast = TRUE))
                              xmat <- xmat[-1, ]
                              ## print(.x)
                              ## print(xmat)
                              x <- log(as.vector(as.vector(.x$value)))
                              ## print(x)
                              ## print(as.vector(x %*% t(xmat)))
                              ## print(exp(as.vector(x %*% t(xmat))))
                              data.frame(
                                      Dist = .x$Dist[-1],
                                      Values = exp(as.vector(x %*% t(xmat)))
                              )
                      }, .progress = TRUE
              ))

      eff_summ <- eff |> dplyr::select(-data) |> 
                unnest(c(data1)) |>
                ungroup() |>
                group_by(A_SECTOR, SHELF, Dist) |> 
        posterior::summarise_draws(
          median,
          HDInterval::hdi,
          Pl = ~ mean(.x < 1),
          Pg = ~ mean(.x > 1)
        )


      eff_summ |>
        ggplot(aes(x = median, y = Dist)) +
        geom_vline(xintercept = 1, linetype = "dashed") +
        geom_pointrange(aes(xmin = lower, xmax = upper)) +
        scale_x_continuous("Effect (Before - After) on a fold scale",
          trans = "log2",
          breaks = c(0.5, 1, 2)
          ## breaks = scales::log_breaks(n = 8, base = 2)
          ## breaks = scales::pretty_breaks(n = 8)
          ## breaks = scales::extended_breaks(n = 8)
        ) +
        facet_grid(A_SECTOR ~ SHELF, scales = "free")
    }
  }
}



#####


## SecShelf, Dist.time and Disturbances (zi full RE) 
{
  ## Prepare data
  {
    data <- full_data

    ## Focus on only the necessary variables
    data <- data |>
      dplyr::select(
        n.points, total.points, Dist.time, s, c, d, b, u,
        AIMS_REEF_NAME, Site, Transect, SecShelf, Dist.number,
        SecShelf
      ) |>
      filter(!SecShelf %in% c("CG M", "PC M", "PC O")) |>
              mutate(across(c(s, c, d, b, u), \(x) ifelse(Dist.time == "Before", 0, x))) 
  }
  ## glmmTMB
  {
    ## Fit model
    {
      mod_glmmTMB <- glmmTMB(cbind(n.points, total.points-n.points) ~ 1+ (s + c + d + b + u) +
                               (1 | Dist.time:SecShelf:(s + c + d + b + u)) +
                               (1|AIMS_REEF_NAME) +
                               (1|Site) +
                               (1|Transect),
        ziformula = ~1 + (1|AIMS_REEF_NAME) +
                               (1|Site) +
                               (1|Transect),
        data = data,
        family = "binomial", 
        REML = TRUE
      )

      summary(mod_glmmTMB)
      ranef(mod_glmmTMB) |> str()
      ranef(mod_glmmTMB)[[1]][[2]]
      newdata <- crossing(Dist.time = data$Dist.time) |>
        crossing(s = c(0,1), d = c(0, 1), c = c(0, 1), b = c(0, 1), u = c(0, 1),
          SecShelf = as.character(unique(data$SecShelf)[1:10])) |>
        mutate(AIMS_REEF_NAME = NA, Site = NA, Transect = NA, Dist.number =  NA)
      newdata <- newdata |>
        filter((Dist.time == "Before" & s == 0 & d == 0 & c == 0 & b == 0 & u == 0) |
                 (Dist.time == "After" & s == 1 & d == 0 & c == 0 & b == 0 & u == 0) |
                 (Dist.time == "After" & s == 0 & d == 1 & c == 0 & b == 0 & u == 0) |
                 (Dist.time == "After" & s == 0 & d == 0 & c == 1 & b == 0 & u == 0) |
                 (Dist.time == "After" & s == 0 & d == 0 & c == 0 & b == 1 & u == 0) |
                 (Dist.time == "After" & s == 0 & d == 0 & c == 0 & b == 0 & u == 1) 
        )
      p <- predict(mod_glmmTMB, newdata = newdata, se.fit = TRUE)
      newdata <- newdata |>
        bind_cols(fit = p$fit, se = p$se.fit) |>
        mutate(Pred = plogis(fit), lower = plogis(fit - 2 * se), upper = plogis(fit + 2 * se))
      newdata <- newdata |>
              mutate(Dist = case_when(
                      s == 1 ~ "s",
                      c == 1 ~ "c",
                      d == 1 ~ "d",
                      b == 1 ~ "b",
                      u == 1 ~ "u",
                      .default = "Before"
              )) |>
        separate(SecShelf, into = c("A_SECTOR", "SHELF"), sep = " ")
      newdata
      newdata |> ggplot(aes(y = Pred, x = Dist, colour = factor(Dist.time))) +
        geom_pointrange(aes(ymin = lower, ymax = upper), position = position_dodge(width = 0.5)) +
        facet_grid(A_SECTOR ~ SHELF, scales = "free")
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
      mod_glmmTMB |>
        emmeans(~ c, type = "response") |>
        as.data.frame() |>
        separate(SecShelf, into = c("A_SECTOR", "SHELF"), remove = FALSE) |> 
        ggplot(aes(y = prob, x = A_SECTOR, colour = Dist.time)) +
        geom_pointrange(aes(ymin = asymp.LCL, ymax = asymp.UCL), position = position_dodge(width = 0.5)) +
        facet_grid(~SHELF)
    }
    ## Contrasts
    {
      mod_glmmTMB |>
              emmeans(~ Dist.time | SecShelf, type = "response") |>
              contrast(method = list(Dist.time = c(-1, 1))) |>
              summary(infer = TRUE) |>
              as.data.frame() |>
              separate(SecShelf, into = c("A_SECTOR", "SHELF"), remove = FALSE) |>
              ggplot(aes(x = odds.ratio, y = A_SECTOR)) +
              geom_vline(xintercept = 1, linetype = "dashed") +
              geom_pointrange(aes(xmin = asymp.LCL, xmax = asymp.UCL)) +
        scale_x_continuous("Effect (Before - After) on a fold scale", trans = "log2", breaks = scales::extended_breaks(n = 8)) +
        facet_grid(~SHELF)
    }
    
    ## Fit model
    {
      data_pred <- data |>
        bind_rows(newdata |>
                    mutate(SecShelf =  paste(A_SECTOR, SHELF)) |>
                    dplyr::select(Dist.time, SecShelf, s, b, c, d, u)) |> 
        mutate(SS = paste(Dist.time, SecShelf, s, b, c, d, u))
      i_newdata <- (nrow(data) + 1):nrow(data_pred)
      mod <- inla(n.points ~ Dist.time + (s + b + c + d + u) +
                    f(model = "iid", SS) +
                    f(model = "iid", AIMS_REEF_NAME) +
                    f(model = "iid", Site) +
                    f(model = "iid", Transect),
        data = data_pred,
        Ntrials = data_pred$total.points,
        ## lincomb =  lincomb,
        family = "zeroinflatedbinomial1", #"binomial" 
        control.predictor = list(link = 1, compute = TRUE),
        control.compute = list(config = TRUE,
          dic = TRUE, waic = TRUE, cpo = TRUE
        )
      )
      summary(mod)
      ## autoplot(mod)
      newdata_pred <- newdata |>
        bind_cols(mod$summary.fitted.values[i_newdata, ])
      newdata_pred 
      newdata_pred |> ggplot(aes(y = `0.5quant`, x = Dist, colour = factor(Dist.time))) +
        geom_pointrange(aes(ymin = `0.025quant`, ymax = `0.975quant`), position = position_dodge(width = 0.5)) +
        facet_grid(A_SECTOR ~ SHELF, scales = "free")


      
  draws <- inla.posterior.sample(n=1000, result = mod)
  cellmeans <- newdata |> bind_cols(sapply(draws, function(x) x$latent[i_newdata]))
cellmeans <- cellmeans |>
        pivot_longer(cols = matches("^[\\.]{3}[0-9]*"), names_to = ".draws") |>
        dplyr::select(Dist.time, A_SECTOR, SHELF, Dist, .draws, value) |>
        posterior::as_draws() |>
        mutate(.draw = .draws) |>
        dplyr::select(-.draws) |>
        mutate(value = plogis(value)) |> 
        group_by(Dist.time, A_SECTOR, SHELF, Dist) |>
        summarise_draws(median, HDInterval::hdi)

      cellmeans |> ggplot(aes(y = median, x = Dist, colour = factor(Dist.time))) +
        geom_pointrange(aes(ymin = lower, ymax = upper), position = position_dodge(width = 0.5)) +
        facet_grid(A_SECTOR ~ SHELF, scales = "free")
      
  contents <- mod$misc$configs$contents
      preds <- posterior_predict.inla(mod, newdata = data_pred[i_newdata, ])
    }
  }
  ## brms
  {
    form <- bf(
      n.points | trials(total.points) ~ 1+ (s + c + d + b + u) + 
        (1 | Dist.time:SecShelf:(s + c + d + b + u)) +
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
      silent =  0#,
      ## refresh = 100
    )
    summary(mod_brm)
    save(mod_brm, file = "../data/modelled/mod_brm_7.1.RData")
    load(file = "../data/modelled/mod_brm_7.1.RData")
  }
  ## INLA
  {
    ## Prepare data
    {
      data <- full_data

      ## Focus on only the necessary variables
      data <- data |>
        dplyr::select(
          n.points, total.points, Dist.time, s, c, d, b, u,
          AIMS_REEF_NAME, Site, Transect, SecShelf, Dist.number
        ) |>
        mutate(across(c(s, c, d, b, u), \(x) ifelse(Dist.time == "Before", 0, x)))
      ## Cellmeans
      newdata <- data.frame(Dist.time = unique(data$Dist.time)) |>
        crossing(
          SecShelf = data$SecShelf,
          s = data$s, c = data$c, d = data$d, b = data$b, u = data$u
        )

      ## Restrict this to only where the sum of the rows is one so
      ## that our newdata is one row per disturbance (for Before and
      ## After)
      newdata <- newdata |>
              rowwise() |>
              filter(sum(c_across(where(is.numeric))) == 1)
      ## Further compress this to just a single Before (no disturbances)
      ## and single disturbance types for each After
      newdata <- newdata |>
        filter(Dist.time == "After") |>
        bind_rows(data.frame(
          Dist.time = "Before",
          s = 0, c = 0, d = 0, b = 0, u = 0
        ) |> crossing(SecShelf = newdata$SecShelf))

      data_pred <- data |>
              bind_rows(newdata) |>
              mutate(SecShelf = forcats::fct_relevel(SecShelf, "CA M"))
      i_newdata <- (nrow(data) + 1):nrow(data_pred)

      ## Pre-defined contrasts
      ## Compare Afters (for each disturbance) to Before (no disturbance)
      ## for each sector/shelf
      Xmat <- model.matrix(~ SecShelf*Dist.time * (s + b + c + d + u), data = newdata)

      nd <- newdata |>
        mutate(Dist = case_when(s == 1 ~ "s",
          c == 1 ~ "c",
          d == 1 ~ "d",
          b == 1 ~ "b",
          u == 1 ~ "u")) |>
        dplyr::select(-any_of(c("s", "c", "d", "b", "u"))) |>
        bind_cols(Xmat) |>
        filter(!is.na(Dist)) |>
        mutate(`(Intercept)` = 0) |>
        dplyr::select(-Dist.time, -Dist, -SecShelf) |>
       as.matrix()

      lincomb <- inla.make.lincombs(as.data.frame(nd))

      nd1 <- newdata |>
        mutate(Dist = case_when(s == 1 ~ "s",
          c == 1 ~ "c",
          d == 1 ~ "d",
          b == 1 ~ "b",
          u == 1 ~ "u")) |>
        dplyr::select(-any_of(c("s", "c", "d", "b", "u"))) |> 
        filter(!is.na(Dist)) 
    }
    ## Fit model
    {
      mod <- inla(n.points ~ SecShelf * Dist.time * (s + b + c + d + u) +
                    f(model = "iid", Dist.number) +
                    f(model = "iid", AIMS_REEF_NAME) +
                    f(model = "iid", Site) +
                    f(model = "iid", Transect),
        data = data_pred,
        Ntrials = data_pred$total.points,
        ## lincomb =  lincomb,
        family = "zeroinflatedbinomial1", #"binomial" 
        control.predictor = list(link = 1, compute = TRUE),
        control.compute = list(config = TRUE,
          dic = TRUE, waic = TRUE, cpo = TRUE
        )
      )
      summary(mod)
      ## autoplot(mod)
    }
    ## Diagnostics
    {
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

      mod.cpo <- data.frame(cpo = mod$cpo$cpo,
        pit = mod$cpo$pit,
        failure = mod$cpo$failure) |>
        filter(failure == 0) |>
        mutate(row = 1:n())

      mod$cpo$failure

      g1 <- mod.cpo |>
        ggplot(aes(y = cpo, x = row)) +
        geom_point()
      g2 <- mod.cpo |>
        ggplot(aes(x = cpo)) +
        geom_histogram()
      g3 <- mod.cpo |>
        ggplot(aes(y = pit, x = row)) +
        geom_point()
      g4 <- mod.cpo |>
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
    ## DHARMa version 1
    {
      preds <- posterior_predict.inla(mod, newdata = data)
      mod_resids <- createDHARMa(simulatedResponse = t(preds),
        observedResponse = data$n.points,
        fittedPredictedResponse = apply(preds, 2, mean),
        integerResponse = TRUE)
      mod_resids |> plot()
      mod_resids |> testDispersion()
      mod_resids |> testZeroInflation()
    }
    ## Partial plots
    {
      newdata_pred <- newdata |>
        bind_cols(mod$summary.fitted.values[i_newdata, ])
      newdata_pred 

      newdata_pred <- newdata_pred |>
        mutate(Dist = case_when(
          s == 1 ~ "s",
          c == 1 ~ "c",
          d == 1 ~ "d",
          b == 1 ~ "b",
          u == 1 ~ "u"
        ))
      newdata_pred <- newdata_pred |>
        separate(SecShelf, into = c("A_SECTOR", "SHELF"), remove = FALSE)
      ## newdata_before <- data.frame(Dist = unique(newdata_pred$Dist)) |>
      ##   crossing(
      ##     A_SECTOR = newdata_pred$A_SECTOR,
      ##     SHELF = newdata_pred$SHELF
      ##   ) |>
      ##   bind_cols(newdata_pred |> filter(Dist.time == "Before") |> dplyr::select(-SecShelf, -A_SECTOR, -SHELF, -s, -c, -d, -b, -u, -Dist))
      
        newdata_pred |>
                ## bind_rows(newdata_before) |> 
              separate(SecShelf, into = c("A_SECTOR", "SHELF"), remove = FALSE) |> 
              ggplot(aes(y = `0.5quant`, x = factor(Dist), color = Dist.time)) +
              geom_pointrange(
                      aes(
                              ymin = `0.025quant`,
                              ymax = `0.975quant`
                      ),
                      position = position_dodge(width = 0.5)
              ) +
      facet_grid(A_SECTOR ~ SHELF)
    }
    ## Contrasts
    {
      nd.pred <- nd1 |>
        bind_cols(mod$summary.lincomb.derived) |>
        mutate(across(where(is.numeric), exp))
      nd.pred |> ggplot(aes(x = `0.5quant`, y = Dist)) +
        geom_vline(xintercept = 1, linetype = "dashed") +
        geom_pointrange(aes(xmin = `0.025quant`, xmax = `0.975quant`)) +
        scale_x_continuous("Effect (Before - After) on a fold scale",
          trans = "log2", breaks = scales::extended_breaks(n = 8)
        )
    }
  }
}
