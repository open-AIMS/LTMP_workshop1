## Retrieve the data
load(file = "../data/processed/data_q2.R")


## Strip out all unnecessary variables

full_data <- data
data <- full_data


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
    ## Cellmeans
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
    }
    ## Effects
    {
      eff_summ_raw <- cellmeans_summ_raw |>
        group_by(type, A_SECTOR, SHELF) |>
        summarise(Mean = Mean[Dist.time == "After"] / Mean[Dist.time == "Before"])
      save(eff_summ_raw, file = "../data/modelled/eff_summ_raw_3.1.RData")
    }
  }
  ## glmmTMB
  {
    ## Fit model
    {
      if (rerun_models) {
        mod_glmmTMB <- glmmTMB(
          cbind(n.points, total.points - n.points) ~ Dist.time * SecShelf +
            (1 | AIMS_REEF_NAME) +
            (1 | Site) +
            (1 | Transect),
          ziformula = ~ 1 + (1 | AIMS_REEF_NAME) + (1 | Site) + (1 | Transect),
          data = data,
          family = "binomial",
          REML = TRUE
        )

        summary(mod_glmmTMB)
        save(mod_glmmTMB, file = "../data/modelled/mod_glmmTMB_3.1.RData")
      }
    }
    ## DHARMa
    {
      load(file = "../data/modelled/mod_glmmTMB_3.1.RData")
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
      if (rerun_models) {
        form <- bf(
          n.points | trials(total.points) ~ Dist.time * SecShelf +
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
        save(mod_brm, file = "../data/modelled/mod_brm_3.1.RData")
      }
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
      if (rerun_models) {
        mod <- inla(
          n.points ~ Dist.time * SecShelf +
            f(model = "iid", AIMS_REEF_NAME) +
            f(model = "iid", Site) +
            f(model = "iid", Transect),
          data = data_pred,
          Ntrials = data_pred$total.points,
          lincomb = lincomb,
          family = "zeroinflatedbinomial2", # "binomial",
          control.predictor = list(link = 1, compute = TRUE),
          control.compute = list(
            config = TRUE,
            dic = TRUE, waic = TRUE, cpo = TRUE,
            return.marginals.predictor = TRUE
          )
        )
        ## summary(mod)
        ## autoplot(mod)
        save(mod, file = "../data/modelled/mod_3.1.RData")
      }
    }
    ## diagnostics
    {
      load(file = "../data/modelled/mod_3.1.RData")
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
      load(file = "../data/modelled/mod_3.1.RData")
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
      load(file = "../data/modelled/mod_3.1.RData")
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
      cellmeans <- newdata |> cbind(sapply(draws, function(x) x$latent[i_newdata]))
      cellmeans_inla <- cellmeans |>
              pivot_longer(cols = matches("^[0-9]*$"), names_to = ".draws") |>
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
