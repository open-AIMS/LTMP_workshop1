## Retrieve the data
load(file = "../data/processed/data_q2.R")


## Strip out all unnecessary variables

full_data <- data
data <- full_data

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

## Start with an intercept only model
{
  ## glmmTMB
  {
    
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
  ## brms
  {
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
    load(file = "../data/modelled/mod_brm_1.1.RData")
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
      ## summary(mod)
      ## autoplot(mod)
    }

    ## Diagnostics
    {
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
      wrap_elements(~testUniformity(mod_resids)) +
        wrap_elements(~plotResiduals(mod_resids)) +
        wrap_elements(~testDispersion(mod_resids)) 

      testDispersion(mod_resids)
      testZeroInflation(mod_resids)
    }
  }
}

## Switch to zero-inflated model (intercept only)
{
  ## Data preparation
  {
    data <- full_data
  }
  ## glmmTMB
  {
    
    mod_glmmTMB <- glmmTMB(cbind(n.points, total.points-n.points) ~ 1 +
                             (1|AIMS_REEF_NAME) +
                             (1|Site) +
                             (1|Transect),
      ziformula =  ~ 1,
      data = data,
      family = "binomial", 
      REML = TRUE
    )

    summary(mod_glmmTMB)
    resids <- simulateResiduals(mod_glmmTMB, plot = TRUE)
    testDispersion(resids)
    testZeroInflation(resids)
    save(mod_glmmTMB, file = "../data/modelled/mod_glmmTMB_1.2.RData")
  }
  ## brms
  {
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
      silent =  0,
    )
    summary(mod_brm)
    save(mod_brm, file = "../data/modelled/mod_brm_1.2.RData")
    load(file = "../data/modelled/mod_brm_1.2.RData")
  }
}

## Switch to zero-inflated model (RE)
{
  ## glmmTMB
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
    resids <- simulateResiduals(mod_glmmTMB, plot = TRUE)
    testDispersion(resids)
    testZeroInflation(resids)
    save(mod_glmmTMB, file = "../data/modelled/mod_glmmTMB_1.3.RData")
  }
}

## Add Before/After
{
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
      mod_glmmTMB |>
        emmeans(~Dist.time, type = "response") |>
        as.data.frame() |> 
        ggplot(aes(y = prob, x = Dist.time)) +
        geom_pointrange(aes(ymin = asymp.LCL, ymax = asymp.UCL))
    }
    ## Contrasts
    {
      mod_glmmTMB |>
              emmeans(~Dist.time, type = "response") |>
              contrast(method = list(Dist.time = c(-1, 1))) |>
              summary(infer = TRUE) |> 
        as.data.frame() |> 
        ggplot(aes(x = odds.ratio, y = contrast)) +
              geom_vline(xintercept = 1, linetype = "dashed") +
        geom_pointrange(aes(xmin = asymp.LCL, xmax = asymp.UCL)) +
        scale_x_continuous("Effect (Before - After) on a fold scale", trans = "log2", breaks = scales::extended_breaks(n = 8))
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
      newdata_pred <- newdata |>
        bind_cols(mod$summary.fitted.values[i_newdata, ])
      newdata_pred 

      newdata_pred |> ggplot(aes(y = `0.5quant`, x = Dist.time)) +
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
      mod_glmmTMB |>
        emmeans(~ Dist.time * SecShelf, type = "response") |>
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
      mod_glmmTMB |>
        emmeans(~c, type = "response") |>
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
  ## glmmTMB
  {
    ## Fit model
    {
      mod_glmmTMB <- glmmTMB(
        cbind(n.points, total.points - n.points) ~ Dist.time + (s + c + d + b + u) +
          (1 | AIMS_REEF_NAME) +
          (1 | Site) +
          (1 | Transect),
        ziformula = ~ 1 + (1 | AIMS_REEF_NAME) + (1 | Site) + (1 | Transect),
        data = data,
        family = "binomial",
        REML = TRUE,
        control = glmmTMBControl(
          optimizer = "optim",
          optArgs = "BFGS"
          ## optArgs = "Nelder-Mead"
        )
      )

      summary(mod_glmmTMB)
      save(mod_glmmTMB, file = "../data/modelled/mod_glmmTMB_4.2.RData")
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
      newdata_pred |>
        ggplot(aes(y = `0.5quant`, x = factor(Dist), color = Dist.time)) +
        geom_pointrange(
          aes(
            ymin = `0.025quant`,
            ymax = `0.975quant`
          ),
          position = position_dodge(width = 0.5)
        )
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
  }
}



## SecShel, Dist.time and Disturbances (zi Intercept only) 
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
      mutate(across(c(s, c, d, b, u), \(x) ifelse(Dist.time == "Before", 0, x))) |>
      mutate(SecShelfP = paste(SecShelf, s, c, d, b, u))
  }
  ## glmmTMB
  {
    ## Fit model
    {
      mod_glmmTMB <- glmmTMB(cbind(n.points, total.points-n.points) ~ Dist.time + SecShelfP,# +
                               ## (1|AIMS_REEF_NAME) +
                               ## (1|Site) +
                               ## (1|Transect),
        ## ziformula = ~1,
        data = data,
        family = "binomial", 
        REML = TRUE
      )

      summary(mod_glmmTMB)
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
