
load(file = "../data/processed/data_q1.RData")
full_data <- data

## Start with an intercept only model
{
  ## INLA
  {
    ## Data preparation
    {
      ## Focus on only the necessary variables
      data <- full_data |>
        dplyr::select(
          AIMS_REEF_NAME, Site, Transect, n.points,
          total.points
        )
      data.pred <- data
    }
    ## Fit model
    {
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
      wrap_elements(~testUniformity(mod_resids)) +
        wrap_elements(~plotResiduals(mod_resids)) +
        wrap_elements(~testDispersion(mod_resids)) 

      mod_resids |> testDispersion()
      mod_resids |> testZeroInflation()
      }
  }
  ## glmmTMB
  {
    library(glmmTMB)
    mod_glmmTMB <- glmmTMB(cbind(n.points, total.points-n.points) ~ 1 +
                             (1|AIMS_REEF_NAME) +
                             (1|Site) +
                             (1|Transect),
      data = data,
      family = "binomial", 
      REML = TRUE
    )
    resids <- simulateResiduals(mod_glmmTMB, plot = TRUE)
    resids |> testDispersion()
    resids |> testZeroInflation()
  }
}

## add REPORT_YEAR
{
  ## INLA
  {
    ## Prepare data
    {
      data <- full_data

      data <- data |>
              dplyr::select(
                AIMS_REEF_NAME, Site, Transect, n.points,
                total.points, cREPORT_YEAR
              )
      ## Cellmeans
      newdata <- data.frame(cREPORT_YEAR = sort(unique(data$cREPORT_YEAR)))

      data_pred <- data |> bind_rows(newdata)
      i_newdata <- (nrow(data) + 1):nrow(data_pred)
    }
    ## Fit model
    {
      mod <- inla(n.points ~ cREPORT_YEAR +
                    f(model = "iid", AIMS_REEF_NAME) +
                    f(model = "iid", Site) +
                    f(model = "iid", Transect),
        data = data_pred,
        Ntrials = data_pred$total.points,
        family = "zeroinflatedbinomial2", #"binomial", 
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
    ## Partial plots version 1
    {
      newdata_pred <- newdata |>
        bind_cols(mod$summary.fitted.values[i_newdata, ])
      newdata_pred 

      newdata_pred |>
        ggplot(aes(y = `0.5quant`,
          x = as.numeric(as.character(cREPORT_YEAR)))) +
        geom_line() +
                geom_pointrange(aes(
                        ymin = `0.025quant`,
                        ymax = `0.975quant`
                )) +
        scale_x_continuous("")
    }
    ## Partial plots version 2
    {
      newdata_fitted <- posterior_fitted.inla(mod, newdata = newdata)
      newdata_fitted <- newdata_fitted |>
              group_by(cREPORT_YEAR) |>
              summarise_draws(
                      median,
                      HDInterval::hdi
              )
      newdata_fitted 

      newdata_fitted |>
        ggplot(aes(y = median,
          x = as.numeric(as.character(cREPORT_YEAR)))) +
        geom_line() +
                geom_pointrange(aes(
                        ymin = lower,
                        ymax = upper
                )) +
        scale_x_continuous("")
    }
  }
  ## glmmTMB
  {
    ## Fit model
    {
      library(glmmTMB)
      mod_glmmTMB <- glmmTMB(cbind(n.points, total.points-n.points) ~ cREPORT_YEAR +
                               (1|AIMS_REEF_NAME) +
                               (1|Site) +
                               (1|Transect),
        ziformula = ~1,
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
        emmeans(~cREPORT_YEAR, type = "response") |>
        as.data.frame() |>
        ggplot(aes(y = prob, x = as.numeric(as.character(cREPORT_YEAR)))) +
        geom_line() +
        geom_pointrange(aes(ymin = asymp.LCL, ymax = asymp.UCL))
    }
  }
}


## add REPORT_YEAR and Sector/Shelf
{
  ## INLA
  {
    ## Prepare data
    {
      data <- full_data

      data <- data |>
        dplyr::select(AIMS_REEF_NAME, Site, Transect, n.points,
          total.points, cREPORT_YEAR, A_SECTOR, SHELF) |>
        mutate(
          SecShelfYr = paste(A_SECTOR, SHELF, cREPORT_YEAR),
          SecShelf = paste(A_SECTOR, SHELF)
        )
      ## what combinations are missing
      data |>
        group_by(SecShelfYr) |>
        summarise(Points = sum(n.points)) |>
        arrange(desc(Points))
      ## 492 missing - lets exclude those

      excludes <- data |>
        group_by(SecShelfYr) |>
        summarise(Points = sum(n.points)) |>
        filter(Points == 0) |>
        mutate(EXCLUDE = TRUE)
      data <- data |>
        left_join(excludes) |>
        filter(is.na(EXCLUDE)) |>
        droplevels() |>
        mutate(
          SecShelfYr = forcats::fct_relevel(SecShelfYr, "IN M 2016")
          ## SecShelfYr = forcats::fct_relevel(SecShelfYr, "SW M 2014")
        ) |>
        dplyr::select(-Points, -EXCLUDE)

      ## Cellmeans
      newdata <- crossing(SecShelfYr = data$SecShelfYr) |>
        separate(SecShelfYr,
          into = c("A_SECTOR", "SHELF", "cREPORT_YEAR"),
          remove = FALSE
        )

      data_pred <- data |> bind_rows(newdata)
      i_newdata <- (nrow(data) + 1):nrow(data_pred)
    }
    ## Fit model
    {
      mod <- inla(n.points ~ SecShelfYr +
                    f(model = "iid", AIMS_REEF_NAME) +
                    f(model = "iid", Site) +
                    f(model = "iid", Transect),
        data = data_pred,
        Ntrials = data_pred$total.points,
        family = "zeroinflatedbinomial1", #"zeroinflatedbetabinomial2", #"binomial", 
        control.predictor = list(link = 1, compute = TRUE),
        control.compute = list(config = TRUE,
          dic = TRUE, waic = TRUE, cpo = TRUE
        )#,
        ## control.fixed = control.fixed
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
      sum(mod$cpo$cop, na.omit = TRUE)
      -mean(log(na.omit(mod$cpo$cpo)))

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

      x <- list(
        y = na.omit(data$n.points / data$total.points),
        yrep = (preds)
      )
      class(x) <- "foo"
      pp_check.foo(x, type = "overlaid")
    }
    ## Partial plots
    {
      newdata_pred <- newdata |>
        bind_cols(mod$summary.fitted.values[i_newdata, ])
      newdata_pred

      ## Put the zeros back
      excludes <- excludes |>
        mutate(`0.5quant` = 0, `0.025quant` = 0, `0.975quant` = 0) |>
        dplyr::select(-Points, -EXCLUDE) |>
        separate(SecShelfYr,
          into = c("A_SECTOR", "SHELF", "cREPORT_YEAR"),
          remove = FALSE
        )

      newdata_pred <- newdata_pred |>
        full_join(excludes) |>
        mutate(A_SECTOR = factor(A_SECTOR,
          levels = c("CG", "PC", "CL", "CA", "IN", "TO", "WH", "PO", "SW", "CB")
        ))

      newdata_pred |>
        ggplot(aes(y = `0.5quant`,
          x = as.numeric(as.character(cREPORT_YEAR)))) +
        geom_line() +
        geom_pointrange(aes(ymin = `0.025quant`, ymax = `0.975quant`)) +
        facet_grid(A_SECTOR ~ SHELF, scales = "free")
    }
    ## Partial plots 2
    {
      newdata_fitted <- mod |>
        posterior_fitted.inla(newdata)

      newdata_fitted <- newdata_fitted |>
        group_by(A_SECTOR, SHELF, cREPORT_YEAR) |>
        dplyr::select(-SecShelfYr) |>
        posterior::summarise_draws(median,
          HDInterval::hdi)

      excludes <- excludes |>
        mutate(median = 0, lower = 0, upper = 0) |>
        dplyr::select(-any_of(c("Points", "EXCLUDE"))) |>
        separate(SecShelfYr,
          into = c("A_SECTOR", "SHELF", "cREPORT_YEAR"),
          remove = TRUE
        )

      newdata_fitted <- newdata_fitted |>
        full_join(excludes) |>
        mutate(A_SECTOR = factor(A_SECTOR,
          levels = c("CG", "PC", "CL", "CA", "IN", "TO", "WH", "PO", "SW", "CB")
        ))

      newdata_fitted |>
        ggplot(aes(y = median, x = as.numeric(as.character(cREPORT_YEAR)))) +
        geom_line() +
        geom_pointrange(aes(ymin = lower, ymax = upper)) +
        facet_grid(A_SECTOR ~ SHELF, scales = "free")
    }
  }
  ## glmmTMB
  {
    ## Fit model
    {
      library(glmmTMB)
      mod_glmmTMB <- glmmTMB(cbind(n.points, total.points-n.points) ~ SecShelfYr +
                               (1|AIMS_REEF_NAME) +
                               (1|Site) +
                               (1|Transect),
        ## ziformula = ~SecShelfYr,
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
      excludes <- excludes |>
        mutate(prob = 0, asymp.LCL = 0, asymp.UCL = 0) |>
        dplyr::select(-Points, -EXCLUDE) |>
        separate(SecShelfYr, into = c("A_SECTOR", "SHELF", "cREPORT_YEAR"), remove = FALSE) 
      mod_glmmTMB |>
        emmeans(~SecShelfYr, type = "response") |>
        as.data.frame() |>
        separate(SecShelfYr, into = c("A_SECTOR", "SHELF", "cREPORT_YEAR"), remove = FALSE) |>
        full_join(excludes) |> 
        mutate(A_SECTOR = factor(A_SECTOR, levels = c("CG", "PC", "CL", "CA", "IN", "TO", "WH", "PO", "SW", "CB"))) |> 
        ggplot(aes(y = prob, x = as.numeric(as.character(cREPORT_YEAR)))) +
        geom_line() +
        geom_pointrange(aes(ymin = asymp.LCL, ymax = asymp.UCL)) +
        facet_grid(A_SECTOR~SHELF, scales = "free")
    }
  }
}
