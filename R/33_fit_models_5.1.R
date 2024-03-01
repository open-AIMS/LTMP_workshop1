## Retrieve the data
load(file = "../data/processed/data_q2.R")


## Strip out all unnecessary variables

full_data <- data
data <- full_data

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
    ## Cellmeans
    {
      cellmeans_summ_raw <- data |>
        filter((s == 0 & d == 0 & c == 0 & b == 0 & u == 0) |
                 (s == 1 & d == 0 & c == 0 & b == 0 & u == 0) |
                 (s == 0 & d == 1 & c == 0 & b == 0 & u == 0) |
                 (s == 0 & d == 0 & c == 1 & b == 0 & u == 0) |
                 (s == 0 & d == 0 & c == 0 & b == 1 & u == 0) |
                 (s == 0 & d == 0 & c == 0 & b == 0 & u == 1) 
        ) |> 
        mutate(Dist = case_when(s == 1 ~ "s",
          c == 1 ~ "c",
          d == 1 ~ "d",
          b == 1 ~ "b",
          u == 1 ~ "u",
          .default = "Before")) |>
        dplyr::select(-any_of(c("s", "c", "d", "b", "u"))) |> 
        group_by(Dist) |> 
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
    }
    ## Effects
    {
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
  }
  ## glmmTMB
  {
    ## Fit model
    {
      if (rerun) {
        mod_glmmTMB <- glmmTMB(cbind(n.points, total.points-n.points) ~ s+c+d+b+u +
                                 (1|AIMS_REEF_NAME) +
                                 (1|Site) +
                                 (1|Transect),
          ziformula = ~ 1 + (1|Site) + (1|Transect), #(1|AIMS_REEF_NAME) + (1|Site) + (1|Transect),
          data = data,
          family = "binomial", 
          REML = TRUE
        )

        summary(mod_glmmTMB)
        save(mod_glmmTMB, file = "../data/modelled/mod_glmmTMB_5.1.RData")
      }
    }
    ## DHARMa
    {
      load(file = "../data/modelled/mod_glmmTMB_5.1.RData")
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

      newdata <- crossing(s = c(0,1), c = c(0, 1), d = c(0, 1),  b = c(0, 1), u = c(0, 1)) |>
        filter((s == 0 & d == 0 & c == 0 & b == 0 & u == 0) |
                 (s == 1 & d == 0 & c == 0 & b == 0 & u == 0) |
                 (s == 0 & d == 1 & c == 0 & b == 0 & u == 0) |
                 (s == 0 & d == 0 & c == 1 & b == 0 & u == 0) |
                 (s == 0 & d == 0 & c == 0 & b == 1 & u == 0) |
                 (s == 0 & d == 0 & c == 0 & b == 0 & u == 1) 
        ) |> 
        mutate(Dist = case_when(
                s == 1 ~ "s",
                c == 1 ~ "c",
                d == 1 ~ "d",
                b == 1 ~ "b",
                u == 1 ~ "u",
                .default = "Before"
        ), Dist = factor(Dist, levels = c("Before", "s", "c", "d", "b", "u"))) |>
        arrange(Dist) |> 
        dplyr::select(-any_of(c("s", "c", "d", "b", "u")))
      Xmat <- cbind(1, contr.treatment(6, contrasts = TRUE))
      cellmeans <- emmobj(
        fixef(mod_glmmTMB)[[1]],
        vcov(mod_glmmTMB)[[1]],
        levels = 1:6,
        linfct = Xmat
      ) |>
        as.data.frame() |>
        mutate(across(c(estimate, asymp.LCL, asymp.UCL), plogis))

      cellmeans_summ_glmmTMB <- newdata |> bind_cols(cellmeans)
      save(cellmeans_summ_glmmTMB, file = "../data/modelled/cellmeans_summ_glmmTMB_5.1.RData")
    }
    ## Contrasts
    {
      XXmat <- cbind(0, contr.treatment(6))[-1, ]
      eff <- emmobj(
        bhat = fixef(mod_glmmTMB)[[1]],
        V = vcov(mod_glmmTMB)[[1]],
        levels = 2:6,
        linfct = XXmat
      ) |> 
        as.data.frame() |>
        mutate(across(c(estimate, asymp.LCL, asymp.UCL), exp))
      eff_summ_glmmTMB <- newdata |>
              filter(Dist != "Before") |>
        droplevels() |>
        bind_cols(eff)
      save(eff_summ_glmmTMB, file = "../data/modelled/eff_summ_glmmTMB_5.1.RData")
    }
    ## Partial plots - old
    {
      if (1 == 2) {
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
    }
    ## Contrasts - old
    {
      if (i == 2) {
        eff_summ_glmmTMB <- 
          cellmeans_summ_glmmTMB |>
          mutate(Values = Pred) |> 
          before_vs_afters()
        save(eff_summ_glmmTMB, file = "../data/modelled/eff_summ_glmmTMB_5.1.RData")
      }
    }
  }
  ## brms
  {
    ## Fit model
    {
      if (rerun) {
        form <- bf(
          n.points | trials(total.points) ~ s + c + d + b + u +
            (1 | AIMS_REEF_NAME) +
            (1 | Site) +
            (1 | Transect),
          zi = ~ 1 + (1 | Site) + (1 | Transect),
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
        dplyr::rename(median = estimate, lower = asymp.LCL, upper = asymp.UCL),
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
        dplyr::rename(median = estimate, lower = asymp.LCL, upper = asymp.UCL),
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
