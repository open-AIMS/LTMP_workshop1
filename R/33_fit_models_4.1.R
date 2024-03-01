## Retrieve the data
load(file = "../data/processed/data_q2.R")


## Strip out all unnecessary variables

full_data <- data
data <- full_data

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
      save(cellmeans_summ_raw, file = "../data/modelled/cellmeans_summ_raw_4.1.RData")
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
      save(eff_summ_raw, file = "../data/modelled/eff_summ_raw_4.1.RData")
    }
  }
  ## glmmTMB
  {
    ## Fit model
    {
      if (rerun_models) {
        mod_glmmTMB <- glmmTMB(
          cbind(n.points, total.points - n.points) ~ s + c + d + b + u +
            (1 | AIMS_REEF_NAME) +
            (1 | Site) +
            (1 | Transect),
          ziformula = ~1,
          data = data,
          family = "binomial",
          REML = TRUE
        )

        summary(mod_glmmTMB)
        save(mod_glmmTMB, file = "../data/modelled/mod_glmmTMB_4.1.RData")
      }
    }
    ## DHARMa
    {
      load(file = "../data/modelled/mod_glmmTMB_4.1.RData")
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
      save(cellmeans_summ_glmmTMB, file = "../data/modelled/cellmeans_summ_glmmTMB_4.1.RData")
            ## a <-
            ##         mod_glmmTMB |>
            ##         emmeans(~ s + c + d + b + u, type = "response") |>
            ##         contrast(method = list(c(1, rep(0, 31))))

            ## plogis(fixef(mod_glmmTMB)[[1]] %*% t(Xmat))

            ## emmobj(a@bhat, a@V, levels = 1:6, linfct = Xmat) |>
            ##         as.data.frame() |>
            ##         mutate(across(c(estimate, asymp.LCL, asymp.UCL), plogis))
            ## qdrg(~ s + d + c + b + u, newdata, a@bhat, a@V)) |> summary()
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
      save(eff_summ_glmmTMB, file = "../data/modelled/eff_summ_glmmTMB_4.1.RData")
    }
    ## Partial plot - old method
    {
      if (1 != 2) {
        newdata <- crossing(s = c(0,1), d = c(0, 1), c = c(0, 1), b = c(0, 1), u = c(0, 1)) |>
          mutate(AIMS_REEF_NAME = NA, Site = NA, Transect = NA, Dist.number =  NA)
        newdata <- newdata |>
          filter((s == 0 & d == 0 & c == 0 & b == 0 & u == 0) |
                   (s == 1 & d == 0 & c == 0 & b == 0 & u == 0) |
                   (s == 0 & d == 1 & c == 0 & b == 0 & u == 0) |
                   (s == 0 & d == 0 & c == 1 & b == 0 & u == 0) |
                   (s == 0 & d == 0 & c == 0 & b == 1 & u == 0) |
                   (s == 0 & d == 0 & c == 0 & b == 0 & u == 1) 
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

        cellmeans_summ_glmmTMB <- newdata |> dplyr::select(Dist, Pred, lower, upper)
        ## save(cellmeans_summ_glmmTMB, file = "../data/modelled/cellmeans_summ_glmmTMB_4.1.RData")
        
        ## newdata |>
        ##   ggplot(aes(y = Pred, x = Dist, colour = factor(Dist.time))) +
        ##   geom_pointrange(aes(ymin = lower, ymax = upper),
        ##     position = position_dodge(width = 0.5)) 
      }
    }
    ## Contrasts - old method
    {
      if (1 != 2) {
        eff_summ_glmmTMB <- 
          cellmeans_summ_glmmTMB |>
          mutate(Values = Pred) |> 
          before_vs_afters()
        ## save(eff_summ_glmmTMB, file = "../data/modelled/eff_summ_glmmTMB_4.1.RData")
      }
    }
  }
  ## brms
  {
    ## Fit model
    {
      if (rerun_models) {
        form <- bf(
          n.points | trials(total.points) ~ s + c + d + b + u + 
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
    }
    ## Partial plot
    {
      load(file = "../data/modelled/mod_brm_4.1.RData")

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

      fit <- mod_brm |>
        as.data.frame() |>
        dplyr::select(matches("^b_[^z].*", perl = TRUE))

      cellmeans_brm <- cbind(newdata, t(as.matrix(fit) %*% t(Xmat))) |>
        pivot_longer(cols = matches("^[0-9]*$"), names_to = ".draws")
      cellmeans_summ_brm <- cellmeans_brm |> 
        dplyr::select(-.draws) |>
        group_by(Dist) |>
        mutate(value =  plogis(value)) |>
        summarise_draws(median, HDInterval::hdi)

      save(cellmeans_summ_brm, file = "../data/modelled/cellmeans_summ_brm_4.1.RData")
      ## cellmeans <- a |> plogis() |> summarise_draws(median, HDInterval::hdi)

      
      ## newdata <- crossing(s = c(0,1), c = c(0, 1), d = c(0, 1),  b = c(0, 1), u = c(0, 1)) |>
      ##   filter((s == 0 & d == 0 & c == 0 & b == 0 & u == 0) |
      ##            (s == 1 & d == 0 & c == 0 & b == 0 & u == 0) |
      ##            (s == 0 & d == 1 & c == 0 & b == 0 & u == 0) |
      ##            (s == 0 & d == 0 & c == 1 & b == 0 & u == 0) |
      ##            (s == 0 & d == 0 & c == 0 & b == 1 & u == 0) |
      ##            (s == 0 & d == 0 & c == 0 & b == 0 & u == 1) 
      ##   ) |> 
      ##   mutate(Dist = case_when(
      ##           s == 1 ~ "s",
      ##           c == 1 ~ "c",
      ##           d == 1 ~ "d",
      ##           b == 1 ~ "b",
      ##           u == 1 ~ "u",
      ##           .default = "Before"
      ##   ), Dist = factor(Dist, levels = c("Before", "s", "c", "d", "b", "u"))) |>
      ##   arrange(Dist) |> 
      ##   dplyr::select(-any_of(c("s", "c", "d", "b", "u")))
      ## Xmat <- cbind(1, contr.treatment(6, contrasts = TRUE))
      ## cellmeans <- emmobj(
      ##   post.beta = p,
      ##   levels = 1:6,
      ##   linfct = Xmat
      ## ) |>
      ##   as.data.frame() |>
      ##   mutate(across(c(estimate, asymp.LCL, asymp.UCL), plogis))

      ## cellmeans_summ_glmmTMB <- newdata |> bind_cols(cellmeans)

      
      ## newdata <- crossing(Dist.time = data$Dist.time) |>
      ##   crossing(s = c(0,1), d = c(0, 1), c = c(0, 1), b = c(0, 1), u = c(0, 1)) |>
      ##   mutate(AIMS_REEF_NAME = NA, Site = NA, Transect = NA, Dist.number =  NA)
      ## newdata <- newdata |>
      ##   filter((Dist.time == "Before" & s == 0 & d == 0 & c == 0 & b == 0 & u == 0) |
      ##            (Dist.time == "After" & s == 1 & d == 0 & c == 0 & b == 0 & u == 0) |
      ##            (Dist.time == "After" & s == 0 & d == 1 & c == 0 & b == 0 & u == 0) |
      ##            (Dist.time == "After" & s == 0 & d == 0 & c == 1 & b == 0 & u == 0) |
      ##            (Dist.time == "After" & s == 0 & d == 0 & c == 0 & b == 1 & u == 0) |
      ##            (Dist.time == "After" & s == 0 & d == 0 & c == 0 & b == 0 & u == 1)) |>
      ##   mutate(total.points = 1)

      ## p <- t(brms::posterior_epred(mod_brm, newdata = newdata, se.fit = TRUE))
      ## cellmeans_brm <-
      ##   newdata |>
      ##   cbind(p) |>
      ##   pivot_longer(cols = matches("^[0-9]*$"), names_to = ".draws") |>
      ##   mutate(Dist = case_when(
      ##     s == 1 ~ "s",
      ##     c == 1 ~ "c",
      ##     d == 1 ~ "d",
      ##     b == 1 ~ "b",
      ##     u == 1 ~ "u",
      ##     .default = "Before"
      ##   )) |>
      ##   dplyr::select(-any_of(c(
      ##     "s", "d", "c", "b", "u", "AIMS_REEF_NAME",
      ##     "Dist.time", "Site", "Transect", "Dist.number", "total.points"
      ##   )))

      ## cellmeans_summ_brm <- cellmeans_brm |>
      ##   dplyr::select(-.draws) |> 
      ##   group_by(Dist) |>
      ##   summarise_draws(median, HDInterval::hdi) |> 
      ##   dplyr::select(Dist, median, lower, upper)
      ## save(cellmeans_summ_brm, file = "../data/modelled/cellmeans_summ_brm_4.1.RData")
      
    }
    ## Contrasts
    {
      
      XXmat <- cbind(0, contr.treatment(6))[-1, ]

      eff_brm <- cbind(newdata[-1,], t(as.matrix(fit) %*% t(XXmat))) |>
        pivot_longer(cols = matches("^[0-9]*$"), names_to = ".draws")
      eff_summ_brm <- eff_brm |> 
        dplyr::select(-.draws) |>
        group_by(Dist) |>
        mutate(value =  exp(value)) |>
        summarise_draws(median, HDInterval::hdi)

      ## eff_brm <- cellmeans_brm |>
      ##   mutate(Values = value) |>
      ##   nest(.by = .draws) |>
      ##   mutate(eff = map(
      ##     .x = data,
      ##     .f = ~ before_vs_afters(.x)
      ##   )) |>
      ##   dplyr::select(-data) |>
      ##   unnest(c(eff))

      ## eff_summ_brm <- eff_brm |> 
      ##   dplyr::select(-.draws) |>
      ##   group_by(Dist) |> 
      ##   summarise_draws(
      ##     median,
      ##     HDInterval::hdi
      ##   )
      
      save(eff_summ_brm, file = "../data/modelled/eff_summ_brm_4.1.RData")
    }
  }
  ## INLA
  {
    ## Prepare data
    {
      data <- full_data

      data <- data |>
        dplyr::select(n.points, total.points, Dist.time, s, c, d, b, u,
          AIMS_REEF_NAME, Site, Transect, Dist.time) |> 
        mutate(across(c(s, c, d, b, u), \(x) ifelse(Dist.time == "Before", 0, x))) 
      ## data <- data |>
      ##         filter(!SecShelf %in% c("CG M", "PC M", "PC O")) |>
      ##         droplevels()
      ## Cellmeans
      newdata <- crossing(s = c(0,1), c = c(0, 1), d = c(0, 1), b = c(0, 1), u = c(0, 1)) |>
        mutate(AIMS_REEF_NAME = NA, Site = NA, Transect = NA, Dist.number =  NA)
      newdata <- newdata |>
        filter((s == 0 & d == 0 & c == 0 & b == 0 & u == 0) |
                 (s == 1 & d == 0 & c == 0 & b == 0 & u == 0) |
                 (s == 0 & d == 1 & c == 0 & b == 0 & u == 0) |
                 (s == 0 & d == 0 & c == 1 & b == 0 & u == 0) |
                 (s == 0 & d == 0 & c == 0 & b == 1 & u == 0) |
                 (s == 0 & d == 0 & c == 0 & b == 0 & u == 1)) |>
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
      if (rerun_models) {
        mod <- inla(
          n.points ~ s + c + d + b + u +
            f(model = "iid", AIMS_REEF_NAME) +
            f(model = "iid", Site) +
            f(model = "iid", Transect),
          data = data_pred,
          Ntrials = data_pred$total.points,
          family = "zeroinflatedbinomial1", # "binomial",
          control.predictor = list(link = 1, compute = TRUE),
          control.compute = list(
            config = TRUE,
            dic = TRUE, waic = TRUE, cpo = TRUE,
            return.marginals.predictor = TRUE
          )
        )
        summary(mod)
        ## autoplot(mod)
        save(mod, file = "../data/modelled/mod_4.1.RData")
      }
    }
    ## Partial plots - version 1
    {
      load(file = "../data/modelled/mod_4.1.RData")
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

if (1 != 2) {
  
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
    ## Cellmeans
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
      save(eff_summ_raw, file = "../data/modelled/eff_summ_raw_4.1.RData")
    }
  }
  ## glmmTMB
  {
    ## Fit model
    {
      if (rerun_models) {
        mod_glmmTMB <- glmmTMB(
          cbind(n.points, total.points - n.points) ~ Dist.time + (s + c + d + b + u) +
            (1 | AIMS_REEF_NAME) +
            (1 | Site) +
            (1 | Transect),
          ziformula = ~1,
          data = data,
          family = "binomial",
          REML = TRUE
        )

        summary(mod_glmmTMB)
        save(mod_glmmTMB, file = "../data/modelled/mod_glmmTMB_4.1.RData")
      }
    }
    ## DHARMa
    {
      load(file = "../data/modelled/mod_glmmTMB_4.1.RData")
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
      if (rerun_models) {
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
        dplyr::select(Dist, median, lower, upper)
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
      if (rerun_models) {
        mod <- inla(
          n.points ~ Dist.time + (s + c + d + b + u) +
            f(model = "iid", AIMS_REEF_NAME) +
            f(model = "iid", Site) +
            f(model = "iid", Transect),
          data = data_pred,
          Ntrials = data_pred$total.points,
          family = "zeroinflatedbinomial2", # "binomial",
          control.predictor = list(link = 1, compute = TRUE),
          control.compute = list(
            config = TRUE,
            dic = TRUE, waic = TRUE, cpo = TRUE,
            return.marginals.predictor = TRUE
          )
        )
        summary(mod)
        ## autoplot(mod)
        save(mod, file = "../data/modelled/mod_4.1.RData")
      }
    }
    ## Partial plots - version 1
    {
      load(file = "../data/modelled/mod_4.1.RData")
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

}
