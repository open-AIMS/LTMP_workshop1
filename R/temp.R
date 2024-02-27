






data |>
        filter(Dist.time %in% c("Before", "After")) |>
        group_by(Dist.number) |>
        count()

dat <-
  data |>
  filter(Dist.time %in% c("Before", "After"), !is.na(Dist.number)) |>
  mutate(
    Site = paste(AIMS_REEF_NAME, SITE_NO),
    Transect = paste(Site, TRANSECT_NO)) |> 
  droplevels() |>
  mutate(Dist.time = factor(Dist.time, levels = c("Before", "After"))) |>
  ## filter(!A_SECTOR %in% c("IN", "PC")) |>
  ## droplevels() |>
  mutate(SecShelf = factor(paste(A_SECTOR, SHELF))) |>
  mutate(presence = as.integer(n.points > 0))

## Exploratory data analysis

dat |>
        group_by(Dist.time, SecShelf) |>
        summarise(cover = mean(n.points / total.points)) |>
        ggplot(aes(y = cover, x = SecShelf, colour = Dist.time)) +
        geom_point()

dat |>
        ggplot(aes(x = n.points)) +
  geom_histogram() +
  facet_grid(SecShelf~Dist.time)

newdata <- data.frame(Dist.time = factor(c("Before", "After"),
  levels = c("Before", "After"))) |>
  crossing(A_SECTOR = dat$A_SECTOR, SHELF = dat$SHELF) |>
  mutate(SecShelf = factor(paste(A_SECTOR, SHELF)))

dat.pred <- dat |> bind_rows(newdata)
i.newdata <- (nrow(dat) +1):nrow(dat.pred)

## mod <- inla(cover ~ Dist.number*Dist.time * A_SECTOR*SHELF +
## ## mod <- inla(cover ~ Dist.time * A_SECTOR +
##               f(model = "iid", AIMS_REEF_NAME) +
##               f(model = "iid", Site) +
##               f(model = "iid", Transect),
##         data = dat.pred,
##         ## control.predictor = list(link = 1, compute = TRUE),
##         control.compute = list(config = TRUE,
##                     dic = TRUE, waic = TRUE, cpo = TRUE
##                     )
## )

## summary(mod)

## newdata |> bind_cols(mod$summary.fitted.values[i.newdata, ])

Xmat <- model.matrix(~Dist.time*SecShelf, data=newdata)
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

## mod <- inla(n.points ~ Dist.number + Dist.time * SecShelf +
mod <- inla(presence ~ Dist.number + Dist.time * SecShelf +
## mod <- inla(cover ~ Dist.time * A_SECTOR +
              f(model = "iid", REPORT_YEAR) +
              f(model = "iid", AIMS_REEF_NAME) +
              f(model = "iid", Site) +
              f(model = "iid", Transect),
  data = dat.pred,
  #Ntrials = dat.pred$total.points,
lincomb = lincomb,
  family = "binomial", #"betabinomial", #"binomial" "zeroinflatedbinomial" 
        control.predictor = list(link = 1, compute = TRUE),
        control.compute = list(config = TRUE,
                    dic = TRUE, waic = TRUE, cpo = TRUE
                    )
)

summary(mod)

newdata.pred <- newdata |> bind_cols(mod$summary.fitted.values[i.newdata, ])
newdata.pred 

newdata.pred |> ggplot(aes(y = `0.5quant`, x = A_SECTOR, colour = Dist.time)) +
newdata.pred |> ggplot(aes(y = `0.5quant`, x = SecShelf, colour = Dist.time)) +
  geom_point() 
  ## geom_pointrange(aes(ymin = `0.025quant`, ymax = `0.975quant`)) +
  #facet_grid(~SHELF)

nd.pred <- nd1 |>
        bind_cols(mod$summary.lincomb.derived) |>
        mutate(across(where(is.numeric), exp))
nd.pred |> ggplot(aes(x = `0.5quant`, y = A_SECTOR)) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  geom_pointrange(aes(xmin = `0.025quant`, xmax = `0.975quant`)) +
  facet_grid(SHELF~.)

###

newdata <- data.frame(Dist.time = factor(c("Before", "After"),
  levels = c("Before", "After"))) |>
  crossing(A_SECTOR = dat$A_SECTOR, SHELF = dat$SHELF) |>
          ## s =  c(0,1), b = c(0,1), c = c(0,1)) |>
          mutate(SecShelf = factor(paste(A_SECTOR, SHELF))) |>
          crossing(disturb = factor(c("s", "b", "c"))) |>
          mutate(
                  s = ifelse(disturb == "s", 1, 0),
                  c = ifelse(disturb == "c", 1, 0),
                  b = ifelse(disturb == "b", 1, 0)
          ) |>
  dplyr::select(-disturb)
    

dat.pred <- dat |>
        bind_rows(newdata) |>
        ## mutate(s = factor(s), b = factor(b), c = factor(c))
        mutate(s = ifelse(s == 0, 0, 1), b = ifelse(b == 0, 0, 1), c = ifelse(c == 0, 0, 1))
i.newdata <- (nrow(dat) +1):nrow(dat.pred)
mod <- inla(n.points ~ Dist.number*Dist.time * SecShelf + (s+b+c) +
## mod <- inla(cover ~ Dist.time * A_SECTOR +
              f(model = "iid", REPORT_YEAR) +
              f(model = "iid", AIMS_REEF_NAME) +
              f(model = "iid", Site) +
              f(model = "iid", Transect),
  data = dat.pred,
  Ntrials = dat.pred$total.points,
lincomb = lincomb,
  family = "binomial", #"betabinomial", #"binomial" "zeroinflatedbinomial" 
        control.predictor = list(link = 1, compute = TRUE),
        control.compute = list(config = TRUE,
                    dic = TRUE, waic = TRUE, cpo = TRUE
                    )
)

summary(mod)

newdata.pred <- newdata |>
        bind_cols(mod$summary.fitted.values[i.newdata, ]) |>
        pivot_longer(cols = c(s, b, c), names_to = "disturb") |>
        filter(!mean == 0, value != 0, `0.975quant` < 0.5)
newdata.pred |>
        as.data.frame() |>
        head()

newdata.pred |> ggplot(aes(y = `0.5quant`, x = A_SECTOR, colour = Dist.time)) +
  geom_pointrange(aes(ymin = `0.025quant`, ymax = `0.975quant`)) +
  facet_grid(disturb~SHELF)

nd.pred <- nd1 |>
        bind_cols(mod$summary.lincomb.derived) |>
        mutate(across(where(is.numeric), exp))
nd.pred |> ggplot(aes(x = `0.5quant`, y = A_SECTOR)) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  geom_pointrange(aes(xmin = `0.025quant`, xmax = `0.975quant`)) +
  facet_grid(SHELF~.)





library(brms)
mod_brm <- brm(bf(n.points | trials(total.points) ~ SecShelfYr +
                         (1|AIMS_REEF_NAME) +
                         (1|Site) +
                         (1|Transect)),
  data = data,
  family = "binomial",
  iter =  1000,
  warmup =  500,
  chains =  3, cores =  3,
  thin = 5,
  backend = "cmdstanr"
)






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

data |>
        group_by(Dist.time, s, b, c, d, u) |>
        summarise(A = mean(n.points / total.points))
      
      mod <- inla(I(n.points/total.points) ~ s + b + c + d + u +
                    f(model = "iid", Dist.number) +
                    f(model = "iid", AIMS_REEF_NAME) +
                    f(model = "iid", Site) +
                    f(model = "iid", Transect),
        data = data_pred,
        ## lincomb =  lincomb,
        family = "gaussian", #"binomial" 
        control.predictor = list(link = 1, compute = TRUE),
        control.compute = list(config = TRUE,
          dic = TRUE, waic = TRUE, cpo = TRUE
        )
      )
      summary(mod)
      mod <- inla(n.points ~ s + b + c + d + u +
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

      library(glmmTMB)
      mod_g <- glmmTMB(I(n.points/total.points) ~ s + b + c + d + u +
                               (1|AIMS_REEF_NAME) +
                               (1|Site) +
                               (1|Transect),
        ziformula = ~1,
        data = data,
        family = "binomial", 
        REML = TRUE
      )

      summary(mod_g)


      mod_b <- brm(bf(n.points | trials(total.points) ~  s + b + c + d + u +
                        (1|AIMS_REEF_NAME) +
                        (1|Site) +
                        (1|Transect)
      ),
      data = data,
      family = "binomial", 
      iter =  1000, warmup =  500,
      chains = 3, cores = 3,
      thin = 5,
      backend = "cmdstanr",
      control = list(adapt_delta = 0.99)
      )
      summary(mod_b)
      mod_b <- brm(bf(n.points | trials(total.points) ~  s + b + c + d + u +
                        (1|AIMS_REEF_NAME) +
                        (1|Site) +
                        (1|Transect),
        zi =  ~ 1 +(1|AIMS_REEF_NAME) +
                        (1|Site) +
                        (1|Transect)
      ),
      data = data,
      family = "zero_inflated_binomial", 
      iter =  1000, warmup =  500,
      chains = 3, cores = 3,
      thin = 5,
      backend = "cmdstanr",
      control = list(adapt_delta = 0.99)
      )
      summary(mod_b)

mod_b |>
        emmeans(~ s + b + c+ d + u, type = "response") |>
        as.data.frame() |>
        ggplot(aes(y = prob, x = s, colour = factor(b))) +
        geom_point()
      
      
      mod_g <- glmmTMB(cbind(n.points, total.points - n.points) ~ s + b + c + d + u +
                               (1|AIMS_REEF_NAME) +
                               (1|Site) +
                               (1|Transect),
        ziformula = ~ 1, #s + b + c + d + u,
        data = data,
        family = "binomial", 
        REML = TRUE
      )

      summary(mod_g)
      AIC(mod_g)

      resids <- simulateResiduals(mod_g, plot = TRUE)
      wrap_elements(~testUniformity(resids)) +
        wrap_elements(~plotResiduals(resids)) +
        wrap_elements(~testDispersion(resids)) 
      testDispersion(resids)
      testZeroInflation(resids)
      
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
