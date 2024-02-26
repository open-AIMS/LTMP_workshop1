
## Retrieve the data
load(file = "../data/processed/data_q1.RData")

## Explore the spatio-temporal domain using simple means and CI's
## Note in particular, the CI's are completely invalid
data |>
  group_by(A_SECTOR, SHELF, cREPORT_YEAR) |>
  summarise(across(cover, list(mu = mean, med = median,
    sd = sd, n = length))) |>
  mutate(
    lower = cover_mu - 2 * cover_sd / sqrt(cover_n),
    upper = cover_mu + 2 * cover_sd / sqrt(cover_n)
  ) |> 
  mutate(A_SECTOR = factor(A_SECTOR, levels = c(
    "CG", "PC", "CL", "CA", "IN", "TO", "WH", "PO", "SW", "CB"
  ))) |>
  ggplot(aes(y = cover_mu, x = as.numeric(as.character(cREPORT_YEAR)))) +
  geom_line() +
  geom_pointrange(aes(ymin = lower, ymax = upper)) +
  facet_grid(A_SECTOR ~ SHELF, scales = "free")

## Conclusions:
## - substantial combinations with all 0's
## - this could be highly problematic

## Explore which Sector/Shelf/Year combinations have only 0 counts
data |>
  group_by(A_SECTOR, SHELF, cREPORT_YEAR) |>
  summarise(Points = sum(n.points)) |>
  filter(Points == 0)

## 66 Sector/Shelf/Years with all 0's It is likely we will need to
## exclude these from the analyses and put them back in at the end -
## they will destabalise the model without very specific and strong
## priors

## these are likely to be very problematic

## Lets explore distributions over the space
data |>
  group_by(A_SECTOR, SHELF) |>
  ggplot(aes(x = n.points)) +
  geom_histogram() +
  facet_grid(A_SECTOR ~ SHELF)

## Could be zero-inflated - more so in some sector/shelf

data |>
  group_by(cREPORT_YEAR) |>
  ggplot(aes(x = n.points)) +
  geom_histogram() +
  facet_wrap(~cREPORT_YEAR)

## More zero-inflation
