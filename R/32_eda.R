## Retrieve the data
load(file = "../data/processed/data_q2.R")


data |> dim()

## How many combinations have only 0's
data |>
  group_by(Dist.time, A_SECTOR, SHELF, Site, Transect) |>
  summarise(Points = sum(n.points)) |>
  filter(Points == 0)

## Raw means to get a sense of the data
data |>
  group_by(Dist.time, A_SECTOR, SHELF) |>
  summarise(cover = mean(n.points / total.points),
    cover_sd = sd(n.points / total.points),
    lower = cover - cover_sd,
    upper = cover + cover_sd) |>
  ggplot(aes(y = cover, x = A_SECTOR, colour = Dist.time)) +
  geom_pointrange(aes(ymin = lower, ymax = upper), position = position_dodge(width = 0.4)) +
  facet_grid(~SHELF)

data |>
  group_by(Dist.time, A_SECTOR, SHELF) |>
  summarise(cover = mean(n.points / total.points),
    cover_sd = sd(n.points / total.points),
    lower = cover - cover_sd,
    upper = cover + cover_sd) |>
  ggplot(aes(y = cover, x = Dist.time, colour = Dist.time)) +
  geom_hline(yintercept = 0) +
  geom_pointrange(aes(ymin = lower, ymax = upper),
    position = position_dodge(width = 0.4)) +
  facet_grid(A_SECTOR ~ SHELF, scales = "free")

data |>
  ggplot(aes(x = n.points)) +
  geom_histogram()

data |>
  ggplot(aes(x = n.points)) +
  geom_histogram() +
  facet_grid(SecShelf ~ Dist.time)

## Zero-inflation...
