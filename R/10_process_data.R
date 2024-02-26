
## Read in the primary data
data <- get(load("../data/primary/seriatopora.RData"))

data |> glimpse()

## We restrict to just the Seriatopora (G_SER)
data <- data |>
        filter(COMP_2021 == "G_SER") |>
        droplevels()

save(data, file = "../data/primary/data.RData")

## read in the disturbance lookup
lookup <- read_csv("../data/primary/dist.lookup 5.csv")

save(lookup, file = "../data/primary/lookup.RData")
