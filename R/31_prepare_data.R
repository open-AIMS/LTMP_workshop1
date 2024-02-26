
## load the data
load(file = "../data/primary/data.RData")
load(file = "../data/primary/lookup.RData")

## join the lookup to the data
data |> colnames()
lookup |> colnames()

## What fields do they have in common
colnames(data)[colnames(data) %in% colnames(lookup)]

## But some fields have changed names!!

data |>
  left_join(lookup,
    by = c("A_SECTOR", "FULLREEF_ID", "REPORT_YEAR" = "year")
  )

## Pay attention to that error message

data |>
  slice(4716) |>
  dplyr::select(P_CODE, A_SECTOR, SHELF, FULLREEF_ID, REPORT_YEAR)

## Ah - multiple disturbances in a single year
lookup |>
  filter(A_SECTOR == "CA", FULLREEF_ID == "16064S", year == 2012) |>
  as.data.frame()

## lets collapse these together and dummy code the disturbance types
lookup |> dim()

## We just want Before and After, and each needs to have the same DISTURBANCE_TYPE

## 1. Remove cases where:
##    - DISTURBANCE_TYPE is NA
##    - Dist.number is 0 (sequence without a known Before)
##    - Dist.type is "Pre"
##    - Dist.time is "Recovery"
## 2. Remove sequences (Dist.number) within REEF that do not have a
##    Before and an After
## 3. Create a DIST_TYPE that is the collapsed concatenation of the
##    non "n" and non NA DISTURBANCE_TYPE values within each sequence
## 4. Use pivot wider to create dummy codes for each disturbance type
## lookup1 <- lookup
## lookup <- lookup1

lookup |>
        filter(REEF == "Pompey Reef No.1") |>
        as.data.frame()
## - Dist.number 1, does not have a Before
## - Dist.time == "Before", Dist.number == 2 has DISTURBANCE_TYPE of "n"
## - Before (n), During (s), After (u)
lookup |>
        filter(REEF == "Arlington Reef") |>
        as.data.frame()
## - Dist.time == "Before" for Dist.number == 1 is a "n"
## - Dist.time == "After" for Dist.number == 1 has two rows
## - Dist.number == 2, Before, During and After all have different DISTURBANCE_TYPE (n, d, b)

lookup <-
  lookup |>
  ## remove cases without disturbances
  filter(!is.na(DISTURBANCE_TYPE),
    Dist.number > 0,
    Dist.time != "Recovery",
    Dist.type != "Pre") |>
  ## within each REEF and Dist.number, assign a DISTURBANCE_TYPE for Before to be the same as After
  ## collapse disturbances
  group_by(FULLREEF_ID, REEF, Dist.number) |>
  filter(sum(str_detect(Dist.time, "Before|After")) == 2) |>
  nest() |>
  # filter(REEF == "Agincourt Reef No.1") |>
  mutate(data1 = purrr::map(
    .x = data,
    .f = ~ {
      # print(.x)
      Before <- .x |>
        filter(Dist.time == "Before") |>
        dplyr::select(-DISTURBANCE_TYPE)
      After <- .x |>
        filter(Dist.time == "After") |>
        dplyr::select(-DISTURBANCE_TYPE)
      Dists <- .x |>
        pull(DISTURBANCE_TYPE) |>
        unique()
      Dists <- Dists[Dists != "n"]
      Before <- Before |> reframe(Before, DISTURBANCE_TYPE = Dists)
      After <- After |> reframe(After, DISTURBANCE_TYPE = Dists)
      bind_rows(Before, After)
      # print("hi")
      # print(as.data.frame(a))
      ## asdf
    }
  )) |>
  dplyr::select(-data) |> 
  unnest(c(data1)) |>
  pivot_wider(id_cols = everything(),
    values_fn = \(x) sum(!is.na(x)),
    values_fill = 0,
    names_from = DISTURBANCE_TYPE, values_from = DISTURBANCE_TYPE) 

## Try the join again

data <-
  data |>
  left_join(lookup,
    by = c("A_SECTOR", "FULLREEF_ID", "REPORT_YEAR" = "year")
  )

## Final prep is to add some derived variables

data <-
  data |>
  mutate(
    Site = paste(AIMS_REEF_NAME, SITE_NO),
    Transect = paste(Site, TRANSECT_NO)
  ) |>
  filter(Dist.time %in% c("Before", "After")) |>
  filter(!is.na(Dist.number)) |>
  droplevels() |> 
  mutate(Dist.time = factor(Dist.time, levels = c("Before", "After"))) |>
  mutate(SecShelf = factor(paste(A_SECTOR, SHELF)))

save(data, file = "../data/processed/data_q2.R")

data |>
        group_by(Dist.number) |>
        count()
