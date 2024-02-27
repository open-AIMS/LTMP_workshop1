## load the data
load(file = "../data/primary/data.RData")

## ---- Q1_process_data
data <-
  data |>
  mutate(
    Site = paste(AIMS_REEF_NAME, SITE_NO),
    Transect = paste(Site, TRANSECT_NO),
    SecShelf = factor(paste(A_SECTOR, SHELF))
  )
## ----end

## Save the data
save(data, file = "../data/processed/data_q1.RData")
