library(here)
library(dplyr)
spp_hist <- read.csv(here("data/spp_hist.csv"), skip = 1)
spp_hist %>% arrange(desc(Grand.Total)) %>%
  head()
str(spp_hist)
spp_hist <- spp_hist[-nrow(spp_hist),]
hist(spp_hist$Grand.Total)

head(spp_hist)

