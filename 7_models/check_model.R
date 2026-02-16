library(brms)
library(StanHeaders)
db_mod_bundance <- readR
DS("C:/Users/PC/Dropbox/CO_DBdata/db_mod_abundance.rds")
db1 <- as.data.frame(db_mod_bundance)

plot(db_mod_bundance)
