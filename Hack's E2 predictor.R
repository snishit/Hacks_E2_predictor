rm(list=ls())
gc()

# Epigenetic prediction of 17β-estradiol
# Hack LM, Nishitani S, Knight AK, Kilaru V, Maddox SA, Seligowski AV, Jovanovic T, Ressler KJ, Smith AK*, Michopoulos V
# Epigenetic prediction of 17β-estradiol and relationship to trauma-related outcomes in women
# Compr Psychoneuroendocrinol, 2021, 9;6:100045.
# PMID: 35757356 PMCID: PMC9216622 DOI: 10.1016/j.cpnec.2021.100045

# GSE176394
# The .csv file "GSE176394_Novakovic_GAHT_beta_values.csv" used in this demonstration is here.
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE176394

library(data.table)
meth <- fread("GSE176394_Novakovic_GAHT_beta_values.csv", data.table = F)
meth <- meth[,-1]
row.names(meth) <- meth$probe_id
meth2 <- meth

library(wateRmelon)
meth2 <- meth2[,-1]
mage <- agep(meth2, coeff=NULL, method="horvath")
#write.csv(mage,"mage.csv")

sel_CpG <- read.csv("Hack_35probes.csv",header=T)
meth <- meth[meth$probe_id %in% sel_CpG$TargetID,]
datMethUsed <- t(meth[,-1])
datMethUsed <- as.data.frame(datMethUsed)
dat <- cbind(datMethUsed,mage$horvath.age)
names(dat)[36] <- "age"
write.csv(dat, "dat.csv", row.names = F)

# load RF model (Hack LM et al.,2021)
RF <- readRDS("FINAL_RF_70.30Split.rds")

library(randomForest)
RF_results <- predict(object = RF, dat[,c(1:36)])
RF_results <- as.data.frame(RF_results)
row.names(RF_results) <- row.names(dat)
x <- exp(RF_results$RF_results)
RF_results <- cbind(RF_results, x)
colnames(RF_results)[c(1:2)] <- c("LogE2","E2 (pg/mL)")	
write.csv(RF_results, "E2.csv")
