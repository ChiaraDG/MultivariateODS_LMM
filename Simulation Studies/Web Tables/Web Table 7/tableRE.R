
#########################################################################################################
##################### Relative Efficiency Table #########################################################
#########################################################################################################

##################################################################################
# Scenario (a)
results  <- read.csv("./OffsetImputation/BivariateLMM/MIResults/SampleMoreMI.csv")
var.true <- results %>% group_by(Sampling) %>% select(matches("est.")) %>% 
  summarise_all(var) 
df.res   <- matrix(0, nrow = 5, ncol = 12)
# ODS + ACML
df.res <-  var.true[var.true$Sampling == "CompleteCaseMLE", -1]/
  var.true[var.true$Sampling == "ACML Slope", -1]
# BDS + ACML
df.res <-  rbind(df.res, var.true[var.true$Sampling == "CompleteCaseMLE", -1]/
  var.true[var.true$Sampling == "ACML SLP BLUP", -1])
# RS + MI
df.res <-  rbind(df.res, var.true[var.true$Sampling == "CompleteCaseMLE", -1]/
                   var.true[var.true$Sampling == "IIM RS", -1])
# ODS + MI
df.res <-  rbind(df.res, var.true[var.true$Sampling == "CompleteCaseMLE", -1]/
                   var.true[var.true$Sampling == "IIM ODS", -1])
# BDS + MI
df.res <-  rbind(df.res, var.true[var.true$Sampling == "CompleteCaseMLE", -1]/
                   var.true[var.true$Sampling == "IIM BDS", -1])


#### put results together ####
df.res <- cbind(Sampling = c("ODS + ACML", "BDS + ACML", "RS + MI", "ODS + MI", "BDS + MI"),
                  df.res)
##################################################################################

##################################################################################
# Scenario (b)
results <- read.csv("./OffsetImputation/BivariateLMM/MIResults/SampleMoreReallyHighConfGrpMI.csv")
var.true <- results %>% group_by(Sampling) %>% select(matches("est.")) %>% 
  summarise_all(var) 
df.res   <- matrix(0, nrow = 5, ncol = 12)
# ODS + ACML
df.res <-  var.true[var.true$Sampling == "CompleteCaseMLE", -1]/
  var.true[var.true$Sampling == "ACML Slope", -1]
# BDS + ACML
df.res <-  rbind(df.res, var.true[var.true$Sampling == "CompleteCaseMLE", -1]/
                   var.true[var.true$Sampling == "ACML SLP BLUP", -1])
# RS + MI
df.res <-  rbind(df.res, var.true[var.true$Sampling == "CompleteCaseMLE", -1]/
                   var.true[var.true$Sampling == "IIM RS", -1])
# ODS + MI
df.res <-  rbind(df.res, var.true[var.true$Sampling == "CompleteCaseMLE", -1]/
                   var.true[var.true$Sampling == "IIM ODS", -1])
# BDS + MI
df.res <-  rbind(df.res, var.true[var.true$Sampling == "CompleteCaseMLE", -1]/
                   var.true[var.true$Sampling == "IIM BDS", -1])
# put results together
df.res <- cbind(Sampling = c("ODS + ACML", "BDS + ACML", "RS + MI", "ODS + MI", "BDS + MI"),
                df.res)
##################################################################################


##################################################################################
# Scenario (c)
results <- read.csv("./OffsetImputation/BivariateLMM/MIResults/SampleMoreExtremeDesignMI.csv")
var.true <- results %>% group_by(Sampling) %>% select(matches("est.")) %>% 
  summarise_all(var) 
df.res   <- matrix(0, nrow = 5, ncol = 12)
# ODS + ACML
df.res <-  var.true[var.true$Sampling == "CompleteCaseMLE", -1]/
  var.true[var.true$Sampling == "ACML Slope", -1]
# BDS + ACML
df.res <-  rbind(df.res, var.true[var.true$Sampling == "CompleteCaseMLE", -1]/
                   var.true[var.true$Sampling == "ACML SLP BLUP", -1])
# RS + MI
df.res <-  rbind(df.res, var.true[var.true$Sampling == "CompleteCaseMLE", -1]/
                   var.true[var.true$Sampling == "IIM RS", -1])
# ODS + MI
df.res <-  rbind(df.res, var.true[var.true$Sampling == "CompleteCaseMLE", -1]/
                   var.true[var.true$Sampling == "IIM ODS", -1])
# BDS + MI
df.res <-  rbind(df.res, var.true[var.true$Sampling == "CompleteCaseMLE", -1]/
                   var.true[var.true$Sampling == "IIM BDS", -1])
# put results together
df.res <- cbind(Sampling = c("ODS + ACML", "BDS + ACML", "RS + MI", "ODS + MI", "BDS + MI"),
                df.res)
##################################################################################

##################################################################################
# Scenario (d)
results <- read.csv("./OffsetImputation/BivariateLMM/MIResults/SampleMoreHighConfTimeMI.csv")
var.true <- results %>% group_by(Sampling) %>% select(matches("est.")) %>% 
  summarise_all(var) 
df.res   <- matrix(0, nrow = 5, ncol = 12)
# ODS + ACML
df.res <-  var.true[var.true$Sampling == "CompleteCaseMLE", -1]/
  var.true[var.true$Sampling == "ACML Slope", -1]
# BDS + ACML
df.res <-  rbind(df.res, var.true[var.true$Sampling == "CompleteCaseMLE", -1]/
                   var.true[var.true$Sampling == "ACML SLP BLUP", -1])
# RS + MI
df.res <-  rbind(df.res, var.true[var.true$Sampling == "CompleteCaseMLE", -1]/
                   var.true[var.true$Sampling == "IIM RS", -1])
# ODS + MI
df.res <-  rbind(df.res, var.true[var.true$Sampling == "CompleteCaseMLE", -1]/
                   var.true[var.true$Sampling == "IIM ODS", -1])
# BDS + MI
df.res <-  rbind(df.res, var.true[var.true$Sampling == "CompleteCaseMLE", -1]/
                   var.true[var.true$Sampling == "IIM BDS", -1])
# put results together
df.res <- cbind(Sampling = c("ODS + ACML", "BDS + ACML", "RS + MI", "ODS + MI", "BDS + MI"),
                df.res)
