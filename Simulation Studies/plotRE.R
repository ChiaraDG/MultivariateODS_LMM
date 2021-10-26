#########################################################################################################
##################### Relative Efficiency Plot ##########################################################
#########################################################################################################

# Code to generate figure 2 after running the simulation for all the scenarios in the paper

library(tidyverse)
library(latex2exp)
library(ggpubr)
library(viridis)

# Main Scenario
results  <- read.csv("./OffsetImputation/BivariateLMM/MIResults/SampleMoreMI.csv")
var.true <- results %>% group_by(Sampling) %>% select(matches("est.")) %>% summarise_all(var) %>% 
  select(Sampling, est.tg.Y1, est.tg.Y2) 
# compute relative efficiency for the first and the second oucome
var.true[,"releff.Y1"] <-  as.numeric(var.true[var.true$Sampling == "CompleteCaseMLE", "est.tg.Y1"])/
  var.true[,"est.tg.Y1"]
var.true[,"releff.Y1"] <- round(var.true[,"releff.Y1"], 2)
var.true[,"releff.Y2"] <-  as.numeric(var.true[var.true$Sampling == "CompleteCaseMLE", "est.tg.Y2"])/
  var.true[,"est.tg.Y2"]
var.true[,"releff.Y2"] <- round(var.true[,"releff.Y2"], 2)
# get the variabled of relative efficiency we need
var.true  <- var.true %>% select(Sampling, releff.Y1, releff.Y2)
releff.main <- data.frame(Scenario = "(a)", var.true)

# delta = -2
results <- read.csv("./OffsetImputation/BivariateLMM/MIResults/SampleMoreReallyHighConfGrpMI.csv")
var.true <- results %>% group_by(Sampling) %>% select(matches("est.")) %>% summarise_all(var) %>% 
  select(Sampling, est.tg.Y1, est.tg.Y2) 
# compute relative efficiency for the first and the second oucome
var.true[,"releff.Y1"] <-  as.numeric(var.true[var.true$Sampling == "CompleteCaseMLE", "est.tg.Y1"])/
  var.true[,"est.tg.Y1"]
var.true[,"releff.Y1"] <- round(var.true[,"releff.Y1"], 2)
var.true[,"releff.Y2"] <-  as.numeric(var.true[var.true$Sampling == "CompleteCaseMLE", "est.tg.Y2"])/
  var.true[,"est.tg.Y2"]
var.true[,"releff.Y2"] <- round(var.true[,"releff.Y2"], 2)
# get the variabled of relative efficiency we need
var.true  <- var.true %>% select(Sampling, releff.Y1, releff.Y2)
releff.reallyhighgrpconf <- data.frame(Scenario = "(b)", var.true)

# Extreme design
results <- read.csv("./OffsetImputation/BivariateLMM/MIResults/SampleMoreExtremeDesignMI.csv")
var.true <- results %>% group_by(Sampling) %>% select(matches("est.")) %>% summarise_all(var) %>% 
  select(Sampling, est.tg.Y1, est.tg.Y2) 
# compute relative efficiency for the first and the second oucome
var.true[,"releff.Y1"] <-  as.numeric(var.true[var.true$Sampling == "CompleteCaseMLE", "est.tg.Y1"])/
  var.true[,"est.tg.Y1"]
var.true[,"releff.Y1"] <- round(var.true[,"releff.Y1"], 2)
var.true[,"releff.Y2"] <-  as.numeric(var.true[var.true$Sampling == "CompleteCaseMLE", "est.tg.Y2"])/
  var.true[,"est.tg.Y2"]
var.true[,"releff.Y2"] <- round(var.true[,"releff.Y2"], 2)
# get the variabled of relative efficiency we need
var.true  <- var.true %>% select(Sampling, releff.Y1, releff.Y2)
releff.extreme <- data.frame(Scenario = "(c)", var.true)

# high time*conf
results <- read.csv("./OffsetImputation/BivariateLMM/MIResults/SampleMoreHighConfTimeMI.csv")
var.true <- results %>% group_by(Sampling) %>% select(matches("est.")) %>% summarise_all(var) %>% 
  select(Sampling, est.tg.Y1, est.tg.Y2) 
# compute relative efficiency for the first and the second oucome
var.true[,"releff.Y1"] <-  as.numeric(var.true[var.true$Sampling == "CompleteCaseMLE", "est.tg.Y1"])/
  var.true[,"est.tg.Y1"]
var.true[,"releff.Y1"] <- round(var.true[,"releff.Y1"], 2)
var.true[,"releff.Y2"] <-  as.numeric(var.true[var.true$Sampling == "CompleteCaseMLE", "est.tg.Y2"])/
  var.true[,"est.tg.Y2"]
var.true[,"releff.Y2"] <- round(var.true[,"releff.Y2"], 2)
# get the variabled of relative efficiency we need
var.true  <- var.true %>% select(Sampling, releff.Y1, releff.Y2)
releff.conftime <- data.frame(Scenario = "(d)", var.true)

# put scenarios together
tmp          <- rbind(releff.main, releff.reallyhighgrpconf, releff.extreme, releff.conftime)
# for the higher conf*time, the naming of bds acl and ods acl are swapped
tmp$Sampling <- c(rep(c("ODS + ACML", "BDS + ACML", "RS + ML",
                        "BDS + MI", "ODS + MI", "RS + MI"), 3), 
                      c("BDS + ACML", "ODS + ACML", "RS + ML", "BDS + MI", "ODS + MI", "RS + MI"))
# No need for us to plot the random (ML)
tmp          <- tmp %>% filter(Sampling != "RS + ML") 
tmp          <- tmp %>% mutate(Sampling = 
                                 fct_relevel(Sampling,
                                             "ODS + ACML", "BDS + ACML", "RS + MI", 
                                             "ODS + MI", "BDS + MI"))

p1 <- ggplot(tmp, aes(x = Scenario, y = releff.Y1, fill = Sampling)) + 
  geom_bar(position="dodge", stat="identity") + 
  labs(y = "Relative Efficiency", x = "Setting",
       title = unname(TeX(c("$\\beta_{1st}$")))) + 
  theme_bw() +
  geom_hline(yintercept = 1, linetype = "dashed") +
  theme(plot.title = element_text(hjust = 0.5), legend.title = element_blank())  +  
  #scale_fill_manual(values=Blue2OrangeRed14Steps)
  scale_fill_viridis(discrete=TRUE)


p2 <-  ggplot(tmp, aes(x = Scenario, y = releff.Y2, fill = Sampling)) + 
  geom_bar(position="dodge", stat="identity") + 
  labs(y = "Relative Efficiency", x = "Setting",
       title = unname(TeX(c("$\\beta_{2st}$")))) + 
  theme_bw() +
  geom_hline(yintercept = 1, linetype = "dashed") +
  #scale_fill_manual(values=Blue2OrangeRed14Steps) +
  scale_fill_viridis(discrete=TRUE) +
  theme(plot.title = element_text(hjust = 0.5), legend.title = element_blank())

ggarrange(p1, p2, ncol=1, nrow=2, common.legend = TRUE, legend="bottom")
ggsave("releffhist.eps", width = 8, height = 10, dpi = 300)
