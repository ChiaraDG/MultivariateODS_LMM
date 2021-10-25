# load the data
LHSres <- read.csv("LHSstudyFEVFVC.csv")

# get the values for the time variant covariates
# TIME
LHSres$lb.time.Y1 <- LHSres$est.time.Y1 - qnorm(0.975)*sqrt(LHSres$var.time.Y1)
LHSres$ub.time.Y1 <- LHSres$est.time.Y1 + qnorm(0.975)*sqrt(LHSres$var.time.Y1)
LHSres$lb.time.Y2 <- LHSres$est.time.Y2 - qnorm(0.975)*sqrt(LHSres$var.time.Y2)
LHSres$ub.time.Y2 <- LHSres$est.time.Y2 + qnorm(0.975)*sqrt(LHSres$var.time.Y2)
# GROUP BY TIME
LHSres$lb.tg.Y1 <- LHSres$est.tg.Y1 - qnorm(0.975)*sqrt(LHSres$var.tg.Y1)
LHSres$ub.tg.Y1 <- LHSres$est.tg.Y1 + qnorm(0.975)*sqrt(LHSres$var.tg.Y1)
LHSres$lb.tg.Y2 <- LHSres$est.tg.Y2 - qnorm(0.975)*sqrt(LHSres$var.tg.Y2)
LHSres$ub.tg.Y2 <- LHSres$est.tg.Y2 + qnorm(0.975)*sqrt(LHSres$var.tg.Y2)
# TIME* FEMALE
LHSres$lb.tf.Y1 <- LHSres$est.tf.Y1 - qnorm(0.975)*sqrt(LHSres$var.tf.Y1)
LHSres$ub.tf.Y1 <- LHSres$est.tf.Y1 + qnorm(0.975)*sqrt(LHSres$var.tf.Y1)
LHSres$lb.tf.Y2 <- LHSres$est.tf.Y2 - qnorm(0.975)*sqrt(LHSres$var.tf.Y2)
LHSres$ub.tf.Y2 <- LHSres$est.tf.Y2 + qnorm(0.975)*sqrt(LHSres$var.tf.Y2)
# TIME* AGE
LHSres$lb.ta.Y1 <- LHSres$est.ta.Y1 - qnorm(0.975)*sqrt(LHSres$var.ta.Y1)
LHSres$ub.ta.Y1 <- LHSres$est.ta.Y1 + qnorm(0.975)*sqrt(LHSres$var.ta.Y1)
LHSres$lb.ta.Y2 <- LHSres$est.ta.Y2 - qnorm(0.975)*sqrt(LHSres$var.ta.Y2)
LHSres$ub.ta.Y2 <- LHSres$est.ta.Y2 + qnorm(0.975)*sqrt(LHSres$var.ta.Y2)
# GROUP
LHSres$lb.grp.Y1 <- LHSres$est.grp.Y1 - qnorm(0.975)*sqrt(LHSres$var.grp.Y1)
LHSres$ub.grp.Y1 <- LHSres$est.grp.Y1 + qnorm(0.975)*sqrt(LHSres$var.grp.Y1)
LHSres$lb.grp.Y2 <- LHSres$est.grp.Y2 - qnorm(0.975)*sqrt(LHSres$var.grp.Y2)
LHSres$ub.grp.Y2 <- LHSres$est.grp.Y2 + qnorm(0.975)*sqrt(LHSres$var.grp.Y2)

# time coefficients
time.coef <- LHSres %>% group_by(Sampling) %>% 
  dplyr::select(est.grp.Y1, est.grp.Y2,
         est.time.Y1, est.time.Y2, est.tf.Y1, est.tf.Y2, est.ta.Y1, est.ta.Y2, 
         est.tg.Y1, est.tg.Y2) %>%
  summarise_all(mean)
time.sd <- LHSres %>% group_by(Sampling) %>% 
  dplyr::select(var.grp.Y1, var.grp.Y2,
    var.time.Y1, var.time.Y2, var.tf.Y1, var.tf.Y2, var.ta.Y1, var.ta.Y2, var.tg.Y1, var.tg.Y2) %>%
  summarise_all(function(x){sqrt(mean(x))})
time.lb   <- LHSres %>% group_by(Sampling) %>% 
  dplyr::select(lb.grp.Y1, lb.grp.Y2,
    lb.time.Y1, lb.time.Y2, lb.tf.Y1, lb.tf.Y2, lb.ta.Y1, lb.ta.Y2, lb.tg.Y1, lb.tg.Y2) %>% 
  summarise_all(mean)
time.ub   <- LHSres %>% group_by(Sampling) %>% 
  dplyr::select(ub.grp.Y1, ub.grp.Y2,
    ub.time.Y1, ub.time.Y2, ub.tf.Y1, ub.tf.Y2, ub.ta.Y1, ub.ta.Y2, ub.tg.Y1, ub.tg.Y2) %>% 
  summarise_all(mean)

time.tot  <- merge(time.coef, time.sd, by = "Sampling")
time.tot  <- merge(time.tot, time.lb, by = "Sampling")
time.tot  <- merge(time.tot, time.ub, by = "Sampling")

# add full cohort results
names(df) <- NULL
fc.time <- c(df["SNP",c(1,6)], df["Time",c(1,6)], df["Time*Female",c(1,6)], 
             df["Time*Age",c(1,6)], df["Time*SNP", c(1,6)],
             df["SNP", c(2, 7)], df["Time",c(2,7)], df["Time*Female",c(2,7)], df["Time*Age",c(2,7)], 
             df["Time*SNP", c(2,7)],
             df["SNP", c(3, 8)], df["Time",c(3,8)], df["Time*Female",c(3,8)], df["Time*Age",c(3,8)], 
             df["Time*SNP", c(3,8)],
             df["SNP", c(4, 9)], df["Time",c(4,9)], df["Time*Female",c(4,9)], df["Time*Age",c(4,9)], 
             df["Time*SNP", c(4,9)])
fc.time               <- as.vector(t(unlist(fc.time)))
fc.time               <- c(NA, fc.time)
names(fc.time)        <- names(time.tot)
time.tot              <- rbind(time.tot, fc.time)
time.tot$Sampling     <- c("ODS + ACML", 
                           "BDS + ACML", 
                           "RS + ML", "RS + MI", 
                           "BDS Slp + MI", 
                           "ODS Slp + MI",  
                           "FC")

time.tot <- time.tot %>% 
  mutate(Sampling = fct_relevel(Sampling,
                                "BDS + MI",
                                "ODS + MI", 
                                "RS + MI",
                                "BDS + ACML", 
                                "ODS + ACML", 
                                "RS + ML", "FC"))


time.tot.Y1 <- time.tot %>% 
  dplyr::select(Sampling, est.grp.Y1, var.grp.Y1,
         est.time.Y1, var.time.Y1, est.tg.Y1, var.tg.Y1,
         est.ta.Y1, var.ta.Y1, est.tf.Y1, var.tf.Y1)
time.tot.Y2 <- time.tot %>% 
  dplyr::select(Sampling, est.grp.Y2, var.grp.Y2,
         est.time.Y2, var.time.Y2, est.tg.Y2, var.tg.Y2,
         est.ta.Y2, var.ta.Y2, est.tf.Y2, var.tf.Y2)


kable(time.tot.Y1, digits = 3,
      col.names = c("Sampling", "grp", "se(grp)", "time", "se(time)",
                    "time*grp", "se(time*grp)", "time*age", "se(time*age)",
                    "time*sex", "se(time*sex)")) %>% add_header_above(c(" " = 1, "FEV" = 10)) %>% 
  kable_styling()

kable(time.tot.Y2, digits = 3,
      col.names = c("Sampling", "grp", "se(grp)", "time", "se(time)",
                    "time*grp", "se(time*grp)", "time*age", "se(time*age)",
                    "time*sex", "se(time*sex)")) %>% 
  add_header_above(c(" " = 1, "FVC" = 10)) %>% 
  kable_styling()


p3 <- ggplot(time.tot, aes(x = Sampling, y = est.tg.Y1)) + geom_point() +
  geom_errorbar(aes(ymin = lb.tg.Y1, ymax = ub.tg.Y1), width=.1) + 
  ylim(c(-0.25, 0.25)) +
  labs(title = "FEV", x = "",  y = TeX("Study Year $\\times$ Presence of the T-allele")) + 
  coord_flip()  + theme_bw()  +
  annotate("text", x = 7.1, y = 0.19, label= "-0.078 (95% CI: -0.123	-0.032)", size = 2.8) + 
  annotate("text", x = 6.1, y = 0.19, label= "-0.079 (95% CI: -0.161, 0.003)", size = 2.8) + 
  annotate("text", x = 5.1, y = 0.19, label= "-0.084 (95% CI: -0.160, -0.008)", size = 2.8) + 
  annotate("text", x = 4.1, y = 0.19, label= "-0.073 (95% CI: -0.147, 0.001)", size = 2.8) + 
  annotate("text", x = 3.1, y = 0.19, label= "-0.074 (95% CI: -0.148, 0.000)", size = 2.8) + 
  annotate("text", x = 2.1, y = 0.19, label= "-0.080 (95% CI: -0.142, -0.019)", size = 2.8) +
  annotate("text", x = 1.1, y = 0.19, label= "-0.070 (95% CI: -0.126, -0.008)", size = 2.8) +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_hline(yintercept = 0, linetype="dotted")


p4 <- ggplot(time.tot, aes(x = Sampling, y = est.tg.Y2)) + geom_point() +
  geom_errorbar(aes(ymin = lb.tg.Y2, ymax = ub.tg.Y2), width=.1) + 
  ylim(c(-0.25, 0.25)) +
  labs(title = "FVC", x = "",  y = TeX("Study Year $\\times$ Presence of the T-allele")) + 
  coord_flip() + theme_bw() + 
  annotate("text", x = 7.1, y = 0.19, label= "-0.074 (95% CI: -0.129, -0.019)", size = 2.8) + 
  annotate("text", x = 6.1, y = 0.19, label= "-0.077 (95% CI: -0.177, 0.023)", size = 2.8) + 
  annotate("text", x = 5.1, y = 0.19, label= "-0.083 (95% CI: -0.175, 0.009)", size = 2.8) + 
  annotate("text", x = 4.1, y = 0.19, label= "-0.081 (95% CI: -0.175, 0.013)", size = 2.8) + 
  annotate("text", x = 3.1, y = 0.19, label= "-0.090 (95% CI: -0.180, 0.000)", size = 2.8) + 
  annotate("text", x = 2.1, y = 0.19, label= "-0.081 (95% CI: -0.155, -0.007)", size = 2.8) +
  annotate("text", x = 1.1, y = 0.19, label= "-0.060 (95% CI: -0.127, 0.018)", size = 2.8) +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_hline(yintercept = 0, linetype="dotted")


p7 <- ggplot(time.tot, aes(x = Sampling, y = est.tf.Y1)) + geom_point() +
  geom_errorbar(aes(ymin = lb.tf.Y1, ymax = ub.tf.Y1), width=.1) + 
  ylim(c(-0.15, 0.30)) +
  labs(title = "FEV", x = "", 
       y = TeX("Study Year $\\times$ Female")) + 
  coord_flip() + theme_bw() + 
  annotate("text", x = 7.1, y = 0.24, label= "0.083 (95% CI: 0.035, 0.130)", size = 2.8) + 
  annotate("text", x = 6.1, y = 0.24, label= "0.080 (95% CI: -0.004, 0.164)", size = 2.8) + 
  annotate("text", x = 5.1, y = 0.24, label= "0.092 (95% CI: 0.011, 0.172)", size = 2.8) + 
  annotate("text", x = 4.1, y = 0.24, label= "0.088 (95% CI: 0.008, 0.168)", size = 2.8) + 
  annotate("text", x = 3.1, y = 0.24, label= "0.082 (95% CI: 0.035, 0.130)", size = 2.8) + 
  annotate("text", x = 2.1, y = 0.24, label= "0.083 (95% CI: 0.035, 0.130)", size = 2.8) + 
  annotate("text", x = 1.1, y = 0.24, label= "0.084 (95% CI: 0.036, 0.131)", size = 2.8) +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_hline(yintercept = 0, linetype="dotted")

p8 <- ggplot(time.tot, aes(x = Sampling, y = est.tf.Y2)) + geom_point() +
  geom_errorbar(aes(ymin = lb.tf.Y2, ymax = ub.tf.Y2), width=.1) + 
  ylim(c(-0.15, 0.30)) +
  labs(title = "FVC", x = "", 
        y = TeX("Study Year $\\times$ Female")) + coord_flip() + 
  theme_bw() + 
  annotate("text", x = 7.1, y = 0.24, label= "-0.021 (95% CI: -0.078, 0.036)", size = 2.8) + 
  annotate("text", x = 6.1, y = 0.24, label= "-0.021 (95% CI: -0.123, 0.081)", size = 2.8) + 
  annotate("text", x = 5.1, y = 0.24, label= "-0.003 (95% CI: -0.099, 0.093)", size = 2.8) + 
  annotate("text", x = 4.1, y = 0.24, label= "-0.022 (95% CI: -0.122, 0.078)", size = 2.8) + 
  annotate("text", x = 3.1, y = 0.24, label= "-0.021 (95% CI: -0.078, 0.036)", size = 2.8) + 
  annotate("text", x = 2.1, y = 0.24, label= "-0.020 (95% CI: -0.078, 0.036)", size = 2.8) + 
  annotate("text", x = 1.1, y = 0.24, label= "-0.020 (95% CI: -0.078, 0.036)", size = 2.8) +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_hline(yintercept = 0, linetype="dotted")


ggpubr::ggarrange(p3, p4, p7, p8, ncol = 2, nrow = 2)
# How to save the LHS plot
ggsave("lhsCI.eps", width = 14, height = 8, dpi = 300)
