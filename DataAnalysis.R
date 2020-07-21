###########################################
####Parameter-As-Outcome Model (PAOM)######
###########################################


###Fitting psychometric functions and extracting means and standard deviations
Parameters = PsychometricFunctions$par
Parameters2 = Parameters %>%
  filter(parn == "p1") %>%
  select(ID,ConditionOfInterest,Mean=par, StandardValues)
Parameters2$SD = Parameters$par[Parameters$parn == "p2"]
Parameters = Parameters2

###performing ANOVA
###computes pvalues for ANOVA
require(lmerTest) 

ANOVA_Mean = aov(Mean ~ as.factor(ConditionOfInterest)*StandardValues,Parameters)
summary(ANOVA_Mean)

ANOVA_SD = aov(SD ~ as.factor(ConditionOfInterest)*StandardValues,Parameters)
summary(ANOVA_SD)


LM = lm(Mean/StandardValues ~ ConditionOfInterest*StandardValues,
   data = Parameters)
summary(LM)

LMM = lmer(Mean/StandardValues ~ ConditionOfInterest*StandardValues + (1 | ID),
         data = Parameters)
summary(LMM)


Parameters$ConditionOfInterest[Parameters$ConditionOfInterest == 1] = "Condition of Interest"
Parameters$ConditionOfInterest[Parameters$ConditionOfInterest == 0] = "Baseline"

Plot_LMM_Mean = ggplot(Parameters,aes(StandardValues,Mean/StandardValues,color = ID)) +
  geom_point(size = 4) +
  facet_grid(.~ConditionOfInterest) +
  scale_color_manual(name = "",
                     values = colorRampPalette(c(BlauUB,Yellow, Red))(5)) +
  geom_smooth(method='lm', se = FALSE) +
  xlab("Standard Values (m/s)") +
  ylab("Normalized Mean (m/s)")

Plot_LM_Mean = ggplot(Parameters,aes(StandardValues,Mean/StandardValues)) +
  geom_point(size = 4) +
  facet_grid(.~ConditionOfInterest) +
  geom_smooth(method='lm',color = "black", se = FALSE) +
  xlab("Standard Values (m/s)") +
  ylab("Normalized Mean (m/s)")
plot_grid(Plot_LM_Mean,Plot_LMM_Mean, labels = "AUTO")
ggsave("(Figure 4) Means.jpg", w = 12, h = 5)



GLMM_RandomIntercepts_JND = glmer(cbind(Yes, Total - Yes) ~ ConditionOfInterest*Difference + (1| ID), 
             family = binomial(link = "probit"),
             data = Psychometric)
GLMM2_RandomInterceptsAndSlopes_JND = glmer(cbind(Yes, Total - Yes) ~ ConditionOfInterest*Difference + (Difference| ID), 
             family = binomial(link = "probit"),
             data = Psychometric)
GLMM3_RandomInterceptsAndTwoSlopes_JND = glmer(cbind(Yes, Total - Yes) ~ ConditionOfInterest*Difference + (Difference+ConditionOfInterest| ID), 
                                            family = binomial(link = "probit"),
                                            data = Psychometric)
summary(GLMM_RandomIntercepts_JND)
summary(GLMM2_RandomInterceptsAndSlopes_JND)
summary(GLMM3_RandomInterceptsAndTwoSlopes_JND)

GLMM_RandomIntercepts_Null_JND = glmer(cbind(Yes, Total - Yes) ~ ConditionOfInterest + Difference + (1| ID), 
                              family = binomial(link = "probit"),
                              data = Psychometric)
GLMM2_RandomInterceptsAndSlopes_Null_JND = glmer(cbind(Yes, Total - Yes) ~ ConditionOfInterest + Difference + (Difference| ID), 
                                        family = binomial(link = "probit"),
                                        data = Psychometric)
GLMM3_RandomInterceptsTwoAndSlopes_Null_JND = glmer(cbind(Yes, Total - Yes) ~ ConditionOfInterest + Difference + (Difference + ConditionOfInterest| ID), 
                                                 family = binomial(link = "probit"),
                                                 data = Psychometric)

anova(GLMM_RandomIntercepts_JND,GLMM_RandomIntercepts_Null_JND)
anova(GLMM2_RandomInterceptsAndSlopes_JND,GLMM2_RandomInterceptsAndSlopes_Null_JND)
anova(GLMM3_RandomInterceptsAndTwoSlopes_JND,GLMM3_RandomInterceptsTwoAndSlopes_Null_JND)





GLMM_RandomIntercepts_PSE = glmer(cbind(Yes, Total - Yes) ~ ConditionOfInterest*Difference + (1| ID), 
                              family = binomial(link = "probit"),
                              data = Psychometric)
GLMM2_RandomInterceptsAndSlopes_PSE = glmer(cbind(Yes, Total - Yes) ~ ConditionOfInterest + Difference + (Difference| ID), 
                                        family = binomial(link = "probit"),
                                        data = Psychometric)
GLMM3_RandomInterceptsTwoAndSlopes_PSE = glmer(cbind(Yes, Total - Yes) ~ ConditionOfInterest + Difference + (Difference + ConditionOfInterest| ID), 
                                            family = binomial(link = "probit"),
                                            data = Psychometric)

GLMM_RandomIntercepts_Null_PSE = glmer(cbind(Yes, Total - Yes) ~ Difference + (1| ID), 
                                   family = binomial(link = "probit"),
                                   data = Psychometric)
GLMM2_RandomInterceptsAndSlopes_Null_PSE = glmer(cbind(Yes, Total - Yes) ~ Difference + (Difference| ID), 
                                             family = binomial(link = "probit"),
                                             data = Psychometric)
GLMM3_RandomInterceptsAndTwoSlopes_Null_PSE = glmer(cbind(Yes, Total - Yes) ~ Difference + (Difference + ConditionOfInterest| ID), 
                                                 family = binomial(link = "probit"),
                                                 data = Psychometric)

anova(GLMM_RandomIntercepts_PSE,GLMM_RandomIntercepts_Null_PSE)
anova(GLMM2_RandomInterceptsAndSlopes_PSE,GLMM2_RandomInterceptsAndSlopes_Null_PSE)
anova(GLMM3_RandomInterceptsTwoAndSlopes_PSE,GLMM3_RandomInterceptsAndTwoSlopes_Null_PSE)


####Testing assumptions is hard. The DHARMa package helps a lot with that: 
require(DHARMa)


Sim1 = simulateResiduals(GLMM_RandomIntercepts_JND)
plot(Sim1)

Sim2 = simulateResiduals(GLMM2_RandomInterceptsAndSlopes_JND)
plot(Sim2)

Sim3 = simulateResiduals(GLMM3_RandomInterceptsAndTwoSlopes_JND)
plot(Sim3)


Sim4 = simulateResiduals(GLMM_RandomIntercepts_PSE)
plot(Sim4)

Sim5 = simulateResiduals(GLMM2_RandomInterceptsAndSlopes_PSE)
plot(Sim5)

Sim6 = simulateResiduals(GLMM3_RandomInterceptsTwoAndSlopes_PSE)
plot(Sim6)


GLMM4_RandomInterceptsTimesTwo_JND = glmer(cbind(Yes, Total - Yes) ~ ConditionOfInterest*Difference + 
                                            (1| ID) + 
                                            (1| StandardValues), 
                                  family = binomial(link = "probit"),
                                  data = Psychometric)
GLMM5_RandomInterceptsAndSlopesTimesTwo_JND = glmer(cbind(Yes, Total - Yes) ~ ConditionOfInterest*Difference + 
                                                      (Difference| ID) + 
                                                      (Difference| StandardValues), 
                                            family = binomial(link = "probit"),
                                            data = Psychometric)
GLMM6_RandomInterceptsAndTwoSlopesTimesTwo_JND = glmer(cbind(Yes, Total - Yes) ~ ConditionOfInterest*Difference + 
                                                         (Difference+ConditionOfInterest| ID) +
                                                         (Difference+ConditionOfInterest| StandardValues), 
                                               family = binomial(link = "probit"),
                                               data = Psychometric)

Sim7 = simulateResiduals(GLMM4_RandomInterceptsTimesTwo_JND)
plot(Sim7)

Sim8 = simulateResiduals(GLMM5_RandomInterceptsAndSlopesTimesTwo_JND)
plot(Sim8)

Sim9 = simulateResiduals(GLMM6_RandomInterceptsAndTwoSlopesTimesTwo_JND)
plot(Sim9)

hello = testResiduals(Sim9)
hello$uniformity$p.value
hello$dispersion$p.value

testDispersion(Sim9)
?simulateResiduals