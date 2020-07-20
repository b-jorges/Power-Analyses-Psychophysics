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
             family = binomial(link = "logit"),
             data = Psychometric)
GLMM2_RandomInterceptsAndSlopes_JND = glmer(cbind(Yes, Total - Yes) ~ ConditionOfInterest*Difference + (Difference| ID), 
             family = binomial(link = "logit"),
             data = Psychometric)
summary(GLMM_RandomIntercepts)
summary(GLMM2_RandomInterceptsAndSlopes)

GLMM_RandomIntercepts_Null_JND = glmer(cbind(Yes, Total - Yes) ~ ConditionOfInterest + Difference + (1| ID), 
                              family = binomial(link = "probit"),
                              data = Psychometric)
GLMM2_RandomInterceptsAndSlopes_Null_JND = glmer(cbind(Yes, Total - Yes) ~ ConditionOfInterest + Difference + (Difference| ID), 
                                        family = binomial(link = "probit"),
                                        data = Psychometric)

anova(GLMM_RandomIntercepts,GLMM_RandomIntercepts_Null)
anova(GLMM2_RandomInterceptsAndSlopes,GLMM2_RandomInterceptsAndSlopes_Null)





GLMM_RandomIntercepts_PSE = glmer(cbind(Yes, Total - Yes) ~ ConditionOfInterest + Difference + (1| ID), 
                              family = binomial(link = "probit"),
                              data = Psychometric)
GLMM2_RandomInterceptsAndSlopes_PSE = glmer(cbind(Yes, Total - Yes) ~ ConditionOfInterest + Difference + (Difference| ID), 
                                        family = binomial(link = "probit"),
                                        data = Psychometric)
GLMM3_RandomInterceptsAndSlopes_PSE = glmer(cbind(Yes, Total - Yes) ~ ConditionOfInterest*Difference + (Difference + ConditionOfInterest| ID), 
                                            family = binomial(link = "probit"),
                                            data = Psychometric)

GLMM_RandomIntercepts_Null_PSE = glmer(cbind(Yes, Total - Yes) ~ Difference + (1| ID), 
                                   family = binomial(link = "probit"),
                                   data = Psychometric)
GLMM2_RandomInterceptsAndSlopes_Null_PSE = glmer(cbind(Yes, Total - Yes) ~ Difference + (Difference| ID), 
                                             family = binomial(link = "probit"),
                                             data = Psychometric)

anova(GLMM_RandomIntercepts_PSE,GLMM2_RandomInterceptsAndSlopes_PSE)
anova(GLMM_RandomIntercepts_JND,GLMM2_RandomInterceptsAndSlopes_JND)


require(DHARMa)


summary(GLMM_RandomIntercepts_PSE)

summary(GLMM_RandomIntercepts_PSE)
Sim = simulateResiduals(GLMM_RandomIntercepts_PSE)
plot(Sim)

Sim = simulateResiduals(GLMM2_RandomInterceptsAndSlopes_JND)
plot(Sim)

Sim = simulateResiduals(GLMM_RandomIntercepts_PSE)
plot(Sim)

Sim = simulateResiduals(GLMM3_RandomInterceptsAndSlopes_PSE)
plot(Sim)

plot(predict(GLMM3_RandomInterceptsAndSlopes_PSE),resid(GLMM3_RandomInterceptsAndSlopes_PSE))
plot(predict(GLMM3_RandomInterceptsAndSlopes_PSE),residuals(GLMM3_RandomInterceptsAndSlopes_PSE,"response"))


ggplot(data.frame(residuals(GLMM3_RandomInterceptsAndSlopes_PSE,"response")),aes(residuals.GLMM3_RandomInterceptsAndSlopes_PSE...response..)) +
  geom_density()

colnames(data.frame(residuals(GLMM3_RandomInterceptsAndSlopes_PSE,"response")))

plot(residuals(GLMM3_RandomInterceptsAndSlopes_PSE,"response"))

GLMM3_RandomInterceptsAndSlopes_PSE


GLMM_RandomIntercepts_PSE

Sim = simulateResiduals(GLMM2_RandomInterceptsAndSlopes_Null_PSE)
plot(Sim)
plot(GLMM3_RandomInterceptsAndSlopes_PSE)

devtools::install_github("goodekat/ggResidpanel")
require(ggResidpanel)
resid_panel(GLMM2_RandomInterceptsAndSlopes_PSE)
resid_panel(GLMM3_RandomInterceptsAndSlopes_PSE)
