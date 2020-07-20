###Pull the whole repository
require(dplyr)
require(tidyverse)
require(lme4)
require(ggplot2)
require(cowplot)
theme_set(theme_cowplot())
require(quickpsy)
require(brms)
require(rstan)
#require(lmerTest)
#require(DHARMa)

Where_Am_I <- function(path=T){
  if (path == T){
    dirname(rstudioapi::getSourceEditorContext()$path)
  }
  else {
    rstudioapi::getSourceEditorContext()$path
  }
}

binomial_smooth <- function(...) {
  geom_smooth(method = "glm", method.args = list(family = "binomial"), ...)}

setwd(Where_Am_I())

source("Utilities/parabolic.r")
source("Utilities/functions.r")
source("Utilities/colourschemes.r")
source("Utilities/PowerFunctions.r")

#optimize for fitting of Bayesian Linear Mixed Models (packages "rstan", "bmrs")
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
Sys.setenv(LOCAL_CPPFLAGS = '-march=corei7')


#set.seed(9121)


ID = paste0("s",1:15)
ConditionOfInterest = c(0,1)
StandardValues = c(5,8)
reps = 1:100
PSE_Difference = -0.1
JND_Difference = 0.25
Multiplicator_PSE_Standard = 0
Multiplicator_SD_Standard = 0.15
Type_ResponseFunction = "Normal"
SD_ResponseFunction = 0.1
Mean_Variability_Between = 0.1
SD_Variability_Between = 0.1

Dataframe = data.frame()

Psychometric = expand.grid(ID=ID, ConditionOfInterest=ConditionOfInterest, StandardValues=StandardValues, reps = reps)

Psychometric = Psychometric %>%
  group_by(ID) %>%#
  mutate(PSE_Factor_ID = rnorm(1,1,Mean_Variability_Between), #how much variability is in the means of the psychometric functions between subjects?
         SD_Factor_ID = rnorm(1,1,SD_Variability_Between)) #how much variability is in the standard deviations of the psychometric functions between subjects?

Psychometric = Psychometric %>%
  mutate(
    Mean_Standard = StandardValues+StandardValues*Multiplicator_PSE_Standard,
    SD_Standard = StandardValues*Multiplicator_SD_Standard,
    Mean = (Mean_Standard + (ConditionOfInterest==1)*Mean_Standard*PSE_Difference),
    SD = abs(SD_Standard + (ConditionOfInterest==1)*SD_Standard*JND_Difference))

Psychometric = Psychometric %>%
  mutate(
    Mean = Mean*PSE_Factor_ID,
    SD = SD*SD_Factor_ID)

if (Type_ResponseFunction == "normal"){
  
  Psychometric = Psychometric %>%
    mutate(
      staircase_factor = pnorm(length(reps),1,SD_ResponseFunction*(1+ConditionOfInterest*JND_Difference)))
  
} else if (Type_ResponseFunction == "Cauchy"){
  
  Psychometric = Psychometric %>%
    mutate(
      staircase_factor = rcauchy(length(reps),1,SD_ResponseFunction*(1+ConditionOfInterest*JND_Difference)))
  
  #} else if (Type_ResponseFunction == "uniform"){
  #  
  #  Psychometric = Psychometric %>%
  #  mutate(
  #      staircase_factor = seq(SD_ResponseFunction[1],SD_ResponseFunction[2],(SD_ResponseFunction[2]-SD_ResponseFunc#tion[1]/6)))
  
} else{
  
  print("distribution not valid")
  
}

Psychometric = Psychometric %>%
  mutate(
    staircase_factor = rcauchy(length(reps),1,SD_ResponseFunction), 
    Presented_TestStimulusStrength = Mean*staircase_factor,
    Difference = Presented_TestStimulusStrength - StandardValues)

Psychometric = Psychometric %>%
  mutate(
    AnswerProbability = pnorm(Presented_TestStimulusStrength,Mean,SD),
    
    ##get binary answers ("Test was stronger" yes/no) from probabilities for each trial
    Answer = as.numeric(rbernoulli(length(AnswerProbability),AnswerProbability))
  )

###prepare for glmer() - needs sum of YES/Total per stimulus strength and condition
Psychometric = Psychometric %>%
  filter(abs(staircase_factor-1) < 0.75) %>%
  group_by(ID,ConditionOfInterest,StandardValues,Difference) %>%
  mutate(Yes = sum(Answer==1),
         Total = length(ConditionOfInterest))



GLMM = glmer(cbind(Yes, Total - Yes) ~ ConditionOfInterest*Difference + (ConditionOfInterest+Difference| ID) + (ConditionOfInterest+Difference| StandardValues), 
             family = binomial(link = "logit"), 
             data = Psychometric,
             nAGQ = 0,
             control = glmerControl(optimizer = "nloptwrap"))

GLMM2 = glmer(cbind(Yes, Total - Yes) ~ ConditionOfInterest*Difference + (Difference| ID) + (Difference| StandardValues), 
              family = binomial(link = "logit"), 
              data = Psychometric,
              nAGQ = 0,
              control = glmerControl(optimizer = "nloptwrap"))

GLMM3 = glmer(cbind(Yes, Total - Yes) ~ ConditionOfInterest*Difference + (ConditionOfInterest| ID) + (ConditionOfInterest| StandardValues), 
              family = binomial(link = "logit"), 
              data = Psychometric,
              nAGQ = 0,
              control = glmerControl(optimizer = "nloptwrap"))

GLMM4 = glmer(cbind(Yes, Total - Yes) ~ ConditionOfInterest*Difference + (1| ID) + (1| StandardValues), 
              family = binomial(link = "logit"), 
              data = Psychometric,
              nAGQ = 0,
              control = glmerControl(optimizer = "nloptwrap"))

AICs = c(summary(GLMM)$AIC[[1]],summary(GLMM2)$AIC[[1]],summary(GLMM3)$AIC[[1]],summary(GLMM4)$AIC[[1]])

PvalueAccuracy_1 = summary(GLMM)$coefficients[14]
PvalueAccuracy_2 = summary(GLMM2)$coefficients[14]
PvalueAccuracy_3 = summary(GLMM3)$coefficients[14]
PvalueAccuracy_4 = summary(GLMM4)$coefficients[14]

PvaluePrecision_1 = summary(GLMM)$coefficients[16]
PvaluePrecision_2 = summary(GLMM2)$coefficients[16]
PvaluePrecision_3 = summary(GLMM3)$coefficients[16]
PvaluePrecision_4 = summary(GLMM4)$coefficients[16]

pValues_GLMM_Accuracy = c(PvalueAccuracy_1,PvalueAccuracy_2,PvalueAccuracy_3,PvalueAccuracy_4)
pValues_GLMM_Precision = c(PvaluePrecision_1,PvaluePrecision_2,PvaluePrecision_3,PvaluePrecision_4)


Summ1$`Pr(>Chisq)`

Parameters = quickpsy(Psychometric,Difference,Answer,grouping = .(ID,ConditionOfInterest,StandardValues), bootstrap = "none")$par
Parameters2 = data.frame(Parameters) %>%
  filter(parn == "p1") %>%
  dplyr::select(ID,ConditionOfInterest,Mean=par, StandardValues)
Parameters2$SD = Parameters$par[Parameters$parn == "p2"]
Parameters = Parameters2


LMM4_Mean = lmer(Mean ~ ConditionOfInterest + (1| ID) + (1| StandardValues),
            data = Psychometric)

AICs = c(AICs,
         summary(LMM1)$AIC[[1]],summary(LMM2)$AIC[[1]],summary(LMM3)$AIC[[1]],summary(LMM4)$AIC[[1]])
Kurtosis = c(0,0,0,0,
             descdist(data.frame(residuals(LMM1))$residuals.LMM1_Mean., discrete = FALSE)$kurtosis,
             descdist(data.frame(residuals(LMM2))$residuals.LMM2_Mean., discrete = FALSE)$kurtosis,
             descdist(data.frame(residuals(LMM3))$residuals.LMM3_Mean., discrete = FALSE)$kurtosis,
             descdist(data.frame(residuals(LMM4))$residuals.LMM4_Mean., discrete = FALSE)$kurtosis)

Skewness = c(0,0,0,0,
             descdist(data.frame(residuals(LMM1))$residuals.LMM1_Mean., discrete = FALSE)$skewness,
             descdist(data.frame(residuals(LMM2))$residuals.LMM2_Mean., discrete = FALSE)$skewness,
             descdist(data.frame(residuals(LMM3))$residuals.LMM3_Mean., discrete = FALSE)$skewness,
             descdist(data.frame(residuals(LMM4))$residuals.LMM4_Mean., discrete = FALSE)$skewness)

LMM1_Accuracy = summary(LMM1)

summary(LMM1)$coefficients[14]

Designation = c("GLMM1","GLMM2","GLMM3","GLMM4",
                "LMM1","LMM2","LMM3","LMM4")

rbind(Dataframe,
      data.frame(AICs = AICs,
                 Kurtosis = Kurtosis, 
                 Skewness = Skewness,
                 Designation = Designation,
                 pValues_GLMM_Accuracy = pValues_GLMM_Accuracy,
                 pValues_GLMM_Precision = pValues_GLMM_Precision))