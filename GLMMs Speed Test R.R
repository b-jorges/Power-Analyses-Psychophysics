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
require(lmerTest)

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


ConditionOfInterest = c(0,1)
StandardValues = c(5,8)
Range_reps = c(30,40,50,60)
Range_PSE_Difference = 0.1
Range_JND_Difference = 0.2
Multiplicator_PSE_Standard = 0
Multiplicator_SD_Standard = 0.108
SD_ResponseFunction = 0.1
Mean_Variability_Between = 0.1
SD_Variability_Between = 0.1
nIterations = 10
Range_Participants = c(10,12,14,16,18,20)


TotalNumber = length(Range_reps)*length(Range_PSE_Difference)*length(Range_JND_Difference)*length(Range_Participants)
CurrentRunthrough = 0
n = c()
TimelyPowerfulDataframe = c()

for (reps in (Range_reps)){
  for (n in (Range_Participants)){

    CurrentRunthrough = CurrentRunthrough + 1
    print(paste0("This is runthrough N° ", CurrentRunthrough, " out of ", TotalNumber))
    
    for (j in 1:nIterations){
      Psychometric = SimulatePsychometricFunction_Staircase(ID = paste0("s",1:n), 
                                  ConditionOfInterest, 
                                  StandardValues, 
                                  1:reps, 
                                  PSE_Difference, 
                                  JND_Difference, 
                                  Multiplicator_PSE_Standard, 
                                  Multiplicator_SD_Standard, 
                                  SD_ResponseFunction, 
                                  Mean_Variability_Between, 
                                  SD_Variability_Between)
      
      TimeStartTrial = Sys.time() #get time at beginning of trial
      GLMM = glmer(cbind(Yes, Total - Yes) ~ ConditionOfInterest*Difference + (ConditionOfInterest+Difference| ID) + (ConditionOfInterest+Difference| StandardValues), 
                   family = binomial(link = "logit"), 
                   data = Psychometric,
                   nAGQ = 0)
      Duration_NelderMead_nAGQ0 = Sys.time() - TimeStartTrial #get duration of fitting
      Pvalues_NelderMead_nAGQ0 = summary(GLMM)$coefficients[c(14,16)]
      AIC_NelderMead_nAGQ0 = AIC(GLMM)
      
      TimeStartTrial = Sys.time() #get time at beginning of trial
      GLMM = glmer(cbind(Yes, Total - Yes) ~ ConditionOfInterest*Difference + (ConditionOfInterest+Difference| ID) + (ConditionOfInterest+Difference| StandardValues), 
                   family = binomial(link = "logit"), 
                   data = Psychometric,
                   nAGQ = 1)
      Duration_NelderMead_nAGQ1 = Sys.time() - TimeStartTrial #get duration of fitting
      Pvalues_NelderMead_nAGQ1 = summary(GLMM)$coefficients[c(14,16)]
      AIC_NelderMead_nAGQ1 = AIC(GLMM)
      
      TimeStartTrial = Sys.time() #get time at beginning of trial
      GLMM = glmer(cbind(Yes, Total - Yes) ~ ConditionOfInterest*Difference + (ConditionOfInterest+Difference| ID) + (ConditionOfInterest+Difference| StandardValues), 
                   family = binomial(link = "logit"), 
                   data = Psychometric,
                   nAGQ = 0,
                   glmerControl(optimizer = "bobyqa"))
      Duration_Bobyqa_nAGQ0 = Sys.time() - TimeStartTrial #get duration of fitting
      Pvalues_Bobyqa_nAGQ0 = summary(GLMM)$coefficients[c(14,16)]      
      AIC_Bobyqa_nAGQ0 = AIC(GLMM)
      
      TimeStartTrial = Sys.time() #get time at beginning of trial
      GLMM = glmer(cbind(Yes, Total - Yes) ~ ConditionOfInterest*Difference + (ConditionOfInterest+Difference| ID) + (ConditionOfInterest+Difference| StandardValues), 
                   family = binomial(link = "logit"), 
                   data = Psychometric,
                   nAGQ = 1,
                   glmerControl(optimizer = "bobyqa"))
      Duration_Bobyqa_nAGQ1 = Sys.time() - TimeStartTrial #get duration of fitting
      Pvalues_Bobyqa_nAGQ1 = summary(GLMM)$coefficients[c(14,16)]
      AIC_Bobyqa_nAGQ1 = AIC(GLMM)
      
      TimelyPowerfulDataframe = rbind(TimelyDataframe,rbind(reps,n, j,
      Duration_NelderMead_nAGQ0, Pvalues_NelderMead_nAGQ0,AIC_NelderMead_nAGQ0,
      Duration_NelderMead_nAGQ1,Pvalues_NelderMead_nAGQ1,AIC_NelderMead_nAGQ1,
      Duration_Bobyqa_nAGQ0,Pvalues_Bobyqa_nAGQ0,AIC_Bobyqa_nAGQ0,
      Duration_Bobyqa_nAGQ1,Pvalues_Bobyqa_nAGQ1,AIC_Bobyqa_nAGQ1))
    }
  }
}

colnames(TimelyPowerfulDataframe) = c("reps", "n", "iteration", 
                                      "Duration_NelderMead_nAGQ0", "Pvalues_NelderMead_nAGQ0", "AIC_NelderMead_nAGQ0",
                                      "Duration_NelderMead_nAGQ1", "Pvalues_NelderMead_nAGQ1", "AIC_NelderMead_nAGQ1",
                                      "Duration_Bobyqa_nAGQ0", "Pvalues_Bobyqa_nAGQ0", "AIC_Bobyqa_nAGQ0",
                                      "Duration_Bobyqa_nAGQ1", "Pvalues_Bobyqa_nAGQ1", "AIC_Bobyqa_nAGQ1")
write.csv(TimelyPowerfulDataframe,"DurationsR.csv")

