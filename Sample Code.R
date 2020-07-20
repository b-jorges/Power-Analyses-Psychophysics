require(dplyr)
require(tidyverse)
require(lme4)
require(DHARMa)


set.seed(1)


ID = paste0("S0",1:10)
ConditionOfInterest = c(0,1)
StandardValues = c(5,6,7,8)
reps = 1:100

PSE_Difference = 0.5
JND_Difference = 0.5

Multiplicator_PSE_Standard = 0
Multiplicator_SD_Standard = 0.15
Type_ResponseFunction = "Normal"
SD_ResponseFunction = 0.1
Mean_Variability_Between = 0.1
SD_Variability_Between = 0.1


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

###prepare for glmer()
Psychometric = Psychometric %>%
  filter(abs(staircase_factor-1) < 0.75) %>%
  group_by(ID,ConditionOfInterest,StandardValues,Difference) %>%
  mutate(Yes = sum(Answer==1),
         Total = length(ConditionOfInterest))




GLMM_RandomIntercepts_PSE = glmer(cbind(Yes, Total - Yes) ~ ConditionOfInterest + Difference + (1| ID), 
                                  family = binomial(link = "probit"),
                                  data = Psychometric)

Sim = simulateResiduals(GLMM_RandomIntercepts_PSE)
plot(Sim)
