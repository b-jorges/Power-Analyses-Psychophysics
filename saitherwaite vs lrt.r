ConditionOfInterest = c(0,1)
StandardValues = c(5,8)
Range_reps = c(50)
PSE_Difference = 0.0
JND_Difference = 0.0
Multiplicator_PSE_Standard = 0
Multiplicator_SD_Standard = 0.108
SD_ResponseFunction = 0.1
Mean_Variability_Between = 0.00001
SD_Variability_Between = 0.00001
nIterations = 50
Range_Participants = c(20)

TotalNumber = length(Range_reps)*length(Range_PSE_Difference)*length(Range_JND_Difference)*length(Range_Participants)
CurrentRunthrough = 0
n = c()
TimelyRow = c()
TimelyPowerfulDataframe_LTRvsSatherwaite = c()

for (reps in (Range_reps)){
  for (n in (Range_Participants)){
    TimeStartTrial = Sys.time()
    
    CurrentRunthrough = CurrentRunthrough + 1
    print(paste0("This is runthrough N° ", CurrentRunthrough))
    
    Pvalues_Precision_LTR = c()
    Pvalues_Precision_Satherwaite = c()
    
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
      
      GLMM = glmer(cbind(Yes, Total - Yes) ~ ConditionOfInterest*Difference + (ConditionOfInterest+Difference| ID) + (ConditionOfInterest+Difference| StandardValues), 
                   family = binomial(link = "probit"), 
                   data = Psychometric,
                   nAGQ = 0)
      
      GLMM1 = glmer(cbind(Yes, Total - Yes) ~ ConditionOfInterest + Difference + (ConditionOfInterest+Difference| ID) + (ConditionOfInterest+Difference| StandardValues), 
                    family = binomial(link = "probit"), 
                    data = Psychometric,
                    nAGQ = 0)
      
      p = anova(GLMM,GLMM1)$`Pr(>Chisq)`[2]
      p2 = summary(GLMM)$coefficients[16]
      
      print(paste0("LTR: ",p))
      print(paste0("Satherwaite: ",p2))
      
      Pvalues_Precision_LTR = c(Pvalues_Precision_LTR,p)
      Pvalues_Precision_Satherwaite = c(Pvalues_Precision_Satherwaite,p2)
    }
    
    print(Sys.time() - TimeStartTrial)
    
    TimelyRow = c(mean(Pvalues_Precision_LTR < 0.05),mean(Pvalues_Precision_Satherwaite < 0.05),reps, n, nIterations, Sys.time() - TimeStartTrial)
    
    TimelyPowerfulDataframe_LTRvsSatherwaite = rbind(TimelyPowerfulDataframe_LTRvsSatherwaite,TimelyRow)
  }
}


TimelyPowerfulDataframe_LTRvsSatherwaite
TimelyPowerfulDataframe_LTRvsSatherwaite3
TimelyPowerfulDataframe_LTRvsSatherwaite2

halo = summary(GLMM)

