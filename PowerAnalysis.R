set.seed(1)

RangeNs = c(10,12,14,16,18,20)
RangeRepetitions = c(40,70,100)

ConditionOfInterest = c(0,1)
StandardValues = c(5,6,7,8)
PSE_Difference = 0.025
JND_Difference = 0.05
Multiplicator_PSE_Standard = 0
Multiplicator_SD_Standard = 0.15
Type_ResponseFunction = "normal"
SD_ResponseFunction = 0.1
Mean_Variability_Between = 0.2
SD_Variability_Between = 0.2

nIterations = 200

PowerfulDataframe = data.frame()

for (nParticipants in RangeNs){
  for (reps in RangeRepetitions){
    
    TimeStartTrial = Sys.time() #get time at beginning of trial
    
    for(i in 1:nIterations){
      
    Psychometric = SimulatePsychometricData(nParticipants,
                                            ConditionOfInterest,
                                            StandardValues,
                                            reps,
                                            PSE_Difference,
                                            JND_Difference,
                                            Multiplicator_PSE_Standard,
                                            Multiplicator_SD_Standard,
                                            Type_ResponseFunction,
                                            SD_ResponseFunction,
                                            Mean_Variability_Between,
                                            SD_Variability_Between)
    

    GLMM = glmer(cbind(Yes, Total - Yes) ~ ConditionOfInterest*Difference + (ConditionOfInterest+Difference| ID) + (ConditionOfInterest+Difference| StandardValues), 
                 family = binomial(link = "logit"), 
                 data = Psychometric,
                 nAGQ = 1,
                 glmerControl(optimizer = "nloptwrap"))

    summary(GLMM)$coefficients[14]
    summary(GLMM)$coefficients[16]
    
    PowerfulDataframe = rbind(PowerfulDataframe,c(nParticipants=nParticipants,
                                                  reps=reps, 
                                                  pvalue_PSE = summary(GLMM)$coefficients[14],
                                                  pvalue_JND = summary(GLMM)$coefficients[16],
                                                  iteration = i))
    colnames(PowerfulDataframe) = c("nParticipants","reps","pvalue_PSE","pvalue_JND","iteration")
    
    }
    print(paste0(nIterations, " iterations took ", round(Sys.time() - TimeStartTrial), " seconds. The power for the current run through (",nParticipants," Participants, ", reps, " Repetitions) is ",mean(PowerfulDataframe$pvalue_PSE[PowerfulDataframe$nParticipants == nParticipants & PowerfulDataframe$reps == reps] < 0.05), " (PSE) and ", mean(PowerfulDataframe$pvalue_JND[PowerfulDataframe$nParticipants == nParticipants & PowerfulDataframe$reps == reps] < 0.05)," (JND)."))
  }
}

alpha = 0.05

PowerfulDataframe = PowerfulDataframe %>% group_by(nParticipants,reps) %>% 
  mutate(Power_PSE = mean(pvalue_PSE < alpha),
         Power_JND = mean(pvalue_JND < alpha))

PowerfulDataframe %>% group_by(nParticipants,reps) %>% 
  slice(1)

ggplot(PowerfulDataframe, aes(reps,Power_PSE, color = as.factor(nParticipants))) +
  geom_line() +
  xlab("Number of trials per staircase") +
  ylab("Power") +
  scale_color_discrete(name = "Number of\nParticipants")