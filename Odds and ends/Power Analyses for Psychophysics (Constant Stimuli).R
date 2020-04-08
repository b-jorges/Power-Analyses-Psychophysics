SimulatePsychometricFunction_ConstantStimuli = function(ID, Condition1, Condition2, reps, PSE_Difference, JND_Difference, Mean_Standard, SD_Standard){
  Psychometric = expand.grid(ID=ID, Condition1=Condition1, Condition2=Condition2, reps = reps)
  
  Psychometric = Psychometric %>%
    group_by(ID) %>%#
    mutate(PSE_Factor_ID = rnorm(1,1,0.1), #how much variability is in the means of the psychometric functions between subjects?
           SD_Factor_ID = rnorm(1,1,0.1)) #how much variability is in the standard deviations of the psychometric functions between subjects?
  
  #Next, we are going to simulating PSEs and JNDs for each condition
  Psychometric = Psychometric %>%
    mutate(
      #PSE is calculated as the strength of the comparison stimulus; in the test condition, we add the 
      Mean = (Condition2 + (Condition==1)*PSE_Difference)*PSE_Factor_ID,
      SD = (SD_Standard + (Condition==1)*JND_Difference)*SD_Factor_ID,
      
      #now we draw the comparison stimulus strengths according to 
      staircase_factor = rcauchy(length(reps),1,SD_ResponseFuntion), #here we draw a multiplier for each trial 
      Presented_TestStimulusStrength = Mean*staircase_factor, ###translates from values around 1 to values around PSE
      
      #We assume that the responses can be adequately captured by a Cummulative Gaussian Distribution ... 
      #with the Means and SDs we simulated above.
      #This line gives us the response probability for each trial:
      AnswerProbability = pnorm(Presented_TestStimulusStrength,Mean,SD),
      
      ####Get difference between Standard stimulus strength and Test stimulus strength
      Difference = Presented_TestStimulusStrength - Mean,
      
      ##get binary answers ("Test was stronger" yes/no) from probabilities for each trial
      Answer = as.numeric(rbernoulli(length(AnswerProbability),AnswerProbability))
    )
  
  
  ###prepare for glmer() - needs sum of YES/Total per stimulus strength and condition
  Psychometric = Psychometric %>%
    group_by(ID,Motion,velH,Difference) %>%
    mutate(Yes = sum(Answer==1),
           Total = length(Motion))
  
  Psychometric #output above data frame
}


Analyze_Pychometric_Precision = function(Psychometric){
  
  TimeBeginning = Sys.time()
  
  ###precision = slope for Difference (the steeper the slope, the more sensitive) ... 
  ###so if there is an interaction between Motion Condition and Difference, that means that it changes the slope
  mod1 = glmer(cbind(Yes, Total - Yes) ~ as.factor(Condition1)*Difference + (Difference  | ID) + (Difference  | Condition2), 
               family = binomial(link = "probit"), 
               data = Psychometric)
  mod2 = glmer(cbind(Yes, Total - Yes) ~ as.factor(Condition1) + Difference + (Difference  | ID) + (Difference  | Condition2),
               family = binomial(link = "probit"), 
               data = Psychometric)
  
  print(TimeBeginning - Sys.time())
  
  p = anova(mod1,mod2)$`Pr(>Chisq)`[2] ##Model 1 beats model 2
  
  print(p)
  p
}

Analyze_Pychometric_Accuracy = function(Psychometric){
  
  TimeBeginning = Sys.time()
  
  mod1 = glmer(cbind(Yes, Total - Yes) ~ Condition1 + (Difference  | ID)  + (Difference  | Condition2),
               family = binomial(link = "probit"), 
               data = Psychometric)
  
  mod2 = glmer(cbind(Yes, Total - Yes) ~ (Difference  | ID) + (Difference  | Condition2),
               family = binomial(link = "probit"), 
               data = Psychometric)
  
  print(TimeBeginning - Sys.time())
  
  p = anova(mod1,mod2)$`Pr(>Chisq)`[2] ##Model 1 beats model 2
  
  print(p)
  p
  
}

NumbersOfSubjects = c(10,12,14,16,18,20) #for how many subjects are we simulating the power?

PowerPerN_Precision = c()
for (i in NumbersOfSubjects){
  ID = paste0("s",1:i)
  
  Power_Precision = c()
  
  nIterations = 500 #how many data frames should we simulate per number of subjects? more data frames makes for a more reliable estimate of the power
  
  out <- replicate(nIterations, { #this function performs the number of iterations indicated above and stores the result for each in "out"
    Analyze_Pychometric_Precision(SimulatePsychometricFunction_Staircase(ID=ID, Motion=Motion, velH=velH, reps=reps, PSE_Diff = 1/8, JND_Diff = 0.025))})
  hist(out) ###Distribution of p values
  
  Power_Precision = mean(out < 0.05) ###Power is the times the difference between the two models is significant
  
  PowerPerN_Precision = c(PowerPerN_Precision,Power_Precision) ###This is a vector with the power for each n for the JNDs
  
  paste0("For ", i, " subjects: ", Power_Precision) #gives an estimate of the power for each n
}

PowerPerN_Accuracy = c()
for (i in NumbersOfSubjects){
  ID = paste0("s",1:i)
  
  Power_Accuracy = c()
  
  nIterations = 500 #how many data frames should we simulate per number of subjects? more data frames makes for a more reliable estimate of the power
  
  out2 <- replicate(nIterations, { #this function performs the number of iterations indicated above and stores the result for each in "out2"
    Analyze_Pychometric_Accuracy(SimulatePsychometricFunction_Staircase(ID=ID, Motion=Motion, velH=velH, reps=reps, PSE_Diff = 1/8, JND_Diff = 0.025))})
  
  hist(out2) ###Distribution of p values
  
  Power_Accuracy = mean(out2 < 0.05) ###Power is the times the difference between the two models is significant
  
  PowerPerN_Accuracy = c(PowerPerN_Accuracy,Power_Accuracy) ###This is a vector with the power for each n for the PSEs
  
  paste0("For ", i, " subjects: ", Power_Accuracy) #prints estimate of the power for each n
}