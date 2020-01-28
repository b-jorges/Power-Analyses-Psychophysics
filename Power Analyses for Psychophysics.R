###Pull the whole repository
require(dplyr)
require(tidyverse)
require(lme4)

theme_set(theme_cowplot())
set.seed(912)

#####This is the function from which I choose which comparison values we put into the modelling part: 
#super high curtosis because we will get many values around the PSE
DataFrame3 = data.frame(x = rcauchy(1000,1,0.02),
                        y = rnorm(1000,1,0.1))
ggplot(DataFrame3, aes(x)) + ####just to get an idea of what this function looks like 
  geom_density() +
  coord_cartesian(xlim=c(0.5,1.5))

ID = paste0("s",1:16) #For how many subjects do you want to simulate the power?
Condition1 = c(0,1) #Values of a categorical variable of interest (e. g. high contrast and low contrast)
Condition2 = c(-8,-6.6, 6.6,8) #intensities of the comparison stimuli
reps = seq(1,55,1) #how many trials do we expect for each condition?
Mean_Difference = 1/8 #by how much does the PSE differ between test and comparison conditions?
Mean_Difference = 1/3 #by how much does the JND differ between test and comparison conditions?
Mean_Standard = Condition2 #Mean of the cummulative Gaussian in the standard condition
SD_Standard = Condition2*0.1 #SD of the cummulative Gaussian in the standard condition
SD_ResponseFuntion = 0.04 #Standard deviation of the distribution we draw the answers from (normal distribution), or its scale (Cauchy distribution)

####This function serves to simulate a full dataset of responses, which then can serve to fit psychometric functions
SimulatePsychometricFunction_Staircase = function(ID, Condition1, Condition2, reps, PSE_Difference, JND_Difference, Mean_Standard, SD_Standard, SD_ResponseFuntion){
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
      
      #Now we draw stimulus strengths. For the typical staircase, they are not uniformly distributed, but rather concentrated around the PSE
      #A normal distribution with a low standard deviation might be appropriate, or a normal function with an extremely high kurtosis, 
      #another option is the Cauchy function;
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


ggplot(Psychometric,aes(Difference,Answer,color=Condition1)) +
  binomial_smooth() +
  facet_grid(velH~ID) +
  xlab("Ratio Target/Comparison") +
  ylab("Probability Target Bigger") +
  coord_cartesian(xlim = c(-5,-1))


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


######power for very low estimates of effect size: difference of 0.1 m/s in PSEs, and JNDs 1/4 higher when self-motion is simulated
PowerPerN_Precision = c()

NumbersOfSubjects = c(10,12,14,16,18,20)

for (i in NumbersOfSubjects){
  ID = paste0("s",1:i)
  Power_Precision = c()
  nIterations = 500
  out <- replicate(nIterations, {
    Analyze_Pychometric_Precision(SimulatePsychometricFunction(ID=ID, Motion=Motion, velH=velH, reps=reps, PSE_Diff = 1/8, JND_Diff = 0.025))})
  hist(out) ###Distribution of p values
  Power_Precision = mean(out < 0.05) ###Power is the times the difference between the two models is significant
  PowerPerN_Precision = c(PowerPerN_Precision,Power_Precision) ###This is the power for the JNDs
  paste0("For ", i, " subjects: ", Power_Precision)
}

PowerPerN_Accuracy = c()
for (i in NumbersOfSubjects){
  ID = paste0("s",1:i)
  Power_Accuracy = c()
  nIterations = 500
  out2 <- replicate(nIterations, {
    Analyze_Pychometric_Accuracy(SimulatePsychometricFunction(ID=ID, Motion=Motion, velH=velH, reps=reps, PSE_Diff = 1/8, JND_Diff = 0.025))})
  hist(out2) ###Distribution of p values
  Power_Accuracy = mean(out2 < 0.05) ###Power is the times the difference between the two models is significant
  PowerPerN_Accuracy = c(PowerPerN_Accuracy,Power_Accuracy) ###This is the power for the PSEs
  paste0("For ", i, " subjects: ", Power_Accuracy)
}