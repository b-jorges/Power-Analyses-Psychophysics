###Pull the whole repository
require(dplyr)
require(tidyverse)
require(lme4)
require(ggplot2)
require(cowplot)
theme_set(theme_cowplot())

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
set.seed(912)

ID = paste0("s",1:10) #For how many subjects do you want to simulate the power?
ConditionOfInterest = c(0,1) #Values of a categorical variable of interest (e. g. high contrast and low contrast)
SecondaryCondition = c(6.6, 8, 10) #intensities of the comparison stimuli
reps = seq(1,55,1) #how many trials do we expect for each condition?
PSE_Difference = 1/8 #by how much does the PSE differ between test and comparison conditions, in units of the comparison PSE
JND_Difference = 1/3 #by how much does the JND differ between test and comparison conditions, in units of the comparison JND
PSE_Standard = SecondaryCondition #Mean of the cummulative Gaussian in the standard condition
SD_Standard = SecondaryCondition*0.1 #SD of the cummulative Gaussian in the standard condition, see "XXXXX.r" on how to get standard deviations from Weber Fractions
SD_ResponseFuntion = 0.06 #Standard deviation of the distribution we draw the presented stimulus strengths from (normal distribution), or its scale (Cauchy distribution)

####This function serves to simulate a full dataset of responses, which then can serve to fit psychometric functions
SimulatePsychometricFunction_Staircase = function(ID, ConditionOfInterest, SecondaryCondition, reps, PSE_Difference, JND_Difference, Mean_Standard, SD_Standard, SD_ResponseFuntion){
  Psychometric = expand.grid(ID=ID, ConditionOfInterest=ConditionOfInterest, SecondaryCondition=SecondaryCondition, reps = reps)
  
  Psychometric = Psychometric %>%
    group_by(ID) %>%#
    mutate(PSE_Factor_ID = rnorm(1,1,0.1), #how much variability is in the means of the psychometric functions between subjects?
           SD_Factor_ID = rnorm(1,1,0.1)) #how much variability is in the standard deviations of the psychometric functions between subjects?
 
  #Next, we are going to simulating PSEs and JNDs for each condition
  Psychometric = Psychometric %>%
    mutate(
      #PSE is calculated as the strength of the comparison stimulus; in the test condition, we add the 
      Mean = (SecondaryCondition + (ConditionOfInterest==1)*SecondaryCondition*PSE_Difference)*PSE_Factor_ID,
      SD = abs((SD_Standard + (ConditionOfInterest==1)*SD_Standard*JND_Difference)*SD_Factor_ID),
      
      #Now we draw stimulus strengths. For the typical staircase, they are not uniformly distributed, but rather concentrated around the PSE
      #A normal distribution with a low standard deviation might be appropriate, or a normal function with an extremely high kurtosis, 
      #another option is the Cauchy function;
      staircase_factor = rcauchy(length(reps),1,SD_ResponseFuntion), #here we draw a multiplier for each trial 
      Presented_TestStimulusStrength = Mean*staircase_factor, ###translates from values around 1 to values around PSE

      ####Get difference between Standard stimulus strength and Test stimulus strength
      Difference = Presented_TestStimulusStrength - SecondaryCondition,
      
      #We assume that the responses can be adequately captured by a Cummulative Gaussian Distribution ... 
      #with the Means and SDs we simulated above.
      #This line gives us the response probability for each trial:
      AnswerProbability = pnorm(Presented_TestStimulusStrength,Mean,SD),
      
      ##get binary answers ("Test was stronger" yes/no) from probabilities for each trial
      Answer = as.numeric(rbernoulli(length(AnswerProbability),AnswerProbability))
    )

  ###prepare for glmer() - needs sum of YES/Total per stimulus strength and condition
  Psychometric = Psychometric %>%
    filter(abs(Difference) < 0.5*abs(Mean)) %>%
    group_by(ID,ConditionOfInterest,SecondaryCondition,Difference) %>%
    mutate(Yes = sum(Answer==1),
           Total = length(ConditionOfInterest))
  
  Psychometric #output above data frame
}

Psychometric = SimulatePsychometricFunction_Staircase(ID, ConditionOfInterest, SecondaryCondition, reps, PSE_Difference, JND_Difference, Mean_Standard, SD_Standard, SD_ResponseFuntion)

ggplot(Psychometric, aes(Difference, Answer, color = as.factor(ConditionOfInterest))) +
  binomial_smooth() +
  facet_grid(SecondaryCondition~ID) +
  geom_vline(xintercept = 0, color = "grey") +
  geom_hline(yintercept = 0.5, color = "grey") +
  xlab("Difference between Comparison and Test") +
  ylab("Probability to choose Test")

Analyze_Pychometric_Precision = function(Psychometric){

  TimeBeginning = Sys.time()
  
  ###precision = slope for Difference (the steeper the slope, the more sensitive) ... 
  ###so if there is an interaction between Motion Condition and Difference, that means that it changes the slope
  mod1 = glmer(cbind(Yes, Total - Yes) ~ as.factor(ConditionOfInterest)*Difference + (Difference  | ID) + (Difference  | SecondaryCondition), 
               family = binomial(link = "probit"), 
               data = Psychometric,
               nAGQ = 0,
               control = glmerControl(optimizer = "nloptwrap"))
  mod2 = glmer(cbind(Yes, Total - Yes) ~ as.factor(ConditionOfInterest) + Difference + (Difference  | ID) + (Difference  | SecondaryCondition),
               family = binomial(link = "probit"), 
               data = Psychometric,
               nAGQ = 0,
               control = glmerControl(optimizer = "nloptwrap"))

  p = anova(mod1,mod2)$`Pr(>Chisq)`[2] ##Model 1 beats model 2

  print(TimeBeginning - Sys.time())  
  print(p)
  
  p
}


Analyze_Pychometric_Accuracy = function(Psychometric){

  TimeBeginning = Sys.time()
  
  mod1 = glmer(cbind(Yes, Total - Yes) ~ ConditionOfInterest + (Difference  | ID)  + (Difference  | SecondaryCondition),
               family = binomial(link = "probit"), 
               data = Psychometric,
               nAGQ = 0,
               control = glmerControl(optimizer = "nloptwrap"))
  
  mod2 = glmer(cbind(Yes, Total - Yes) ~ (Difference  | ID) + (Difference  | SecondaryCondition),
               family = binomial(link = "probit"), 
               data = Psychometric,
               nAGQ = 0,
               control = glmerControl(optimizer = "nloptwrap"))

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
    Analyze_Pychometric_Precision(SimulatePsychometricFunction_Staircase(ID, ConditionOfInterest, SecondaryCondition, reps, PSE_Difference, JND_Difference, Mean_Standard, SD_Standard, SD_ResponseFuntion))})
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
    Analyze_Pychometric_Accuracy(SimulatePsychometricFunction_Staircase(ID, ConditionOfInterest, SecondaryCondition, reps, PSE_Difference, JND_Difference, Mean_Standard, SD_Standard, SD_ResponseFuntion))})

  hist(out2) ###Distribution of p values

  Power_Accuracy = mean(out2 < 0.05) ###Power is the times the difference between the two models is significant
  
  PowerPerN_Accuracy = c(PowerPerN_Accuracy,Power_Accuracy) ###This is a vector with the power for each n for the PSEs
  
  paste0("For ", i, " subjects: ", Power_Accuracy) #prints estimate of the power for each n
}
