###Pull the whole repository
require(dplyr)
require(tidyverse)
require(lme4)
require(ggplot2)
require(cowplot)
theme_set(theme_cowplot())
require(quickpsy)


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
set.seed(9121)


ID = paste0("s",1:5)
ConditionOfInterest = c(0,1)
StandardValues = c(2,4)
reps = seq(1,55,1)
PSE_Difference = 0.0
JND_Difference = 0.3
Multiplicator_PSE_Standard = 0
Multiplicator_SD_Standard = 0.108
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

ResponseDistributions = data.frame(
  Value=c(rcauchy(1650,1,0.05),
          rnorm(1650,1,0.1),
          rep(c(0.7,0.85,1,1.15,1.3),1650/5)),
  label = c(rep("Cauchy",1650),
            rep("Normal",1650),
            rep("Uniform",1650))
) %>% filter(abs(Value-1) < 0.5)

ggplot(ResponseDistributions %>% filter(label %in% c("Cauchy","Normal")), aes(Value,color = label)) +
  geom_density(size=2) +
#  coord_cartesian(xlim=c(0.5,1.5)) +
  xlab("Stimulus Intensity") +
  ylab("Density") +
  scale_color_manual(name = "Distribution\nType",
                     values = c(Red,BlauUB),
                     labels = c("Cauchy","Gaussian"))
ggsave("Figure1 Distributions.jpg", w = 6, h = 4)

if (Type_ResponseFunction == "normal"){
  
  Psychometric = Psychometric %>%
    mutate(
      staircase_factor = pnorm(length(reps),1,SD_ResponseFunction))
  
} else if (Type_ResponseFunction == "Cauchy"){
  
  Psychometric = Psychometric %>%
    mutate(
      staircase_factor = rcauchy(length(reps),1,SD_ResponseFunction))
  
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

PsychometricFunctions = quickpsy(Psychometric,Difference,Answer,grouping = .(ConditionOfInterest,ID,StandardValues), bootstrap = "none")

plot(PsychometricFunctions) +
  scale_color_manual(name = "",
                     values = c(Red,BlauUB),
                     labels = c("Control","Experimental")) +
  xlab("Difference between Comparison and Test") +
  ylab("Probability to choose Test") +
  geom_vline(linetype = 2, xintercept = 0, color = "grey") +
  geom_hline(linetype = 2, yintercept = 0.5, color = "grey")
ggsave("Figure02.jpg", w = 10, h = 5)


GLMM_Accuracy = glmer(cbind(Yes, Total - Yes) ~ ConditionOfInterest + (Difference  | ID)  + (Difference  | StandardValues),
                      family = binomial(link = "probit"), 
                      data = Psychometric,
                      nAGQ = 0,
                      control = glmerControl(optimizer = "nloptwrap"))

require(lmerTest)
summary(GLMM_Accuracy)$coefficients


GLMM_Precision = glmer(cbind(Yes, Total - Yes) ~ ConditionOfInterest*Difference + (1| ID) + (1| StandardValues), 
                       family = binomial(link = "probit"), 
                       data = Psychometric,
                       nAGQ = 0,
                       control = glmerControl(optimizer = "nloptwrap"))

require(lmerTest)
summary(GLMM_Precision)$coefficients


SimulatePsychometricFunction_Staircase = function(ID, ConditionOfInterest, StandardValues, reps, PSE_Difference, JND_Difference, 
                                                  Multiplicator_PSE_Standard, Multiplicator_SD_Standard, SD_ResponseFunction, Mean_Variability_Between = 0.1, SD_Variability_Between = 0.1){
  Psychometric = expand.grid(ID=ID, ConditionOfInterest=ConditionOfInterest, StandardValues=StandardValues, reps = reps)
  
  Psychometric = Psychometric %>%
    group_by(ID) %>%#
    mutate(PSE_Factor_ID = rnorm(1,1,Mean_Variability_Between),
           SD_Factor_ID = rnorm(1,1,SD_Variability_Between))
  
  Psychometric = Psychometric %>%
    mutate(
      Mean_Standard = StandardValues+StandardValues*Multiplicator_PSE_Standard,
      SD_Standard = StandardValues*Multiplicator_SD_Standard,
      Mean = (Mean_Standard + (ConditionOfInterest==1)*StandardValues*PSE_Difference)*PSE_Factor_ID,
      SD = abs((SD_Standard + (ConditionOfInterest==1)*SD_Standard*JND_Difference)*SD_Factor_ID),
      staircase_factor = rcauchy(length(reps),1,SD_ResponseFunction), 
      Presented_TestStimulusStrength = Mean*staircase_factor,
      Difference = Presented_TestStimulusStrength - StandardValues,
      AnswerProbability = pnorm(Presented_TestStimulusStrength,Mean,SD),
      Answer = as.numeric(rbernoulli(length(AnswerProbability),AnswerProbability))
    )
  
  Psychometric = Psychometric %>%
    filter(abs(staircase_factor-1) < 0.75) %>%
    group_by(ID,ConditionOfInterest,StandardValues,Difference) %>%
    mutate(Yes = sum(Answer==1),
           Total = length(ConditionOfInterest))
  
  Psychometric
}

Analyze_Pychometric_Accuracy_GLMM = function(Psychometric){
  
  TimeBeginning = Sys.time()
  
  GLMM_Accuracy = glmer(cbind(Yes, Total - Yes) ~ ConditionOfInterest + (Difference  | ID)  + (Difference  | StandardValues),
                        family = binomial(link = "probit"), 
                        data = Psychometric,
                        nAGQ = 0,
                        control = glmerControl(optimizer = "nloptwrap"))
  
  p = summary(GLMM_Accuracy)$coefficients[8]
  
  #print(TimeBeginning - Sys.time()) ###This is two show how long each iteration takes
  #print(p)
  
  p
}

Analyze_Pychometric_Precision_GLMM = function(Psychometric){
  
  TimeBeginning = Sys.time()
  
  GLMM_Precision = glmer(cbind(Yes, Total - Yes) ~ ConditionOfInterest*Difference + (Difference  | ID) + (Difference  | StandardValues), 
                         family = binomial(link = "probit"), 
                         data = Psychometric,
                         nAGQ = 0,
                         control = glmerControl(optimizer = "nloptwrap"))
  
  p = summary(GLMM_Precision)$coefficients[16]
  
  
  #print(p)
  
  p
}


####################################################################################
##################Getting power for a couple for GLMMs##############################
####################################################################################

NumbersOfSubjects = c(10,12,14,16,18,20) #for how many subjects are we simulating the power?

nIterations = 50
pvalue = 0.05

Power = data.frame()
for (i in NumbersOfSubjects){
  
  ID = paste0("s",1:i)
  TimeBeginning = Sys.time()
  Dataframe_Temp = c()
  
  for (j in 1:nIterations){
    
    Dataframe = SimulatePsychometricFunction_Staircase(ID, ConditionOfInterest, StandardValues, reps, PSE_Difference, JND_Difference, Multiplicator_PSE_Standard, Multiplicator_SD_Standard, SD_ResponseFunction, Mean_Variability_Between, SD_Variability_Between)
    
    Parameters = GetParametersOfPsychometricFunction(Dataframe)
    
    p = c(Analyze_Pychometric_Accuracy_GLMM(Dataframe),
          Analyze_Pychometric_Precision_GLMM(Dataframe))
    
    Dataframe_Temp = rbind(Dataframe_Temp,p)
    
    if ((j/25) %in% 1:(nIterations/25)){ #print the number of the current iteration every 50 trials
      (print(j))
    }
  }
  
  Power = rbind(Power,
                data.frame(value = c(mean(Dataframe_Temp[,1] < pvalue),
                                     mean(Dataframe_Temp[,2] < pvalue)),
                           label = c("Accuracy GLMM",
                                     "Precision GLMM"),
                           nSubjects = i))
  
  print(paste0("This iteration has taken ", Sys.time() - TimeBeginning))  ###This is two show how long each iteration takes
  print(paste0("Accuracy GLMM for ", i, " subjects: ", mean(Dataframe_Temp[,1] < pvalue))) #outputs an estimate of the power for each n
  print(paste0("Precision GLMM for ", i, " subjects: ", mean(Dataframe_Temp[,2] < pvalue))) #outputs an estimate of the power for each n
}

ggplot(Power,aes(nSubjects,value, color = label)) +
  geom_line(size = 2) +
  xlab("N° of Subjects") +
  ylab("Power") +
  scale_color_manual(values = c(BlauUB,Red)) +
  geom_hline(linetype = 2, yintercept = 0.8) +
  geom_hline(linetype = 3, yintercept = 0.95) +
  scale_x_continuous(breaks=c(10,12,14,16,18,20))

####################################################################################
##################Comparing Two-Level approach and GLMMs############################
####################################################################################
Parameters = quickpsy(Psychometric,Difference,Answer,grouping = .(ID,ConditionOfInterest,StandardValues), bootstrap = "none")$par
plot(quickpsy(Psychometric,Difference,Answer,grouping = .(ConditionOfInterest,ID,StandardValues), bootstrap = "none"))
Parameters2 = Parameters %>%
  filter(parn == "p1") %>%
  select(ID,ConditionOfInterest,Mean=par, StandardValues)
Parameters2$SD = Parameters$par[Parameters$parn == "p2"]
Parameters = Parameters2

ANOVA_Mean = aov(Mean ~ as.factor(ConditionOfInterest)*StandardValues,Parameters)
summary(ANOVA_Mean)

ANOVA_SD = aov(SD ~ as.factor(ConditionOfInterest)*StandardValues,Parameters)
summary(ANOVA_SD)

GetParametersOfPsychometricFunction = function(Psychometric){
  Parameters = quickpsy(Psychometric,Difference,Answer,grouping = .(ID,ConditionOfInterest,StandardValues), bootstrap = "none")$par
  
  Parameters2 = Parameters %>%
    filter(parn == "p1") %>%
    select(ID,ConditionOfInterest,Mean=par, StandardValues)
  Parameters2$SD = Parameters$par[Parameters$parn == "p2"]
  Parameters2
}

Analyze_Pychometric_Accuracy_2Level = function(Parameters){
  
  ANOVA_Mean = aov(Mean ~ as.factor(ConditionOfInterest)*StandardValues,Parameters)
  Coefficients = summary(ANOVA_Mean)[[1]]
  Coefficients$`Pr(>F)`[1]
}
Analyze_Pychometric_Precision_2Level = function(Parameters){
  
  ANOVA_SD = aov(SD ~ as.factor(ConditionOfInterest)*StandardValues,Parameters)
  Coefficients = summary(ANOVA_SD)[[1]]
  Coefficients$`Pr(>F)`[1]
}

Analyze_Pychometric_Accuracy_2Level_LMM = function(Parameters){
  
  mod = lmer(Mean ~ ConditionOfInterest + (1  | ID) + (1  | StandardValues), 
        data = Parameters)
  summary(mod)$coefficients[10]
}

Analyze_Pychometric_Precision_2Level_LMM = function(Parameters){
  
  mod = lmer(SD ~ ConditionOfInterest + (1  | ID) + (1  | StandardValues), 
             data = Parameters)
  summary(mod)$coefficients[10]
}

Power = data.frame()
nIterations = 200
pvalue = 0.05
NumbersOfSubjects = c(10,12,14,16,18,20)
#NumbersOfSubjects = c(2,3,4,5,6)

Mean_Variability_Between = 0.1
SD_Variability_Between = 0.1

ComparePowers = function(ConditionOfInterest, StandardValues, reps, PSE_Difference, JND_Difference, 
                         Multiplicator_PSE_Standard, Multiplicator_SD_Standard, SD_ResponseFunction, Mean_Variability_Between = 0.1, SD_Variability_Between = 0.1,
                         NumbersOfSubjects){
  for (i in NumbersOfSubjects){
    
    ID = paste0("s",1:i)
    TimeBeginning = Sys.time()
    Dataframe_Temp = c()
    
    for (j in 1:nIterations){
      
      Dataframe = SimulatePsychometricFunction_Staircase(ID, ConditionOfInterest, StandardValues, reps, PSE_Difference, JND_Difference, 
                                                         Multiplicator_PSE_Standard, Multiplicator_SD_Standard, SD_ResponseFunction, Mean_Variability_Between = 0.1, SD_Variability_Between = 0.1)
      
      Parameters = GetParametersOfPsychometricFunction(Dataframe)
      
      p = c(Analyze_Pychometric_Accuracy_GLMM(Dataframe),
            Analyze_Pychometric_Precision_GLMM(Dataframe),
            Analyze_Pychometric_Accuracy_2Level(Parameters),    
            Analyze_Pychometric_Precision_2Level(Parameters),
            Analyze_Pychometric_Accuracy_2Level_LMM(Parameters),
            Analyze_Pychometric_Precision_2Level_LMM(Parameters))
      
      
      Dataframe_Temp = rbind(Dataframe_Temp,p)
      
      if ((j/25) %in% 1:40){
        (print(j))
      }
    }
    
    Power = rbind(Power,
                  data.frame(value = c(mean(Dataframe_Temp[,1] < pvalue),
                                       mean(Dataframe_Temp[,2] < pvalue),
                                       mean(Dataframe_Temp[,3] < pvalue),
                                       mean(Dataframe_Temp[,4] < pvalue),
                                       mean(Dataframe_Temp[,5] < pvalue),
                                       mean(Dataframe_Temp[,6] < pvalue)),
                             label = c("Accuracy GLMM",
                                       "Precision GLMM",
                                       "Accuracy Two-Level",
                                       "Precision Two-Level",
                                       "Accuracy Two-Level LMM",
                                       "Precision Two-Level LMM"),
                             reps = reps[length(reps)],
                             PSE_Difference = PSE_Difference,
                             JND_Difference = JND_Difference,
                             StandardValues = paste0(StandardValues[1],StandardValues[length(StandardValues)]),
                             nStandardValues = length(StandardValues),
                             TrialsPerSubject = length(StandardValues)*length(reps)*
                             SD_ResponseFunction = SD_ResponseFunction,
                             nSubjects = i))
    
    print(paste0("This iteration has taken ", Sys.time() - TimeBeginning))  ###This is two show how long each iteration takes
    print(paste0("Accuracy GLMM for ", i, " subjects: ", mean(Dataframe_Temp[,1] < pvalue))) #outputs an estimate of the power for each n
    print(paste0("Precision GLMM for ", i, " subjects: ", mean(Dataframe_Temp[,2] < pvalue))) #outputs an estimate of the power for each n
    print(paste0("Accuracy 2Level for ", i, " subjects: ", mean(Dataframe_Temp[,3] < pvalue))) #outputs an estimate of the power for each n
    print(paste0("Precision 2Level for ", i, " subjects: ", mean(Dataframe_Temp[,4] < pvalue))) #outputs an estimate of the power for each n
    print(paste0("Accuracy 2Level LMM for ", i, " subjects: ", mean(Dataframe_Temp[,5] < pvalue))) #outputs an estimate of the power for each n
    print(paste0("Precision 2Level LMM for ", i, " subjects: ", mean(Dataframe_Temp[,6] < pvalue))) #outputs an estimate of the power for each n
  }

Power
}

Powers1 = ComparePowers(ConditionOfInterest, StandardValues, reps = 1:55, PSE_Difference=-0.05, JND_Difference=-0.05, 
                        Multiplicator_PSE_Standard, Multiplicator_SD_Standard, SD_ResponseFunction, 
                        Mean_Variability_Between = 0.1, SD_Variability_Between = 0.1, 
                        NumbersOfSubjects)
write.csv(Powers1,"Powers1.csv") #small negative differences

Powers2 = ComparePowers(ConditionOfInterest, StandardValues, reps = 1:25, PSE_Difference = 0.1, JND_Difference = 0.2, 
                        Multiplicator_PSE_Standard, Multiplicator_SD_Standard, SD_ResponseFunction, 
                        Mean_Variability_Between = 0.1, SD_Variability_Between = 0.1, 
                        NumbersOfSubjects)
write.csv(Powers2,"Powers2.csv") #not super big differences, fewer trials

Powers3 = ComparePowers(ConditionOfInterest, StandardValues = c(2,4), reps = 1:55, PSE_Difference = 0.1, JND_Difference = 0.3, 
                        Multiplicator_PSE_Standard = 0.1, Multiplicator_SD_Standard, SD_ResponseFunction, 
                        Mean_Variability_Between, SD_Variability_Between, 
                        NumbersOfSubjects)
write.csv(Powers3,"Powers3.csv") #only two standard variables

Powers4 = ComparePowers(ConditionOfInterest, StandardValues = c(2,4), reps = 1:25, PSE_Difference = 0.2, JND_Difference = 0.3, 
                        Multiplicator_PSE_Standard, Multiplicator_SD_Standard, SD_ResponseFunction = 0.1, 
                        Mean_Variability_Between, SD_Variability_Between, 
                        NumbersOfSubjects)
write.csv(Powers4,"Powers4.csv") #only two standard and fewer trials

Powers5 = ComparePowers(ConditionOfInterest, StandardValues, reps = 1:25, PSE_Difference=0, JND_Difference=0.5, 
                        Multiplicator_PSE_Standard, Multiplicator_SD_Standard, SD_ResponseFunction, 
                        Mean_Variability_Between = 0.1, SD_Variability_Between = 0.1, 
                        NumbersOfSubjects)
write.csv(Powers5,"Powers5.csv") #huge difference in JND, no difference in PSE

Powers6 = ComparePowers(ConditionOfInterest, StandardValues, reps = 1:25, PSE_Difference = 0.2, JND_Difference = 2, 
                        Multiplicator_PSE_Standard, Multiplicator_SD_Standard, SD_ResponseFunction, 
                        Mean_Variability_Between = 0.1, SD_Variability_Between = 0.1, 
                        NumbersOfSubjects)
write.csv(Powers6,"Powers6.csv") #huge ratio between PSE and JND

Powers7 = ComparePowers(ConditionOfInterest, StandardValues, reps = 1:55, PSE_Difference = 0.1, JND_Difference = 0.3, 
                        Multiplicator_PSE_Standard = 0.1, Multiplicator_SD_Standard, SD_ResponseFunction = 0.2, 
                        Mean_Variability_Between, SD_Variability_Between, 
                        NumbersOfSubjects)
write.csv(Powers3,"Powers7.csv") #broader response function

Powers8 = ComparePowers(ConditionOfInterest, StandardValues, reps = 1:25, PSE_Difference = 0.1, JND_Difference = 0.3, 
                        Multiplicator_PSE_Standard, Multiplicator_SD_Standard, SD_ResponseFunction = 0.2, 
                        Mean_Variability_Between, SD_Variability_Between, 
                        NumbersOfSubjects)
write.csv(Powers4,"Powers8.csv") #broader response function, fewer trials

ggplot(Power,aes(nSubjects,value, color = label)) +
  geom_line(size = 2) +
  xlab("N° of Subjects") +
  ylab("Power") +
  scale_color_manual(values = c(BlauUB,Red,LightBlauUB,LightRed)) +
  geom_hline(linetype = 2, yintercept = 0.8) +
  geom_hline(linetype = 3, yintercept = 0.95) +
  scale_x_continuous(breaks=c(10,12,14,16,18,20))
