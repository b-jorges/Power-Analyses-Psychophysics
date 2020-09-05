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
require(DHARMa)

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


ID = paste0("s",1:5)
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

GLMM1.5 = glmer(cbind(Yes, Total - Yes) ~ ConditionOfInterest*Difference + (ConditionOfInterest+Difference| ID) + (ConditionOfInterest+Difference| StandardValues), 
             family = binomial(link = "probit"), 
             data = Psychometric,
             nAGQ = 0,
             control = glmerControl(optimizer = "nloptwrap"))

GLMM2 = glmer(cbind(Yes, Total - Yes) ~ ConditionOfInterest*Difference + (Difference| ID) + (Difference| StandardValues), 
             family = binomial(link = "logit"), 
             data = Psychometric,
             nAGQ = 0,
             control = glmerControl(optimizer = "nloptwrap"))

GLMM3 = glmer(cbind(Yes, Total - Yes) ~ ConditionOfInterest*Difference + (1| ID) + (1| StandardValues), 
              family = binomial(link = "logit"), 
              data = Psychometric,
              nAGQ = 0,
              control = glmerControl(optimizer = "nloptwrap"))


haha = anova(GLMM,GLMM2)

haha$`Pr(>Chisq)`[2]

SummaryGLMM2 = summary(GLMM2)
SummaryGLMM = summary(GLMM)

SummaryGLMM$AICtab
SummaryGLMM2$AICtab

plot(GLMM)
plot(GLMM2)

hallo = anova(GLMM,GLMM2)
hallo


Psychometric %>% 
  group_by(ConditionOfInterest, StandardValues) %>% 
  slice(1)




PsychometricFunctions = quickpsy(Psychometric,Difference,Answer,grouping = .(ConditionOfInterest,ID,StandardValues), bootstrap = "none")
PsychometricFunctions$par

plot(PsychometricFunctions) +
  scale_color_manual(name = "",
                     values = c(Red,BlauUB),
                     labels = c("Control","Experimental")) +
  xlab("Difference between Comparison and Test") +
  ylab("Probability to choose Test") +
  geom_vline(linetype = 2, xintercept = 0, color = "grey") +
  geom_hline(linetype = 2, yintercept = 0.5, color = "grey")
ggsave("Figure02.jpg", w = 10, h = 5)

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
Parameters2 = Parameters %>%
  filter(parn == "p1") %>%
  select(ID,ConditionOfInterest,Mean=par, StandardValues)
Parameters2$SD = Parameters$par[Parameters$parn == "p2"]
Parameters = Parameters2

require(lmerTest)

ANOVA_Mean = aov(Mean ~ as.factor(ConditionOfInterest)*StandardValues,Parameters)
Pvalue_Mean_ANOVA = summary(ANOVA_Mean)[[1]][["Pr(>F)"]][1]
summary(GLMM)
summary(GLMM2)

LMM_Mean_Big = lmer(Mean ~ ConditionOfInterest*Difference + (ConditionOfInterest+Difference| ID) + (ConditionOfInterest+Difference| StandardValues),
              data = Psychometric)
LMM_Mean_Small = lmer(Mean ~ ConditionOfInterest*Difference + (Difference| ID) + (Difference| StandardValues),
                data = Psychometric)
LMM_Mean_VerySmall = lmer(Mean ~ ConditionOfInterest*Difference + (1| ID) + (1| StandardValues),
                      data = Psychometric)
Hello$residuals.GLMM.
Hello = data.frame(residuals(GLMM))


ggplot(data.frame(residuals(GLMM)),aes(residuals.GLMM.)) +
  geom_density()

ggplot(data.frame(residuals(GLMM2)),aes(residuals.GLMM2.)) +
  geom_density()

ggplot(data.frame(residuals(GLMM3)),aes(residuals.GLMM3.)) +
  geom_density()

ggplot(data.frame(residuals(LMM_Mean_Big)),aes(residuals.LMM_Mean_Big.)) +
  geom_density()
ggplot(data.frame(residuals(LMM_Mean_Small)),aes(residuals.LMM_Mean_Small.)) +
  geom_density()
ggplot(data.frame(residuals(LMM_Mean_VerySmall)),aes(residuals.LMM_Mean_VerySmall.)) +
  geom_density()

library(fitdistrplus)
Hello = descdist(data.frame(residuals(LMM_Mean_VerySmall))$residuals.LMM_Mean_VerySmall., discrete = FALSE)
Hello2 = descdist(data.frame(residuals(LMM_Mean_Small))$residuals.LMM_Mean_Small., discrete = FALSE)
Hello3 = descdist(data.frame(residuals(LMM_Mean_Big))$residuals.LMM_Mean_Big., discrete = FALSE)

plot(fitdist(data.frame(residuals(LMM_Mean_VerySmall))$residuals.LMM_Mean_VerySmall., "norm"))
plot(fitdist(data.frame(residuals(LMM_Mean_Small))$residuals.LMM_Mean_Small., "norm"))
plot(fitdist(data.frame(residuals(LMM_Mean_Big))$residuals.LMM_Mean_Big., "norm"))

Hello$kurtosis
Hello2$kurtosis
Hello3$kurtosis
Hello$skewness
Hello2$skewness
Hello3$skewness


ggplot(Hello,aes(residuals.GLMM.)) +
  geom_density()


ANOVA_SD = aov(SD ~ as.factor(ConditionOfInterest)*StandardValues,Parameters)
Pvalue_SD_ANOVA = summary(ANOVA_SD)[[1]][["Pr(>F)"]][1]


Power = data.frame()
nIterations = 200
pvalue = 0.05
NumbersOfSubjects = c(10,12,14,16,18,20)
#NumbersOfSubjects = c(2,3,4,5,6)





                           


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
write.csv(Powers7,"Powers7.csv") #broader response function

Powers8 = ComparePowers(ConditionOfInterest, StandardValues, reps = 1:25, PSE_Difference = 0.1, JND_Difference = 0.3, 
                        Multiplicator_PSE_Standard, Multiplicator_SD_Standard, SD_ResponseFunction = 0.2, 
                        Mean_Variability_Between, SD_Variability_Between, 
                        NumbersOfSubjects)
write.csv(Powers8,"Powers8.csv") #broader response function, fewer trials

Powers9 = ComparePowers(ConditionOfInterest, StandardValues, reps = 1:55, PSE_Difference = 0, JND_Difference = 0.3, 
                        Multiplicator_PSE_Standard = 0.1, Multiplicator_SD_Standard, SD_ResponseFunction = 0.1, 
                        Mean_Variability_Between = 0.2, SD_Variability_Between = 0.2, 
                        NumbersOfSubjects)
write.csv(Powers9,"Powers9.csv") #broader response function

Powers10 = ComparePowers(ConditionOfInterest, StandardValues, reps = 1:25, PSE_Difference = 0.15, JND_Difference = 0.3, 
                        Multiplicator_PSE_Standard, Multiplicator_SD_Standard, SD_ResponseFunction = 0.2, 
                        Mean_Variability_Between = 0.2, SD_Variability_Between = 0.2, 
                        NumbersOfSubjects)
write.csv(Powers10,"Powers10.csv") #broader response function, fewer trials

Powers11 = ComparePowers(ConditionOfInterest, StandardValues, reps = 1:55, PSE_Difference = 0, JND_Difference = 0.3, 
                        Multiplicator_PSE_Standard = 0.1, Multiplicator_SD_Standard, SD_ResponseFunction = 0.2, 
                        Mean_Variability_Between, SD_Variability_Between, 
                        NumbersOfSubjects)
write.csv(Powers11,"Powers11.csv") #broader response function

Powers12 = ComparePowers(ConditionOfInterest, StandardValues, reps = 1:25, PSE_Difference = 0, JND_Difference = 0.3, 
                        Multiplicator_PSE_Standard, Multiplicator_SD_Standard, SD_ResponseFunction = 0.2, 
                        Mean_Variability_Between, SD_Variability_Between, 
                        NumbersOfSubjects)
write.csv(Powers12,"Powers12.csv") #broader response function, fewer trials

Powers13 = ComparePowers(ConditionOfInterest, StandardValues, reps = 1:55, PSE_Difference = 0, JND_Difference = 0.2, 
                        Multiplicator_PSE_Standard = 0.1, Multiplicator_SD_Standard, SD_ResponseFunction = 0.1, 
                        Mean_Variability_Between = 0.2, SD_Variability_Between = 0.2, 
                        NumbersOfSubjects)
write.csv(Powers13,"Powers13.csv") #broader response function

Powers14 = ComparePowers(ConditionOfInterest, StandardValues, reps = 1:25, PSE_Difference = 0, JND_Difference = 0.2, 
                         Multiplicator_PSE_Standard, Multiplicator_SD_Standard, SD_ResponseFunction = 0.1, 
                         Mean_Variability_Between = 0.1, SD_Variability_Between = 0.1, 
                         NumbersOfSubjects)
write.csv(Powers14,"Powers14.csv") #broader response function, fewer trials

Powers15 = ComparePowers(ConditionOfInterest, StandardValues, reps = 1:55, PSE_Difference = -0.05, JND_Difference = 0.3, 
                         Multiplicator_PSE_Standard = 0.1, Multiplicator_SD_Standard, SD_ResponseFunction = 0.2, 
                         Mean_Variability_Between, SD_Variability_Between, 
                         NumbersOfSubjects)
write.csv(Powers15,"Powers15.csv") #broader response function

Powers16 = ComparePowers(ConditionOfInterest, StandardValues, reps = 1:25, PSE_Difference = -0.05, JND_Difference = 0.3, 
                         Multiplicator_PSE_Standard, Multiplicator_SD_Standard, SD_ResponseFunction = 0.2, 
                         Mean_Variability_Between, SD_Variability_Between, 
                         NumbersOfSubjects)
write.csv(Powers16,"Powers16.csv") #broader response function, fewer trials

Powers17 = ComparePowers(ConditionOfInterest, StandardValues, reps = 1:55, PSE_Difference = -0.05, JND_Difference = 0.2, 
                         Multiplicator_PSE_Standard = 0.1, Multiplicator_SD_Standard, SD_ResponseFunction = 0.1, 
                         Mean_Variability_Between = 0.2, SD_Variability_Between = 0.2, 
                         NumbersOfSubjects)
write.csv(Powers17,"Powers17.csv") #broader response function

Powers18 = ComparePowers(ConditionOfInterest, StandardValues, reps = 1:25, PSE_Difference = -0.05, JND_Difference = 0.2, 
                         Multiplicator_PSE_Standard, Multiplicator_SD_Standard, SD_ResponseFunction = 0.1, 
                         Mean_Variability_Between = 0.1, SD_Variability_Between = 0.1, 
                         NumbersOfSubjects)
write.csv(Powers18,"Powers18.csv") #broader response function, fewer trials


Powers19 = ComparePowers(ConditionOfInterest, StandardValues, reps = 1:55, PSE_Difference = 0.1, JND_Difference = 0.05, 
                         Multiplicator_PSE_Standard, Multiplicator_SD_Standard, SD_ResponseFunction = 0.1, 
                         Mean_Variability_Between, SD_Variability_Between, 
                         NumbersOfSubjects)
write.csv(Powers19,"Powers19.csv") #broader response function

Powers20 = ComparePowers(ConditionOfInterest, StandardValues, reps = 1:25, PSE_Difference = 0.1, JND_Difference = 0.05, 
                         Multiplicator_PSE_Standard, Multiplicator_SD_Standard, SD_ResponseFunction = 0.1, 
                         Mean_Variability_Between, SD_Variability_Between, 
                         NumbersOfSubjects)
write.csv(Powers20,"Powers20.csv") #broader response function, fewer trials

Powers21 = ComparePowers(ConditionOfInterest, StandardValues, reps = 1:55, PSE_Difference = 0.1, JND_Difference = 0.2, 
                         Multiplicator_PSE_Standard, Multiplicator_SD_Standard, SD_ResponseFunction = 0.1, 
                         Mean_Variability_Between = 0.2, SD_Variability_Between = 0.2, 
                         NumbersOfSubjects)
write.csv(Powers21,"Powers21.csv") #broader response function

Powers22 = ComparePowers(ConditionOfInterest, StandardValues, reps = 1:25, PSE_Difference = 0.1, JND_Difference = 0.2, 
                         Multiplicator_PSE_Standard, Multiplicator_SD_Standard, SD_ResponseFunction = 0.1, 
                         Mean_Variability_Between = 0.1, SD_Variability_Between = 0.1, 
                         NumbersOfSubjects)
write.csv(Powers22,"Powers22.csv") #broader response function, fewer trials

Powers23 = ComparePowers(ConditionOfInterest, StandardValues, reps = 1:25, PSE_Difference = 0.05, JND_Difference = 0.05, 
                         Multiplicator_PSE_Standard, Multiplicator_SD_Standard, SD_ResponseFunction = 0.1, 
                         Mean_Variability_Between=0.1, SD_Variability_Between = 0.1, 
                         NumbersOfSubjects)
write.csv(Powers23,"Powers23.csv") #broader response function

Powers24 = ComparePowers(ConditionOfInterest, StandardValues, reps = 1:25, PSE_Difference = 0.05, JND_Difference = 0.05, 
                         Multiplicator_PSE_Standard, Multiplicator_SD_Standard, SD_ResponseFunction = 0.1, 
                         Mean_Variability_Between=0.1, SD_Variability_Between=0.1, 
                         NumbersOfSubjects)
write.csv(Powers24,"Powers24.csv") #broader response function, fewer trials

Powers25 = ComparePowers(ConditionOfInterest, StandardValues, reps = 1:25, PSE_Difference = -0.1, JND_Difference = 0.1, 
                         Multiplicator_PSE_Standard, Multiplicator_SD_Standard, SD_ResponseFunction = 0.1, 
                         Mean_Variability_Between = 0.1, SD_Variability_Between = 0.1, 
                         NumbersOfSubjects)
write.csv(Powers25,"Powers25.csv") #broader response function

Powers26 = ComparePowers(ConditionOfInterest, StandardValues, reps = 1:25, PSE_Difference = -0.1, JND_Difference = 0.1, 
                         Multiplicator_PSE_Standard, Multiplicator_SD_Standard, SD_ResponseFunction = 0.1, 
                         Mean_Variability_Between = 0.1, SD_Variability_Between = 0.1, 
                         NumbersOfSubjects)
write.csv(Powers26,"Powers26.csv") #broader response function, fewer trials

PowerFrame <- rbind(read.csv(header = T, file = paste0(Where_Am_I(),"/Data/Powers1.csv")),
                    read.csv(header = T, file = paste0(Where_Am_I(),"/Data/Powers2.csv")),
                    read.csv(header = T, file = paste0(Where_Am_I(),"/Data/Powers3.csv")),
                    read.csv(header = T, file = paste0(Where_Am_I(),"/Data/Powers4.csv")),
                    read.csv(header = T, file = paste0(Where_Am_I(),"/Data/Powers5.csv")),
                    read.csv(header = T, file = paste0(Where_Am_I(),"/Data/Powers6.csv")),
                    read.csv(header = T, file = paste0(Where_Am_I(),"/Data/Powers7.csv")),
                    read.csv(header = T, file = paste0(Where_Am_I(),"/Data/Powers8.csv")),
                    read.csv(header = T, file = paste0(Where_Am_I(),"/Data/Powers9.csv")),
                    read.csv(header = T, file = paste0(Where_Am_I(),"/Data/Powers10.csv")),
                    read.csv(header = T, file = paste0(Where_Am_I(),"/Data/Powers11.csv")),
                    read.csv(header = T, file = paste0(Where_Am_I(),"/Data/Powers12.csv")),
                    read.csv(header = T, file = paste0(Where_Am_I(),"/Data/Powers13.csv")),
                    read.csv(header = T, file = paste0(Where_Am_I(),"/Data/Powers14.csv")),
                    read.csv(header = T, file = paste0(Where_Am_I(),"/Data/Powers15.csv")),
                    read.csv(header = T, file = paste0(Where_Am_I(),"/Data/Powers16.csv")),
                    read.csv(header = T, file = paste0(Where_Am_I(),"/Data/Powers17.csv")),
                    read.csv(header = T, file = paste0(Where_Am_I(),"/Data/Powers19.csv")),
                    read.csv(header = T, file = paste0(Where_Am_I(),"/Data/Powers20.csv")),
                    read.csv(header = T, file = paste0(Where_Am_I(),"/Data/Powers21.csv")),
                    read.csv(header = T, file = paste0(Where_Am_I(),"/Data/Powers22.csv")),
                    read.csv(header = T, file = paste0(Where_Am_I(),"/Data/Powers23.csv")),
                    read.csv(header = T, file = paste0(Where_Am_I(),"/Data/Powers24.csv")),
                    read.csv(header = T, file = paste0(Where_Am_I(),"/Data/Powers25.csv")),
                    read.csv(header = T, file = paste0(Where_Am_I(),"/Data/Powers26.csv")))


ggplot(PowerFrame %>% filter(label != c("Accuracy Two-Level LMM","Precision Two-Level LMM", 
                                        "Accuracy GLMM Only Intercepts", "Precision GLMM Only Intercepts") &
                               reps == 25),
       aes(nSubjects,value, color = label)) +
  geom_line(size = 2) +
  facet_grid(JND_Difference~PSE_Difference) +
  xlab("N° of Subjects") +
  ylab("Power") +
  geom_hline(linetype = 2, yintercept = 0.8) +
  geom_hline(linetype = 3, yintercept = 0.95) +
  scale_x_continuous(breaks=c(10,12,14,16,18,20))






########################################################################
##############compare power for GLMM and Two-Level approach#############
########################################################################
Dataframe_wide_Big = rbind(read.csv(header = T, file = paste0(Where_Am_I(),"/Data/Powers_BiggerModel.csv"))) %>%
                  select(power_Accuracy,power_Precision,power_Accuracy_Twolevel,power_Precision_Twolevel,power_Accuracy_ANOVA,
                         power_Precision_ANOVA,n, PSE_Difference, JND_Difference, reps)
Dataframe_Powers_Big = data.frame(PSE_Difference = rep(Dataframe_wide_Big$PSE_Difference,6),
                              JND_Difference = rep(Dataframe_wide_Big$JND_Difference,6),
                              n = rep(Dataframe_wide_Big$n,6),
                              reps = rep(Dataframe_wide_Big$reps,6),
                              power = c(Dataframe_wide_Big$power_Accuracy,
                                        Dataframe_wide_Big$power_Precision,
                                        Dataframe_wide_Big$power_Accuracy_Twolevel,
                                        Dataframe_wide_Big$power_Precision_Twolevel,
                                        Dataframe_wide_Big$power_Accuracy_ANOVA,
                                        Dataframe_wide_Big$power_Precision_ANOVA),
                              label = c(rep("Accuracy_GLMM",length(Dataframe_wide_Big$reps)),
                                        rep("Precision_GLMM",length(Dataframe_wide_Big$reps)),
                                        rep("Accuracy_Twolevel",length(Dataframe_wide_Big$reps)),
                                        rep("Precision_Twolevel",length(Dataframe_wide_Big$reps)),
                                        rep("Accuracy_Twolevel_ANOVA",length(Dataframe_wide_Big$reps)),
                                        rep("Precision_Twolevel_ANOVA",length(Dataframe_wide_Big$reps))))

Dataframe_wide_Small = rbind(read.csv(header = T, file = paste0(Where_Am_I(),"/Data/Powers_SmallerModel.csv"))) %>%
                  select(power_Accuracy,power_Precision,power_Accuracy_Twolevel,power_Precision_Twolevel,power_Accuracy_ANOVA,
                         power_Precision_ANOVA,n,PSE_Difference, JND_Difference, reps)
Dataframe_Powers_Small = data.frame(PSE_Difference = rep(Dataframe_wide_Small$PSE_Difference,6),
                                  JND_Difference = rep(Dataframe_wide_Small$JND_Difference,6),
                                  n = rep(Dataframe_wide_Small$n,6),
                                  reps = rep(Dataframe_wide_Small$reps,6),
                                  power = c(Dataframe_wide_Small$power_Accuracy,
                                            Dataframe_wide_Small$power_Precision,
                                            Dataframe_wide_Small$power_Accuracy_Twolevel,
                                            Dataframe_wide_Small$power_Precision_Twolevel,
                                            Dataframe_wide_Small$power_Accuracy_ANOVA,
                                            Dataframe_wide_Small$power_Precision_ANOVA),
                                  label = c(rep("Accuracy_GLMM",length(Dataframe_wide_Small$reps)),
                                            rep("Precision_GLMM",length(Dataframe_wide_Small$reps)),
                                            rep("Accuracy_Twolevel",length(Dataframe_wide_Small$reps)),
                                            rep("Precision_Twolevel",length(Dataframe_wide_Small$reps)),
                                            rep("Accuracy_Twolevel_ANOVA",length(Dataframe_wide_Small$reps)),
                                            rep("Precision_Twolevel_ANOVA",length(Dataframe_wide_Small$reps))))


Powers1 = ggplot(Dataframe_Powers_Big %>% filter(reps == 40), aes(n,power,color = label)) +
  geom_line(size = 2) +
  facet_grid(JND_Difference~PSE_Difference) +
  xlab("N° of Participants") +
  ylab("Power") +
  geom_hline(linetype = 2, yintercept = 0.8) +
  geom_hline(linetype = 1, yintercept = 0.05) +
  geom_hline(linetype = 3, yintercept = 0.95) +
  scale_x_continuous(breaks=c(10,15,20)) +
  scale_y_continuous(breaks=c(0.25,0.75)) +
  scale_color_manual(values = c(BlauUB,LightBlauUB,Red,LightRed,Yellow,LightYellow), 
                     name = "") +
  ggtitle("A. 40 Repetitions")

Powers2 = ggplot(Dataframe_Powers_Big %>% filter(reps == 70), aes(n,power,color = label)) +
  geom_line(size = 2) +
  facet_grid(JND_Difference~PSE_Difference) +
  xlab("N° of Participants") +
  ylab("Power") +
  geom_hline(linetype = 2, yintercept = 0.8) +
  geom_hline(linetype = 1, yintercept = 0.05) +
  geom_hline(linetype = 3, yintercept = 0.95) +
  scale_x_continuous(breaks=c(10,15,20)) +
  scale_y_continuous(breaks=c(0.25,0.75)) +
  scale_color_manual(values = c(BlauUB,LightBlauUB,Red,LightRed,Yellow,LightYellow), 
                     name = "") +
  ggtitle("A. 70 Repetitions")

Powers3 = ggplot(Dataframe_Powers_Big %>% filter(reps == 100), aes(n,power,color = label)) +
  geom_line(size = 2) +
  facet_grid(JND_Difference~PSE_Difference) +
  xlab("N° of Participants") +
  ylab("Power") +
  geom_hline(linetype = 2, yintercept = 0.8) +
  geom_hline(linetype = 1, yintercept = 0.05) +
  geom_hline(linetype = 3, yintercept = 0.95) +
  scale_x_continuous(breaks=c(10,15,20)) +
  scale_y_continuous(breaks=c(0.25,0.75)) +
  scale_color_manual(values = c(BlauUB,LightBlauUB,Red,LightRed,Yellow,LightYellow), 
                     name = "") +
  ggtitle("A. 100 Repetitions")

plot_shared_legend(Powers1,Powers2, Powers3)
ggsave("Figures/PowersBigModel.jpg", w = 12, h = 8)


Powers4 = ggplot(Dataframe_Powers_Small %>% filter(reps == 40), aes(n,power,color = label)) +
  geom_line(size = 2) +
  facet_grid(JND_Difference~PSE_Difference) +
  xlab("N° of Participants") +
  ylab("Power") +
  geom_hline(linetype = 2, yintercept = 0.8) +
  geom_hline(linetype = 1, yintercept = 0.05) +
  geom_hline(linetype = 3, yintercept = 0.95) +
  scale_x_continuous(breaks=c(10,15,20)) +
  scale_y_continuous(breaks=c(0.25,0.75)) +
  scale_color_manual(values = c(BlauUB,LightBlauUB,Red,LightRed,Yellow,LightYellow), 
                     name = "") +
  ggtitle("A. 40 Repetitions")

Powers5 = ggplot(Dataframe_Powers_Small %>% filter(reps == 70), aes(n,power,color = label)) +
  geom_line(size = 2) +
  facet_grid(JND_Difference~PSE_Difference) +
  xlab("N° of Participants") +
  ylab("Power") +
  geom_hline(linetype = 2, yintercept = 0.8) +
  geom_hline(linetype = 1, yintercept = 0.05) +
  geom_hline(linetype = 3, yintercept = 0.95) +
  scale_x_continuous(breaks=c(10,15,20)) +
  scale_y_continuous(breaks=c(0.25,0.75)) +
  scale_color_manual(values = c(BlauUB,LightBlauUB,Red,LightRed,Yellow,LightYellow), 
                     name = "") +
  ggtitle("A. 70 Repetitions")
?brm
Powers6 = ggplot(Dataframe_Powers_Small %>% filter(reps == 100), aes(n,power,color = label)) +
  geom_line(size = 2) +
  facet_grid(JND_Difference~PSE_Difference) +
  xlab("N° of Participants") +
  ylab("Power") +
  geom_hline(linetype = 2, yintercept = 0.8) +
  geom_hline(linetype = 1, yintercept = 0.05) +
  geom_hline(linetype = 3, yintercept = 0.95) +
  scale_x_continuous(breaks=c(10,15,20)) +
  scale_y_continuous(breaks=c(0.25,0.75)) +
  scale_color_manual(values = c(BlauUB,LightBlauUB,Red,LightRed,Yellow,LightYellow), 
                     name = "") +
  ggtitle("A. 100 Repetitions")

plot_shared_legend(Powers4,Powers5, Powers6)
ggsave("Figures/Powers_SmallModel.jpg", w = 12, h = 8)
########################################################################
########################################################################
########################################################################


########################################################################
####################Compare Optimizer Configurations####################
########################################################################
Dataframe_pvalues1 = rbind(read.csv(header = T, file = paste0(Where_Am_I(),"/Data/ComparisonMethodsSmallEffectBigModel.csv")))
Dataframe_pvalues1$Size = "Bigger Model"
#Dataframe_pvalues2 = rbind(read.csv(header = T, file = paste0(Where_Am_I(),"/Data/ComparisonMethods_SmallModel.csv")))
#Dataframe_pvalues2$Size = "Smaller Model"
#colnames(Dataframe_pvalues2)[1] <- "iteration"



Dataframe_pvalues = rbind(Dataframe_pvalues1,Dataframe_pvalues2)

Dataframe_pvalues = Dataframe_pvalues %>%
  mutate(Optimizer = case_when(
    label == "JuliaAIC_NeldMeader_AGP0" ~ "Julia: Nelder-Mead, fast",
    label == "JuliaAIC_bobyqa_AGP0" ~ "Julia: BOBYQA, fast",
    label == "JuliaAIC_NeldMeader_AGP1" ~ "Julia: Nelder-Mead, slow",
    label == "JuliaAIC_bobyqa_AGP1" ~ "Julia: BOBYQA, slow",
    label == "NelderMead_nAGQ0" ~ "R: Nelder-Mead, fast",
    label == "NelderMead_nAGQ1" ~ "R: Nelder-Mead, slow",
    label == "Bobyqa_nAGQ0" ~ "R: BOBYQA, fast",
    label == "Bobyqa_nAGQ1" ~ "R: BOBYQA, slow",
    label == "nloptwrap_nAGQ0" ~ "R: nloptwrap, fast",
    label == "nloptwrap_nAGQ1" ~ "R: nloptwrap, slow",
    label == "JuliaLRT" ~ "Julia: LRT",
    label == "R: LRT" ~ "R: LRT")
    )%>%
  group_by(n,reps,label,PSE_Difference,JND_Difference,Size) %>%
  mutate(nTrials = case_when(
           reps == 30 ~ "30 repetitions",
           reps == 40 ~ "40 repetitions",
           reps == 50 ~ "50 repetitions",
           reps == 60 ~ "60 repetitions"),
         Effect = case_when(
           PSE_Difference == 0.025 ~ "Small Effect",
           PSE_Difference == 0 ~ "No Effect")
  )
         
Dataframe_pvalues = Dataframe_pvalues %>%         
  group_by(reps,n,iteration,Effect,Size) %>% 
  mutate(AIC_Ratio = AIC/AIC[2],
         Duration_LRT_Julia = Duration[9] + Duration[12],
         Duration_LRT_R = Duration[5] + Duration[7]) %>% 
  group_by(reps,n,Optimizer,Size) %>% 
  mutate(Median_AIC_Ratio = median(AIC_Ratio))

Dataframe_pvalues$Duration[Dataframe_pvalues$label == "JuliaLRT"] = 
        Dataframe_pvalues$Duration_LRT_Julia[Dataframe_pvalues$label == "JuliaLRT"]
Dataframe_pvalues$Duration[Dataframe_pvalues$label == "R: LRT"] = 
  Dataframe_pvalues$Duration_LRT_R[Dataframe_pvalues$label == "R: LRT"]

Dataframe_pvalues = Dataframe_pvalues %>%         
  group_by(reps,n,label,Effect,Size) %>% 
    mutate(MedianDuration = median(Duration),
           SEDuration = SE(Duration),
           SE_Duration_n_reps_label = SE(Duration))

#######Timing
TimingPlot1_Bigger = ggplot(Dataframe_pvalues %>% 
         filter(Size == "Bigger Model") ,
       aes(n,MedianDuration, color = Optimizer)) +
  geom_point(size=2) +
  geom_line(size=1) +
  facet_grid(Effect~nTrials) +
  ylab("Median Duration (s)") +
  scale_color_manual(values = rainbow(12), name = "Method") +
  scale_x_continuous(breaks = c(10,15,20)) +
  ggtitle("A. All Configurations")


TimingPlot2_Bigger = ggplot(Dataframe_pvalues %>% 
         filter(Optimizer  %in% c("Julia: BOBYQA, fast",
                                  "R: nloptwrap, fast",
                                  "Julia: LRT",
                                  "R: LRT",
                                  "R: nloptwrap, slow",
                                  "Julia: BOBYQA, slow") &
                Size == "Bigger Model"),
       aes(n,MedianDuration, color = Optimizer)) +
  geom_point(size=2) +
  geom_line(size=1) +
  facet_grid(Effect~nTrials) +
  ylab("Median Duration (s)") +
  scale_color_manual(values = rainbow(6), name = "Method") +
  scale_x_continuous(breaks = c(10,15,20)) +
  ggtitle("B. Fastest Configurations")

plot_shared_legend(TimingPlot1_Bigger,TimingPlot2_Bigger)
ggsave("Figures/Different Durations Bigger.jpg",w=12,h=6)

TimingPlot1_Smaller = ggplot(Dataframe_pvalues %>% 
                              filter(Size == "Smaller Model") ,
                            aes(n,MedianDuration, color = Optimizer)) +
  geom_point(size=2) +
  geom_line(size=1) +
  facet_grid(Effect~nTrials) +
  ylab("Median Duration (s)") +
  scale_color_manual(values = rainbow(12), name = "Method") +
  scale_x_continuous(breaks = c(10,15,20)) +
  ggtitle("A. All Configurations")


TimingPlot2_Smaller = ggplot(Dataframe_pvalues %>% 
                              filter(Optimizer  %in% c("Julia: BOBYQA, fast",
                                                       "R: nloptwrap, fast",
                                                       "Julia: LRT",
                                                       "R: LRT",
                                                       "R: nloptwrap, slow",
                                                       "Julia: BOBYQA, slow") &
                                       Size == "Smaller Model"),
                            aes(n,MedianDuration, color = Optimizer)) +
  geom_point(size=2) +
  geom_line(size=1) +
  facet_grid(Effect~nTrials) +
  ylab("Median Duration (s)") +
  scale_color_manual(values = rainbow(6), name = "Method") +
  scale_x_continuous(breaks = c(10,15,20)) +
  ggtitle("B. Fastest Configurations")

plot_shared_legend(TimingPlot1_Smaller,TimingPlot2_Smaller)
ggsave("Figures/Different Durations Smaller",w=12,h=6)


#######False Positives
Dataframe_pvalues$Bin_Accuracy = 0
Dataframe_pvalues$Bin_Interaction = 0
for (i in (1:length(Dataframe_pvalues$iteration))){
  print(i)  
  Bins = seq(0.025,0.975,0.025)
  Dataframe_pvalues$Bin_Accuracy[i] = Bins[which.min(abs(Bins-Dataframe_pvalues$Pvalues_Accuracy[i]))]+0.0125
  Dataframe_pvalues$Bin_Interaction[i] = Bins[which.min(abs(Bins-Dataframe_pvalues$Pvalues_Interaction[i]))]+0.0125

}

Dataframe_pvalues = Dataframe_pvalues %>%
  group_by(Bin_Accuracy,Optimizer,PSE_Difference,Size) %>%
  mutate(BinCountAccuracy = length(Bin_Accuracy))%>%
  group_by(Bin_Interaction,Optimizer,PSE_Difference,Size) %>%
  mutate(BinCountInteraction = length(Bin_Interaction))

PlotAccuracy_Bigger = ggplot(Dataframe_pvalues %>% 
         filter(Optimizer != "Julia: LRT" & Optimizer != "R: LRT" &
                  Effect == "No Effect" &
                  Size == "Bigger Model"),
       aes(Bin_Accuracy-0.0125,Optimizer, fill = as.factor(BinCountAccuracy))) +
  geom_tile() +
  xlab("") +
  scale_fill_manual(values=colorRampPalette(c(BlauUB,Red,Yellow))(27)) +
  theme(legend.position = "") +
  ylab("") +
  ggtitle("A. Accuracy, No Effect")

PlotInteraction_Bigger = ggplot(Dataframe_pvalues %>% filter(Effect == "No Effect" &
                                                        Size == "Bigger Model"),
                aes(Bin_Interaction-0.025,Optimizer, fill = as.factor(BinCountInteraction))) +
  geom_tile() +
  xlab("") +
  scale_fill_manual(values=colorRampPalette(c(BlauUB,Red, Yellow))(32)) +
  theme(legend.position = "") +
  ylab("") +
  ggtitle("B. Interaction, No Effect")


PlotAccuracy2_Bigger = ggplot(Dataframe_pvalues %>% 
                        filter(Effect == "Small Effect" &
                                 Optimizer != "Julia: LRT" &
                                 Optimizer != "R: LRT" &
                                 Size == "Bigger Model"),
                      aes(Bin_Accuracy-0.025,Optimizer, fill = as.factor(BinCountAccuracy))) +
  geom_tile() +
  xlab("") +
  ylab("") +
  scale_fill_manual(values=colorRampPalette(c(BlauUB,Red,Yellow))(40)) +
  theme(legend.position = "") +
  ggtitle("C. Accuracy, Small Effect")

PlotInteraction2_Bigger = ggplot(Dataframe_pvalues %>% 
                           filter(Effect == "Small Effect" &
                                    Size == "Bigger Model"),
                         aes(Bin_Interaction-0.025,Optimizer, fill = as.factor(BinCountInteraction))) +
  geom_tile() +
  xlab("") +
  ylab("") +
  scale_fill_manual(values=colorRampPalette(c(BlauUB,Red,Yellow))(41)) +
  theme(legend.position = "") +
  ggtitle("D. Interaction, Small Effect")
plot_grid(PlotAccuracy_Bigger,PlotInteraction_Bigger,PlotAccuracy2_Bigger,PlotInteraction2_Bigger, nrow = 2)
ggsave("Figures/False Positives Bigger Model.jpg",w=12,h=8)



PlotAccuracy_Small = ggplot(Dataframe_pvalues %>% 
                        filter(Optimizer != "Julia: LRT" & Optimizer != "R: LRT" &
                                 Effect == "No Effect" &
                                 Size == "Smaller Model"),
                      aes(Bin_Accuracy-0.025,Optimizer, fill = as.factor(BinCountAccuracy))) +
  geom_tile() +
  xlab("") +
  scale_fill_manual(values=colorRampPalette(c(BlauUB,Red,Yellow))(25)) +
  theme(legend.position = "") +
  ylab("") +
  ggtitle("A. Accuracy, No Effect")

PlotInteraction_Small = ggplot(Dataframe_pvalues %>% filter(Effect == "No Effect" &
                                                        Size == "Smaller Model"),
                         aes(Bin_Interaction-0.025,Optimizer, fill = as.factor(BinCountInteraction))) +
  geom_tile() +
  xlab("") +
  scale_fill_manual(values=colorRampPalette(c(BlauUB,Red, Yellow))(24)) +
  theme(legend.position = "") +
  ylab("") +
  ggtitle("B. Interaction, No Effect")


PlotAccuracy2_Small = ggplot(Dataframe_pvalues %>% 
                         filter(Effect == "Small Effect" &
                                  Optimizer != "Julia: LRT" &
                                  Optimizer != "R: LRT" &
                                  Size == "Smaller Model"),
                       aes(Bin_Accuracy-0.025,Optimizer, fill = as.factor(BinCountAccuracy))) +
  geom_tile() +
  xlab("") +
  ylab("") +
  scale_fill_manual(values=colorRampPalette(c(BlauUB,Red,Yellow))(16)) +
  theme(legend.position = "") +
  ggtitle("C. Accuracy, Small Effect")

PlotInteraction2_Small = ggplot(Dataframe_pvalues %>% 
                            filter(Effect == "Small Effect" &
                                     Size == "Smaller Model"),
                          aes(Bin_Interaction-0.025,Optimizer, fill = as.factor(BinCountInteraction))) +
  geom_tile() +
  xlab("") +
  ylab("") +
  scale_fill_manual(values=colorRampPalette(c(BlauUB,Red,Yellow))(30)) +
  theme(legend.position = "") +
  ggtitle("D. Interaction, Small Effect")
plot_grid(PlotAccuracy_Small,PlotInteraction_Small,PlotAccuracy2_Small,PlotInteraction2_Small, nrow = 2)
ggsave("Figures/False Positives Smaller Model.jpg",w=12,h=8)

############AICs
ggplot(Dataframe_pvalues %>%
         filter(Optimizer != "Julia: LRT" & 
                  Optimizer != "R: LRT" &
                  Size == "Bigger Model"),
       aes(n,Median_AIC_Ratio, color = Optimizer)) +
  geom_point(size=2) +
  geom_line(size=1) +
  facet_grid(Effect~nTrials) +
  ylab("AIC DIfference") +
  scale_color_manual(values = rainbow(10), name = "Method") +
  scale_x_continuous(breaks = c(10,15,20))
ggsave("Figures/AIC differences.jpg",w=12,h=6)

