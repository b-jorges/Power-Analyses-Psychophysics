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

pnorm(10.73,10,10*0.108)

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
                       family = binomial(link = "probit"), 
                       data = Psychometric,
                       nAGQ = 0,
                       control = glmerControl(optimizer = "nloptwrap"))

GLMM2 = glmer(cbind(Yes, Total - Yes) ~ ConditionOfInterest*Difference + (Difference| ID) + (Difference| StandardValues), 
             family = binomial(link = "probit"), 
             data = Psychometric,
             nAGQ = 0,
             control = glmerControl(optimizer = "nloptwrap"))


summary(GLMM2)
summary(GLMM)

hallo = anova(GLMM,GLMM2)
hallo


Psychometric %>% 
  group_by(ConditionOfInterest, StandardValues) %>% 
  slice(1)




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



Power = data.frame()
nIterations = 200
pvalue = 0.05
NumbersOfSubjects = c(10,12,14,16,18,20)
#NumbersOfSubjects = c(2,3,4,5,6)

Mean_Variability_Between = 0.1
SD_Variability_Between = 0.1



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
Dataframe_wide = rbind(read.csv(header = T, file = paste0(Where_Am_I(),"/Data/SamplePowers1_40_reps.csv")),
                  read.csv(header = T, file = paste0(Where_Am_I(),"/Data/SamplePowers2_60_reps.csv"))) %>%
                  select(power_Accuracy,power_Precision,power_Accuracy_Twolevel,power_Precision_Twolevel,n, PSE_Difference, JND_Difference, reps)
Dataframe_Powers = data.frame(PSE_Difference = rep(Dataframe_wide$PSE_Difference,4),
                              JND_Difference = rep(Dataframe_wide$JND_Difference,4),
                              n = rep(Dataframe_wide$n,4),
                              reps = rep(Dataframe_wide$reps,4),
                              power = c(Dataframe_wide$power_Accuracy,
                                        Dataframe_wide$power_Precision,
                                        Dataframe_wide$power_Accuracy_Twolevel,
                                        Dataframe_wide$power_Precision_Twolevel),
                              label = c(rep("Accuracy_GLMM",length(Dataframe_wide$reps)),
                                        rep("Precision_GLMM",length(Dataframe_wide$reps)),
                                        rep("Accuracy_Twolevel",length(Dataframe_wide$reps)),
                                        rep("Precision_Twolevel",length(Dataframe_wide$reps))))


ggplot(Dataframe_Powers %>% filter(reps == 60), aes(n,power,color = label)) +
  geom_line(size = 2) +
  facet_grid(JND_Difference~PSE_Difference) +
  xlab("N° of Subjects") +
  ylab("Power") +
  geom_hline(linetype = 2, yintercept = 0.8) +
  geom_hline(linetype = 1, yintercept = 0.05) +
  geom_hline(linetype = 3, yintercept = 0.95) +
  scale_x_continuous(breaks=c(10,12,14,16,18,20)) +
  scale_color_manual(values = c(BlauUB,LightBlauUB,Red,LightRed), 
                     name = "")
########################################################################
########################################################################
########################################################################



########################################################################
##############compare Speed for R and Julia#############
########################################################################
Dataframe = read.csv(header = T, file = paste0(Where_Am_I(),"/Data/Durations_Julia.csv"))
Dataframe_Julia = gather(Dataframe,analysis,duration,
                         c(DurationGLMM_NeldMeader_AGP0,DurationGLMM_bobyqa_AGP0),factor_key = TRUE)
Julia_AIC = gather(Dataframe,analysis2,AIC,
                         c(AIC_NeldMeader_AGP0,AIC_bobyqa_AGP0),factor_key = TRUE)$AIC
Dataframe_Julia$AIC = Julia_AIC
Dataframe_Julia$duration = Dataframe_Julia$duration/1000
Dataframe_Julia = rbind(Dataframe_Julia %>% select(nIteration,n,reps,analysis,duration,AIC))
colnames(Dataframe_Julia) = c("iteration","n","reps","analysis","duration", "AIC")
Dataframe_Julia$Program = "Julia"

Dataframe_R = read.csv(header = T, file = paste0(Where_Am_I(),"/Data/Durations_R.csv"))
Dataframe_R1 = gather(Dataframe_R,analysis,duration,
                     c(Duration_NelderMead_nAGQ0,Duration_NelderMead_nAGQ1,Duration_Bobyqa_nAGQ0,
                       Duration_Bobyqa_nAGQ1,Duration_nloptwrap_nAGQ0,Duration_nloptwrap_nAGQ1),factor_key = TRUE)
AIC_R2 = gather(Dataframe_R,analysis,AIC,
                     c(AIC_NelderMead_nAGQ0,AIC_NelderMead_nAGQ1,AIC_Bobyqa_nAGQ0,
                       AIC_Bobyqa_nAGQ1,AIC_nloptwrap_nAGQ0,AIC_nloptwrap_nAGQ1),factor_key = TRUE)$AIC
Dataframe_R1$AIC = AIC_R2
Dataframe_R = Dataframe_R1 %>% select(iteration, n, reps, analysis, duration,AIC)
colnames(Dataframe_R) = c("iteration","n","reps","analysis", "duration","AIC")
Dataframe_R$Program = "R"

Dataframe = rbind(Dataframe_Julia,Dataframe_R) %>%
  group_by(n,reps,Program,analysis) %>%
  mutate(MeanDuration = median(duration),
         MeanAIC = median(AIC))

Dataframe = Dataframe %>%
  group_by(n,reps,Program) %>%
  mutate(MedianAIC_n_reps = median(AIC)) %>%
  ungroup() %>%
  mutate(AICDifferences = MeanAIC-MedianAIC_n_reps,
         Optimizer = case_when(
           analysis == "DurationGLMM_NeldMeader_AGP0" ~ "Julia: Nelder-Mead, nAGQ 0",
           analysis == "DurationGLMM_bobyqa_AGP0" ~ "Julia: bobyqa, nAGQ 0",
           analysis == "Duration_NelderMead_nAGQ0" ~ "R: Nelder-Mead, nAGQ 0",
           analysis == "Duration_NelderMead_nAGQ1" ~ "R: Nelder-Mead, nAGQ 1",
           analysis == "Duration_Bobyqa_nAGQ0" ~ "R: bobyqa, nAGQ 0",
           analysis == "Duration_Bobyqa_nAGQ1" ~ "R: bobyqa, nAGQ 1",
           analysis == "Duration_nloptwrap_nAGQ0" ~ "R: nloptwrap, nAGQ 0",
           analysis == "Duration_nloptwrap_nAGQ1" ~ "R: nloptwrap, nAGQ 1"
          )
         )

ggplot(Dataframe,aes(n, MeanDuration, color = analysis)) +
  geom_point() +
  geom_line() +
  facet_grid(.~reps) +
  coord_cartesian(ylim = c(0,10)) +
  scale_color_manual(values = colorRampPalette(c(BlauUB,Yellow,Red,"green"))(8))

ggplot(Dataframe,aes(n, AICDifferences, color = Optimizer)) +
  geom_point() +
  geom_line() +
  facet_grid(.~reps) +
  scale_color_manual(values = colorRampPalette(c(BlauUB,Yellow,Red,"green"))(8)) +
  coord_cartesian(ylim = c(-0.5,0.5))


########################################################################
##############compare AICs for different algorithms#####################
########################################################################
Dataframe_AICs = read.csv(header = T, file = paste0(Where_Am_I(),"/Data/AICs.csv"))

Dataframe_AICs = Dataframe_AICs %>%
  mutate(Optimizer = case_when(
    label == "JuliaAIC_NeldMeader_AGP0" ~ "Julia: Nelder-Mead, nAGQ 0",
    label == "JuliaAIC_bobyqa_AGP0" ~ "Julia: bobyqa, nAGQ 0",
    label == "NelderMead_nAGQ0" ~ "R: Nelder-Mead, nAGQ 0",
    label == "NelderMead_nAGQ1" ~ "R: Nelder-Mead, nAGQ 1",
    label == "Bobyqa_nAGQ0" ~ "R: bobyqa, nAGQ 0",
    label == "Bobyqa_nAGQ1" ~ "R: bobyqa, nAGQ 1",
    label == "nloptwrap_nAGQ0" ~ "R: nloptwrap, nAGQ 0",
    label == "nloptwrap_nAGQ1" ~ "R: nloptwrap, nAGQ 1")
    )%>%
  group_by(n,reps) %>%
  mutate(MedianAIC_n_reps = median(AIC),
         Median_pvalue_Accuracy_n_reps = median(Pvalues_Accuracy),
         Median_pvalue_Interaction_n_reps = median(Pvalues_Interaction)) %>%
  group_by(n,reps,label) %>%
  mutate(MedianDuration = median(Duration),
         SE_Duration_n_reps_label = SE(Duration),
         MedianAIC_Difference = median(AIC)-MedianAIC_n_reps,
         Median_Pvalue_Accuracy_Difference = median(Pvalues_Accuracy) - Median_pvalue_Accuracy_n_reps,
         Median_Pvalue_Interaction_Difference = median(Pvalues_Interaction) - Median_pvalue_Interaction_n_reps)


ggplot(Dataframe_AICs,aes(n, MedianDuration, color = Optimizer)) +
  geom_point() +
  geom_line() +
  geom_errorbar(aes(ymin = MedianDuration-SE_Duration_n_reps_label, 
                    ymax = MedianDuration+SE_Duration_n_reps_label), width = 0.7) +
  facet_grid(.~reps) +
  ylab("Median Duration (s)") +
  scale_x_continuous(breaks =  c(10,12,14,16,18,20), name = "N° Participants") +
  scale_color_manual(values = colorRampPalette(c(BlauUB,Yellow))(8), name = "") +
  theme(legend.position = c(0.05,0.8))
ggsave("Figures/Duration for each Optimizer.jpg", w = 10, h = 5)

ggplot(Dataframe_AICs,aes(n, MedianAIC_Difference, color = Optimizer)) +
  geom_point(size=2) +
  geom_line() +
  facet_grid(.~reps) +
  ylab("Median AIC difference from median across optimizers") +
  scale_x_continuous(breaks =  c(10,12,14,16,18,20), name = "N° Participants") +
  scale_color_manual(values = colorRampPalette(c(BlauUB,Yellow))(8))
ggsave("Figures/AICs for each Optimizer.jpg", w = 10, h = 5)


########################################################################
###########################compare p values#############################
########################################################################
Dataframe_pvalues = read.csv(header = T, file = paste0(Where_Am_I(),"/Data/Pvalues_Julia.csv"))

Dataframe_pvalues = Dataframe_AICs %>%
  mutate(Optimizer = case_when(
    label == "JuliaAIC_NeldMeader_AGP0" ~ "Julia: Nelder-Mead, nAGQ 0",
    label == "JuliaAIC_bobyqa_AGP0" ~ "Julia: bobyqa, nAGQ 0",
    label == "NelderMead_nAGQ0" ~ "R: Nelder-Mead, nAGQ 0",
    label == "NelderMead_nAGQ1" ~ "R: Nelder-Mead, nAGQ 1",
    label == "Bobyqa_nAGQ0" ~ "R: bobyqa, nAGQ 0",
    label == "Bobyqa_nAGQ1" ~ "R: bobyqa, nAGQ 1",
    label == "nloptwrap_nAGQ0" ~ "R: nloptwrap, nAGQ 0",
    label == "nloptwrap_nAGQ1" ~ "R: nloptwrap, nAGQ 1")
  )%>%
  group_by(n,reps) %>%
  mutate(MedianAIC_n_reps = median(AIC),
         Median_pvalue_Accuracy_n_reps = median(Pvalues_Accuracy),
         Median_pvalue_Interaction_n_reps = median(Pvalues_Interaction)) %>%
  group_by(n,reps,label) %>%
  mutate(MedianDuration = median(Duration),
         SE_Duration_n_reps_label = SE(Duration),
         MedianAIC_Difference = median(AIC)-MedianAIC_n_reps,
         Median_Pvalue_Accuracy_Difference = median(Pvalues_Accuracy) - Median_pvalue_Accuracy_n_reps,
         Median_Pvalue_Interaction_Difference = median(Pvalues_Interaction) - Median_pvalue_Interaction_n_reps)


ggplot(Dataframe_pvalues,aes(Pvalues_Accuracy, color = Optimizer)) +
  geom_density(size = 2) +
  coord_cartesian(ylim = c(0,5)) +
  facet_grid(n~reps)
ggsave("Figures/Duration for each Optimizer.jpg", w = 10, h = 5)


ggplot(Dataframe_pvalues,aes(round(Pvalues_Interaction,2), color = Optimizer, fill = Optimizer)) +
  geom_histogram(bins = 20) +
  facet_grid(.~Optimizer)


Dataframe_pvalues$Bin_Accuracy = 0
Dataframe_pvalues$Bin_Interaction = 0
for (i in (1:length(Dataframe_pvalues$iteration))){
  
  Bins = seq(0.025,0.975,0.05)
  Dataframe_pvalues$Bin_Accuracy[i] = Bins[which.min(abs(Bins-Dataframe_pvalues$Pvalues_Accuracy[i]))]+0.025
  Dataframe_pvalues$Bin_Interaction[i] = Bins[which.min(abs(Bins-Dataframe_pvalues$Pvalues_Interaction[i]))]+0.025
  print(i)
    
}

Dataframe_pvalues = Dataframe_pvalues %>%
  group_by(Bin_Accuracy,Optimizer) %>%
  mutate(BinCountAccuracy = length(Bin_Accuracy))%>%
  group_by(Bin_Interaction,Optimizer) %>%
  mutate(BinCountInteraction = length(Bin_Interaction))

ggplot(Dataframe_pvalues,aes(Bin_Accuracy,Optimizer, fill = BinCountAccuracy)) +
  geom_tile() +
  xlab("")

ggplot(Dataframe_pvalues,aes(Bin_Interaction,Optimizer, fill = BinCountInteraction)) +
  geom_tile() +
  xlab("")
