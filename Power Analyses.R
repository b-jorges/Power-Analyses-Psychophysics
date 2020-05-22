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

haha = anova(GLMM,GLMM2)
haha$`Pr(>Chisq)`[2]
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
Dataframe_wide = rbind(read.csv(header = T, file = paste0(Where_Am_I(),"/Data/Final Sample Powers.csv"))) %>%
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


Powers1 = ggplot(Dataframe_Powers %>% filter(reps == 40), aes(n,power,color = label)) +
  geom_line(size = 2) +
  facet_grid(JND_Difference~PSE_Difference) +
  xlab("N° of Subjects") +
  ylab("Power") +
  geom_hline(linetype = 2, yintercept = 0.8) +
  geom_hline(linetype = 1, yintercept = 0.05) +
  geom_hline(linetype = 3, yintercept = 0.95) +
  scale_x_continuous(breaks=c(10,15,20)) +
  scale_y_continuous(breaks=c(0.25,0.75)) +
  scale_color_manual(values = c(BlauUB,LightBlauUB,Red,LightRed), 
                     name = "") +
  ggtitle("A. 40 Repetitions")

Powers2 = ggplot(Dataframe_Powers %>% filter(reps == 60), aes(n,power,color = label)) +
  geom_line(size = 2) +
  facet_grid(JND_Difference~PSE_Difference) +
  xlab("N° of Subjects") +
  ylab("Power") +
  geom_hline(linetype = 2, yintercept = 0.8) +
  geom_hline(linetype = 1, yintercept = 0.05) +
  geom_hline(linetype = 3, yintercept = 0.95) +
  scale_x_continuous(breaks=c(10,15,20)) +
  scale_y_continuous(breaks=c(0.25,0.75)) +
  scale_color_manual(values = c(BlauUB,LightBlauUB,Red,LightRed), 
                     name = "") +
  ggtitle("A. 60 Repetitions")

Powers3 = ggplot(Dataframe_Powers %>% filter(reps == 100), aes(n,power,color = label)) +
  geom_line(size = 2) +
  facet_grid(JND_Difference~PSE_Difference) +
  xlab("N° of Subjects") +
  ylab("Power") +
  geom_hline(linetype = 2, yintercept = 0.8) +
  geom_hline(linetype = 1, yintercept = 0.05) +
  geom_hline(linetype = 3, yintercept = 0.95) +
  scale_x_continuous(breaks=c(10,15,20)) +
  scale_y_continuous(breaks=c(0.25,0.75)) +
  scale_color_manual(values = c(BlauUB,LightBlauUB,Red,LightRed), 
                     name = "") +
  ggtitle("A. 100 Repetitions")

plot_shared_legend(Powers1,Powers2)
ggsave("Figures/Powers.jpg", w = 12, h = 8)
########################################################################
########################################################################
########################################################################





########################################################################
####################Compare Optimizer Configurations####################
########################################################################
Dataframe_pvalues = rbind(read.csv(header = T, file = paste0(Where_Am_I(),"/Data/Pvalues_Julia2.csv")))


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
    label == "JuliaLRT" ~ "Julia: BOBYQA, fast, LRT",
    label == "R: LRT" ~ "R: LRT")
    )%>%
  group_by(n,reps,label,PSE_Difference,JND_Difference) %>%
  mutate(MedianDuration = median(Duration),
         SEDuration = SE(Duration),
         SE_Duration_n_reps_label = SE(Duration),
         nTrials = case_when(
           reps == 30 ~ "30 repetitions",
           reps == 40 ~ "40 repetitions",
           reps == 50 ~ "50 repetitions",
           reps == 60 ~ "60 repetitions"),
         Effect = case_when(
           PSE_Difference == 0.025 ~ "Small Effect",
           PSE_Difference == 0 ~ "No Effect")
  )
         
Dataframe_pvalues = Dataframe_pvalues %>%         
  group_by(reps,n,iteration,Effect) %>% 
  mutate(AIC_Ratio = AIC/AIC[2]) %>% 
  group_by(reps,n,Optimizer) %>% 
  mutate(Median_AIC_Ratio = median(AIC_Ratio))


#######Timing
TimingPlot1 = ggplot(Dataframe_pvalues %>% 
         filter(Optimizer != "Julia: BOBYQA, fast, LRT" & Optimizer != "R: LRT"),
       aes(n,MedianDuration, color = Optimizer)) +
  geom_point(size=2) +
  geom_line(size=1) +
  facet_grid(Effect~nTrials) +
  ylab("Median Duration (s)") +
  scale_color_manual(values = rainbow(10), name = "Method") +
  scale_x_continuous(breaks = c(10,15,20)) +
  ggtitle("A. All Configurations")


TimingPlot2 = ggplot(Dataframe_pvalues %>% 
         filter(Optimizer  %in% c("Julia: BOBYQA, fast",
                                  "Julia: BOBYQA, slow",
                                  "R: nloptwrap, slow",
                                  "R: nloptwrap, fast")),
       aes(n,MedianDuration, color = Optimizer)) +
  geom_point(size=2) +
  geom_line(size=1) +
  facet_grid(Effect~nTrials) +
  ylab("Median Duration (s)") +
  scale_color_manual(values = rainbow(10)[c(1,2,9,10)], name = "Method") +
  scale_x_continuous(breaks = c(10,15,20)) +
  ggtitle("B. Fastest Confirgurations")

plot_shared_legend(TimingPlot1,TimingPlot2)
ggsave("Figures/Different Durations.jpg",w=12,h=6)


#######False Positives
Dataframe_pvalues$Bin_Accuracy = 0
Dataframe_pvalues$Bin_Interaction = 0
for (i in (1:length(Dataframe_pvalues$iteration))){
  print(i)  
  Bins = seq(0.025,0.975,0.05)
  Dataframe_pvalues$Bin_Accuracy[i] = Bins[which.min(abs(Bins-Dataframe_pvalues$Pvalues_Accuracy[i]))]+0.025
  Dataframe_pvalues$Bin_Interaction[i] = Bins[which.min(abs(Bins-Dataframe_pvalues$Pvalues_Interaction[i]))]+0.025

}

Dataframe_pvalues = Dataframe_pvalues %>%
  group_by(Bin_Accuracy,Optimizer,PSE_Difference) %>%
  mutate(BinCountAccuracy = length(Bin_Accuracy))%>%
  group_by(Bin_Interaction,Optimizer,PSE_Difference) %>%
  mutate(BinCountInteraction = length(Bin_Interaction))

PlotAccuracy = ggplot(Dataframe_pvalues %>% 
         filter(Optimizer != "Julia: BOBYQA, fast, LRT" & Optimizer != "R: LRT" &
                  Effect == "No Effect"),
       aes(Bin_Accuracy-0.025,Optimizer, fill = as.factor(BinCountAccuracy))) +
  geom_tile() +
  xlab("") +
  scale_fill_manual(values=colorRampPalette(c(BlauUB,Red,Yellow))(37)) +
  theme(legend.position = "") +
  ylab("") +
  ggtitle("A. Accuracy, No Effect")

PlotInteraction = ggplot(Dataframe_pvalues %>% filter(Effect == "No Effect"),
                aes(Bin_Interaction-0.025,Optimizer, fill = as.factor(BinCountInteraction))) +
  geom_tile() +
  xlab("") +
  scale_fill_manual(values=colorRampPalette(c(BlauUB,Red, Yellow))(45)) +
  theme(legend.position = "") +
  ylab("") +
  ggtitle("B. Interaction, No Effect")

PlotAccuracy2 = ggplot(Dataframe_pvalues %>% 
                        filter(Effect == "Small Effect" &
                                 Optimizer != "Julia: BOBYQA, fast, LRT" &
                                 Optimizer != "R: LRT"),
                      aes(Bin_Accuracy-0.025,Optimizer, fill = as.factor(BinCountAccuracy))) +
  geom_tile() +
  xlab("") +
  ylab("") +
  scale_fill_manual(values=colorRampPalette(c(BlauUB,Red,Yellow))(38)) +
  theme(legend.position = "") +
  ggtitle("C. Accuracy, Small Effect")

PlotInteraction2 = ggplot(Dataframe_pvalues %>% 
                           filter(Effect == "Small Effect"),
                         aes(Bin_Interaction-0.025,Optimizer, fill = as.factor(BinCountInteraction))) +
  geom_tile() +
  xlab("") +
  ylab("") +
  scale_fill_manual(values=colorRampPalette(c(BlauUB,Red,Yellow))(46)) +
  theme(legend.position = "") +
  ggtitle("D. Interaction, Small Effect")
plot_grid(PlotAccuracy,PlotInteraction,PlotAccuracy2,PlotInteraction2, nrow = 2)
ggsave("Figures/False Positives.jpg",w=12,h=8)


############AICs
ggplot(Dataframe_pvalues %>%
         filter(Optimizer != "Julia: BOBYQA, fast, LRT" & 
                  Optimizer != "R: LRT"),
       aes(n,Median_AIC_Ratio, color = Optimizer)) +
  geom_point(size=2) +
  geom_line(size=1) +
  facet_grid(Effect~nTrials) +
  ylab("AIC DIfference") +
  scale_color_manual(values = rainbow(10), name = "Method") +
  scale_x_continuous(breaks = c(10,15,20))
ggsave("Figures/AIC differences.jpg",w=12,h=6)
