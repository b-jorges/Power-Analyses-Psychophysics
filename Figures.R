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
source("Utilities/PowerFunctions.r")
set.seed(9121)


##################################################################
#############FIGURE 1: one psychometric function #################
##################################################################
Psychometric = SimulatePsychometricFunction_Staircase(ID = paste("s",1), 
                                       ConditionOfInterest = c(0,1), 
                                       StandardValues = c(8), 
                                       reps = 1:1000, 
                                       PSE_Difference = 0, 
                                       JND_Difference = 0, 
                                       Multiplicator_PSE_Standard = 0, 
                                       Multiplicator_SD_Standard = 0.05, 
                                       SD_ResponseFunction = 0.1, 
                                       Mean_Variability_Between = 0, 
                                       SD_Variability_Between = 0)

PsychometricFunctions = quickpsy(Psychometric,Difference,Answer, bootstrap = "none")

ggplot(PsychometricFunctions$curves, aes(x,y)) +
  geom_line(size = 2) +
  geom_point(data = Psychometric, aes(Difference, Answer), alpha = 0.2) +
  xlab("Difference between Comparison and Test") +
  ylab("Probability to choose Test") +
  coord_cartesian(xlim = c(PsychometricFunctions$curves$x[which.min(abs(PsychometricFunctions$curves$y-0.5))]-2.5,
                           PsychometricFunctions$curves$x[which.min(abs(PsychometricFunctions$curves$y-0.5))]+2.5)) +
  geom_vline(linetype = 2, xintercept = PsychometricFunctions$curves$x[which.min(abs(PsychometricFunctions$curves$y-0.5))], color = "grey") +
  geom_vline(linetype = 1, xintercept = PsychometricFunctions$curves$x[which.min(abs(PsychometricFunctions$curves$y-0.75))], color = "grey") +
  geom_vline(linetype = 1, xintercept = PsychometricFunctions$curves$x[which.min(abs(PsychometricFunctions$curves$y-0.25))], color = "grey") +
  geom_hline(linetype = 2, yintercept = 0.5, color = "grey") +
  geom_hline(linetype = 1, yintercept = 0.25, color = "grey") +
  geom_hline(linetype = 1, yintercept = 0.75, color = "grey") +
  geom_segment(aes(x=PsychometricFunctions$curves$x[which.min(abs(PsychometricFunctions$curves$y-0.25))],
                   y=0.25,xend=PsychometricFunctions$curves$x[which.min(abs(PsychometricFunctions$curves$y-0.5))],
                   yend = 0.25),
               color = Yellow, 
               size = 2) +
  geom_segment(aes(x=PsychometricFunctions$curves$x[which.min(abs(PsychometricFunctions$curves$y-0.5))],
                   y=0.75,
                   xend=PsychometricFunctions$curves$x[which.min(abs(PsychometricFunctions$curves$y-0.75))],
                   yend = 0.75),
               color = Yellow, 
               size = 2) +
  geom_segment(aes(x=PsychometricFunctions$curves$x[which.min(abs(PsychometricFunctions$curves$y-0.5))],
                   y=0,
                   xend=PsychometricFunctions$curves$x[which.min(abs(PsychometricFunctions$curves$y-0.5))],
                   yend = 0.5),
               color = Red, 
               size = 2)
ggsave("Figure01_PsychometricFunction.jpg", w = 4, h = 4)
##################################################################
#################################FIGURE 1 END#####################
##################################################################


##################################################################
#############FIGURE 2: response distributions ####################
##################################################################

ResponseDistributions = data.frame(
  Value=c(rcauchy(1650,1,0.05),
          rnorm(1650,1,0.05),
          rep(c(1),1)),
  label = c(rep("Cauchy",1650),
            rep("Normal",1650),
            rep("Uniform",1))
) %>% filter(abs(Value-1) < 0.5)

ggplot(ResponseDistributions %>% filter(label %in% c("Cauchy","Normal", "Uniform")), aes(Value,color = label)) +
  geom_density(size=2) +
  geom_segment(aes(x=0.7,y=2,xend=1.3,yend = 2),color = Yellow, size = 2) +
  geom_segment(aes(x=0.7,y=0,xend=0.7,yend = 2),color = Yellow, size = 2) +
  geom_segment(aes(x=1.3,y=0,xend=1.3,yend = 2),color = Yellow, size = 2) +
  coord_cartesian(xlim=c(0.65,1.35)) +
  xlab("Stimulus Intensity") +
  ylab("Density") +
  scale_color_manual(name = "Distribution\nType",
                     values = c(Red,BlauUB, Yellow),
                     labels = c("Cauchy","Gaussian", "Uniform"))
ggsave("Figure1 Distributions.jpg", w = 6, h = 4)

##################################################################
#################################FIGURE 2 END#####################
##################################################################


##################################################################
#############FIGURE 3: five psychometric functions #################
##################################################################
Psychometric = SimulatePsychometricFunction_Staircase(ID = paste0("s",1:5), 
                                                      ConditionOfInterest = c(0,1), 
                                                      StandardValues = c(5,8), 
                                                      reps = 1:55, 
                                                      PSE_Difference = 0.2, 
                                                      JND_Difference = 0.4, 
                                                      Multiplicator_PSE_Standard = 0, 
                                                      Multiplicator_SD_Standard = 0.108, 
                                                      SD_ResponseFunction = 0.1, 
                                                      Mean_Variability_Between = 0.1, 
                                                      SD_Variability_Between = 0.1)

PsychometricFunctions = quickpsy(Psychometric,Difference,Answer,grouping = .(ConditionOfInterest,ID,StandardValues), bootstrap = "none")

plot(PsychometricFunctions) +
  scale_color_manual(name = "",
                     values = c(Red,BlauUB),
                     labels = c("Control","Experimental")) +
  xlab("Difference between Comparison and Test") +
  ylab("Probability to choose Test") +
  geom_vline(linetype = 2, xintercept = 0, color = "grey") +
  geom_hline(linetype = 2, yintercept = 0.5, color = "grey")
ggsave("Figure3 Psychometric Functions.jpg", w = 10, h = 5)
##################################################################
#################################FIGURE 3 END#####################
##################################################################



##################################################################
#############FIGURE 3: LMM explainer #################
##################################################################

require(dplyr)
require(ggplot2)
require(cowplot)
require(lme4)
require(lmerTest)
theme_set(theme_cowplot())
ID = paste0("S0",1:3)
Condition = c("Baseline","Manipulation")
trial = 1:100
LMM_Example = expand.grid(ID=ID,Condition=Condition,trial=trial)
LMM_Example = LMM_Example %>% 
  mutate(Effect_ID = case_when(
            ID == "S01" ~ 1,
            ID == "S02" ~ -0.25,
            ID == "S03" ~ 0),
         Effect_Manipulation = case_when(
            Condition == "Baseline" ~ 0,
            Condition == "Manipulation" ~ 1
         ),
         Stimulus = runif(length(ID)),
         Response = rnorm(length(ID),Effect_ID +Stimulus*Effect_Manipulation,0.3))

Plot_LMM = ggplot(LMM_Example,aes(Stimulus,Response,color = ID)) +
  geom_point(size = 2) +
  facet_grid(.~Condition) +
  scale_color_manual(name = "",
                     values = c(Red, Yellow, BlauUB))

Plot_LM = ggplot(LMM_Example,aes(Stimulus,Response)) +
  geom_point(size = 2) +
  facet_grid(.~Condition)

LM_Fit = lm(Response ~ Stimulus*Condition,
   data = LMM_Example)
summary(LM_Fit)

LMM_Fit = lmer(Response ~ Stimulus*Condition + (1 | ID),
   data = LMM_Example)
summary(LMM_Fit)

plot_grid(Plot_LM,Plot_LMM)
ggsave("(Figure 4) LMM.jpg", w = 12, h = 5)
##################################################################
#################################FIGURE 3 END#####################
##################################################################