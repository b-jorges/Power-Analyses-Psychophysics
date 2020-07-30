require(dplyr)
require(tidyverse)
require(lme4)
require(purrr)

set.seed(1)

nParticipants = 10
ConditionOfInterest = c(0,1)
StandardValues = c(5,6,7,8)
reps = 1:100
PSE_Difference = -0.1
JND_Difference = 0.25
Multiplicator_PSE_Standard = 0
Multiplicator_SD_Standard = 0.15
Type_ResponseFunction = "normal"
SD_ResponseFunction = 0.20
Mean_Variability_Between = 0.15
SD_Variability_Between = 0.15
Dataframe = c()

FitCumGaussian = function(par,Difference,Prediction){
  (mean(pnorm(Difference,par[1],par[2])-Prediction)^2)
}

for (i in 1:50){
  
  Beginning = Sys.time()
  Dataframe1 = SimulatePsychometricData(nParticipants,
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
  
  Model1 = glm(cbind(Yes, Total - Yes) ~ ConditionOfInterest*Difference, 
               family = binomial(link = "probit"),
               data = Dataframe1,
               )
  Dataframe1$Prediction_Model1 = predict(Model1, type = "response", newdata = Dataframe1)
  
  Model2 = glmer(cbind(Yes, Total - Yes) ~ ConditionOfInterest*Difference + 
                   (1| ID), 
                 family = binomial(link = "probit"),
                 data = Dataframe1,
                 nAGQ = 0,
                 glmerControl(optimizer = "nloptwrap"))
  Dataframe1$Prediction_Model2 = predict(Model2, type = "response", newdata = Dataframe1)
  
  Model3 = glmer(cbind(Yes, Total - Yes) ~ ConditionOfInterest*Difference + 
                   (1 + Difference| ID), 
                 family = binomial(link = "probit"),
                 data = Dataframe1,
                 nAGQ = 0,
                 glmerControl(optimizer = "nloptwrap"))
  Dataframe1$Prediction_Model3 = predict(Model3, type = "response", newdata = Dataframe1)
  
  Model4 = glmer(cbind(Yes, Total - Yes) ~ ConditionOfInterest*Difference + 
                   (1 + ConditionOfInterest| ID), 
                 family = binomial(link = "probit"),
                 data = Dataframe1,
                 nAGQ = 0,
                 glmerControl(optimizer = "nloptwrap"))
  Dataframe1$Prediction_Model4 = predict(Model4, type = "response", newdata = Dataframe1)
  
  Model5 = glmer(cbind(Yes, Total - Yes) ~ ConditionOfInterest*Difference + 
                   (1 + Difference + ConditionOfInterest| ID), 
                 family = binomial(link = "probit"),
                 data = Dataframe1,
                 nAGQ = 0,
                 glmerControl(optimizer = "nloptwrap"))
  Dataframe1$Prediction_Model5 = predict(Model5, type = "response", newdata = Dataframe1)
  
  Model6 = glmer(cbind(Yes, Total - Yes) ~ ConditionOfInterest*Difference + 
                   (1|StandardValues), 
                 family = binomial(link = "probit"),
                 data = Dataframe1,
                 nAGQ = 0,
                 glmerControl(optimizer = "nloptwrap"))
  Dataframe1$Prediction_Model6 = predict(Model6, type = "response", newdata = Dataframe1)
  
  Model7 = glmer(cbind(Yes, Total - Yes) ~ ConditionOfInterest*Difference + 
                   (1| ID) +
                   (1|StandardValues), 
                 family = binomial(link = "probit"),
                 data = Dataframe1,
                 nAGQ = 0,
                 glmerControl(optimizer = "nloptwrap"))
  Dataframe1$Prediction_Model7 = predict(Model7, type = "response", newdata = Dataframe1)
  
  Model8 = glmer(cbind(Yes, Total - Yes) ~ ConditionOfInterest*Difference + 
                   (1 + Difference| ID) +
                   (1|StandardValues), 
                 family = binomial(link = "probit"),
                 data = Dataframe1,
                 nAGQ = 0,
                 glmerControl(optimizer = "nloptwrap"))
  Dataframe1$Prediction_Model8 = predict(Model8, type = "response", newdata = Dataframe1)
  
  Model9 = glmer(cbind(Yes, Total - Yes) ~ ConditionOfInterest*Difference + 
                   (1 + ConditionOfInterest| ID) +
                   (1|StandardValues), 
                 family = binomial(link = "probit"),
                 data = Dataframe1,
                 nAGQ = 0,
                 glmerControl(optimizer = "nloptwrap"))
  Dataframe1$Prediction_Model9 = predict(Model9, type = "response", newdata = Dataframe1)
  
  Model10 = glmer(cbind(Yes, Total - Yes) ~ ConditionOfInterest*Difference + 
                    (1 + Difference + ConditionOfInterest| ID) +
                    (1|StandardValues), 
                  family = binomial(link = "probit"),
                  data = Dataframe1,
                  nAGQ = 0,
                  glmerControl(optimizer = "nloptwrap"))
  Dataframe1$Prediction_Model10 = predict(Model10, type = "response", newdata = Dataframe1)
  
  Model11 = glmer(cbind(Yes, Total - Yes) ~ ConditionOfInterest*Difference +
                    (1 + Difference |StandardValues), 
                  family = binomial(link = "probit"),
                  data = Dataframe1,
                  nAGQ = 0,
                  glmerControl(optimizer = "nloptwrap"))
  Dataframe1$Prediction_Model11 = predict(Model11, type = "response", newdata = Dataframe1)
  
  Model12 = glmer(cbind(Yes, Total - Yes) ~ ConditionOfInterest*Difference + 
                    (1| ID) +
                    (1 + Difference |StandardValues), 
                  family = binomial(link = "probit"),
                  data = Dataframe1,
                  nAGQ = 0,
                  glmerControl(optimizer = "nloptwrap"))
  Dataframe1$Prediction_Model12 = predict(Model12, type = "response", newdata = Dataframe1)
  
  Model13 = glmer(cbind(Yes, Total - Yes) ~ ConditionOfInterest*Difference + 
                    (1 + Difference| ID) +
                    (1 + Difference |StandardValues), 
                  family = binomial(link = "probit"),
                  data = Dataframe1,
                  nAGQ = 0,
                  glmerControl(optimizer = "nloptwrap"))
  Dataframe1$Prediction_Model13 = predict(Model13, type = "response", newdata = Dataframe1)
  
  Model14 = glmer(cbind(Yes, Total - Yes) ~ ConditionOfInterest*Difference + 
                    (1 + ConditionOfInterest| ID) +
                    (1 + Difference |StandardValues), 
                  family = binomial(link = "probit"),
                  data = Dataframe1,
                  nAGQ = 0,
                  glmerControl(optimizer = "nloptwrap"))
  Dataframe1$Prediction_Model14 = predict(Model14, type = "response", newdata = Dataframe1)
  
  Model15 = glmer(cbind(Yes, Total - Yes) ~ ConditionOfInterest*Difference + 
                    (1 + Difference + ConditionOfInterest| ID) +
                    (1 + Difference |StandardValues), 
                  family = binomial(link = "probit"),
                  data = Dataframe1,
                  nAGQ = 0,
                  glmerControl(optimizer = "nloptwrap"))
  Dataframe1$Prediction_Model15 = predict(Model15, type = "response", newdata = Dataframe1)
  
  Model16 = glmer(cbind(Yes, Total - Yes) ~ ConditionOfInterest*Difference + 
                    (1 + ConditionOfInterest|StandardValues), 
                  family = binomial(link = "probit"),
                  data = Dataframe1,
                  nAGQ = 0,
                  glmerControl(optimizer = "nloptwrap"))
  Dataframe1$Prediction_Model16 = predict(Model16, type = "response", newdata = Dataframe1)
  
  Model17 = glmer(cbind(Yes, Total - Yes) ~ ConditionOfInterest*Difference + 
                    (1| ID) +
                    (1 + ConditionOfInterest|StandardValues), 
                  family = binomial(link = "probit"),
                  data = Dataframe1,
                  nAGQ = 0,
                  glmerControl(optimizer = "nloptwrap"))
  Dataframe1$Prediction_Model17 = predict(Model17, type = "response", newdata = Dataframe1)
  
  Model18 = glmer(cbind(Yes, Total - Yes) ~ ConditionOfInterest*Difference + 
                    (1 + Difference| ID) +
                    (1 + ConditionOfInterest|StandardValues), 
                  family = binomial(link = "probit"),
                  data = Dataframe1,
                  nAGQ = 0,
                  glmerControl(optimizer = "nloptwrap"))
  Dataframe1$Prediction_Model18 = predict(Model18, type = "response", newdata = Dataframe1)
  
  Model19 = glmer(cbind(Yes, Total - Yes) ~ ConditionOfInterest*Difference + 
                    (1 + ConditionOfInterest| ID) +
                    (1 + ConditionOfInterest|StandardValues), 
                  family = binomial(link = "probit"),
                  data = Dataframe1,
                  nAGQ = 0,
                  glmerControl(optimizer = "nloptwrap"))
  Dataframe1$Prediction_Model19 = predict(Model19, type = "response", newdata = Dataframe1)
  
  Model20 = glmer(cbind(Yes, Total - Yes) ~ ConditionOfInterest*Difference + 
                    (1 + Difference + ConditionOfInterest| ID) +
                    (1 + ConditionOfInterest|StandardValues), 
                  family = binomial(link = "probit"),
                  data = Dataframe1,
                  nAGQ = 0,
                  glmerControl(optimizer = "nloptwrap"))
  Dataframe1$Prediction_Model20 = predict(Model20, type = "response", newdata = Dataframe1)
  
  Model21 = glmer(cbind(Yes, Total - Yes) ~ ConditionOfInterest*Difference +
                    (1 + Difference + ConditionOfInterest|StandardValues), 
                  family = binomial(link = "probit"),
                  data = Dataframe1,
                  nAGQ = 0,
                  glmerControl(optimizer = "nloptwrap"))
  Dataframe1$Prediction_Model21 = predict(Model21, type = "response", newdata = Dataframe1)
  
  Model22 = glmer(cbind(Yes, Total - Yes) ~ ConditionOfInterest*Difference + 
                    (1| ID) +
                    (1 + Difference + ConditionOfInterest|StandardValues), 
                  family = binomial(link = "probit"),
                  data = Dataframe1,
                  nAGQ = 0,
                  glmerControl(optimizer = "nloptwrap"))
  Dataframe1$Prediction_Model22 = predict(Model22, type = "response", newdata = Dataframe1)
  
  Model23 = glmer(cbind(Yes, Total - Yes) ~ ConditionOfInterest*Difference + 
                    (1 + Difference| ID) +
                    (1 + Difference + ConditionOfInterest|StandardValues), 
                  family = binomial(link = "probit"),
                  data = Dataframe1,
                  nAGQ = 0,
                  glmerControl(optimizer = "nloptwrap"))
  Dataframe1$Prediction_Model23 = predict(Model23, type = "response", newdata = Dataframe1)
  
  Model24 = glmer(cbind(Yes, Total - Yes) ~ ConditionOfInterest*Difference + 
                    (1 + ConditionOfInterest| ID) +
                    (1 + Difference + ConditionOfInterest|StandardValues), 
                  family = binomial(link = "probit"),
                  data = Dataframe1,
                  nAGQ = 0,
                  glmerControl(optimizer = "nloptwrap"))
  Dataframe1$Prediction_Model24 = predict(Model24, type = "response", newdata = Dataframe1)
  
  Model25 = glmer(cbind(Yes, Total - Yes) ~ ConditionOfInterest*Difference + 
                    (1 + Difference + ConditionOfInterest| ID) +
                    (1 + Difference + ConditionOfInterest|StandardValues), 
                  family = binomial(link = "probit"),
                  data = Dataframe1,
                  nAGQ = 0,
                  glmerControl(optimizer = "nloptwrap"))
  Dataframe1$Prediction_Model25 = predict(Model25, type = "response", newdata = Dataframe1)
  
  Dataframe1 = Dataframe1 %>%
    group_by(ConditionOfInterest,ID,StandardValues) %>%
    mutate(Mean_Model1 = optim(par = c(0,0.2),FitCumGaussian,par,Difference,Prediction_Model1)$par[1],
           SD_Model1 = optim(par = c(0,0.2),FitCumGaussian,par,Difference,Prediction_Model1)$par[2],
           Mean_Model2 = optim(par = c(0,0.2),FitCumGaussian,par,Difference,Prediction_Model2)$par[1],
           SD_Model2 = optim(par = c(0,0.2),FitCumGaussian,par,Difference,Prediction_Model2)$par[2],
           Mean_Model3 = optim(par = c(0,0.2),FitCumGaussian,par,Difference,Prediction_Model3)$par[1],
           SD_Model3 = optim(par = c(0,0.2),FitCumGaussian,par,Difference,Prediction_Model3)$par[2],
           Mean_Model4 = optim(par = c(0,0.2),FitCumGaussian,par,Difference,Prediction_Model4)$par[1],
           SD_Model4 = optim(par = c(0,0.2),FitCumGaussian,par,Difference,Prediction_Model4)$par[2],
           Mean_Model5 = optim(par = c(0,0.2),FitCumGaussian,par,Difference,Prediction_Model5)$par[1],
           SD_Model5 = optim(par = c(0,0.2),FitCumGaussian,par,Difference,Prediction_Model5)$par[2],
           Mean_Model6 = optim(par = c(0,0.2),FitCumGaussian,par,Difference,Prediction_Model6)$par[1],
           SD_Model6 = optim(par = c(0,0.2),FitCumGaussian,par,Difference,Prediction_Model6)$par[2],
           Mean_Model7 = optim(par = c(0,0.2),FitCumGaussian,par,Difference,Prediction_Model7)$par[1],
           SD_Model7 = optim(par = c(0,0.2),FitCumGaussian,par,Difference,Prediction_Model7)$par[2],
           Mean_Model8 = optim(par = c(0,0.2),FitCumGaussian,par,Difference,Prediction_Model8)$par[1],
           SD_Model8 = optim(par = c(0,0.2),FitCumGaussian,par,Difference,Prediction_Model8)$par[2],
           Mean_Model9 = optim(par = c(0,0.2),FitCumGaussian,par,Difference,Prediction_Model9)$par[1],
           SD_Model9 = optim(par = c(0,0.2),FitCumGaussian,par,Difference,Prediction_Model9)$par[2],
           Mean_Model10 = optim(par = c(0,0.2),FitCumGaussian,par,Difference,Prediction_Model10)$par[1],
           SD_Model10 = optim(par = c(0,0.2),FitCumGaussian,par,Difference,Prediction_Model10)$par[2],
           Mean_Model11 = optim(par = c(0,0.2),FitCumGaussian,par,Difference,Prediction_Model11)$par[1],
           SD_Model11 = optim(par = c(0,0.2),FitCumGaussian,par,Difference,Prediction_Model11)$par[2],
           Mean_Model12 = optim(par = c(0,0.2),FitCumGaussian,par,Difference,Prediction_Model12)$par[1],
           SD_Model12 = optim(par = c(0,0.2),FitCumGaussian,par,Difference,Prediction_Model12)$par[2],
           Mean_Model13 = optim(par = c(0,0.2),FitCumGaussian,par,Difference,Prediction_Model13)$par[1],
           SD_Model13 = optim(par = c(0,0.2),FitCumGaussian,par,Difference,Prediction_Model13)$par[2],
           Mean_Model14 = optim(par = c(0,0.2),FitCumGaussian,par,Difference,Prediction_Model14)$par[1],
           SD_Model14 = optim(par = c(0,0.2),FitCumGaussian,par,Difference,Prediction_Model14)$par[2],
           Mean_Model15 = optim(par = c(0,0.2),FitCumGaussian,par,Difference,Prediction_Model15)$par[1],
           SD_Model15 = optim(par = c(0,0.2),FitCumGaussian,par,Difference,Prediction_Model15)$par[2],
           Mean_Model16 = optim(par = c(0,0.2),FitCumGaussian,par,Difference,Prediction_Model16)$par[1],
           SD_Model16 = optim(par = c(0,0.2),FitCumGaussian,par,Difference,Prediction_Model16)$par[2],
           Mean_Model17 = optim(par = c(0,0.2),FitCumGaussian,par,Difference,Prediction_Model17)$par[1],
           SD_Model17 = optim(par = c(0,0.2),FitCumGaussian,par,Difference,Prediction_Model17)$par[2],
           Mean_Model18 = optim(par = c(0,0.2),FitCumGaussian,par,Difference,Prediction_Model18)$par[1],
           SD_Model18 = optim(par = c(0,0.2),FitCumGaussian,par,Difference,Prediction_Model18)$par[2],
           Mean_Model19 = optim(par = c(0,0.2),FitCumGaussian,par,Difference,Prediction_Model19)$par[1],
           SD_Model19 = optim(par = c(0,0.2),FitCumGaussian,par,Difference,Prediction_Model19)$par[2],
           Mean_Model20 = optim(par = c(0,0.2),FitCumGaussian,par,Difference,Prediction_Model20)$par[1],
           SD_Model20 = optim(par = c(0,0.2),FitCumGaussian,par,Difference,Prediction_Model20)$par[2],
           Mean_Model21 = optim(par = c(0,0.2),FitCumGaussian,par,Difference,Prediction_Model21)$par[1],
           SD_Model21 = optim(par = c(0,0.2),FitCumGaussian,par,Difference,Prediction_Model21)$par[2],
           Mean_Model22 = optim(par = c(0,0.2),FitCumGaussian,par,Difference,Prediction_Model22)$par[1],
           SD_Model22 = optim(par = c(0,0.2),FitCumGaussian,par,Difference,Prediction_Model22)$par[2],
           Mean_Model23 = optim(par = c(0,0.2),FitCumGaussian,par,Difference,Prediction_Model23)$par[1],
           SD_Model23 = optim(par = c(0,0.2),FitCumGaussian,par,Difference,Prediction_Model23)$par[2],
           Mean_Model24 = optim(par = c(0,0.2),FitCumGaussian,par,Difference,Prediction_Model24)$par[1],
           SD_Model24 = optim(par = c(0,0.2),FitCumGaussian,par,Difference,Prediction_Model24)$par[2],
           Mean_Model25 = optim(par = c(0,0.2),FitCumGaussian,par,Difference,Prediction_Model25)$par[1],
           SD_Model25 = optim(par = c(0,0.2),FitCumGaussian,par,Difference,Prediction_Model25)$par[2])
  
  Dataframe = rbind(Dataframe,Dataframe1 %>% 
                      group_by(ConditionOfInterest,ID,StandardValues) %>%
                      slice(1))
  
  print(Sys.time() - Beginning)
}


Dataframe2 = data.frame(Model = rep(paste0("Model",1:25),each=nParticipants*length(ConditionOfInterest)*length(StandardValues)*50),
                        StandardValues = rep(Dataframe$StandardValues,25),
                        ConditionOfInterest = rep(Dataframe$ConditionOfInterest,25),
                        ID = rep(Dataframe$ID,25),
                        Mean = c(Dataframe$Mean_Model1,Dataframe$Mean_Model2,Dataframe$Mean_Model3,Dataframe$Mean_Model4,
                                Dataframe$Mean_Model5,Dataframe$Mean_Model6,Dataframe$Mean_Model7,Dataframe$Mean_Model8,
                                Dataframe$Mean_Model9,Dataframe$Mean_Model10,Dataframe$Mean_Model11,Dataframe$Mean_Model12,
                                Dataframe$Mean_Model13,Dataframe$Mean_Model14,Dataframe$Mean_Model15,Dataframe$Mean_Model16,
                                Dataframe$Mean_Model17,Dataframe$Mean_Model18,Dataframe$Mean_Model19,Dataframe$Mean_Model20,
                                Dataframe$Mean_Model21,Dataframe$Mean_Model22,Dataframe$Mean_Model23,Dataframe$Mean_Model24,
                                Dataframe$Mean_Model25),
                        SD = c(Dataframe$SD_Model1,Dataframe$SD_Model2,Dataframe$SD_Model3,Dataframe$SD_Model4,
                                Dataframe$SD_Model5,Dataframe$SD_Model6,Dataframe$SD_Model7,Dataframe$SD_Model8,
                                Dataframe$SD_Model9,Dataframe$SD_Model10,Dataframe$SD_Model11,Dataframe$SD_Model12,
                                Dataframe$SD_Model13,Dataframe$SD_Model14,Dataframe$SD_Model15,Dataframe$SD_Model16,
                                Dataframe$SD_Model17,Dataframe$SD_Model18,Dataframe$SD_Model19,Dataframe$SD_Model20,
                                Dataframe$SD_Model21,Dataframe$SD_Model22,Dataframe$SD_Model23,Dataframe$SD_Model24,
                                Dataframe$SD_Model25))
Dataframe2 = Dataframe2 %>% 
  mutate(ActualPSEs = case_when(
    ConditionOfInterest == 1 ~ -0.1*StandardValues,
    ConditionOfInterest == 0 ~ 0))

ggplot(Dataframe2,aes(Model,Mean)) +
  geom_boxplot() + 
  facet_grid(ConditionOfInterest~StandardValues) +
  geom_hline(aes(yintercept = Actual_PSEs))

ggplot(Dataframe2,aes(Model,SD)) +
  geom_boxplot() + 
  facet_grid(ConditionOfInterest~StandardValues)