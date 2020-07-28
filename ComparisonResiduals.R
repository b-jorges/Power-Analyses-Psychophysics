require(dplyr)
require(tidyverse)
require(lme4)
require(purrr)

set.seed(1)

nParticipants = 10
ID = paste0("S0",1:5)
ConditionOfInterest = c(0,1)
StandardValues = c(5,6,7,8)
reps = 1:100
PSE_Difference = -0.1
JND_Difference = 0.25
Multiplicator_PSE_Standard = 0
Multiplicator_SD_Standard = 0.15
Type_ResponseFunction = "normal"
SD_ResponseFunction = 0.1
Mean_Variability_Between = 0.15
SD_Variability_Between = 0.15

for(i in 1:5){
  
  start_time <- Sys.time()

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
               data = Dataframe1)
Residuals1 = simulateResiduals(Model1)

Model2 = glmer(cbind(Yes, Total - Yes) ~ ConditionOfInterest*Difference + 
                 (1| ID), 
               family = binomial(link = "probit"),
               data = Dataframe1)
Residuals2 = simulateResiduals(Model2)

Model3 = glmer(cbind(Yes, Total - Yes) ~ ConditionOfInterest*Difference + 
                 (1 + Difference| ID), 
               family = binomial(link = "probit"),
               data = Dataframe1)
Residuals3 = simulateResiduals(Model3)

Model4 = glmer(cbind(Yes, Total - Yes) ~ ConditionOfInterest*Difference + 
                 (1 + ConditionOfInterest| ID), 
               family = binomial(link = "probit"),
               data = Dataframe1)
Residuals4 = simulateResiduals(Model4)

Model5 = glmer(cbind(Yes, Total - Yes) ~ ConditionOfInterest*Difference + 
                 (1 + Difference + ConditionOfInterest| ID), 
               family = binomial(link = "probit"),
               data = Dataframe1)
Residuals5 = simulateResiduals(Model5)

Model6 = glmer(cbind(Yes, Total - Yes) ~ ConditionOfInterest*Difference + 
                 (1|StandardValues), 
               family = binomial(link = "probit"),
               data = Dataframe1)
Residuals6 = simulateResiduals(Model6)

Model7 = glmer(cbind(Yes, Total - Yes) ~ ConditionOfInterest*Difference + 
                 (1| ID) +
                 (1|StandardValues), 
               family = binomial(link = "probit"),
               data = Dataframe1)
Residuals7 = simulateResiduals(Model7)

Model8 = glmer(cbind(Yes, Total - Yes) ~ ConditionOfInterest*Difference + 
                 (1 + Difference| ID) +
                 (1|StandardValues), 
               family = binomial(link = "probit"),
               data = Dataframe1)
Residuals8 = simulateResiduals(Model8)

Model9 = glmer(cbind(Yes, Total - Yes) ~ ConditionOfInterest*Difference + 
                 (1 + ConditionOfInterest| ID) +
                 (1|StandardValues), 
               family = binomial(link = "probit"),
               data = Dataframe1)
Residuals9 = simulateResiduals(Model9)

Model10 = glmer(cbind(Yes, Total - Yes) ~ ConditionOfInterest*Difference + 
                 (1 + Difference + ConditionOfInterest| ID) +
                 (1|StandardValues), 
               family = binomial(link = "probit"),
               data = Dataframe1)
Residuals10 = simulateResiduals(Model10)

Model11 = glmer(cbind(Yes, Total - Yes) ~ ConditionOfInterest*Difference +
                 (1 + Difference |StandardValues), 
               family = binomial(link = "probit"),
               data = Dataframe1)
Residuals11 = simulateResiduals(Model11)

Model12 = glmer(cbind(Yes, Total - Yes) ~ ConditionOfInterest*Difference + 
                 (1| ID) +
                 (1 + Difference |StandardValues), 
               family = binomial(link = "probit"),
               data = Dataframe1)
Residuals12 = simulateResiduals(Model12)

Model13 = glmer(cbind(Yes, Total - Yes) ~ ConditionOfInterest*Difference + 
                 (1 + Difference| ID) +
                 (1 + Difference |StandardValues), 
               family = binomial(link = "probit"),
               data = Dataframe1)
Residuals13 = simulateResiduals(Model13)

Model14 = glmer(cbind(Yes, Total - Yes) ~ ConditionOfInterest*Difference + 
                 (1 + ConditionOfInterest| ID) +
                 (1 + Difference |StandardValues), 
               family = binomial(link = "probit"),
               data = Dataframe1)
Residuals14 = simulateResiduals(Model14)

Model15 = glmer(cbind(Yes, Total - Yes) ~ ConditionOfInterest*Difference + 
                 (1 + Difference + ConditionOfInterest| ID) +
                 (1 + Difference |StandardValues), 
               family = binomial(link = "probit"),
               data = Dataframe1)
Residuals15 = simulateResiduals(Model15)

Model16 = glmer(cbind(Yes, Total - Yes) ~ ConditionOfInterest*Difference + 
                 (1 + ConditionOfInterest|StandardValues), 
               family = binomial(link = "probit"),
               data = Dataframe1)
Residuals16 = simulateResiduals(Model16)

Model17 = glmer(cbind(Yes, Total - Yes) ~ ConditionOfInterest*Difference + 
                 (1| ID) +
                 (1 + ConditionOfInterest|StandardValues), 
               family = binomial(link = "probit"),
               data = Dataframe1)
Residuals17 = simulateResiduals(Model17)

Model18 = glmer(cbind(Yes, Total - Yes) ~ ConditionOfInterest*Difference + 
                 (1 + Difference| ID) +
                 (1 + ConditionOfInterest|StandardValues), 
               family = binomial(link = "probit"),
               data = Dataframe1)
Residuals18 = simulateResiduals(Model18)

Model19 = glmer(cbind(Yes, Total - Yes) ~ ConditionOfInterest*Difference + 
                 (1 + ConditionOfInterest| ID) +
                 (1 + ConditionOfInterest|StandardValues), 
               family = binomial(link = "probit"),
               data = Dataframe1)
Residuals19 = simulateResiduals(Model19)

Model20 = glmer(cbind(Yes, Total - Yes) ~ ConditionOfInterest*Difference + 
                 (1 + Difference + ConditionOfInterest| ID) +
                 (1 + ConditionOfInterest|StandardValues), 
               family = binomial(link = "probit"),
               data = Dataframe1)
Residuals20 = simulateResiduals(Model20)

Model21 = glmer(cbind(Yes, Total - Yes) ~ ConditionOfInterest*Difference +
                 (1 + Difference + ConditionOfInterest|StandardValues), 
               family = binomial(link = "probit"),
               data = Dataframe1)
Residuals21 = simulateResiduals(Model21)

Model22 = glmer(cbind(Yes, Total - Yes) ~ ConditionOfInterest*Difference + 
                 (1| ID) +
                 (1 + Difference + ConditionOfInterest|StandardValues), 
               family = binomial(link = "probit"),
               data = Dataframe1)
Residuals22 = simulateResiduals(Model22)

Model23 = glmer(cbind(Yes, Total - Yes) ~ ConditionOfInterest*Difference + 
                 (1 + Difference| ID) +
                 (1 + Difference + ConditionOfInterest|StandardValues), 
               family = binomial(link = "probit"),
               data = Dataframe1)
Residuals23 = simulateResiduals(Model23)

Model24 = glmer(cbind(Yes, Total - Yes) ~ ConditionOfInterest*Difference + 
                 (1 + ConditionOfInterest| ID) +
                 (1 + Difference + ConditionOfInterest|StandardValues), 
               family = binomial(link = "probit"),
               data = Dataframe1)
Residuals24 = simulateResiduals(Model24)

Model25 = glmer(cbind(Yes, Total - Yes) ~ ConditionOfInterest*Difference + 
                 (1 + Difference + ConditionOfInterest| ID) +
                 (1 + Difference + ConditionOfInterest|StandardValues), 
               family = binomial(link = "probit"),
               data = Dataframe1)
Residuals25 = simulateResiduals(Model25)

data.frame(Label = c("0/0",	"1/0",	"Diff/0",	"Cond/0",	"Diff+Cond/0",
                     "0/1",	"1/1",	"Diff/1",	"Cond/1",	"Diff+Cond/1",
                     "0/Diff",	"1/Diff",	"Diff/Diff",	"Cond/Diff",	"Diff+Cond/Diff",
                     "0/Cond",	"1/Cond",	"Diff/Cond",	"Cond/Cond",	"Diff+Cond/Cond",
                     "0/Diff+Cond",	"1/Diff+Cond",	"Diff/Diff+Cond",	"Cond/Diff+Cond",	"Diff+Cond/Diff+Cond"),
           SD_ID_Intercept = c(2,sd(ranef(Model2)$ID[[1]]),sd(ranef(Model3)$ID[[1]]),sd(ranef(Model4)$ID[[1]]),
                               sd(ranef(Model5)$ID[[1]]),sd(ranef(Model6)$ID[[1]]),sd(ranef(Model7)$ID[[1]]),sd(ranef(Model8)$ID[[1]]),
                               sd(ranef(Model9)$ID[[1]]),sd(ranef(Model10)$ID[[1]]),sd(ranef(Model11)$ID[[1]]),sd(ranef(Model12)$ID[[1]]),
                               sd(ranef(Model13)$ID[[1]]),sd(ranef(Model14)$ID[[1]]),sd(ranef(Model15)$ID[[1]]),sd(ranef(Model16)$ID[[1]]),
                               sd(ranef(Model17)$ID[[1]]),sd(ranef(Model18)$ID[[1]]),sd(ranef(Model19)$ID[[1]]),sd(ranef(Model20)$ID[[1]]),
                               sd(ranef(Model21)$ID[[1]]),sd(ranef(Model22)$ID[[1]]),sd(ranef(Model23)$ID[[1]]),sd(ranef(Model24)$ID[[1]]),
                               sd(ranef(Model25)$ID[[1]])),
           SD_ID_Intercept = c(2,sd(ranef(Model2)$ID[[2]]),sd(ranef(Model3)$ID[[2]]),sd(ranef(Model4)$ID[[2]]),
                               sd(ranef(Model5)$ID[[2]]),sd(ranef(Model6)$ID[[2]]),sd(ranef(Model7)$ID[[2]]),sd(ranef(Model8)$ID[[2]]),
                               sd(ranef(Model9)$ID[[2]]),sd(ranef(Model10)$ID[[2]]),sd(ranef(Model11)$ID[[2]]),sd(ranef(Model12)$ID[[2]]),
                               sd(ranef(Model13)$ID[[2]]),sd(ranef(Model14)$ID[[2]]),sd(ranef(Model15)$ID[[2]]),sd(ranef(Model16)$ID[[2]]),
                               sd(ranef(Model17)$ID[[2]]),sd(ranef(Model18)$ID[[2]]),sd(ranef(Model19)$ID[[2]]),sd(ranef(Model20)$ID[[2]]),
                               sd(ranef(Model21)$ID[[2]]),sd(ranef(Model22)$ID[[2]]),sd(ranef(Model23)$ID[[2]]),sd(ranef(Model24)$ID[[2]]),
                               sd(ranef(Model25)$ID[[2]])),
           SD_ID_Intercept = c(2,sd(ranef(Model2)$ID[[3]]),sd(ranef(Model3)$ID[[3]]),sd(ranef(Model4)$ID[[3]]),
                               sd(ranef(Model5)$ID[[3]]),sd(ranef(Model6)$ID[[3]]),sd(ranef(Model7)$ID[[3]]),sd(ranef(Model8)$ID[[3]]),
                               sd(ranef(Model9)$ID[[3]]),sd(ranef(Model10)$ID[[3]]),sd(ranef(Model11)$ID[[3]]),sd(ranef(Model12)$ID[[3]]),
                               sd(ranef(Model13)$ID[[3]]),sd(ranef(Model14)$ID[[3]]),sd(ranef(Model15)$ID[[3]]),sd(ranef(Model16)$ID[[3]]),
                               sd(ranef(Model17)$ID[[3]]),sd(ranef(Model18)$ID[[3]]),sd(ranef(Model19)$ID[[3]]),sd(ranef(Model20)$ID[[3]]),
                               sd(ranef(Model21)$ID[[3]]),sd(ranef(Model22)$ID[[3]]),sd(ranef(Model23)$ID[[3]]),sd(ranef(Model24)$ID[[3]]),
                               sd(ranef(Model25)$ID[[3]])),
           SD_SV_Intercept = c(2,sd(ranef(Model2)$StandardValues[[1]]),sd(ranef(Model3)$StandardValues[[1]]),sd(ranef(Model4)$StandardValues[[1]]),
                               sd(ranef(Model5)$StandardValues[[1]]),sd(ranef(Model6)$StandardValues[[1]]),sd(ranef(Model7)$StandardValues[[1]]),sd(ranef(Model8)$StandardValues[[1]]),
                               sd(ranef(Model9)$StandardValues[[1]]),sd(ranef(Model10)$StandardValues[[1]]),sd(ranef(Model11)$StandardValues[[1]]),sd(ranef(Model12)$StandardValues[[1]]),
                               sd(ranef(Model13)$StandardValues[[1]]),sd(ranef(Model14)$StandardValues[[1]]),sd(ranef(Model15)$StandardValues[[1]]),sd(ranef(Model16)$StandardValues[[1]]),
                               sd(ranef(Model17)$StandardValues[[1]]),sd(ranef(Model18)$StandardValues[[1]]),sd(ranef(Model19)$StandardValues[[1]]),sd(ranef(Model20)$StandardValues[[1]]),
                               sd(ranef(Model21)$StandardValues[[1]]),sd(ranef(Model22)$StandardValues[[1]]),sd(ranef(Model23)$StandardValues[[1]]),sd(ranef(Model24)$StandardValues[[1]]),
                               sd(ranef(Model25)$StandardValues[[1]])),
           SD_SV_Intercept = c(2,sd(ranef(Model2)$StandardValues[[2]]),sd(ranef(Model3)$StandardValues[[2]]),sd(ranef(Model4)$StandardValues[[2]]),
                               sd(ranef(Model5)$StandardValues[[2]]),sd(ranef(Model6)$StandardValues[[2]]),sd(ranef(Model7)$StandardValues[[2]]),sd(ranef(Model8)$StandardValues[[2]]),
                               sd(ranef(Model9)$StandardValues[[2]]),sd(ranef(Model10)$StandardValues[[2]]),sd(ranef(Model11)$StandardValues[[2]]),sd(ranef(Model12)$StandardValues[[2]]),
                               sd(ranef(Model13)$StandardValues[[2]]),sd(ranef(Model14)$StandardValues[[2]]),sd(ranef(Model15)$StandardValues[[2]]),sd(ranef(Model16)$StandardValues[[2]]),
                               sd(ranef(Model17)$StandardValues[[2]]),sd(ranef(Model18)$StandardValues[[2]]),sd(ranef(Model19)$StandardValues[[2]]),sd(ranef(Model20)$StandardValues[[2]]),
                               sd(ranef(Model21)$StandardValues[[2]]),sd(ranef(Model22)$StandardValues[[2]]),sd(ranef(Model23)$StandardValues[[2]]),sd(ranef(Model24)$StandardValues[[2]]),
                               sd(ranef(Model25)$StandardValues[[2]])),
           SD_SV_Intercept = c(2,sd(ranef(Model2)$StandardValues[[3]]),sd(ranef(Model3)$StandardValues[[3]]),sd(ranef(Model4)$StandardValues[[3]]),
                               sd(ranef(Model5)$StandardValues[[3]]),sd(ranef(Model6)$StandardValues[[3]]),sd(ranef(Model7)$StandardValues[[3]]),sd(ranef(Model8)$StandardValues[[3]]),
                               sd(ranef(Model9)$StandardValues[[3]]),sd(ranef(Model10)$StandardValues[[3]]),sd(ranef(Model11)$StandardValues[[3]]),sd(ranef(Model12)$StandardValues[[3]]),
                               sd(ranef(Model13)$StandardValues[[3]]),sd(ranef(Model14)$StandardValues[[3]]),sd(ranef(Model15)$StandardValues[[3]]),sd(ranef(Model16)$StandardValues[[3]]),
                               sd(ranef(Model17)$StandardValues[[3]]),sd(ranef(Model18)$StandardValues[[3]]),sd(ranef(Model19)$StandardValues[[3]]),sd(ranef(Model20)$StandardValues[[3]]),
                               sd(ranef(Model21)$StandardValues[[3]]),sd(ranef(Model22)$StandardValues[[3]]),sd(ranef(Model23)$StandardValues[[3]]),sd(ranef(Model24)$StandardValues[[3]]),
                               sd(ranef(Model25)$StandardValues[[3]])),
           
           
           
           
           )
sd(ranef(Model2)$ID[[1]])
sd(ranef(Model5)$ID[[1]])

jpeg(paste0("PlotResiduals",i,".jpeg"), w = 2000, h = 2000)
par(mfrow=c(5,5))
plotResiduals(Residuals1, Dataframe1$Difference, quantreg = T)
plotResiduals(Residuals2, Dataframe1$Difference, quantreg = T)
plotResiduals(Residuals3, Dataframe1$Difference, quantreg = T)
plotResiduals(Residuals4, Dataframe1$Difference, quantreg = T)
plotResiduals(Residuals5, Dataframe1$Difference, quantreg = T)
plotResiduals(Residuals6, Dataframe1$Difference, quantreg = T)
plotResiduals(Residuals7, Dataframe1$Difference, quantreg = T)
plotResiduals(Residuals8, Dataframe1$Difference, quantreg = T)
plotResiduals(Residuals9, Dataframe1$Difference, quantreg = T)
plotResiduals(Residuals10, Dataframe1$Difference, quantreg = T)
plotResiduals(Residuals11, Dataframe1$Difference, quantreg = T)
plotResiduals(Residuals12, Dataframe1$Difference, quantreg = T)
plotResiduals(Residuals13, Dataframe1$Difference, quantreg = T)
plotResiduals(Residuals14, Dataframe1$Difference, quantreg = T)
plotResiduals(Residuals15, Dataframe1$Difference, quantreg = T)
plotResiduals(Residuals16, Dataframe1$Difference, quantreg = T)
plotResiduals(Residuals17, Dataframe1$Difference, quantreg = T)
plotResiduals(Residuals18, Dataframe1$Difference, quantreg = T)
plotResiduals(Residuals19, Dataframe1$Difference, quantreg = T)
plotResiduals(Residuals20, Dataframe1$Difference, quantreg = T)
plotResiduals(Residuals21, Dataframe1$Difference, quantreg = T)
plotResiduals(Residuals22, Dataframe1$Difference, quantreg = T)
plotResiduals(Residuals23, Dataframe1$Difference, quantreg = T)
plotResiduals(Residuals24, Dataframe1$Difference, quantreg = T)
plotResiduals(Residuals25, Dataframe1$Difference, quantreg = T)
dev.off()

Sys.time()-start_time
}
