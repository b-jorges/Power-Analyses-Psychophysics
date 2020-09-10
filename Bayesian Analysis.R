BayesianAnalysis = brm(bf(Yes ~ ConditionOfInterest*Difference + (Difference + ConditionOfInterest | ID) + (Difference + ConditionOfInterest | StandardValues)),
                       data = Psychometric, 
                       family = bernoulli())

BayesianAnalysis2 = brm(bf(Yes ~ ConditionOfInterest*Difference + (Difference | ID) + (Difference | StandardValues)),
                        data = Psychometric, 
                        family = bernoulli())

BayesianAnalysis3 = brm(bf(Yes ~ ConditionOfInterest*Difference + (1 | ID) + (1 | StandardValues)),
                        data = Psychometric, 
                        family = bernoulli())

GLMM_SaveThisOne = GLMM2
summary(GLMM)
coef(GLMM)

summary(BayesianAnalysis)
coef(BayesianAnalysis)
BayesianAnalysis
BayesianAnalysis2
BayesianAnalysis3