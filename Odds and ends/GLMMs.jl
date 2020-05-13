using Pkg
Pkg.activate()
Pkg.instantiate()
using MixedModels
using RCall
using DataFrames, Tables
using Random
using CSV
using RData
using CategoricalArrays
using Statistics
using Dates

R"""
require(dplyr, quietly = TRUE)   # for data wrangling
require(tidyverse, quietly = TRUE)   # for data wrangling
require(lme4)
require(lmerTest)
require(quickpsy)

SimulatePsychometricFunction_Staircase = function(ID, 
                                                  ConditionOfInterest, 
                                                  StandardValues, 
                                                  reps, 
                                                  PSE_Difference, 
                                                  JND_Difference, 
                                                  Multiplicator_PSE_Standard, 
                                                  Multiplicator_SD_Standard, 
                                                  SD_ResponseFunction, 
                                                  Mean_Variability_Between = 0.1, 
                                                  SD_Variability_Between = 0.1){
  Psychometric = expand.grid(ID=ID, ConditionOfInterest=ConditionOfInterest, StandardValues=StandardValues, reps = reps)
  
  Psychometric = Psychometric %>%
    group_by(ID) %>%#
    mutate(PSE_Factor_ID = rnorm(1,1,Mean_Variability_Between),
           SD_Factor_ID = rnorm(1,1,SD_Variability_Between))
  
  Psychometric = Psychometric %>%
    mutate(
      Mean_Standard = StandardValues+StandardValues*Multiplicator_PSE_Standard,
      SD_Standard = StandardValues*Multiplicator_SD_Standard,
      Mean = (Mean_Standard + (ConditionOfInterest==ConditionOfInterest[2])*StandardValues*PSE_Difference)*PSE_Factor_ID,
      SD = abs((SD_Standard + (ConditionOfInterest==ConditionOfInterest[2])*SD_Standard*JND_Difference)*SD_Factor_ID),
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
""";

function SimulateDataframe(n,
        ConditionOfInterest,
        StandardValues,
        reps,
        PSE_Difference,
        JND_Difference,
        Multiplicator_PSE_Standard,
        Multiplicator_SD_Standard,
        SD_ResponseFunction,
        Mean_Variability_Between,
        SD_Variability_Between)
    
    @rput n ConditionOfInterest StandardValues reps PSE_Difference JND_Difference Multiplicator_PSE_Standard Multiplicator_SD_Standard SD_ResponseFunction Mean_Variability_Between SD_Variability_Between

    R"""
    ID = paste0("s",1:n)
        Psychometric = SimulatePsychometricFunction_Staircase(ID,
            ConditionOfInterest,
            StandardValues,
            1:reps,
            PSE_Difference,
            JND_Difference,
            Multiplicator_PSE_Standard,
            Multiplicator_SD_Standard,
            SD_ResponseFunction,
            Mean_Variability_Between,
            SD_Variability_Between)
    
        Parameters = quickpsy(Psychometric,Difference,Answer,grouping = .(ID,ConditionOfInterest,StandardValues), bootstrap = "none")$par
        Parameters2 = Parameters %>%
        filter(parn == "p1") %>%
        select(ID,ConditionOfInterest,Mean=par, StandardValues)
        Parameters2$SD = Parameters$par[Parameters$parn == "p2"]
        FittedPsychometricFunctions = Parameters2
    
        FittedPsychometricFunctions$StandardValues = as.character(FittedPsychometricFunctions$StandardValues)
        Psychometric$StandardValues = as.character(Psychometric$StandardValues)
    
    """
    @rget Psychometric FittedPsychometricFunctions
    
    Psychometric[:StandardValuesAsFactor] = "placeholder"
    
    formula1 = @formula(Answer ~ Difference*ConditionOfInterest + (Difference + ConditionOfInterest |ID) + (Difference + ConditionOfInterest|StandardValues));
    GLMM = fit!(GeneralizedLinearMixedModel(formula1, Psychometric, Binomial()), fast=true)
    
    formula2 = @formula(Mean ~ ConditionOfInterest + (1|ID) + (1|StandardValues));
    TwoLevelMean = fit(MixedModel,formula2, FittedPsychometricFunctions)
    
    formula3 = @formula(SD ~ ConditionOfInterest + (1|ID) + (1|StandardValues));
    TwoLevelSD = fit(MixedModel,formula3, FittedPsychometricFunctions)

    [(coeftable(GLMM)).cols[4][2];(coeftable(GLMM)).cols[4][4];(coeftable(TwoLevelMean)).cols[4][2];(coeftable(TwoLevelSD)).cols[4][2]]
    
end

ConditionOfInterest = [0;1]
StandardValues = [5;8]
Range_reps = [40,60]
Range_PSE_Difference = [-0.1,-0.05,0,0.05,0.1]
Range_JND_Difference = [-0.2,-0.1,0,0.1,0.2]
Multiplicator_PSE_Standard = 0
Multiplicator_SD_Standard = 0.108
SD_ResponseFunction = 0.1
Mean_Variability_Between = 0.1
SD_Variability_Between = 0.1
nIterations = 100
Range_Participants = [10,12,14,16,18,20]

TotalNumber = length(Range_reps)*length(Range_PSE_Difference)*length(Range_JND_Difference)*length(Range_Participants)
CurrentRunthrough = 0
rightnow = Dates.now()

for reps in Range_reps
    for PSE_Difference in Range_PSE_Difference
        for JND_Difference in Range_JND_Difference
            for n in Range_Participants
                
                TimeStartTrial = Dates.now()
                
                Pvalues_Accuracy = []
                Pvalues_Precision = []
                Pvalues_Accuracy_TwoLevel = []
                Pvalues_Precision_TwoLevel = []
                
                for j in 1:nIterations
                Pvalues = SimulateDataframe(n, 
                                      ConditionOfInterest, 
                                      StandardValues, 
                                      reps, 
                                      PSE_Difference, 
                                      JND_Difference, 
                                      Multiplicator_PSE_Standard, 
                                      Multiplicator_SD_Standard, 
                                      SD_ResponseFunction, 
                                      Mean_Variability_Between, 
                                      SD_Variability_Between)
                    Pvalues_Accuracy = [Pvalues_Accuracy;Pvalues[1]]
                    Pvalues_Precision = [Pvalues_Precision;Pvalues[2]]
                    Pvalues_Accuracy_TwoLevel = [Pvalues_Accuracy_TwoLevel;Pvalues[3]]
                    Pvalues_Precision_TwoLevel = [Pvalues_Precision_TwoLevel;Pvalues[4]]
                end
                
                CurrentRunthrough = CurrentRunthrough + 1

                if CurrentRunthrough == 1

                   global PowerfulDataframe = DataFrame(n=n, 
                        ConditionsOfInterest=length(ConditionOfInterest), 
                        StandardValue1=StandardValues[1],
                        StandardValue2=StandardValues[2], reps=reps, 
                        PSE_Difference=PSE_Difference, 
                        JND_Difference=JND_Difference, 
                        Multiplicator_PSE_Standard=Multiplicator_PSE_Standard, 
                        Multiplicator_SD_Standard=Multiplicator_SD_Standard, 
                        SD_ResponseFunction=SD_ResponseFunction, 
                        Mean_Variability_Between=Mean_Variability_Between, 
                        SD_Variability_Between=SD_Variability_Between, 
                        power_Accuracy = mean(Pvalues_Accuracy .< 0.05),  
                        power_Precision = mean(Pvalues_Precision .< 0.05),
                        power_Accuracy_Twolevel = mean(Pvalues_Accuracy_TwoLevel .< 0.05),  
                        power_Precision_Twolevel = mean(Pvalues_Precision_TwoLevel .< 0.05),
                        Duration = ((Dates.now()) - TimeStartTrial))   

                else
                    row = DataFrame(n=n, 
                        ConditionsOfInterest=length(ConditionOfInterest), 
                        StandardValue1=StandardValues[1],StandardValue2=StandardValues[2], 
                        reps=reps, 
                        PSE_Difference=PSE_Difference, 
                        JND_Difference=JND_Difference, 
                        Multiplicator_PSE_Standard=Multiplicator_PSE_Standard, 
                        Multiplicator_SD_Standard=Multiplicator_SD_Standard, 
                        SD_ResponseFunction=SD_ResponseFunction, 
                        Mean_Variability_Between=Mean_Variability_Between, 
                        SD_Variability_Between=SD_Variability_Between, 
                        power_Accuracy = mean(Pvalues_Accuracy .< 0.05),  
                        power_Precision = mean(Pvalues_Precision .< 0.05),
                        power_Accuracy_Twolevel = mean(Pvalues_Accuracy_TwoLevel .< 0.05),  
                        power_Precision_Twolevel = mean(Pvalues_Precision_TwoLevel .< 0.05),
                        Duration=((Dates.now()) - TimeStartTrial))
                    
                    PowerfulDataframe = append!(PowerfulDataframe,row)
                end
                
                print("RUNTHROUGH ", CurrentRunthrough, " out of ", TotalNumber,": ", n, " ", reps, " ", PSE_Difference, " ", JND_Difference, " ", mean(Pvalues_Accuracy .< 0.05), " ", 
                    mean(Pvalues_Precision .< 0.05), " ", PowerfulDataframe[!,:Duration][CurrentRunthrough], " END. ")

            end
            CSV.write(join([reps,"_", PSE_Difference, "_", JND_Difference, ".csv"]),PowerfulDataframe)
        end
    end
end

