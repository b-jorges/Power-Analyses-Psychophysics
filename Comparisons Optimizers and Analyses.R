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

Dataframe_pvalues1 = rbind(read.csv(header = T, file = paste0(Where_Am_I(),"/Data/ComparisonMethodsSmallEffectBigModel.csv")))
Dataframe_pvalues1$Effect = "Small Effect"

Dataframe_pvalues1 = rbind(read.csv(header = T, file = paste0(Where_Am_I(),"/Data/ComparisonMethodsNoEffectBigModel.csv")))
Dataframe_pvalues1$Effect = "No Effect"





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