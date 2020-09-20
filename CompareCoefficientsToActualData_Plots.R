
Dataframe2 = rbind(read.csv(header = T, file = paste0(Where_Am_I(),"/Data/ComparisonMethodsSmallEffectBigModel.csv")))

#####MAIN PLOTS
#model predictions 0.1 0.125
FigurePSEsFromGLMM = ggplot(Dataframe2 %>% filter(PSE_Difference == 0.05 & JND_Difference == 0.125),aes(Model,Mean_Modeled-Mean_Actual)) +
  geom_boxplot() + 
  ylab("PSE Output of GLMM - Actual PSE") +
  geom_hline(aes(yintercept = 0), linetype = 2, size = 1) + 
  coord_cartesian(ylim = c(-3,3)) +
  xlab("") +
  ggtitle("PSE coefficient (PSE Difference = 0.05/JND Difference = 0.125)")

FigureSDsFromGLMM = ggplot(Dataframe2 %>% filter(PSE_Difference == 0.05 & JND_Difference == 0.125),aes(Model,SD_Modeled-SD_Actual)) +
  geom_boxplot() + 
  coord_cartesian(ylim = c(-1,2.5)) +
  geom_hline(aes(yintercept = 0), linetype = 2, size = 1) +
  ylab("SD Output of GLMM - Actual SD") +
  xlab("") +
  ggtitle("JND coefficient (PSE Difference = 0.05/JND Difference = 0.125)")
FigureValuesFromGLMMs = plot_grid(FigurePSEsFromGLMM,FigureSDsFromGLMM, nrow = 2, labels = "AUTO")
ggsave("Figures/FigureValuesFromGLMMs_0.05_0.125.jpeg", w = 12, h = 7.2)
#####

###SE
FigureSEsPSEsFromGLMM = ggplot(Dataframe2 %>% filter(PSE_Difference == 0.05 & JND_Difference == 0.125),
                               aes(Model,SECoI)) +
  geom_boxplot() + 
  ylab("SE for PSE estimate") +
  xlab("") +
  coord_cartesian(ylim = c(0,0.2)) +
  ggtitle("SE of PSE coefficient (PSE Difference = 0.05/JND Difference = 0.125)")
FigureSEsSDsFromGLMM = ggplot(Dataframe2 %>% filter(PSE_Difference == 0.05 & JND_Difference == 0.125),
                              aes(Model,SEInterac)) +
  geom_boxplot() + 
  ylab("SE for SD estimate") +
  xlab("") +
  coord_cartesian(ylim = c(0,0.2)) +
  ggtitle("SE of JND coefficient (PSE Difference = 0.05/JND Difference = 0.125)")
FigureSEsFromGLMMs = plot_grid(FigureSEsPSEsFromGLMM,FigureSEsSDsFromGLMM, nrow = 2, labels = "AUTO")
ggsave(paste0("Figures/FigureSEsFromGLMMs_0.05_0.125.jpeg"), w = 12, h = 7.2)


#####
#p values
Dataframe3 = Dataframe2 %>%
  filter(PSE_Difference == 0.05 & JND_Difference == 0.125) %>% 
  group_by(Model,Repetition) %>% 
  slice(1) %>% 
  group_by(Model) %>% 
  mutate(Power_CoI = sum(PvaluesCoI < 0.05)/length(PvaluesCoI),
         Power_Interac = sum(PvaluesInterac < 0.05)/length(PvaluesInterac))

FigurePvaluesFromGLMM_CoI = ggplot(Dataframe3,aes(Model,Power_CoI)) +
  geom_point(size = 5) + 
  ylab("Power") +
  xlab("") +
  geom_hline(yintercept = 0.05, linetype = 3) +
  geom_hline(yintercept = 0.9, linetype = 2) +
  coord_cartesian(ylim = c(0,1)) +
  ggtitle("Power PSE coefficient (PSE Difference = 0.05/JND Difference = 0.125)")
FigurePvaluesFromGLMM_Interac = ggplot(Dataframe3,aes(Model,Power_Interac)) +
  geom_point(size = 5) + 
  ylab("Power") +
  xlab("") +
  geom_hline(yintercept = 0.05, linetype = 3) +
  geom_hline(yintercept = 0.9, linetype = 2) +
  coord_cartesian(ylim = c(0,1)) +
  ggtitle("Power JND coefficient (PSE Difference = 0.05/JND Difference = 0.125)")
FigureValuesFromGLMMs = plot_grid(FigurePvaluesFromGLMM_CoI,FigurePvaluesFromGLMM_Interac, nrow = 2, labels = "AUTO")
ggsave(paste0("Figures/FigurePValuesFromGLMMs_0.05_0.125.jpeg"), w = 12, h = 7.2)
#####

#####
#AICs
ggplot(Dataframe2 %>% filter(PSE_Difference == 0.05 & JND_Difference == 0.125),aes(Model,AIC_Norm)) +
  geom_boxplot() + 
  xlab("") +
  ylab("AIC - AIC(Model25)") +
  geom_hline(aes(yintercept = 0), linetype = 2, size = 1) +
  ggtitle("AIC (PSE Difference = 0.05/JND Difference = 0.125)")
ggsave(paste0("Figures/Figure AICs_0.05_0.125_Models.jpeg"), w = 12, h = 3.6)
#####



###################################################################################################################
###################################SUPPLEMENTARY PLOTS#############################################################
###################################################################################################################
###################################PSE = -0.05, JND = -0.125 _M0.05_M0.125
#model predictions
FigurePSEsFromGLMM_M0.05_M0.125 = ggplot(Dataframe2 %>% filter(PSE_Difference == -0.05 & JND_Difference == -0.125),aes(Model,Mean_Modeled-Mean_Actual)) +
  geom_boxplot() + 
  ylab("PSE Output of GLMM - Actual PSE") +
  geom_hline(aes(yintercept = 0), linetype = 2, size = 1) + 
  coord_cartesian(ylim = c(-3,3)) +
  xlab("") +
  ggtitle("PSE coefficient (PSE Difference = -0.05/JND Difference = -0.125)")

FigureSDsFromGLMM_M0.05_M0.125 = ggplot(Dataframe2 %>% filter(PSE_Difference == -0.05 & JND_Difference == -0.125),aes(Model,SD_Modeled-SD_Actual)) +
  geom_boxplot() + 
  coord_cartesian(ylim = c(-1,2.5)) +
  geom_hline(aes(yintercept = 0), linetype = 2, size = 1) +
  ylab("SD Output of GLMM - Actual SD") +
  xlab("") +
  ggtitle("JND coefficient (PSE Difference = -0.05/JND Difference = -0.125)")
FigureValuesFromGLMMs_M0.05_M0.125 = plot_grid(FigurePSEsFromGLMM_M0.05_M0.125,FigureSDsFromGLMM_M0.05_M0.125, nrow = 2, labels = "AUTO")
ggsave(paste0("Figures/FigureValuesFromGLMMs_-0.05_-0.125.jpeg"), w = 12, h = 7.2)
#####

###SE
FigureSEsPSEsFromGLMM_M0.05_M0.125 = ggplot(Dataframe2 %>% filter(PSE_Difference == -0.05 & JND_Difference == -0.125),
                                            aes(Model,SECoI)) +
  geom_boxplot() + 
  ylab("SE for PSE estimate") +
  xlab("") +
  coord_cartesian(ylim = c(0,0.2)) +
  ggtitle("SE of PSE coefficient (PSE Difference = -0.05/JND Difference = -0.125)")
FigureSEsSDsFromGLMM_M0.05_M0.125 = ggplot(Dataframe2 %>% filter(PSE_Difference == -0.05 & JND_Difference == -0.125),
                                           aes(Model,SEInterac)) +
  geom_boxplot() + 
  ylab("SE for SD estimate") +
  xlab("") +
  coord_cartesian(ylim = c(0,0.2)) +
  ggtitle("SE of JND coefficient (PSE Difference = -0.05/JND Difference = -0.125)")
FigureSEsFromGLMMs_M0.05_M0.125 = plot_grid(FigureSEsPSEsFromGLMM_M0.05_M0.125,FigureSEsSDsFromGLMM_M0.05_M0.125, nrow = 2, labels = "AUTO")
ggsave(paste0("Figures/FigureSEsFromGLMMs_-0.05_-0.125.jpeg"), w = 12, h = 7.2)


#####
#p values
Dataframe3_M0.05_M0.125 = Dataframe2 %>%
  filter(PSE_Difference == -0.05 & JND_Difference == -0.125) %>% 
  group_by(Model,Repetition) %>% 
  slice(1) %>% 
  group_by(Model) %>% 
  mutate(Power_CoI = sum(PvaluesCoI < 0.05)/length(PvaluesCoI),
         Power_Interac = sum(PvaluesInterac < 0.05)/length(PvaluesInterac))

FigurePvaluesFromGLMM_CoI_M0.05_M0.125 = ggplot(Dataframe3_M0.05_M0.125,aes(Model,Power_CoI)) +
  geom_point(size = 5) + 
  ylab("Power") +
  xlab("") +
  geom_hline(yintercept = 0.05, linetype = 3) +
  geom_hline(yintercept = 0.9, linetype = 2) +
  coord_cartesian(ylim = c(0,1)) +
  ggtitle("Power PSE coefficient (PSE Difference = -0.05/JND Difference = -0.125)")
FigurePvaluesFromGLMM_Interac_M0.05_M0.125 = ggplot(Dataframe3_M0.05_M0.125,aes(Model,Power_Interac)) +
  geom_point(size = 5) + 
  ylab("Power") +
  xlab("") +
  geom_hline(yintercept = 0.05, linetype = 3) +
  geom_hline(yintercept = 0.9, linetype = 2) +
  coord_cartesian(ylim = c(0,1)) +
  ggtitle("Power JND coefficient (PSE Difference = -0.05/JND Difference = -0.125)")
FigureValuesFromGLMMs_M0.05_M0.125 = plot_grid(FigurePvaluesFromGLMM_CoI_M0.05_M0.125,FigurePvaluesFromGLMM_Interac_M0.05_M0.125, nrow = 2, labels = "AUTO")
ggsave(paste0("Figures/FigurePValuesFromGLMMs_-0.05_-0.125.jpeg"), w = 12, h = 7.2)
#####

#####
#AICs
ggplot(Dataframe2 %>% filter(PSE_Difference == -0.05 & JND_Difference == -0.125),aes(Model,AIC_Norm)) +
  geom_boxplot() + 
  xlab("") +
  ylab("AIC - AIC(Model25)") +
  geom_hline(aes(yintercept = 0), linetype = 2, size = 1) +
  ggtitle("AIC (PSE Difference = -0.05/JND Difference = -0.125)")
ggsave(paste0("Figures/Figure AICs_-0.05_-0.125_Models.jpeg"), w = 12, h = 3.6)
#####


###################################PSE = 0, JND = 0.125
FigurePSEsFromGLMM_0_0.125_ = ggplot(Dataframe2 %>% filter(PSE_Difference == 0 & JND_Difference == 0.125),aes(Model,Mean_Modeled-Mean_Actual)) +
  geom_boxplot() + 
  ylab("PSE Output of GLMM - Actual PSE") +
  geom_hline(aes(yintercept = 0), linetype = 2, size = 1) + 
  coord_cartesian(ylim = c(-3,3)) +
  xlab("") +
  ggtitle("PSE coefficient (PSE Difference = 0/JND Difference = 0.125)")

FigureSDsFromGLMM_0_0.125_ = ggplot(Dataframe2 %>% filter(PSE_Difference == 0 & JND_Difference == 0.125),aes(Model,SD_Modeled-SD_Actual)) +
  geom_boxplot() + 
  coord_cartesian(ylim = c(-1,2.5)) +
  geom_hline(aes(yintercept = 0), linetype = 2, size = 1) +
  ylab("SD Output of GLMM - Actual SD") +
  xlab("") +
  ggtitle("JND coefficient (PSE Difference = 0/JND Difference = 0.125)")
plot_grid(FigurePSEsFromGLMM_0_0.125_,FigureSDsFromGLMM_0_0.125_, nrow = 2, labels = "AUTO")
ggsave(paste0("Figures/FigureValuesFromGLMMs_0_0.125.jpeg"), w = 12, h = 7.2)
#####

###SE
FigureSEsPSEsFromGLMM_0_0.125_ = ggplot(Dataframe2 %>% filter(PSE_Difference == 0 & JND_Difference == 0.125),
                                        aes(Model,SECoI)) +
  geom_boxplot() + 
  ylab("SE for PSE estimate") +
  xlab("") +
  coord_cartesian(ylim = c(0,0.2)) +
  ggtitle("SE of PSE coefficient (PSE Difference = 0/JND Difference = 0.125)")
FigureSEsSDsFromGLMM_0_0.125_ = ggplot(Dataframe2 %>% filter(PSE_Difference == 0 & JND_Difference == 0.125),
                                       aes(Model,SEInterac)) +
  geom_boxplot() + 
  ylab("SE for SD estimate") +
  xlab("") +
  coord_cartesian(ylim = c(0,0.2)) +
  ggtitle("SE of JND coefficient (PSE Difference = 0/JND Difference = 0.125)")
FigureSEsFromGLMMs_0_0.125_ = plot_grid(FigureSEsPSEsFromGLMM_0_0.125_,FigureSEsSDsFromGLMM_0_0.125_, nrow = 2, labels = "AUTO")
ggsave(paste0("Figures/FigureSEsFromGLMMs_0_0.125.jpeg"), w = 12, h = 7.2)

#####
#p values
Dataframe4 = Dataframe2 %>%
  filter(PSE_Difference == 0 & JND_Difference == 0.125) %>% 
  group_by(Model,Repetition) %>% 
  slice(1) %>% 
  group_by(Model) %>% 
  mutate(Power_CoI = sum(PvaluesCoI < 0.05)/length(PvaluesCoI),
         Power_Interac = sum(PvaluesInterac < 0.05)/length(PvaluesInterac))

FigurePvaluesFromGLMM_CoI_0_0.125_ = ggplot(Dataframe4,aes(Model,Power_CoI)) +
  geom_point(size = 5) + 
  ylab("Power") +
  xlab("") +
  geom_hline(yintercept = 0.05, linetype = 3) +
  geom_hline(yintercept = 0.9, linetype = 2) +
  coord_cartesian(ylim = c(0,1)) +
  ggtitle("Power PSE coefficient (PSE Difference = 0/JND Difference = 0.125)")
FigurePvaluesFromGLMM_Interac_0_0.125_ = ggplot(Dataframe4,aes(Model,Power_Interac)) +
  geom_point(size = 5) + 
  ylab("Power") +
  xlab("") +
  geom_hline(yintercept = 0.05, linetype = 3) +
  geom_hline(yintercept = 0.9, linetype = 2) +
  coord_cartesian(ylim = c(0,1)) +
  ggtitle("Power JND coefficient (PSE Difference = 0/JND Difference = 0.125)")
plot_grid(FigurePvaluesFromGLMM_CoI_0_0.125_,FigurePvaluesFromGLMM_Interac_0_0.125_, nrow = 2, labels = "AUTO")
ggsave(paste0("Figures/FigurePValuesFromGLMMs_0_0.125.jpeg"), w = 12, h = 7.2)
#####

#####
#AICs
ggplot(Dataframe2 %>% filter(PSE_Difference == 0 & JND_Difference == 0.125),aes(Model,AIC_Norm)) +
  geom_boxplot() + 
  xlab("") +
  ylab("AIC - AIC(Model25)") +
  geom_hline(aes(yintercept = 0), linetype = 2, size = 1) +
  ggtitle("AIC (PSE Difference = 0/JND Difference = 0.125)")
ggsave(paste0("Figures/Figure AICs_0-0.125_Models.jpeg"), w = 12, h = 3.6)
#####


###################################################################################################################
###################################PSE = -0.05, JND = 0.125 _M0.05_0.125
FigurePSEsFromGLMM_M0.05_0.125 = ggplot(Dataframe2 %>% filter(PSE_Difference == -0.05 & JND_Difference == 0.125),aes(Model,Mean_Modeled-Mean_Actual)) +
  geom_boxplot() + 
  ylab("PSE Output of GLMM - Actual PSE") +
  geom_hline(aes(yintercept = 0), linetype = 2, size = 1) + 
  coord_cartesian(ylim = c(-3,3)) +
  xlab("") +
  ggtitle("PSE coefficient (PSE Difference = -0.05/JND Difference = 0.125)")
FigureSDsFromGLMM_M0.05_0.125 = ggplot(Dataframe2 %>% filter(PSE_Difference == -0.05 & JND_Difference == 0.125),aes(Model,SD_Modeled-SD_Actual)) +
  geom_boxplot() + 
  coord_cartesian(ylim = c(-1,2.5)) +
  geom_hline(aes(yintercept = 0), linetype = 2, size = 1) +
  ylab("SD Output of GLMM - Actual SD") +
  xlab("") +
  ggtitle("JND coefficient (PSE Difference = -0.05/JND Difference = 0.125)")
plot_grid(FigurePSEsFromGLMM_M0.05_0.125,FigureSDsFromGLMM_M0.05_0.125, nrow = 2, labels = "AUTO")
ggsave(paste0("Figures/FigureValuesFromGLMMs_-0.05_0.125.jpeg"), w = 12, h = 7.2)
#####


###SE
FigureSEsPSEsFromGLMM_M0.05_0.125 = ggplot(Dataframe2 %>% filter(PSE_Difference == -0.05 & JND_Difference == 0.125),
                                           aes(Model,SECoI)) +
  geom_boxplot() + 
  ylab("SE for PSE estimate") +
  xlab("") +
  coord_cartesian(ylim = c(0,0.2)) +
  ggtitle("SE of PSE coefficient (PSE Difference = -0.05/JND Difference = 0.125)")
FigureSEsSDsFromGLMM_M0.05_0.125 = ggplot(Dataframe2 %>% filter(PSE_Difference == -0.05 & JND_Difference == 0.125),
                                          aes(Model,SEInterac)) +
  geom_boxplot() + 
  ylab("SE for SD estimate") +
  xlab("") +
  coord_cartesian(ylim = c(0,0.2)) +
  ggtitle("SE of JND coefficient (PSE Difference = -0.05/JND Difference = 0.125)")
FigureSEsFromGLMMs_M0.05_0.125 = plot_grid(FigureSEsPSEsFromGLMM_M0.05_0.125,FigureSEsSDsFromGLMM_M0.05_0.125, nrow = 2, labels = "AUTO")
ggsave(paste0("Figures/FigureSEsFromGLMMs_-0.05_0.125.jpeg"), w = 12, h = 7.2)



#####
#p values
Dataframe5 = Dataframe2 %>%
  filter(PSE_Difference == -0.05 & JND_Difference == 0.125) %>% 
  group_by(Model,Repetition) %>% 
  slice(1) %>% 
  group_by(Model) %>% 
  mutate(Power_CoI = sum(PvaluesCoI < 0.05)/length(PvaluesCoI),
         Power_Interac = sum(PvaluesInterac < 0.05)/length(PvaluesInterac))

FigurePvaluesFromGLMM_CoI_M0.05_0.125 = ggplot(Dataframe5,aes(Model,Power_CoI)) +
  geom_point(size = 5) + 
  ylab("Power") +
  xlab("") +
  geom_hline(yintercept = 0.05, linetype = 3) +
  geom_hline(yintercept = 0.9, linetype = 2) +
  coord_cartesian(ylim = c(0,1)) +
  ggtitle("Power PSE coefficient (PSE Difference = -0.05/JND Difference = 0.125)")
FigurePvaluesFromGLMM_Interac_M0.05_0.125 = ggplot(Dataframe5,aes(Model,Power_Interac)) +
  geom_point(size = 5) + 
  ylab("Power") +
  xlab("") +
  geom_hline(yintercept = 0.05, linetype = 3) +
  geom_hline(yintercept = 0.9, linetype = 2) +
  coord_cartesian(ylim = c(0,1)) +
  ggtitle("Power JND coefficient (PSE Difference = -0.05/JND Difference = 0.125)")
plot_grid(FigurePvaluesFromGLMM_CoI_M0.05_0.125,FigurePvaluesFromGLMM_Interac_M0.05_0.125, nrow = 2, labels = "AUTO")
ggsave(paste0("Figures/FigurePValuesFromGLMMs_-0.05_0.125.jpeg"), w = 12, h = 7.2)
#####

#####
#AICs
ggplot(Dataframe2 %>% filter(PSE_Difference == -0.05 & JND_Difference == 0.125),aes(Model,AIC_Norm)) +
  geom_boxplot() + 
  xlab("") +
  ylab("AIC - AIC(Model25)") +
  geom_hline(aes(yintercept = 0), linetype = 2, size = 1) +
  ggtitle("AIC (PSE Difference = -0.05/JND Difference = 0.125)")
ggsave(paste0("Figures/Figure AICs_-0.05_0.125_Models.jpeg"), w = 12, h = 3.6)
#####


###################################################################################################################
###################################PSE = 0.05, JND = 0 _0.05_0
FigurePSEsFromGLMM_0.05_0 = ggplot(Dataframe2 %>% filter(PSE_Difference == 0.05 & JND_Difference == 0),aes(Model,Mean_Modeled-Mean_Actual)) +
  geom_boxplot() + 
  ylab("PSE Output of GLMM - Actual PSE") +
  geom_hline(aes(yintercept = 0), linetype = 2, size = 1) + 
  coord_cartesian(ylim = c(-3,3)) +
  xlab("") +
  ggtitle("PSE coefficient (PSE Difference = 0.05/JND Difference = 0)")
FigureSDsFromGLMM_0.05_0 = ggplot(Dataframe2 %>% filter(PSE_Difference == 0.05 & JND_Difference == 0),aes(Model,SD_Modeled-SD_Actual)) +
  geom_boxplot() + 
  coord_cartesian(ylim = c(-1,2.5)) +
  geom_hline(aes(yintercept = 0), linetype = 2, size = 1) +
  ylab("SD Output of GLMM - Actual SD") +
  xlab("") +
  ggtitle("JND coefficient (PSE Difference = 0.05/JND Difference = 0)")
plot_grid(FigurePSEsFromGLMM_0.05_0,FigureSDsFromGLMM_0.05_0, nrow = 2, labels = "AUTO")
ggsave(paste0("Figures/FigureValuesFromGLMMs_0.05_0.jpeg"), w = 12, h = 7.2)
#####

###SE
FigureSEsPSEsFromGLMM_0.05_0 = ggplot(Dataframe2 %>% filter(PSE_Difference == 0.05 & JND_Difference == 0),
                                      aes(Model,SECoI)) +
  geom_boxplot() + 
  ylab("SE for PSE estimate") +
  xlab("") +
  coord_cartesian(ylim = c(0,0.2)) +
  ggtitle("SE of PSE coefficient (PSE Difference = 0.05/JND Difference = 0)")
FigureSEsSDsFromGLMM_0.05_0 = ggplot(Dataframe2 %>% filter(PSE_Difference == 0.05 & JND_Difference == 0),
                                     aes(Model,SEInterac)) +
  geom_boxplot() + 
  ylab("SE for SD estimate") +
  xlab("") +
  coord_cartesian(ylim = c(0,0.2)) +
  ggtitle("SE of JND coefficient (PSE Difference = 0.05/JND Difference = 0)")
FigureSEsFromGLMMs_0.05_0 = plot_grid(FigureSEsPSEsFromGLMM_0.05_0,FigureSEsSDsFromGLMM_0.05_0, nrow = 2, labels = "AUTO")
ggsave(paste0("Figures/FigureSEsFromGLMMs_0.05_0.jpeg"), w = 12, h = 7.2)


#####
#p values
Dataframe6 = Dataframe2 %>%
  filter(PSE_Difference == 0.05 & JND_Difference == 0) %>% 
  group_by(Model,Repetition) %>% 
  slice(1) %>% 
  group_by(Model) %>% 
  mutate(Power_CoI = sum(PvaluesCoI < 0.05)/length(PvaluesCoI),
         Power_Interac = sum(PvaluesInterac < 0.05)/length(PvaluesInterac))

FigurePvaluesFromGLMM_CoI_0.05_0 = ggplot(Dataframe6,aes(Model,Power_CoI)) +
  geom_point(size = 5) + 
  ylab("Power") +
  xlab("") +
  geom_hline(yintercept = 0.05, linetype = 3) +
  geom_hline(yintercept = 0.9, linetype = 2) +
  coord_cartesian(ylim = c(0,1)) +
  ggtitle("Power PSE coefficient (PSE Difference = 0.05/JND Difference = 0)")
FigurePvaluesFromGLMM_Interac_0.05_0 = ggplot(Dataframe6,aes(Model,Power_Interac)) +
  geom_point(size = 5) + 
  ylab("Power") +
  xlab("") +
  geom_hline(yintercept = 0.05, linetype = 3) +
  geom_hline(yintercept = 0.9, linetype = 2) +
  ggtitle("Power JND coefficient (PSE Difference = 0.05/JND Difference = 0)") +
  facet_grid(Repetitions~.)
plot_grid(FigurePvaluesFromGLMM_CoI_0.05_0,FigurePvaluesFromGLMM_Interac_0.05_0, nrow = 2, labels = "AUTO")
ggsave(paste0("Figures/FigurePValuesFromGLMMs_0.05_0.jpeg"), w = 12, h = 7.2)
#####

#####
#AICs
ggplot(Dataframe2 %>% filter(PSE_Difference == 0.05 & JND_Difference == 0),aes(Model,AIC_Norm)) +
  geom_boxplot() + 
  xlab("") +
  ylab("AIC - AIC(Model25)") +
  geom_hline(aes(yintercept = 0), linetype = 2, size = 1) +
  ggtitle("AIC (PSE Difference = 0.05/JND Difference = 0)")
ggsave(paste0("Figures/Figure AICs_0.05_0_Models.jpeg"), w = 12, h = 3.6)
#####


###################################################################################################################
###################################PSE = 0.05, JND = -0.125##########################################################
#_0.05_M0.125
FigurePSEsFromGLMM_0.05_M0.125 = ggplot(Dataframe2 %>% filter(PSE_Difference == 0.05 & JND_Difference == -0.125),aes(Model,Mean_Modeled-Mean_Actual)) +
  geom_boxplot() + 
  ylab("PSE Output of GLMM - Actual PSE") +
  geom_hline(aes(yintercept = 0), linetype = 2, size = 1) + 
  coord_cartesian(ylim = c(-3,3)) +
  xlab("") +
  ggtitle("PSE coefficient (PSE Difference = 0.05/JND Difference = -0.125)")
FigureSDsFromGLMM_0.05_M0.125 = ggplot(Dataframe2 %>% filter(PSE_Difference == 0.05 & JND_Difference == -0.125),aes(Model,SD_Modeled-SD_Actual)) +
  geom_boxplot() + 
  coord_cartesian(ylim = c(-1,2.5)) +
  geom_hline(aes(yintercept = 0), linetype = 2, size = 1) +
  ylab("SD Output of GLMM - Actual SD") +
  xlab("") +
  ggtitle("JND coefficient (PSE Difference = 0.05/JND Difference = -0.125)")
plot_grid(FigurePSEsFromGLMM_0.05_M0.125,FigureSDsFromGLMM_0.05_M0.125, nrow = 2, labels = "AUTO")
ggsave(paste0("Figures/FigureValuesFromGLMMs_0.05_-0.125.jpeg"), w = 12, h = 7.2)
#####

###SE
FigureSEsPSEsFromGLMM_0.05_M0.125 = ggplot(Dataframe2 %>% filter(PSE_Difference == 0.05 & JND_Difference == -0.125),
                                           aes(Model,SECoI)) +
  geom_boxplot() + 
  ylab("SE for PSE estimate") +
  xlab("") +
  coord_cartesian(ylim = c(0,0.2)) +
  ggtitle("SE of PSE coefficient (PSE Difference = 0.05/JND Difference = -0.125)")
FigureSEsSDsFromGLMM_0.05_M0.125 = ggplot(Dataframe2 %>% filter(PSE_Difference == 0.05 & JND_Difference == -0.125),
                                          aes(Model,SEInterac)) +
  geom_boxplot() + 
  ylab("SE for SD estimate") +
  xlab("") +
  coord_cartesian(ylim = c(0,0.2)) +
  ggtitle("SE of JND coefficient (PSE Difference = 0.05/JND Difference = -0.125)")
FigureSEsFromGLMMs_0.05_M0.125 = plot_grid(FigureSEsPSEsFromGLMM_0.05_M0.125,FigureSEsSDsFromGLMM_0.05_M0.125, nrow = 2, labels = "AUTO")
ggsave(paste0("Figures/FigureSEsFromGLMMs_0.05_-0.125.jpeg"), w = 12, h = 7.2)


#####
#p values
Dataframe7 = Dataframe2 %>%
  filter(PSE_Difference == 0.05 & JND_Difference == -0.125) %>% 
  group_by(Model,Repetition) %>% 
  slice(1) %>% 
  group_by(Model) %>% 
  mutate(Power_CoI = sum(PvaluesCoI < 0.05)/length(PvaluesCoI),
         Power_Interac = sum(PvaluesInterac < 0.05)/length(PvaluesInterac))

FigurePvaluesFromGLMM_CoI_0.05_M0.125 = ggplot(Dataframe7,aes(Model,Power_CoI)) +
  geom_point(size = 5) + 
  ylab("Power") +
  xlab("") +
  geom_hline(yintercept = 0.05, linetype = 3) +
  geom_hline(yintercept = 0.9, linetype = 2) +
  coord_cartesian(ylim = c(0,1)) +
  ggtitle("Power PSE coefficient (PSE Difference = 0.05/JND Difference = -0.125)")
FigurePvaluesFromGLMM_Interac_0.05_M0.125 = ggplot(Dataframe7,aes(Model,Power_Interac)) +
  geom_point(size = 5) + 
  ylab("Power") +
  xlab("") +
  geom_hline(yintercept = 0.05, linetype = 3) +
  geom_hline(yintercept = 0.9, linetype = 2) +
  coord_cartesian(ylim = c(0,1)) +
  ggtitle("Power JND coefficient (PSE Difference = 0.05/JND Difference = -0.125)")
plot_grid(FigurePvaluesFromGLMM_CoI_0.05_M0.125,FigurePvaluesFromGLMM_Interac_0.05_M0.125, nrow = 2, labels = "AUTO")
ggsave(paste0("Figures/FigurePValuesFromGLMMs_0.05_-0.125.jpeg"), w = 12, h = 7.2)
#####

#####
#AICs
ggplot(Dataframe2 %>% filter(PSE_Difference == 0.05 & JND_Difference == -0.125),aes(Model,AIC_Norm)) +
  geom_boxplot() + 
  xlab("") +
  ylab("AIC - AIC(Model25)") +
  geom_hline(aes(yintercept = 0), linetype = 2, size = 1) +
  ggtitle("AIC (PSE Difference = 0.05/JND Difference = -0.125)")
ggsave(paste0("Figures/Figure AICs_0.05_-0.125_Models.jpeg"), w = 12, h = 3.6)
#####


####Power together
Dataframe8 = Dataframe2 %>%
  group_by(Condition_PSEJND,Model,Repetition) %>% 
  slice(1) %>% 
  group_by(Condition_PSEJND,Model) %>% 
  mutate(Power_CoI = sum(PvaluesCoI < 0.05)/length(PvaluesCoI),
         Power_Interac = sum(PvaluesInterac < 0.05)/length(PvaluesInterac),
         Condition_PSEJND2 = 
           case_when(Condition_PSEJND == "-0.125-0.05" ~ "PSE: -0.05, JND: -0.125",
                     Condition_PSEJND == "-0.1250.05" ~ "PSE: 0.05, JND: -0.125",
                     Condition_PSEJND == "0.125-0.05" ~ "PSE: -0.05, JND: 0.125",
                     Condition_PSEJND == "0.1250" ~ "PSE: 0, JND: 0.125",
                     Condition_PSEJND == "0.1250.05" ~ "PSE: 0.05, JND: 0.125",
                     Condition_PSEJND == "00.05" ~ "PSE: 0.05, JND: 0"))

ggplot(Dataframe8,aes(Model,Power_Interac,color = Condition_PSEJND2, shape = Condition_PSEJND2)) +
  geom_point(size = 5) + 
  ylab("Power") +
  xlab("") +
  geom_hline(yintercept = 0.05, linetype = 3) +
  geom_hline(yintercept = 0.9, linetype = 2) +
  ggtitle("Power") +
  scale_color_discrete(name = "") +
  scale_shape_discrete(name = "") +
  theme(legend.position = "bottom")
ggsave(paste0("Figures/Power All.jpeg"), w = 12, h = 3.6)