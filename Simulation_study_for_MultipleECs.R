##############################################################################################################
#-----------     This file includes R codes to conduct the simulation study of leveraging the
#-----------     Gaussian endpoint with unknown variance from multiple (three) external control arms      
##############################################################################################################


############################################################################################
#---------------------------- Loading packages and functions  -----------------------------#
############################################################################################

setwd("~/Code")  #  set the working directory
source("Function_for_MultipleECs.R")

############################################################################################
#----------------------------------- General setting  ------------------------------------#
############################################################################################
wd <- getwd()
setwd("stan")  # Set the working directory to the stan code file
seed =2023
LOC_hetero = seq(2,-2,-0.1)     # Population mean heterogeneity                                              
VAR_hetero = c(1,0.4,2.5)       # Population variance heterogeneity                                         
sd_C = 2.5

############################################################################################
#---------- Generating simulated datasets for the multiple external control arms-----------#
############################################################################################

norm_effect = 0.48

for (i in 1:3){  ## i=1,2,3 indicates Scenario F,G,H
  Creat_dataset_MultipleECs(N_CC = 213, N_CT = 426, N_EC1 = 71, N_EC2 = 71, N_EC3 = 71, 
                            mu_CC = 0, mu_EC1 = 0, mu_EC2 = 0, mu_EC3 = 0,
                            LOC_hetero = LOC_hetero, norm_effect = norm_effect, 
                            sd_C = sd_C, sd_EC1 = sd_C, sd_EC2 = VAR_hetero[i]*sd_C, sd_EC3 = VAR_hetero[i]*sd_C, 
                            N_dataset = 2500, seed = 202305) %>%
    rbindlist() %>%
    saveRDS(paste0(wd,"/Simulation_dataset/data_norm_5_",c("F","G","H")[i],"_multi.rds")) 
}

for (i in 1:3){   ## i=1,2,3 indicates Scenario I,J,K
  Creat_dataset_MultipleECs(N_CC = 213, N_CT = 426, N_EC1 = 71, N_EC2 = 71, N_EC3 = 71, 
                            mu_CC = 0, mu_EC1 = 0, mu_EC2 = sd_C, mu_EC3 = -1*sd_C,
                            LOC_hetero = LOC_hetero, norm_effect = norm_effect, 
                            sd_C = sd_C, sd_EC1 = sd_C, sd_EC2 = VAR_hetero[i]*sd_C, sd_EC3 = VAR_hetero[i]*sd_C, 
                            N_dataset = 2500, seed = 202305) %>%
    rbindlist() %>%
    saveRDS(paste0(wd,"/Simulation_dataset/data_norm_5_",c("I","J","K")[i],"_multi.rds")) 
}

############################################################################################
#----------------------------   Warm-up running (save time)    ----------------------------#
############################################################################################

fit0_n <- Norm_convention(N_CT = 426, Ybar_CT = 0, Yvar_CT = sd_C^2, N_CC = 213, Ybar_CC = 0, Yvar_CC = sd_C^2, seed = seed)

fit1_n_O <- Norm_multi_MPP(N_CC = 213, Ybar_CC = 0, Yvar_CC = sd_C^2, N_EC1 = 71, Ybar_EC1 = 0, Yvar_EC1 = sd_C^2, 
                           N_EC2 = 71, Ybar_EC2 = 1*sd_C, Yvar_EC2 = sd_C^2, N_EC3 = 71, Ybar_EC3 = -1*sd_C, Yvar_EC3 = sd_C^2, type = "MPP1_O", seed = seed)
fit1_n_A <- Norm_multi_MPP(N_CC = 213, Ybar_CC = 0, Yvar_CC = sd_C^2, N_EC1 = 71, Ybar_EC1 = 0, Yvar_EC1 = sd_C^2, 
                           N_EC2 = 71, Ybar_EC2 = 1*sd_C, Yvar_EC2 = sd_C^2, N_EC3 = 71, Ybar_EC3 = -1*sd_C, Yvar_EC3 = sd_C^2, type = "MPP1_A", seed = seed) 
fit2_n_O <- Norm_multi_MPP(N_CC = 213, Ybar_CC = 0, Yvar_CC = sd_C^2, N_EC1 = 71, Ybar_EC1 = 0, Yvar_EC1 = sd_C^2, 
                           N_EC2 = 71, Ybar_EC2 = 1*sd_C, Yvar_EC2 = sd_C^2, N_EC3 = 71, Ybar_EC3 = -1*sd_C, Yvar_EC3 = sd_C^2, type = "MPP2_O", seed = seed) 
fit2_n_A <- Norm_multi_MPP(N_CC = 213, Ybar_CC = 0, Yvar_CC = sd_C^2, N_EC1 = 71, Ybar_EC1 = 0, Yvar_EC1 = sd_C^2, 
                           N_EC2 = 71, Ybar_EC2 = 1*sd_C, Yvar_EC2 = sd_C^2, N_EC3 = 71, Ybar_EC3 = -1*sd_C, Yvar_EC3 = sd_C^2, type = "MPP2_A", seed = seed)
fit3_n_O <- Norm_multi_MBPP(N_CC = 213, Ybar_CC = 0, Yvar_CC = sd_C^2, N_EC1 = 71, Ybar_EC1 = 0, Yvar_EC1 = sd_C^2, 
                            N_EC2 = 71, Ybar_EC2 = 1*sd_C, Yvar_EC2 = sd_C^2, N_EC3 = 71, Ybar_EC3 = -1*sd_C, Yvar_EC3 = sd_C^2, type = "MBPP_O", seed = seed) 
fit3_n_A <- Norm_multi_MBPP(N_CC = 213, Ybar_CC = 0, Yvar_CC = sd_C^2, N_EC1 = 71, Ybar_EC1 = 0, Yvar_EC1 = sd_C^2, 
                            N_EC2 = 71, Ybar_EC2 = 1*sd_C, Yvar_EC2 = sd_C^2, N_EC3 = 71, Ybar_EC3 = -1*sd_C, Yvar_EC3 = sd_C^2, type = "MBPP_A", seed = seed) 
fit4_n_O <- Norm_multi_MCBPP(N_CC = 213, Ybar_CC = 0, Yvar_CC = sd_C^2, N_EC1 = 71, Ybar_EC1 = 0, Yvar_EC1 = sd_C^2, 
                             N_EC2 = 71, Ybar_EC2 = 1*sd_C, Yvar_EC2 = sd_C^2, N_EC3 = 71, Ybar_EC3 = -1*sd_C, Yvar_EC3 = sd_C^2, type = "MCBPP_O", seed = seed) 
fit4_n_A <- Norm_multi_MCBPP(N_CC = 213, Ybar_CC = 0, Yvar_CC = sd_C^2, N_EC1 = 71, Ybar_EC1 = 0, Yvar_EC1 = sd_C^2, 
                             N_EC2 = 71, Ybar_EC2 = 1*sd_C, Yvar_EC2 = sd_C^2, N_EC3 = 71, Ybar_EC3 = -1*sd_C, Yvar_EC3 = sd_C^2, type = "MCBPP_A", seed = seed) 
fit5_n_A <- Norm_multi_MCBPP(N_CC = 213, Ybar_CC = 0, Yvar_CC = sd_C^2, N_EC1 = 71, Ybar_EC1 = 0, Yvar_EC1 = sd_C^2, 
                             N_EC2 = 71, Ybar_EC2 = 1*sd_C, Yvar_EC2 = sd_C^2, N_EC3 = 71, Ybar_EC3 = -1*sd_C, Yvar_EC3 = sd_C^2, type = "rMCBPP_A", seed = seed) 
fit6_n <- Norm_multi_CP_gamma(N_CC = 213, Ybar_CC = 0, Yvar_CC = sd_C^2, N_EC1 = 71, Ybar_EC1 = 0, Yvar_EC1 = sd_C^2, 
                              N_EC2 = 71, Ybar_EC2 = 1*sd_C, Yvar_EC2 = sd_C^2, N_EC3 = 71, Ybar_EC3 = -1*sd_C, Yvar_EC3 = sd_C^2, seed = seed)

gc()

############################################################################################
#----------------------------------- Simulation study -------------------------------------#
############################################################################################

n_cpu = 16 # number of cores

# Do not run the entire code, since it will take hours.
# Recommend to select certain Scenario or certain population mean heterogeneity 

data_name <- paste0("data_norm_5_",c("F","G","H","I","J","K"),"_multi.rds") 

for(i in 1:length(data_name)){  ## select i=1,2,3,5,6 for Scenario F,G,H,I,J,K 
  sfInit(cpus = n_cpu, type = 'SOCK', parallel = T) ### parallel computation
  sfSource(paste0(wd, "/Function_for_MultipleECs.R"))
  sfExportAll()
  readRDS(paste0(wd,"/Simulation_dataset/",data_name[i])) %>%
    # filter(heterogeneity == 0) %>%  ##  select certain population mean heterogeneity instead of all population mean heterogeneities
    sfApply(1, function(y){
      data.table(Norm_multiple_simulation_function(y, seed = seed), heterogeneity = y[21])
    })%>%
    rbindlist() %>%
    data.table() %>%
    saveRDS(file = paste0(wd, "/Sim_rst_for_MultipleECs/simu_5_",c("F","G","H","I","J","K")[i],".rds"))
  sfStop()
  gc()
}

###################################################
#####  Calculate the operating characteristics
###################################################

norm_effect = 0.48

for(i in 1:6){
  
  threshold <- filter(readRDS(file = paste0(wd, "/Sim_rst_for_MultipleECs/simu_5_",c("F","G","H","F","G","H")[i],".rds")), 
                      test=="H0", heterogeneity == 0) %>%                       # Find threshold probability under H0 and no population mean heterogeneity (Scenario F,G,H)
    group_by(prior,strategy) %>%
    summarise(threshold = quantile(post_probability, probs = 0.975))            # Control type I error rate at 0.025
  saveRDS(threshold,paste0(wd, "/Sim_rst_for_MultipleECs/simu_multi_threshlod",i,".rds"))
  
  data <- readRDS(file = paste0(wd, "/Sim_rst_for_MultipleECs/simu_5_",c("F","G","H","I","J","K")[i],".rds"))
  
  rst_H0 <- data %>%
    filter(test=="H0") %>%   
    full_join(threshold, by = c("prior","strategy")) %>%
    group_by(prior,strategy,heterogeneity) %>%
    summarise(ESS = median(ESS), 
              Bias = mean(Delta_mean - 0), 
              MSE = mean((Delta_mean - 0)^2), 
              Type_I_error = mean(post_probability > 0.975),                    # Type I error rate based on 0.975
              Type_I_error_calibrated = mean(post_probability > threshold),     # calibrated Type I error rate based on calibrated threshold probability
              Proportion = median(Proportion)) %>%
    arrange(prior, strategy, heterogeneity) %>%
    as.data.frame() 
  
  rst_H0_no_borrowing <- filter(rst_H0, prior == "No-borrowing") %>%
    dplyr::select(c(Bias,MSE,Proportion,strategy,heterogeneity)) %>%
    set_names(c("Bias_0","MSE_0","Proportion_0","strategy","heterogeneity"))
  
  rst_H0 <- rst_H0 %>%
    full_join(rst_H0_no_borrowing, by = c("heterogeneity","strategy")) %>%
    mutate(relative_bias = abs(Bias - Bias_0), relative_mse = MSE - MSE_0, relative_Proportion = Proportion - Proportion_0) %>%
    saveRDS(file = paste0(wd, "/Sim_rst_for_MultipleECs/simu_5_",c("F","G","H","I","J","K")[i],"_H0.rds"))
  
  rst_H1 <- data %>%
    filter(test=="H1") %>%   
    full_join(threshold, by = c("prior","strategy")) %>%
    group_by(prior,strategy,heterogeneity) %>%
    summarise(ESS = median(ESS), 
              Bias = mean(Delta_mean - norm_effect), 
              MSE = mean((Delta_mean - norm_effect)^2), 
              Power = mean(post_probability > 0.975),                           # Power based on 0.975
              Power_calibrated = mean(post_probability > threshold),            # calibrated Power based on calibrated threshold probability
              Proportion = median(Proportion)) %>%
    arrange(prior, strategy, heterogeneity) %>%
    as.data.frame()
  
  rst_H1_no_borrowing <- filter(rst_H1, prior == "No-borrowing") %>%
    dplyr::select(c(Bias,MSE,Proportion,strategy,heterogeneity)) %>%
    set_names(c("Bias_0","MSE_0","Proportion_0","strategy","heterogeneity"))
  
  rst_H1 <- rst_H1 %>%
    full_join(rst_H1_no_borrowing, by = c("heterogeneity","strategy")) %>%
    mutate(relative_bias = abs(Bias - Bias_0), relative_mse = MSE - MSE_0, relative_Proportion = Proportion - Proportion_0) %>%
    saveRDS(file = paste0(wd, "/Sim_rst_for_MultipleECs/simu_5_",c("F","G","H","I","J","K")[i],"_H1.rds"))
}

###################################################
### Summary of calibrated threshold value: Table S3
###################################################

setwd("~/Code")  #  set the working directory
wd <- getwd()
library(tidyverse)
library(data.table)

Type_prior = c("MPP1","MPP2","MBPP","MCBPP","rMCBPP","CP","Pooled","No-borrowing")

simu_multi_threshlod <- lapply(as.list(1:3),function(x){
  readRDS(paste0(wd, "/Sim_rst_for_MultipleECs/simu_multi_threshlod",x,".rds")) %>% 
    mutate(scenario = paste0("Scenario ",c("F & I","G & J","H & K")[x]))}) %>%
  rbindlist() %>%
  mutate(Deisgn = 5, threshold = round(threshold,4)) %>%
  spread(key = scenario, value = threshold)
Table_S3 <- full_join(data.frame(prior = Type_prior), simu_multi_threshlod ,by = "prior")
saveRDS(Table_S3, file = paste0(wd, "/Sim_rst_for_MultipleECs/Table_S3.rds"))
