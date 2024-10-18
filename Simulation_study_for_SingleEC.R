###############################################################################################################
#--------------   This file includes R codes to conduct the simulation study of leveraging the 
#--------------   Gaussian endpoint with unknown variance from a single external control arm     
################################################################################################################

############################################################################################
#---------------------------- Loading packages and functions  -----------------------------#
############################################################################################

setwd("~/Code")  #  set the working directory
source("Function_for_SingleEC.R")

############################################################################################
#----------------------------------- General setting  ------------------------------------#
############################################################################################
wd <- getwd()
setwd("stan")  # Set the working directory to the stan code file
seed = 2023
LOC_hetero = seq(-2,2,0.1)            # Population mean heterogeneity
VAR_hetero = c(1,1.2,0.8,2.5,0.4)     # Population variance heterogeneity


############################################################################################
#----------- Generating simulated datasets for the single external control arm-------------#
############################################################################################

mu_CC = 0; sd_C = 2.5
N_CT = c(202, 26, 202, 202)
N_CC = c(101, 13, 152, 50)
N_EC = c(101, 13, 50, 152)
norm_effect = c(0.7, 2, 0.7, 0.7)
  
for (i in 1:4){
  for (j in 1:5){
    Creat_dataset_SingleEC(N_CC = N_CC[i], N_CT = N_CT[i], N_EC = N_EC[i], 
                           norm_effect = norm_effect[i], mu_CC = mu_CC, sd_C = sd_C, 
                           LOC_hetero = LOC_hetero, VAR_hetero = VAR_hetero[j],
                           N_dataset = 2500, seed = c(202301,202302,202303,202304)[i]) %>%
      rbindlist() %>%
      saveRDS(paste0(wd,"/Simulation_dataset/data_norm_",i,"_",c("A","B","C","D","E")[j],".rds"))
  }
}

############################################################################################
#----------------------------------- Warm-up running  -------------------------------------#
############################################################################################

mu_CC = 0; sd_C = 2.5

fit0_n <- Norm_convention(N_CT = 202, Ybar_CT = 0, Yvar_CT = sd_C^2, N_CC = 101, Ybar_CC = 0, Yvar_CC = sd_C^2, seed = seed)
fit1_n <- Norm_MPP(N_EC = 101, N_CC = 101, Ybar_EC = 0, Ybar_CC = 0, Yvar_CC = sd_C^2, Yvar_EC = sd_C^2, type = "MPP1", seed = seed)
fit2_n <- Norm_MPP(N_EC = 101, N_CC = 101, Ybar_EC = 0, Ybar_CC = 0, Yvar_CC = sd_C^2, Yvar_EC = sd_C^2, type = "MPP2", seed = seed)
fit3_n <- Norm_MBPP(N_EC = 101, N_CC = 101, Ybar_EC = 0, Ybar_CC = 0, Yvar_CC = sd_C^2, Yvar_EC = sd_C^2, type = "MBPP", seed = seed)
fit4_n <- Norm_MCBPP(N_EC = 101, N_CC = 101, Ybar_EC = 0, Ybar_CC = 0, Yvar_CC = sd_C^2, Yvar_EC = sd_C^2, type = "MCBPP", seed = seed)
fit5_n <- Norm_MCBPP(N_EC = 101, N_CC = 101, Ybar_EC = 0, Ybar_CC = 0, Yvar_CC = sd_C^2, Yvar_EC = sd_C^2, type = "rMCBPP", seed = seed)
fit6_n <- Norm_CP_gamma(N_EC = 101, N_CC = 101, Ybar_EC = 0, Ybar_CC = 0, Yvar_CC = sd_C^2, Yvar_EC = sd_C^2, seed = seed)

############################################################################################
#----------------------------------- Simulation study -------------------------------------#
############################################################################################

n_cpu = 16 # number of cores

# Do not run the entire code, since it will take hours.
# Recommend to select certain Scenario or certain population mean heterogeneity 

for(i in 1:4){ ## select i = 1,2,3,4 for Design 1,2,3,4, respectively
  
  for(j in 1:5){ ## select j = 1,2,3,4,5 for Scenario A,B,C,D,E, respectively
    
    data_name <- paste0("data_norm_",i,"_",c("A","B","C","D","E"),".rds")
    
    sfInit(cpus = n_cpu, type = 'SOCK', parallel = T) ### parallel computation
    sfSource(paste0(wd, "/Function_for_SingleEC.R"))
    sfExportAll()
    readRDS(paste0(wd,"/Simulation_dataset/",data_name[j])) %>%
      # filter(heterogeneity == 0) %>%  ##  select certain population mean heterogeneity instead of all population mean heterogeneity
      sfApply(1, function(y){
        Norm_simulation_function(y, seed = seed) %>%   
          data.table(heterogeneity = y[15])
      })%>%
      rbindlist() %>%
      data.table() %>%
      saveRDS(file = paste0(wd, "/Sim_rst_for_SingleEC/simu_",i,"_",c("A","B","C","D","E")[j],".rds"))
    sfStop()
    gc()
  }
}
gc()


### simulation study of adaptive selection between MPP1 and rMCBPP 
### using different sample variance ratios in Design 1

for(i in 1:2){
  for(j in 1:5){ ## select j = 1,2,3,4,5 for Scenario A,B,C,D,E, respectively
    
    data_name <- paste0("data_norm_",1,"_",c("A","B","C","D","E"),".rds")
    
    sfInit(cpus = n_cpu, type = 'SOCK', parallel = T) ### parallel computation
    sfSource(paste0(wd, "/Function_for_SingleEC.R"))
    sfExportAll()
    readRDS(paste0(wd,"/Simulation_dataset/",data_name[j])) %>%
      # filter(heterogeneity == 0) %>%  ##  select certain population mean heterogeneity instead of all population mean heterogeneity
      sfApply(1, function(y){
        Norm_simulation_adaptive_function(y, tolerance_sd_u = c(1.05,1.2)[i], tolerance_sd_l = c(0.95,0.8)[i], seed = seed) %>% 
          data.table(heterogeneity = y[15], tolerance = c("strict","mild")[i])
      })%>%
      rbindlist() %>%
      data.table() %>%
      saveRDS(file = paste0(wd, "/Sim_rst_for_SingleEC/simu_",i,"_",c("A","B","C","D","E")[j],"_adapt.rds"))
    sfStop()
    gc()
  }
}


###################################################
#####  Calculate the operating characteristics
###################################################

setwd("~/Code")  #  set the working directory
wd <- getwd()
library(tidyverse)
norm_effect = c(0.7, 2, 0.7, 0.7)

for(i in 1:4){ # i: Design 1-4
  
  for(j in 1:5){ # j: different severity of population variance heterogeneity
    
    data <- readRDS((file = paste0(wd, "/Sim_rst_for_SingleEC/simu_",i,"_",c("A","B","C","D","E")[j],".rds")))
    
    threshold <- filter(data, test=="H0", heterogeneity==0) %>%                   # Find threshold value under H0 and no location heterogeneity
      group_by(prior) %>%
      summarise(threshold = quantile(post_probability, probs = 0.975))            # Calibrated threshold probability
    saveRDS(threshold,paste0(wd, "/Sim_rst_for_SingleEC/simu_",i,"_threshlod_",c("A","B","C","D","E")[j],".rds"))
    
    rst_H0 <- data %>%
      filter(test=="H0") %>%   
      full_join(threshold, by = "prior") %>%
      group_by(prior,heterogeneity) %>%
      summarise(ESS = median(ESS), 
                Bias = mean(Delta_mean - 0), 
                MSE = mean((Delta_mean - 0)^2), 
                Type_I_error = mean(post_probability > 0.975),                    # Type I error rate based on 0.975
                Type_I_error_calibrated = mean(post_probability > threshold),     # calibrated Type I error rate based on threshold
                Proportion = median(Proportion)) %>%
      arrange(prior, heterogeneity) %>%
      as.data.frame()
    
    rst_H0_no_borrowing <- filter(rst_H0, prior == "No-borrowing") %>%
      dplyr::select(c(Bias,MSE,Proportion,heterogeneity)) %>%
      set_names(c("Bias_0","MSE_0","Proportion_0","heterogeneity"))
    
    rst_H0 <- rst_H0 %>%
      full_join(rst_H0_no_borrowing, by = "heterogeneity") %>%
      mutate(relative_bias = abs(Bias - Bias_0), relative_mse = MSE - MSE_0, relative_Proportion = Proportion - Proportion_0) %>%
      saveRDS(file = paste0(wd, "/Sim_rst_for_SingleEC/simu_",i,"_H0_",c("A","B","C","D","E")[j],".rds"))
    
    rst_H1 <- data %>%
      filter(test=="H1") %>%   
      full_join(threshold, by = "prior") %>%
      group_by(prior,heterogeneity) %>%
      summarise(ESS = median(ESS), 
                Bias = mean(Delta_mean - norm_effect[i]), 
                MSE = mean((Delta_mean - norm_effect[i])^2), 
                Power = mean(post_probability > 0.975),                           # Power based on 0.975
                Power_calibrated = mean(post_probability > threshold),            # calibrated Power based on threshold
                Proportion = median(Proportion)) %>%
      arrange(prior, heterogeneity) %>%
      as.data.frame() 
    
    rst_H1_no_borrowing <- filter(rst_H1, prior == "No-borrowing") %>%
      dplyr::select(c(Bias,MSE,Proportion,heterogeneity)) %>%
      set_names(c("Bias_0","MSE_0","Proportion_0","heterogeneity"))
    
    rst_H1 <- rst_H1 %>%
      full_join(rst_H1_no_borrowing, by = "heterogeneity") %>%
      mutate(relative_bias = abs(Bias - Bias_0), relative_mse = MSE - MSE_0, relative_Proportion = Proportion - Proportion_0) %>%
      saveRDS(file = paste0(wd, "/Sim_rst_for_SingleEC/simu_",i,"_H1_",c("A","B","C","D","E")[j],".rds"))
    
  }
}


### simulation study of adaptive selection between MPP1 and rMCBPP 
### using different sample variance ratios in Design 1
for(i in 1:2){ # i=1: strict tolerance; i=2: mild tolerance
  for(j in 1:5){ # j: different severity of population variance heterogeneity
    data <- readRDS((file = paste0(wd, "/Sim_rst_for_SingleEC/simu_",i,"_",c("A","B","C","D","E")[j],"_adapt.rds")))
    data$prior[data$prior!="No-borrowing"] = paste0("adaptive_",c("strict","mild")[i])
    
    threshold <- filter(data, test=="H0", heterogeneity==0) %>%                   # Find threshold value under H0 and no location heterogeneity
      group_by(prior) %>%
      summarise(threshold = quantile(post_probability, probs = 0.975))            # Calibrated threshold probability
    saveRDS(threshold,paste0(wd, "/Sim_rst_for_SingleEC/simu_",i,"_threshlod_",c("A","B","C","D","E")[j],"_adapt.rds"))
    
    rst_H0 <- data %>%
      filter(test=="H0") %>%   
      full_join(threshold, by = "prior") %>%
      group_by(prior,heterogeneity) %>%
      summarise(ESS = median(ESS), 
                Bias = mean(Delta_mean - 0), 
                MSE = mean((Delta_mean - 0)^2), 
                Type_I_error = mean(post_probability > 0.975),                    # Type I error rate based on 0.975
                Type_I_error_calibrated = mean(post_probability > threshold),     # calibrated Type I error rate based on threshold
                Proportion = median(Proportion)) %>%
      arrange(prior, heterogeneity) %>%
      as.data.frame()
    
    rst_H0_no_borrowing <- filter(rst_H0, prior == "No-borrowing") %>%
      dplyr::select(c(Bias,MSE,Proportion,heterogeneity)) %>%
      set_names(c("Bias_0","MSE_0","Proportion_0","heterogeneity"))
    
    rst_H0 <- rst_H0 %>%
      full_join(rst_H0_no_borrowing, by = "heterogeneity") %>%
      mutate(relative_bias = abs(Bias - Bias_0), relative_mse = MSE - MSE_0, relative_Proportion = Proportion - Proportion_0) %>%
      saveRDS(file = paste0(wd, "/Sim_rst_for_SingleEC/simu_",i,"_H0_",c("A","B","C","D","E")[j],"_adapt.rds"))
    
    rst_H1 <- data %>%
      filter(test=="H1") %>%   
      full_join(threshold, by = "prior") %>%
      group_by(prior,heterogeneity) %>%
      summarise(ESS = median(ESS), 
                Bias = mean(Delta_mean - 0.7), 
                MSE = mean((Delta_mean - 0.7)^2), 
                Power = mean(post_probability > 0.975),                           # Power based on 0.975
                Power_calibrated = mean(post_probability > threshold),            # calibrated Power based on threshold
                Proportion = median(Proportion)) %>%
      arrange(prior, heterogeneity) %>%
      as.data.frame() 
    
    rst_H1_no_borrowing <- filter(rst_H1, prior == "No-borrowing") %>%
      dplyr::select(c(Bias,MSE,Proportion,heterogeneity)) %>%
      set_names(c("Bias_0","MSE_0","Proportion_0","heterogeneity"))
    
    rst_H1 <- rst_H1 %>%
      full_join(rst_H1_no_borrowing, by = "heterogeneity") %>%
      mutate(relative_bias = abs(Bias - Bias_0), relative_mse = MSE - MSE_0, relative_Proportion = Proportion - Proportion_0) %>%
      saveRDS(file = paste0(wd, "/Sim_rst_for_SingleEC/simu_",i,"_H1_",c("A","B","C","D","E")[j],"_adapt.rds"))
    
  }
}


###################################################
### Summary of calibrated threshold value: Table S2
###################################################

setwd("~/Code")  #  set the working directory
wd <- getwd()
library(tidyverse)
library(data.table)

simu_threshlod <- list()
for(i in 1:4){
  simu_threshlod[[i]] <- lapply(as.list(c("A","B","C","D","E")),function(x){
    readRDS(paste0(wd, "/Sim_rst_for_SingleEC/simu_",i,"_threshlod_",x,".rds")) %>% 
      mutate(scenario = paste0("Scenario ",x))}) %>%
    rbindlist() %>%
    mutate(Design = i, threshold = round(threshold,4)) %>%
    spread(key = scenario, value = threshold)
}
simu_threshlod <- rbindlist(simu_threshlod)

Type_prior = c("MPP1","MPP2","MBPP","MCBPP","rMCBPP","CP","Pooled","No-borrowing")

Table_S2 <- full_join(data.frame(prior = Type_prior), simu_threshlod ,by = "prior") %>% arrange(Design)
saveRDS(Table_S2, file = paste0(wd, "/Sim_rst_for_SingleEC/Table_S2.rds"))
