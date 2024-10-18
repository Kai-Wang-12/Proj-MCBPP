##############################################################################################################
#-----------     This file includes R codes to conduct the numerical study of borrowing the
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
VAR_hetero = c(seq(0.1,2,0.1),seq(2.25,5,0.25))     # Population variance heterogeneity
sd_C <- 2.5;

############################################################################################
#----------------------------------- Warm-up running  -------------------------------------#
############################################################################################

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

gc()

############################################################################################
#----------------------------- Numerical study for Figure 2  ------------------------------#
############################################################################################

n_cpu = 16 # number of cores

#------------- Situation1  ---------------#
sfInit(cpus = n_cpu, type = 'SOCK', parallel = T) ### parallel computation
sfSource(paste0(wd, "/Function_for_MultipleECs.R"))
sfExportAll()
sfLapply(as.list(VAR_hetero), function(x){
  data.frame(Norm_multiple_numerical_function(N = c(426, 213, 71, 71, 71), Ybar = c(0, 0, 0, 0, 0), 
                                              Yvar = c(sd_C^2, sd_C^2, sd_C^2, sd_C^2 * (x)^2, sd_C^2 * (x)^2), seed = seed), 
             heterogeneity = x)
}) %>%
  rbindlist() %>%
  saveRDS(file = paste0(wd, "/Num_rst_for_MultipleECs/num_situation1.rds"))
sfStop()
gc()

#------------- Situation2  ---------------#
sfInit(cpus = n_cpu, type = 'SOCK', parallel = T) ### parallel computation
sfSource(paste0(wd, "/Function_for_MultipleECs.R"))
sfExportAll()
sfLapply(as.list(VAR_hetero), function(x){
  data.frame(Norm_multiple_numerical_function(N = c(426, 213, 71, 71, 71), Ybar = c(0, 0, 0, 1*sd_C, -1*sd_C), 
                                              Yvar = c(sd_C^2, sd_C^2, sd_C^2, sd_C^2 * (x)^2, sd_C^2 * (x)^2), seed = seed), 
             heterogeneity = x)
}) %>%
  rbindlist() %>%
  saveRDS(file = paste0(wd, "/Num_rst_for_MultipleECs/num_situation2.rds"))
sfStop()
gc()

