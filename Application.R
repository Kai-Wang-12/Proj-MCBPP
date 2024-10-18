##################################################################################################################################################################
#-------------------------------------------       This file includes R codes to conduct the application      ---------------------------------------------------#
##################################################################################################################################################################

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
seed =2023

#----------------------------------- Warm-up running  

sd_C = 2.5

fit0_n <- Norm_convention(N_CT = 202, Ybar_CT = 0, Yvar_CT = sd_C^2, N_CC = 101, Ybar_CC = 0, Yvar_CC = sd_C^2, seed = seed)
fit1_n <- Norm_MPP(N_EC = 101, N_CC = 101, Ybar_EC = 0, Ybar_CC = 0, Yvar_CC = sd_C^2, Yvar_EC = sd_C^2, type = "MPP1", seed = seed)
fit2_n <- Norm_MPP(N_EC = 101, N_CC = 101, Ybar_EC = 0, Ybar_CC = 0, Yvar_CC = sd_C^2, Yvar_EC = sd_C^2, type = "MPP2", seed = seed)
fit3_n <- Norm_MBPP(N_EC = 101, N_CC = 101, Ybar_EC = 0, Ybar_CC = 0, Yvar_CC = sd_C^2, Yvar_EC = sd_C^2, type = "MBPP", seed = seed)
fit4_n <- Norm_MCBPP(N_EC = 101, N_CC = 101, Ybar_EC = 0, Ybar_CC = 0, Yvar_CC = sd_C^2, Yvar_EC = sd_C^2, type = "MCBPP", seed = seed)
fit5_n <- Norm_MCBPP(N_EC = 101, N_CC = 101, Ybar_EC = 0, Ybar_CC = 0, Yvar_CC = sd_C^2, Yvar_EC = sd_C^2, type = "rMCBPP", seed = seed)
fit6_n <- Norm_CP_gamma(N_EC = 101, N_CC = 101, Ybar_EC = 0, Ybar_CC = 0, Yvar_CC = sd_C^2, Yvar_EC = sd_C^2, seed = seed)


#######################################################################
#------------------------- Analysis ----------------------------------#
#######################################################################

case1 <- c(280, 139, 271, -2.65, -0.87, -1.03, 3.333^2, 2.833^2, 3.004^2)
rst_1 <- Norm_application_function(case1, seed = 202306)
rst_1
saveRDS(rst_1,file = paste0(wd,"/Application_rst/rst_case1.rds"))

case2 <- c(280, 139, 271, -2.65, -0.87, -1.03, 3.333^2, 2.833^2, 2.833^2)
rst_2 <- Norm_application_function(case2, seed = 202306)
rst_2
saveRDS(rst_2,file = paste0(wd,"/Application_rst/rst_case2.rds"))

#------------- Varying sample SD ratio  ---------------#
SD_hetero = seq(0.65,1.5,0.05) # sample SD ratio
n_cpu = 13 # number of cores
sfInit(cpus = n_cpu, type = 'SOCK', parallel = T) ### parallel computation
sfSource(paste0(wd, "/Function_for_SingleEC.R"))
sfExportAll()
sfLapply(as.list(SD_hetero), function(x){
  rst <- data.frame(Norm_application_function(c(280, 139, 271, -2.65, -0.87, -1.03, 3.333^2, 2.833^2, (2.833*x)^2), seed = 202306), 
                    SD_hetero = x)
}) %>%
  rbindlist() %>%
  saveRDS(file = paste0(wd, "/Application_rst/rst_varying_ratio.rds"))
sfStop()
gc()

#######################################################################
#------------------------- Table 2  ----------------------------------#
#######################################################################

rst_1 <- readRDS(paste0(wd,"/Application_rst/rst_case1.rds"))
rst_2 <- readRDS(paste0(wd,"/Application_rst/rst_case2.rds"))
rst_varying_ratio <- readRDS(paste0(wd, "/Application_rst/rst_varying_ratio.rds"))

rst_case1 <- rst_1 %>%
  group_by(prior) %>%
  mutate(ESS0 = rst_1$ESS[rst_1$prior == "No-borrowing"], Delta_mean0 = rst_1$Delta_mean[rst_1$prior == "No-borrowing"],
         CI_length0 = rst_1$CI_length[rst_1$prior == "No-borrowing"]) %>%
  mutate(relative_ESS = round(ESS - ESS0,4), 
         Est = round(Delta_mean,4),
         relative_Est = round(Delta_mean - Delta_mean0,4),
         relative_CI = round(CI_length - CI_length0,4)) %>%
  as.data.frame() %>%
  dplyr::select(c(prior,relative_ESS,Est,relative_Est,CI_length,relative_CI)) %>%
  mutate(case = 1)

rst_case2 <- rst_2 %>%
  group_by(prior) %>%
  mutate(ESS0 = rst_2$ESS[rst_2$prior == "No-borrowing"], Delta_mean0 = rst_2$Delta_mean[rst_2$prior == "No-borrowing"],
         CI_length0 = rst_2$CI_length[rst_2$prior == "No-borrowing"]) %>%
  mutate(relative_ESS = round(ESS - ESS0,4), 
         Est = round(Delta_mean,4),
         relative_Est = round(Delta_mean - Delta_mean0,4),
         relative_CI = round(CI_length - CI_length0,4)) %>%
  as.data.frame() %>%
  dplyr::select(c(prior,relative_ESS,Est,relative_Est,CI_length,relative_CI)) %>%
  mutate(case = 2)

rst_plot_varying_ratio <- filter(rst_varying_ratio, prior == "No-borrowing") %>% 
  dplyr::select(c(SD_hetero,ESS,Delta_mean,CI_length)) %>% 
  set_names(c("SD_hetero","ESS0","Delta_mean0","CI_length0")) %>%
  full_join(rst_varying_ratio, by = "SD_hetero") %>%
  mutate(relative_ESS = ESS - ESS0, 
         Est = Delta_mean,
         relative_Est = Delta_mean - Delta_mean0, 
         relative_CI = CI_length - CI_length0) %>%
  as.data.frame() %>%
  dplyr::select(c(prior,SD_hetero,relative_ESS,Est,relative_Est,CI_length,relative_CI)) %>%
  saveRDS(file = paste0(wd,"/Application_rst/rst_plot_varying_ratio.rds"))


#########  Table 2
Table_2 <- rbind(rst_case1,rst_case2)
saveRDS(Table_2, file = paste0(wd,"/Application_rst/Table_2.rds"))







