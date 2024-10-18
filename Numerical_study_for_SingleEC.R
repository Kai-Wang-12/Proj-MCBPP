##############################################################################################################
#-------------    This file includes R codes to conduct the numerical study of borrowing the 
#-------------    Gaussian endpoint with unknown variance from a single external control arm        
##############################################################################################################

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
LOC_hetero = seq(-2,2,0.1)                            # Population mean heterogeneity
VAR_hetero = c(seq(0.1,2,0.1),seq(2.25,5,0.25))       # Population variance heterogeneity
sd_CC <- 2.5;  # sample standard error of concurrent controls

############################################################################################
#----------------------------------- Warm-up running  -------------------------------------#
############################################################################################

fit1_n <- Norm_MPP(N_EC = 101, N_CC = 101, Ybar_EC = 0, Ybar_CC = 0, Yvar_CC = sd_CC^2, Yvar_EC = sd_CC^2, type = "MPP1", seed = seed)
fit2_n <- Norm_MPP(N_EC = 101, N_CC = 101, Ybar_EC = 0, Ybar_CC = 0, Yvar_CC = sd_CC^2, Yvar_EC = sd_CC^2, type = "MPP2", seed = seed)
fit3_n <- Norm_MBPP(N_EC = 101, N_CC = 101, Ybar_EC = 0, Ybar_CC = 0, Yvar_CC = sd_CC^2, Yvar_EC = sd_CC^2, type = "MBPP", seed = seed)
fit4_n <- Norm_MCBPP(N_EC = 101, N_CC = 101, Ybar_EC = 0, Ybar_CC = 0, Yvar_CC = sd_CC^2, Yvar_EC = sd_CC^2, type = "MCBPP", seed = seed)
fit5_n <- Norm_MCBPP(N_EC = 101, N_CC = 101, Ybar_EC = 0, Ybar_CC = 0, Yvar_CC = sd_CC^2, Yvar_EC = sd_CC^2, type = "rMCBPP", seed = seed)

gc()

############################################################################################
#-------------------------- Numerical study for Figure 1, S2-4  ---------------------------#
############################################################################################

N1 = c(202, 101, 101)   # Figure1
N2 = c(26, 13, 13)      # FigureS2
N3 = c(202, 152, 50)    # FigureS3
N4 = c(202, 50, 152)    # FigureS4

N = cbind(N1,N2,N3,N4)

n_cpu = 16 # number of cores

#------------- No population mean heterogeneity  ---------------#
sfInit(cpus = n_cpu, type = 'SOCK', parallel = T) ### parallel computation
sfSource(paste0(wd, "/Function_for_SingleEC.R"))
sfExportAll()
sfLapply(as.list(VAR_hetero), function(x){
  rst <- list()
  for(i in 1:4){
    rst[[i]] <- data.frame(Norm_numerical_function(N = N[,i], Ybar = c(0, 0, 0), Yvar = c(sd_CC^2, sd_CC^2, sd_CC^2 * (x)^2), seed = seed), 
               VAR_hetero = x, design = i)
  }
  return(rbindlist(rst))
}) %>%
  rbindlist() %>%
  saveRDS(file = paste0(wd, "/Num_rst_for_SingleEC/num_0.rds"))
sfStop()
gc()


#------------- Different population variance heterogeneity  ---------------#
for(i in 1:5){
  sfInit(cpus = n_cpu, type = 'SOCK', parallel = T) ### parallel computation
  sfSource(paste0(wd, "/Function_for_SingleEC.R"))
  sfExportAll()
  sfLapply(as.list(LOC_hetero), function(x){
    rst <- list()
    for(j in 1:4){
      rst[[j]] <- data.frame(Norm_numerical_function(N = N[,j], Ybar = c(0, 0, x*sd_CC), Yvar = c(sd_CC^2, sd_CC^2, (c(0.4,0.8,1,1.2,2.5)[i]^2)*sd_CC^2), seed = seed), 
                             LOC_hetero = x, design = j)
    }
    return(rbindlist(rst))
  }) %>%
    rbindlist() %>%
    mutate(scenario = c("A","B","C","D","E")[i]) %>%
    saveRDS(file = paste0(wd, "/Num_rst_for_SingleEC/num_",c("A","B","C","D","E")[i],".rds"))
  sfStop()
}

