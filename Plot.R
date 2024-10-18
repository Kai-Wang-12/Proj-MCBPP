
##############################################################################################
#-----------     This file includes R codes to plot the figures in the article   ------------#
##############################################################################################

library(plotfunctions)
library(tidyverse)
library(data.table)

setwd("~/Code")  #  set the working directory
wd = getwd()

source("Plot_function.R")
setwd("Figure")


################################################
### Numerical results: Figure 1 & Figure S1-3
################################################

### load numerical results

data_0 = readRDS(paste0(wd,"/Num_rst_for_SingleEC/num_0.rds"))
data_A = readRDS(paste0(wd,"/Num_rst_for_SingleEC/num_A.rds"))
data_B = readRDS(paste0(wd,"/Num_rst_for_SingleEC/num_B.rds"))
data_C = readRDS(paste0(wd,"/Num_rst_for_SingleEC/num_C.rds"))
data_D = readRDS(paste0(wd,"/Num_rst_for_SingleEC/num_D.rds"))
data_E = readRDS(paste0(wd,"/Num_rst_for_SingleEC/num_E.rds"))

main_text = c(expression(paste("A. MPP1"," (",delta,")")),expression(paste("B. MPP2"," (",delta,")")),
              expression(paste("C. MBPP"," (",delta[1],")")), expression(paste("D. MCBPP"," (",delta[1],"*)")),
              expression(paste("E. rMCBPP"," (",delta[1],"*)")))


plot_num(data = list(filter(rbind(data_A,data_B,data_C,data_D,data_E), design == 1), filter(data_0, design == 1)), 
          Type_prior = c("MPP1","MPP2","MBPP","MCBPP","rMCBPP"),
          col_vector = colorRampPalette(c('blue','red'))(5),
          main = main_text, 
          file = "Figure_1")
plot_num(data = list(filter(rbind(data_A,data_B,data_C,data_D,data_E), design == 2), filter(data_0, design == 2)), 
          Type_prior = c("MPP1","MPP2","MBPP","MCBPP","rMCBPP"),
          col_vector = colorRampPalette(c('blue','red'))(5),
          main = main_text, 
          file = "Figure_S1")
plot_num(data = list(filter(rbind(data_A,data_B,data_C,data_D,data_E), design == 3), filter(data_0, design == 3)), 
          Type_prior = c("MPP1","MPP2","MBPP","MCBPP","rMCBPP"),
          col_vector = colorRampPalette(c('blue','red'))(5),
          main = main_text, 
          file = "Figure_S2")
plot_num(data = list(filter(rbind(data_A,data_B,data_C,data_D,data_E), design == 4), filter(data_0, design == 4)), 
          Type_prior = c("MPP1","MPP2","MBPP","MCBPP","rMCBPP"),
          col_vector = colorRampPalette(c('blue','red'))(5),
          main = main_text, 
          file = "Figure_S3")

################################################
####  Numerical results: Figure S4
################################################

col_vector1 <- alpha(c(3,"violet",4,7),0.9)
col_vector2 <- alpha(c(3,"violet",4,7,2),0.9)

Type_prior <- c("MPP1","MPP2","MBPP","MCBPP","rMCBPP")
legend_lable1 <- c(expression(paste("MPP1"," (",delta,") ")), 
                   expression(paste("MPP2"," (",delta,") ")),
                   expression(paste("MBPP"," (",delta[1],") ")),
                   expression(paste("MCBPP"," (",delta[1],"*) ")))
legend_lable2 <- c(expression(paste("MPP1"," (",delta,") ")), 
                   expression(paste("MPP2"," (",delta,") ")),
                   expression(paste("MBPP"," (",delta[1],") ")),
                   expression(paste("MCBPP"," (",delta[1],"*) ")),
                   expression(paste("rMCBPP"," (",delta[1],"*) ")))

data_1 = readRDS(paste0(wd,"/Num_rst_for_MultipleECs/num_situation1.rds"))
data_2 = readRDS(paste0(wd,"/Num_rst_for_MultipleECs/num_situation2.rds"))

{
  pdf(width = 5*2, height = 5*2, file = "Figure_S4.pdf")
  par(mfcol = c(2,2), mar = c(5,5,4,2))
  for (i in 1:2){
    plot(NA, xlim = c(0.1,5), ylim = c(0,1), xlab = "Variance-specific prior-data conflict", ylab = "Posterior mean", xaxt = "n",
         frame.plot = T, main = c("Situation 1 \n (Overall Strategy)","Situation 1 \n (Arm-based Strategy)")[i], cex.main = 1.75, cex.lab = 1.5, cex.axis=1.5)
    axis(side = 1, at = c(seq(0.1,2,0.1),seq(2.25,5,0.25)), labels = NA, cex.axis=1.5)
    axis(side = 1, cex.axis=1.5, tcl = -0.5)
    for(k in list(c(1,2,4,3),c(1,2,4,5,3))[[i]]){
      points(dplyr::select(filter(data_1, prior == Type_prior[k], strategy == c("Overall","Arm")[i]), c(heterogeneity,post_mean1)), 
             col = list(col_vector1,col_vector2)[[i]][k], pch = c(17,6,5,15,16)[k], cex = 1, lwd = 1)
      lines(dplyr::select(filter(data_1, prior == Type_prior[k], strategy == c("Overall","Arm")[i]), c(heterogeneity,post_mean1)), 
            col = list(col_vector1,col_vector2)[[i]][k], lwd = 1)
      #points(dplyr::select(filter(data_1, prior == Type_prior[k], strategy == c("Overall","Arm")[i]), c(heterogeneity,post_variance1)), col = list(col_vector1,col_vector2)[[i]][k], pch = 0, cex = 0.75, lwd = 1)
      #lines(dplyr::select(filter(data_1, prior == Type_prior[k], strategy == c("Overall","Arm")[i]), c(heterogeneity,post_variance1)), col = list(col_vector1,col_vector2)[[i]][k], lwd = 1)
    }
    abline(v=1,lty="dashed", col="grey")
    legend("bottomright", legend = list(legend_lable1,legend_lable2)[[i]], horiz = F, ncol =2,
           lty = 1, pch = c(17,6,5,15,16), col = list(c(3,"violet",4,7),c(3,"violet",4,7,2))[[i]], cex = 1.25, bty = "n", xpd = NA)
  }
  for (i in 1:2){
    plot(NA, xlim = c(0.1,5), ylim = c(0,1), xlab = "Variance-specific prior-data conflict", ylab = "Posterior mean", xaxt = "n",
         frame.plot = T, main = c("Situation 2 \n (Overall Strategy)","Situation 2 \n (Arm-based Strategy)")[i], cex.main = 1.75, cex.lab = 1.5, cex.axis=1.5)
    axis(side = 1, at = c(seq(0.1,2,0.1),seq(2.25,5,0.25)), labels = NA, cex.axis=1.5)
    axis(side = 1, cex.axis=1.5, tcl = -0.5)
    for(k in list(c(1,2,4,3),c(1,2,4,5,3))[[i]]){
      points(dplyr::select(filter(data_2, prior == Type_prior[k], strategy == c("Overall","Arm")[i]), c(heterogeneity,post_mean1)), 
             col = list(col_vector1,col_vector2)[[i]][k], pch = c(17,6,5,15,16)[k], cex = 1, lwd = 1)
      lines(dplyr::select(filter(data_2, prior == Type_prior[k], strategy == c("Overall","Arm")[i]), c(heterogeneity,post_mean1)), 
            col = list(col_vector1,col_vector2)[[i]][k], lwd = 1)
      #points(dplyr::select(filter(data_2, prior == Type_prior[k], strategy == c("Overall","Arm")[i]), c(heterogeneity,post_variance1)), col = list(col_vector1,col_vector2)[[i]][k], pch = 0, cex = 0.75, lwd = 1)
      #lines(dplyr::select(filter(data_2, prior == Type_prior[k], strategy == c("Overall","Arm")[i]), c(heterogeneity,post_variance1)), col = list(col_vector1,col_vector2)[[i]][k], lwd = 1)
    }
    abline(v=1,lty="dashed", col="grey")
    legend("bottomright", legend = list(legend_lable1,legend_lable2)[[i]], horiz = F, ncol =2,
           lty = 1, pch = c(17,6,5,15,16), col = list(c(3,"violet",4,7),c(3,"violet",4,7,2))[[i]], cex = 1.25, bty = "n", xpd = NA)
  }

  dev.off()
}


################################################
### Design 1: Figure 2,3 & Figure S5,6,20,21
################################################

### load simulation results
data_H0 = lapply(as.list(LETTERS[1:5]), function(z){readRDS(paste0(wd,"/Sim_rst_for_SingleEC/simu_1_H0_",z,".rds")) %>% mutate(scenario = z)}) %>% rbindlist()
data_H1 = lapply(as.list(LETTERS[1:5]), function(z){readRDS(paste0(wd,"/Sim_rst_for_SingleEC/simu_1_H1_",z,".rds")) %>% mutate(scenario = z)}) %>% rbindlist()

plot_simu1(data = list(data_H0,data_H1), select = c(1,3,5), Type_prior = c("MPP1","MPP2","MBPP","MCBPP","rMCBPP","CP","Pooled","No-borrowing"), 
           yrange1 =  c(-0.01,0.2), yrange2 = c(0.3,1), yrange3 = c(-0.5,1), 
           col_vector = c(alpha(c(3,"violet",4,7,2,5),0.9),8,1), file = "Figure_2")
plot_simu2(data = list(data_H0,data_H1), select = c(1,3,5), Type_prior = c("MPP1","MPP2","MBPP","MCBPP","rMCBPP","CP","Pooled","No-borrowing"),
           index = c("relative_bias","relative_mse"),
           yrange = cbind(c(-0.1,0.3),c(-0.04,0.1)),
           col_vector = c(alpha(c(3,"violet",4,7,2,5),0.9),8,1), file = "Figure_3")

plot_simu1(data = list(data_H0,data_H1), select = c(2,4), Type_prior = c("MPP1","MPP2","MBPP","MCBPP","rMCBPP","CP","Pooled","No-borrowing"), 
           yrange1 =  c(-0.01,0.2), yrange2 = c(0.3,1), yrange3 = c(-0.5,1), 
           col_vector = c(alpha(c(3,"violet",4,7,2,5),0.9),8,1), file = "Figure_S5")
plot_simu2(data = list(data_H0,data_H1), select = c(2,4), Type_prior = c("MPP1","MPP2","MBPP","MCBPP","rMCBPP","CP","Pooled","No-borrowing"),
           index = c("relative_bias","relative_mse"),
           yrange = cbind(c(-0.1,0.3),c(-0.04,0.1)),
           col_vector = c(alpha(c(3,"violet",4,7,2,5),0.9),8,1), file = "Figure_S6")

### load simulation results
data_H0_adaptive1 = lapply(as.list(LETTERS[1:5]), function(z){readRDS(paste0(wd,"/Sim_rst_for_SingleEC/simu_1_H0_",z,"_adapt.rds")) %>% mutate(scenario = z)}) %>% rbindlist() 
data_H1_adaptive1 = lapply(as.list(LETTERS[1:5]), function(z){readRDS(paste0(wd,"/Sim_rst_for_SingleEC/simu_1_H1_",z,"_adapt.rds")) %>% mutate(scenario = z)}) %>% rbindlist() 
data_H0_adaptive2 = lapply(as.list(LETTERS[1:5]), function(z){readRDS(paste0(wd,"/Sim_rst_for_SingleEC/simu_2_H0_",z,"_adapt.rds")) %>% mutate(scenario = z)}) %>% rbindlist() 
data_H1_adaptive2 = lapply(as.list(LETTERS[1:5]), function(z){readRDS(paste0(wd,"/Sim_rst_for_SingleEC/simu_2_H1_",z,"_adapt.rds")) %>% mutate(scenario = z)}) %>% rbindlist() 

data_H0_adaptive = rbind(filter(data_H0_adaptive1,prior!="No-borrowing"),filter(data_H0_adaptive2,prior!="No-borrowing"),
                         filter(data_H0, prior %in% c("MPP1","rMCBPP","No-borrowing")))
data_H1_adaptive = rbind(filter(data_H1_adaptive1,prior!="No-borrowing"),filter(data_H1_adaptive2,prior!="No-borrowing"),
                         filter(data_H1, prior %in% c("MPP1","rMCBPP","No-borrowing")))

plot_simu1(data = list(data_H0_adaptive,data_H1_adaptive), Type_prior = c("MPP1","adaptive_strict","adaptive_mild","rMCBPP","No-borrowing"), 
           yrange1 =  c(0,0.18), yrange2 = c(0.4,0.9), yrange3 = c(-0.2,0.6), 
           col_vector = c(alpha(c(3,7,4,2),0.8),1), file = "Figure_S20")
plot_simu2(data = list(data_H0_adaptive,data_H1_adaptive), Type_prior = c("MPP1","adaptive_strict","adaptive_mild","rMCBPP","No-borrowing"),  
           index = c("relative_bias","relative_mse"),
           yrange = cbind(c(-0.05,0.25),c(-0.035,0.06)),
           col_vector = c(alpha(c(3,7,4,2),0.8),1), file = "Figure_S21")

################################################
### Design 2: Figure S7-8
################################################


### load simulation results
data_H0 = lapply(as.list(LETTERS[1:5]), function(z){readRDS(paste0(wd,"/Sim_rst_for_SingleEC/simu_2_H0_",z,".rds")) %>% mutate(scenario = z)}) %>% rbindlist()
data_H1 = lapply(as.list(LETTERS[1:5]), function(z){readRDS(paste0(wd,"/Sim_rst_for_SingleEC/simu_2_H1_",z,".rds")) %>% mutate(scenario = z)}) %>% rbindlist()

plot_simu1(data = list(data_H0,data_H1), Type_prior = c("MPP1","MPP2","MBPP","MCBPP","rMCBPP","CP","Pooled","No-borrowing"), 
           yrange1 = c(-0.02,0.25), yrange2 = c(0,1), yrange3 = c(-0.6,1),
           col_vector = c(alpha(c(3,"violet",7,4,2,5),0.9),8,1), file = "Figure_S7")
plot_simu2(data = list(data_H0,data_H1), Type_prior = c("MPP1","MPP2","MBPP","MCBPP","rMCBPP","CP","Pooled","No-borrowing"), 
           index = c("relative_bias","relative_mse"), 
           yrange = cbind(c(-0.3,1.0),c(-0.5,1.5)),
           col_vector = c(alpha(c(3,"violet",7,4,2,5),0.9),8,1), file = "Figure_S8")

################################################
### Design 3: Figure S9-10
################################################

### load simulation results
data_H0 = lapply(as.list(LETTERS[1:5]), function(z){readRDS(paste0(wd,"/Sim_rst_for_SingleEC/simu_3_H0_",z,".rds")) %>% mutate(scenario = z)}) %>% rbindlist()
data_H1 = lapply(as.list(LETTERS[1:5]), function(z){readRDS(paste0(wd,"/Sim_rst_for_SingleEC/simu_3_H1_",z,".rds")) %>% mutate(scenario = z)}) %>% rbindlist()

plot_simu1(data = list(data_H0,data_H1), Type_prior = c("MPP1","MPP2","MBPP","MCBPP","rMCBPP","CP","Pooled","No-borrowing"), 
           yrange1 = c(-0.01,0.1), yrange2 = c(0.5,1), yrange3 = c(-0.5,1), 
           col_vector = c(alpha(c(3,"violet",7,4,2,5),0.9),8,1), file = "Figure_S9")
plot_simu2(data = list(data_H0,data_H1), Type_prior = c("MPP1","MPP2","MBPP","MCBPP","rMCBPP","CP","Pooled","No-borrowing"), 
           index = c("relative_bias","relative_mse"), 
           yrange = cbind(c(-0.07,0.2),c(-0.015,0.04)),
           col_vector = c(alpha(c(3,"violet",7,4,2,5),0.9),8,1), file = "Figure_S10")

################################################
### Design 4: Figure S11-12
################################################

### load simulation results
data_H0 = lapply(as.list(LETTERS[1:5]), function(z){readRDS(paste0(wd,"/Sim_rst_for_SingleEC/simu_4_H0_",z,".rds")) %>% mutate(scenario = z)}) %>% rbindlist()
data_H1 = lapply(as.list(LETTERS[1:5]), function(z){readRDS(paste0(wd,"/Sim_rst_for_SingleEC/simu_4_H1_",z,".rds")) %>% mutate(scenario = z)}) %>% rbindlist()

plot_simu1(data = list(data_H0,data_H1), Type_prior = c("MPP1","MPP2","MBPP","MCBPP","rMCBPP","CP","Pooled","No-borrowing"), 
           yrange1 = c(-0.05,0.5), yrange2 = c(0,1), yrange3 = c(-0.2,1), 
           col_vector = c(alpha(c(3,"violet",7,4,2,5),0.9),8,1), file = "Figure_S11")
plot_simu2(data = list(data_H0,data_H1), Type_prior = c("MPP1","MPP2","MBPP","MCBPP","rMCBPP","CP","Pooled","No-borrowing"), 
           index = c("relative_bias","relative_mse"), 
           yrange = cbind(c(-0.2,0.6),c(-0.12,0.3)),
           col_vector = c(alpha(c(3,"violet",7,4,2,5),0.9),8,1), file = "Figure_S12")

################################################
### Design 5: Figure 4 
################################################

data_H0_A = readRDS(paste0(wd,"/Sim_rst_for_MultipleECs/simu_5_F_H0.rds")) 
data_H0_B = readRDS(paste0(wd,"/Sim_rst_for_MultipleECs/simu_5_I_H0.rds")) 
data_H1_A = readRDS(paste0(wd,"/Sim_rst_for_MultipleECs/simu_5_F_H1.rds")) 
data_H1_B = readRDS(paste0(wd,"/Sim_rst_for_MultipleECs/simu_5_I_H1.rds")) 

Type_prior1 <- c("MPP1","MPP2","MBPP","MCBPP","CP")
Type_prior2 <- c("MPP1","MPP2","MBPP","MCBPP","rMCBPP")

### Figure 4 
{
  pdf(width = 4.75*2, height = 4.75*2, file = "Figure_4.pdf")
  par(mfrow = c(2,2), mar = c(1.5,1.5,1.5,1.5), omi = c(0.4,0.4,0.25,0.25), xpd = F)
  for (i in 1:2){
    plot(NA, xlim = c(-2,2), ylim = c(-0.3,0.8), xlab = "", xaxt = "n", ylab = "", frame.plot = T, cex.lab = 1.5, cex.axis=1.25)
    axis(side = 1, cex.axis = 1.25)
    axis(side = 1, at = seq(-2,2,0.1), labels = NA)
    mtext(c("Overall strategy","Arm-based strategy")[i], cex = 1.5, side = 3, line = 1, font = 2)
    mtext(c("","Scenario F")[i], cex = 1.5, side = 4, line = 1.25, font = 2)
    mtext(c("Relative proportion of prior ESS","")[i], cex = 1.5, side = 2, line = 2.5)
    for (k in list(c(5,1,2,4,3),c(1,2,4,5,3))[[i]]) {
      points(dplyr::select(filter(data_H1_A, strategy == c("Overall","Arm")[i], prior == list(Type_prior1,Type_prior2)[[i]][k]), c(heterogeneity,relative_Proportion)), 
             col = list(alpha(c(3,"violet",4,7,5),0.9),alpha(c(3,"violet",4,7,2),0.9))[[i]][k], 
             pch = list(c(17,6,5,15,4),c(17,6,5,15,16))[[i]][k], cex = 1.25)
      lines(dplyr::select(filter(data_H1_A, strategy == c("Overall","Arm")[i], prior == list(Type_prior1,Type_prior2)[[i]][k]), c(heterogeneity,relative_Proportion)), 
            col = list(alpha(c(3,"violet",4,7,5),0.9),alpha(c(3,"violet",4,7,2),0.9))[[i]][k], lwd = 1.25)
    }
    points(dplyr::select(filter(data_H1_A, strategy == c("Overall","Arm")[i], prior == "No-borrowing"), c(heterogeneity,relative_Proportion)), col = "black", pch = NA, cex = 1)
    lines(dplyr::select(filter(data_H1_A, strategy == c("Overall","Arm")[i], prior == "No-borrowing"), c(heterogeneity,relative_Proportion)), col = "black", lwd = 1.5)
    abline(v=0, lty = "dashed", col = "grey80")
    legend("topright", legend = list(Type_prior1,Type_prior2)[[i]], horiz = F, 
           lty = 1, pch = list(c(17,6,5,15,4),c(17,6,5,15,16))[[i]], col = list(alpha(c(3,"violet",4,7,5),1),alpha(c(3,"violet",4,7,2),1))[[i]], cex = 1.25, bty = "n", xpd = NA)
    legend("bottomleft", legend = "No borrowing", horiz = F, 
           lty = 1, pch = NA, col = 1, cex = 1.25, bty = "n", xpd = NA)
    
  }
  for (i in 1:2){
    plot(NA, xlim = c(-2,2), ylim = c(-0.25,0.5), xaxt = "n", frame.plot = T, main = "", cex.main = 1.5, cex.lab = 1.5, cex.axis=1.25)
    axis(side = 1, cex.axis = 1.25)
    axis(side = 1, at = seq(-2,2,0.1), labels = NA)
    mtext(c("","Scenario I")[i], cex = 1.5, side = 4, line = 1.5, font = 2)
    mtext(c("Relative proportion of prior ESS","")[i], cex = 1.5, side = 2, line = 2.5)
    mtext("Population mean heterogeneity (D)", cex = 1.5, side = 1, line = 2.5)
    for (k in list(c(5,1,2,4,3),c(1,2,4,5,3))[[i]]) {
      points(dplyr::select(filter(data_H1_B, strategy == c("Overall","Arm")[i], prior == list(Type_prior1,Type_prior2)[[i]][k]), c(heterogeneity,relative_Proportion)), 
             col = list(alpha(c(3,"violet",4,7,5),0.9),alpha(c(3,"violet",4,7,2),0.9))[[i]][k], 
             pch =  list(c(17,6,5,15,4),c(17,6,5,15,16))[[i]][k], cex = 1.25)
      lines(dplyr::select(filter(data_H1_B, strategy == c("Overall","Arm")[i], prior == list(Type_prior1,Type_prior2)[[i]][k]), c(heterogeneity,relative_Proportion)), 
            col = list(alpha(c(3,"violet",4,7,5),0.9),alpha(c(3,"violet",4,7,2),0.9))[[i]][k], lwd = 1.25)
    }
    points(dplyr::select(filter(data_H1_B, strategy == c("Overall","Arm")[i], prior == "No-borrowing"), c(heterogeneity,relative_Proportion)), col = "black", pch = NA, cex = 1)
    lines(dplyr::select(filter(data_H1_B, strategy == c("Overall","Arm")[i], prior == "No-borrowing"), c(heterogeneity,relative_Proportion)), col = "black", lwd = 1.5)
    abline(v=c(-1,0,1), lty = "dashed", col = "grey80")
    legend("topright", legend = list(Type_prior1,Type_prior2)[[i]], horiz = F, 
           lty = 1, pch = list(c(17,6,5,15,4),c(17,6,5,15,16))[[i]], col = list(alpha(c(3,"violet",4,7,5),1),alpha(c(3,"violet",4,7,2),1))[[i]], cex = 1.25, bty = "n", xpd = NA)
    legend("bottomleft", legend = "No borrowing", horiz = F, 
           lty = 1, pch = NA, col = 1, cex = 1.25, bty = "n", xpd = NA)
  }
  dev.off()
}

################################################
### Design 5: Figure S13-19 
################################################

data_H0 = lapply(as.list(LETTERS[6:11]), function(z){readRDS(paste0(wd,"/Sim_rst_for_MultipleECs/simu_5_",z,"_H0.rds")) %>% mutate(scenario = z)}) %>% rbindlist()
data_H1 = lapply(as.list(LETTERS[6:11]), function(z){readRDS(paste0(wd,"/Sim_rst_for_MultipleECs/simu_5_",z,"_H1.rds")) %>% mutate(scenario = z)}) %>% rbindlist()

Type_prior <- c("MPP1","MBPP","rMCBPP","Pooled","MPP2","MCBPP","CP")

plot_simu3(data = data_H1, Type_prior = Type_prior, index = "relative_Proportion", ylab = "Relative proportion of prior ESS", yrange = c(-0.2,1),  file = "Figure_S13")
plot_simu3(data = data_H0, Type_prior = Type_prior, index = "Type_I_error_calibrated", ylab = "Calibrated Type I error rate (H0)", yrange = c(0,0.25),  file = "Figure_S14")
plot_simu3(data = data_H1, Type_prior = Type_prior, index = "Power_calibrated", ylab = "Calibrated power (H1)", yrange = c(0,1),  file = "Figure_S15")
plot_simu3(data = data_H0, Type_prior = Type_prior, index = "relative_bias", ylab = "Relative bias (H0)", yrange = c(-0.01,0.3),  file = "Figure_S16")
plot_simu3(data = data_H1, Type_prior = Type_prior, index = "relative_bias", ylab = "Relative bias (H1)", yrange = c(-0.01,0.3),  file = "Figure_S17")
plot_simu3(data = data_H0, Type_prior = Type_prior, index = "relative_mse", ylab = "Relative MSE (H0)", yrange = c(-0.03,0.08),  file = "Figure_S18")
plot_simu3(data = data_H1, Type_prior = Type_prior, index = "relative_mse", ylab = "Relative MSE (H1)", yrange = c(-0.03,0.08),  file = "Figure_S19")

################################################
### Application results: Figure 5 (varying sample SD ratio)
################################################

data = readRDS(paste0(wd,"/Application_rst/rst_plot_varying_ratio.rds"))
Type_prior = c("MPP1","MPP2","MBPP","MCBPP","rMCBPP","CP","Pooled","No-borrowing")

pdf(width = 4*3, height = 5*1, file = "Figure_5.pdf")
par(mfrow = c(1,3), mar = c(5,6,6,2))
for(i in 1:3){
  plot(NA, xlim = c(0.65,1.5), ylim = list(c(0,300),c(0,0.12),c(-0.35,0))[[i]], xlab = expression(paste("Varying sample SD ratio, ",s[h]/s[c])), 
       ylab = c("Relative prior ESS","Change in the posterior estimation of ATE \n from No borrowing","Change in the width of 95% CI \n from No borrowing")[i], 
       frame.plot = T, cex.lab = 1.75, cex.axis=1.5)
  axis(side = 1, at = round(seq(0.65,1.5,0.05),2), labels = NA, cex.axis=1)
  
  for(k in c(8,7,6,1,2,4,5,3)){
    points(dplyr::select(filter(data, prior == Type_prior[k]), c(SD_hetero,c("relative_ESS","relative_Est","relative_CI")[i])), 
           col = c(alpha(c(3,"violet",4,7,2,5),0.9),8,1)[k], pch = c(17,6,5,15,16,3,4,NA)[k], cex = 1.5, lwd = 1)
    lines(dplyr::select(filter(data, prior == Type_prior[k]), c(SD_hetero,c("relative_ESS","relative_Est","relative_CI")[i])), 
          col = c(alpha(c(3,"violet",4,7,2,5),0.9),8,1)[k], pch = c(17,6,5,15,16,3,4,NA)[k], lwd = 1)
  }
  abline(v=1, lty="dashed",col="grey80")
  legend(c("topright","top","topleft")[i], legend = list(Type_prior[1:3],Type_prior[4:6],Type_prior[7:8])[[i]], ncol = 1, lwd = 1, inset = c(c(-0.1,0,-0.25)[i],-0.25), xpd = T,
         lty = 1, pch = list(c(17,6,5),c(15,16,3),c(4,NA))[[i]], col = list(alpha(c(3,"violet",4),0.75),alpha(c(7,2,5),0.75),c(8,1))[[i]], cex = 1.5, bty = "n")
}
dev.off()
