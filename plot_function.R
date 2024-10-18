################################################
###### plot function
################################################

plot_num <- function(data,        # numerical results
                     Type_prior,  # methods
                     main,        # main of plot 
                     col_vector,  # vector of colors
                     file         # figure name
){
  
  pdf(width = 5*3, height = 5*2, file = paste0(file,".pdf"))
  par(mfcol = c(2,3), mar = c(6,5,4,2))
  
  for (i in 1:length(Type_prior)){
    
    plot(NA, xlim = c(-2,2), ylim = c(0,1), xlab = "", ylab = "Posterior mean", 
         frame.plot = T, main = main[i], font.main = 2, cex.main = 2.75, cex.lab = 2.25, cex.axis=1.75)
    axis(side = 1, at = seq(-2,2,0.1), labels = NA, cex.axis=1.75)
    mtext(expression(paste("Mean-specific prior-data conflict, ",bar(D))), side = 1, line = 4, cex = 1.5)
    
    for(k in 1:5){
      points(dplyr::select(filter(data[[1]], scenario == c("A","B","C","D","E")[k], prior == Type_prior[i]), c(LOC_hetero,post_mean)), 
             col = col_vector[k], pch = 2, cex = 1)
      lines(dplyr::select(filter(data[[1]], scenario == c("A","B","C","D","E")[k], prior == Type_prior[i]), c(LOC_hetero,post_mean)), 
            col = col_vector[k], lwd = 1)
    }
    abline(v=0,lwd=1,col="grey80",lty="dashed")
    
    legend(1,1, legend = c("0.4","0.8","1.0","1.2","2.5"), horiz = F, ncol = 1, title = expression(hat(R)), lwd = 1.5,
           lty = 1, pch = 2, col = col_vector, y.intersp = 1, cex = 2, bty = "n", xpd=NA)
  }
  
  {
    plot(NA, xlim = c(0,5), ylim = c(0,1), xlab = "", ylab = "Posterior mean", 
         frame.plot = T, main = expression(paste("F. All methods under (",bar(D)==0,")")), 
         font.main = 2, cex.main = 2.75, cex.lab = 2.25, cex.axis=1.75)
    axis(side = 1, at = round(c(seq(0.1,2,0.1),seq(2.25,5,0.25)),2), labels = NA, cex.axis=1.75)
    mtext(expression(paste("Variance-specific prior-data conflict, ",hat(R))), side = 1, line = 4, cex = 1.5)
    
    for(k in c(1,2,4,5,3)){
      points(dplyr::select(filter(data[[2]], prior == Type_prior[k]), c(VAR_hetero,post_mean)), col = alpha(c(3,"violet",4,7,2),0.8)[k], 
             pch = c(17,6,5,15,16)[k], cex = 1.5)
      lines(dplyr::select(filter(data[[2]], prior == Type_prior[k]), c(VAR_hetero,post_mean)), col = alpha(c(3,"violet",4,7,2),0.8)[k], lwd = 1)
    }
    abline(v=1,lwd=1,col="grey80",lty="dashed")
    legend("topleft", legend = c(expression(paste("MPP1"," (",delta,")")),expression(paste("MPP2"," (",delta,")")),
                                 expression(paste("rMCBPP"," (",delta[1],"*)"))), horiz = F, ncol = 1, lwd = 1,
           lty = 1, pch = c(17,6,16), col=c(3,"violet",2), cex = 1.75, bty = "n")
    legend("topright", legend = c(expression(paste("MBPP"," (",delta[1],")")), expression(paste("MCBPP"," (",delta[1],"*)"))), horiz = F, ncol = 1, lwd = 1,
           lty = 1, pch = c(5,15), col=c(4,7), cex = 1.75, bty = "n")
    
  }
  
  dev.off()
}
plot_simu1 <- function(data,      # list of simulation results 
                       select = c(1:5), # select the population variance heterogeneity
                       Type_prior,# methods
                       yrange1,   # range of y-axis for calibrated Type I error rate
                       yrange2,   # range of y-axis for calibrated power
                       yrange3,   # range of y-axis for relative proportion of prior ESS
                       col_vector,# vector of colors
                       file       # figure name
){
  
  pdf(width = 5.3*length(select), height = 5.3*3, file = paste0(file,".pdf"))
  
  par(mfrow = c(3,length(select)), mar=c(2,2,2,2), omi = c(0.42,0.4,0.2,0))
  for (i in select){
    plot(NA, xlim = c(-2,2), ylim = yrange1, xlab = "", xaxt = "n", ylab = "", frame.plot = T, cex.axis = 2)
    axis(side = 1, at = seq(-2,2,0.1), labels = NA)
    axis(side = 1, at = seq(-2,2,1), tcl = -0.75, padj = 0.5, cex.axis = 2)
    mtext(c("R = 0.4 ","R = 0.8 ","R = 1.0 ","R = 1.2 ","R = 2.5 ")[i], side = 3, line = 0.5, cex = 2.25, font = 2)
    for(k in c(8,7,6,1,2,4,5,3)){
      points(dplyr::select(filter(data[[1]], scenario == c("E","C","A","B","D")[i], prior == Type_prior[k]), c(heterogeneity,Type_I_error_calibrated)), col = col_vector[k], 
             pch = c(17,6,5,15,16,3,4,NA)[k], cex = 1.75, lwd = 1.25)
      lines(dplyr::select(filter(data[[1]], scenario == c("E","C","A","B","D")[i], prior == Type_prior[k]), c(heterogeneity,Type_I_error_calibrated)), col = col_vector[k], lwd = 1.5)
    }
    abline(v=0,lty="dashed",col="grey60",lwd=1.25)
    #legend("topleft", legend = c("R = 0.4 ","R = 0.8 ","R = 1.0 ","R = 1.2 ","R = 2.5 ")[i], cex = 3, bty = "n", text.font = 2)
    if(i==max(select)){
      if(file!="Figure_S20"){
        legend("topright", legend = Type_prior[1:6], horiz = F, ncol = 1, lwd = 1,
               lty = 1, pch = c(17,6,5,15,16,3,4,NA)[1:6], col=col_vector[1:6], y.intersp = 1, cex = 2, bty = "n")
        legend("bottomleft", legend = Type_prior[7:8], horiz = F, ncol = 1, lwd = 1,
               lty = 1, pch = c(17,6,5,15,16,3,4,NA)[7:8], col=col_vector[7:8], y.intersp = 1, cex = 2, bty = "n")
      } else
      {
        legend("topright", legend = Type_prior[1:4], horiz = F, ncol = 1, lwd = 1,
               lty = 1, pch = c(17,6,5,15,16,3,4,NA)[1:4], col=col_vector[1:6], y.intersp = 1, cex = 2, bty = "n")
        legend("bottomleft", legend = Type_prior[5], horiz = F, ncol = 1, lwd = 1,
               lty = 1, col=1, y.intersp = 1, cex = 2, bty = "n")
      }
    }
    if(i==select[1])mtext("Calibrated Type I error rate (H0)", side = 2, line = 3, cex = 1.75)
  }
  
  
  for (i in select){
    plot(NA, xlim = c(-2,2), ylim = yrange2, #main = c("R = 0.4 ","R = 0.8 ","R = 1.0 ","R = 1.2 ","R = 2.5 ")[i], cex.main = 3.5,
         xlab = "", xaxt = "n", ylab = "", frame.plot = T, cex.axis = 2)
    axis(side = 1, at = seq(-2,2,0.1), labels = NA)
    axis(side = 1, at = seq(-2,2,1), tcl = -0.75, padj = 0.5, cex.axis = 2)
    for(k in c(8,7,6,1,2,4,5,3)){
      points(dplyr::select(filter(data[[2]], scenario == c("E","C","A","B","D")[i], prior == Type_prior[k]), c(heterogeneity,Power_calibrated)), col = col_vector[k], 
             pch = c(17,6,5,15,16,3,4,NA)[k], cex = 1.75, lwd = 1.25)
      lines(dplyr::select(filter(data[[2]], scenario == c("E","C","A","B","D")[i], prior == Type_prior[k]), c(heterogeneity,Power_calibrated)), col = col_vector[k], lwd = 1.5)
    }
    abline(v=0,lty="dashed",col="grey60",lwd=1.25)
    if(i==max(select)){
      if(file!="Figure_S20"){
        legend("topright", legend = Type_prior[1:6], horiz = F, ncol = 1, lwd = 1,
               lty = 1, pch = c(17,6,5,15,16,3,4,NA)[1:6], col=col_vector[1:6], y.intersp = 1, cex = 2, bty = "n")
        legend("bottomleft", legend = Type_prior[7:8], horiz = F, ncol = 1, lwd = 1,
               lty = 1, pch = c(17,6,5,15,16,3,4,NA)[7:8], col=col_vector[7:8], y.intersp = 1, cex = 2, bty = "n")
      } else
      {
        legend("topright", legend = Type_prior[1:4], horiz = F, ncol = 1, lwd = 1,
               lty = 1, pch = c(17,6,5,15,16,3,4,NA)[1:4], col=col_vector[1:6], y.intersp = 1, cex = 2, bty = "n")
        legend("bottomleft", legend = Type_prior[5], horiz = F, ncol = 1, lwd = 1,
               lty = 1, col=1, y.intersp = 1, cex = 2, bty = "n")
      }
    }
    if(i==select[1])mtext("Calibrated power (H1)", side = 2, line = 3, cex = 1.75)
  }
  
  
  for (i in select){
    plot(NA, xlim = c(-2,2), ylim = yrange3, #main = c("R = 0.4 ","R = 0.8 ","R = 1.0 ","R = 1.2 ","R = 2.5 ")[i], cex.main = 3.5, 
         frame.plot = T, xlab = "", xaxt = "n", ylab = "", cex.lab = 2.5, cex.axis = 2, font.lab = 2, font.main = 2)
    axis(side = 1, at = seq(-2,2,0.1), labels = NA)
    axis(side = 1, at = seq(-2,2,1), tcl = -0.75, padj = 0.5, cex.axis = 2)
    for(k in c(8,7,6,1,2,4,5,3)){
      points(dplyr::select(filter(data[[2]], scenario == c("E","C","A","B","D")[i], prior == Type_prior[k]), c(heterogeneity,relative_Proportion)), col = col_vector[k], 
             pch = c(17,6,5,15,16,3,4,NA)[k], cex = 1.75, lwd = 1.25)
      lines(dplyr::select(filter(data[[2]], scenario == c("E","C","A","B","D")[i], prior == Type_prior[k]), c(heterogeneity,relative_Proportion)), col = col_vector[k], lwd = 1.5)
    }
    abline(v=0,lty="dashed",col="grey60",lwd=1.25)
    if(i==max(select)){
      if(file!="Figure_S20"){
        legend("topright", legend = Type_prior[1:4], horiz = F, ncol = 1, lwd = 1,
               lty = 1, pch = c(17,6,5,15,16,3,4,NA)[1:4], col=col_vector[1:4], y.intersp = 1, cex = 2, bty = "n")
        legend("topleft", legend = Type_prior[5:8], horiz = F, ncol = 1, lwd = 1,
               lty = 1, pch = c(17,6,5,15,16,3,4,NA)[5:8], col=col_vector[5:8], y.intersp = 1, cex = 2, bty = "n")
      } else
      {
        legend("topright", legend = Type_prior[1:4], horiz = F, ncol = 1, lwd = 1,
               lty = 1, pch = c(17,6,5,15,16,3,4,NA)[1:4], col=col_vector[1:6], y.intersp = 1, cex = 2, bty = "n")
        legend("bottomleft", legend = Type_prior[5], horiz = F, ncol = 1, lwd = 1,
               lty = 1, col=1, y.intersp = 1, cex = 2, bty = "n")
      }
    }
    if(i==select[1])mtext("Relative proportion of prior ESS", side = 2, line = 3, cex = 1.75)
    mtext("Population mean heterogeneity (D)", side = 1, line = 4, cex = 1.75)
  }
  
  dev.off()
}
plot_simu2 <- function(data,      # list of simulation results
                       select = c(1:5), # select the population variance heterogeneity
                       Type_prior,# methods
                       index = c("relative_bias","relative_mse"),
                       yrange,    # range of y-axis for relative bias and MSE
                       col_vector,# vector of colors
                       file       # figure name
){
  
  pdf(width = 6*4, height = 6*length(select), file = paste0(file,".pdf"))
  par(mfcol = c(length(select),4))
  for(q in 1:length(index)){
    for (j in 1:2){
      if(j==1){par(mar = c(6,5,1,1))}else{par(mar = c(6,3,1,3))}
      for (i in select){
        plot(NA, xlim = c(-2,2), ylim = yrange[,q], xlab = "", ylab = "", xaxt = "n", yaxt = "n",
             frame.plot = T, main = "", cex.lab = 2.5, font.lab = 2)
        axis(side = 1, tcl = -0.75, padj = 0.35, cex.axis = 2.25)
        axis(side = 1, at = seq(-2,2,0.1), labels = NA)
        axis(side = 2, tcl = -0.75, padj = -0.15, cex.axis = 2.25)
        for(k in c(8,7,6,1,2,4,5,3)){
          points(dplyr::select(filter(data[[j]], scenario == c("E","C","A","B","D")[i], prior == Type_prior[k]), c(heterogeneity,index[q])), col = col_vector[k], 
                 pch = c(17,6,5,15,16,3,4,NA)[k], cex = 1.75, lwd = 1.25)
          lines(dplyr::select(filter(data[[j]], scenario == c("E","C","A","B","D")[i], prior == Type_prior[k]), c(heterogeneity,index[q])), col = col_vector[k], lwd = 1.5)
        }
        abline(v=0, lty = "dashed", col = "grey80")
        legend("topright", legend = c("R = 0.4 ","R = 0.8 ","R = 1.0 ","R = 1.2 ","R = 2.5 ")[i], 
               horiz = F, ncol = 1, cex = 3, bty = "n", text.font = 2)
        legend("topleft", legend = c("H0 ","H1")[j], 
               horiz = F, ncol = 1, cex = 3, bty = "n", text.font = 2)
        
        {
          if(file!="Figure_S21"){
            legend("bottomright", legend = Type_prior[1:4], horiz = F, ncol = 1,lwd = 1, 
                   lty = 1, pch = c(17,6,5,15,16,3,44,NA)[1:4], col=col_vector[1:4], y.intersp = 1, cex = 2, bty = "n")
            legend("bottomleft", legend = Type_prior[5:8], horiz = F, ncol = 1, lwd = 1, 
                   lty = 1, pch = c(17,6,5,15,16,3,4,NA)[5:8], col=col_vector[5:8], y.intersp = 1, cex = 2, bty = "n")
          } else
          {
            legend("bottomright", legend = Type_prior[1:3], horiz = F, ncol = 1, lwd = 1, 
                   lty = 1, pch = c(17,6,5,15,16,3,4,NA)[1:3], col=col_vector[1:3], y.intersp = 1, cex = 2, bty = "n")
            legend("bottomleft", legend = Type_prior[4:5], horiz = F, ncol = 1, lwd = 1, 
                   lty = 1, pch = c(17,6,5,15,16,3,4,NA)[c(4,8)], col=col_vector[4:5], y.intersp = 1, cex = 2, bty = "n")
          }
        }
        
        if(j==1){
          mtext(c("Relative posterior bias","Relative posterior MSE")[q], side = 2, at = (sum(yrange[,q]))/2, line = 3, cex = 2)
          mtext("Population mean heterogeneity (D)", side = 1, at = 2.5, line = 3.75, cex = 2)
        }
      }
    }
  }
  
  dev.off()
}
plot_simu3 <- function(data,      # list of simulation results
                       Type_prior,# methods
                       index,     # name of operating characteristic
                       ylab,      # lab of y-axis
                       yrange,    # range of y-axis
                       file       # figure name
){
  
  pdf(width = 5.5*4, height = 5.5*4, file = paste0(file,".pdf"))
  par(mfrow = c(4,4), mar = c(5,5,3,2),xpd=F)
  
  for (i in 1:7){
    plot(NA, xlim = c(-2,2), ylim = yrange, xlab = "Population mean heterogeneity", yaxt="n", xaxt = "n", ylab = ylab, 
         frame.plot = T, main = paste0("A.",Type_prior[i]), cex.main = 3.5, cex.lab = 2.5, cex.axis=1.5)
    axis(side = 1, tcl = -0.75, padj = 0.35, cex.axis = 1.5)
    axis(side = 1, at = seq(-2,2,0.1), labels = NA)
    axis(side = 2, cex.axis=1.5, tcl = -1, padj = -0.25)
    
    for (k in 1:3) {
      points(dplyr::select(filter(data, strategy == "Overall", scenario == c("F","G","H")[k], prior == Type_prior[i]), c(heterogeneity,all_of(index))), col = colorRampPalette(c("blue","purple"))(3)[k], pch = c(0,2,6)[k], cex = 1.25)
      lines(dplyr::select(filter(data, strategy == "Overall", scenario == c("F","G","H")[k], prior == Type_prior[i]), c(heterogeneity,all_of(index))), col = colorRampPalette(c("blue","purple"))(3)[k], lwd = 1.25)
      points(dplyr::select(filter(data, strategy == "Arm", scenario == c("F","G","H")[k], prior == Type_prior[i]), c(heterogeneity,all_of(index))), col = colorRampPalette(c("red","orange"))(3)[k], pch = c(0,2,6)[k], cex = 1.25)
      lines(dplyr::select(filter(data, strategy == "Arm", scenario == c("F","G","H")[k], prior == Type_prior[i]), c(heterogeneity,all_of(index))), col = colorRampPalette(c("red","orange"))(3)[k], lwd = 1.25)
      points(dplyr::select(filter(data, scenario == c("F","G","H")[k], prior == "No-borrowing"), c(heterogeneity,all_of(index))), col = "black", pch = NA, cex = 1)
      lines(dplyr::select(filter(data, scenario == c("F","G","H")[k], prior == "No-borrowing"), c(heterogeneity,all_of(index))), col = "black", lwd = 2)
    }
    abline(v=0,lty = "dashed", col = "grey80", lwd = 2)
  }
  plot(NA, xlim = c(0,1), ylim = c(0,1), yaxt="n", xaxt="n", xlab = NA, ylab = NA, frame.plot = F)
  legend("bottomleft", legend = c("R = 1, 1, 1 (Overall)",
                                  "R = 1, 0.4, 0.4 (Overall)",
                                  "R = 1, 2.5, 2.5 (Overall)",
                                  "R = 1, 1, 1 (Arm-based)",
                                  "R = 1, 0.4, 0.4 (Arm-based)",
                                  "R = 1, 2.5, 2.5 (Arm-based)",
                                  "No-borrowing"), plot = T, horiz = F, ncol = 1, lwd=2,
         lty = 1, pch = c(0,2,6,0,2,6,NA), col=c(c(colorRampPalette(c("blue","purple"))(3),colorRampPalette(c("red","orange"))(3)),1), y.intersp = 1, cex = 2.5, bty = "n", xpd=NA)
  
  
  for (i in 1:7){
    plot(NA, xlim = c(-2,2), ylim = yrange, xlab = "Population mean heterogeneity", yaxt="n", xaxt = "n", ylab = ylab, 
         frame.plot = T, main = paste0("B.",Type_prior[i]), cex.main = 3.5, cex.lab = 2.5, cex.axis=1.5)
    axis(side = 1, tcl = -0.75, padj = 0.35, cex.axis = 1.5)
    axis(side = 1, at = seq(-2,2,0.1), labels = NA)
    axis(side = 2, cex.axis=1.5, tcl = -1, padj = -0.25)
    
    for (k in 1:3) {
      points(dplyr::select(filter(data, strategy == "Overall", scenario == c("I","J","K")[k], prior == Type_prior[i]), c(heterogeneity,all_of(index))), col = colorRampPalette(c("blue","purple"))(3)[k], pch = c(0,2,6)[k], cex = 1.25)
      lines(dplyr::select(filter(data, strategy == "Overall", scenario == c("I","J","K")[k], prior == Type_prior[i]), c(heterogeneity,all_of(index))), col = colorRampPalette(c("blue","purple"))(3)[k], lwd = 1.25)
      points(dplyr::select(filter(data, strategy == "Arm", scenario == c("I","J","K")[k], prior == Type_prior[i]), c(heterogeneity,all_of(index))), col = colorRampPalette(c("red","orange"))(3)[k], pch = c(0,2,6)[k], cex = 1.25)
      lines(dplyr::select(filter(data, strategy == "Arm", scenario == c("I","J","K")[k], prior == Type_prior[i]), c(heterogeneity,all_of(index))), col = colorRampPalette(c("red","orange"))(3)[k], lwd = 1.25)
      points(dplyr::select(filter(data, scenario == c("I","J","K")[k], prior == "No-borrowing"), c(heterogeneity,all_of(index))), col = "black", pch = NA, cex = 1)
      lines(dplyr::select(filter(data, scenario == c("I","J","K")[k], prior == "No-borrowing"), c(heterogeneity,all_of(index))), col = "black", lwd = 2)
    }
    abline(v=c(-1,0,1),lty = "dashed", col = "grey80", lwd = 2)
  }
  plot(NA, xlim = c(0,1), ylim = c(0,1), yaxt="n", xaxt="n", xlab = NA, ylab = NA, frame.plot = F)
  legend("bottomleft", legend = c("R = 1, 1, 1 (Overall)",
                                  "R = 1, 0.4, 0.4 (Overall)",
                                  "R = 1, 2.5, 2.5 (Overall)",
                                  "R = 1, 1, 1 (Arm-based)",
                                  "R = 1, 0.4, 0.4 (Arm-based)",
                                  "R = 1, 2.5, 2.5 (Arm-based)",
                                  "No-borrowing"), plot = T, horiz = F, ncol = 1, lwd=2,
         lty = 1, pch = c(0,2,6,0,2,6,NA), col=c(c(colorRampPalette(c("blue","purple"))(3),colorRampPalette(c("red","orange"))(3)),1), y.intersp = 1, cex = 2.5, bty = "n", xpd=NA)
  
  dev.off()
}
