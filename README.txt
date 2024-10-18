****************************************************************************************
### Description
****************************************************************************************

Title: "Modified Conditional Borrowing-by-part Power Prior for dynamic and parameter-specific information borrowing of the Gaussian endpoint"
Authors: Kai Wang, Han Cao, Chen Yao. 

For questions, comments, or remarks about the code please email: wangkai1995jagger@163.com (Kai Wang).

****************************************************************************************
### Files and subfolders
****************************************************************************************

*1  "Function_for_SingleEC.R" and "Function_for_MultipleECs.R" contain the R functions with respect to borrowing a single external control arm and multiple external control arms (three), respectively.

*2  "Numerical_study_for_SingleEC.R" and "Numerical_study_for_MultipleECs.R" contain the codes to conduct the numerical study in Section 4.4, including borrowing a single external control arm and multiple external control arms (three), respectively.

*3  "Simulation_study_for_SingleEC.R" and "Simulation_study_for_MultipleECs.R" contain the codes to conduct the simulation study in Section 5, including borrowing a single external control arm and multiple external control arms (three), respectively.

*4  "Application.R" contains the codes to reproduce the application results (Table 2) in Section 6 of this article.

*4  "Plot_function.R" and "Plot.R" contain the functions and codes to reproduce the figure presented in this article and the Supporting information.

*5  Folder "stan" contains the stan codes of fully Bayesian dynamic borrowing methods to borrow the Gaussian endpoint with unknown variance from external controls.

*6  Folder "Simulation_dataset" is used to save the simulated datasets.

*7  Folder "Num_rst_for_SingleEC" and "Num_rst_for_MultipleECs" are used to save the results of the numerical study in Section 4.4.
	a) Single external control arm (each dataset takes about 0.3 hours using 15 cores):  
		"num_0.rds", "num_A.rds", "num_B.rds", "num_C.rds", "num_D.rds", "num_E.rds".
	b) Multiple external control arms (each dataset takes about 0.4 hours using 15 cores):
		"num_situation1.rds" and "num_situation2.rds".

*7  Folder "Sim_rst_for_SingleEC" and "Sim_rst_for_MultipleECs" is used to save the intermediate and final results of the simulation study.
	
	a) Single external control arm: 
			i. Intermediate results (each dataset takes about 8-10 hours using 15 cores):  
			(***** Notice: The following intermediate results can be downloaded from: "https://doi.org/10.5281/zenodo.13935144" (Version v2).******)
				"simu_1_A.rds"~"simu_1_E.rds" (Design 1);  
					"simu_1_A_adapt.rds"~"simu_1_E_adapt.rds" (Design 1 & adaptive selection & strict tolerance); 
					"simu_2_A_adapt.rds"~"simu_2_E_adapt.rds" (Design 1 & adaptive selection & mild tolerance);		        
				"simu_2_A.rds"~"simu_2_E.rds" (Design 2); 
				"simu_3_A.rds"~"simu_3_E.rds" (Design 3); 
				"simu_4_A.rds"~"simu_4_E.rds" (Design 4);
	
			ii. Final results of operating characteristics used for plotting are all saved for ready use: 
				"simu_1_H0_A.rds"~"simu_1_H0_E.rds" (Design 1 & H0); "simu_1_H1_A.rds"~"simu_1_H1_E.rds" (Design 1 & H1);
					"simu_1_H0_A_adapt.rds"~"simu_1_H0_E_adapt.rds" (Design 1 & adaptive selection & H0 & strict tolerance); 
					"simu_1_H1_A_adapt.rds"~"simu_1_H1_E_adapt.rds" (Design 1 & adaptive selection & H1 & strict tolerance);
					"simu_2_H0_A_adapt.rds"~"simu_2_H0_E_adapt.rds" (Design 1 & adaptive selection & H0 & mild tolerance); 
					"simu_2_H1_A_adapt.rds"~"simu_2_H1_E_adapt.rds" (Design 1 & adaptive selection & H1 & mild tolerance);
				"simu_2_H0_A.rds"~"simu_2_H0_E.rds" (Design 2 & H0); "simu_2_H1_A.rds"~"simu_2_H1_E.rds" (Design 2 & H1);
				"simu_3_H0_A.rds"~"simu_3_H0_E.rds" (Design 3 & H0); "simu_3_H1_A.rds"~"simu_3_H1_E.rds" (Design 3 & H1); 
				"simu_4_H0_A.rds"~"simu_4_H0_E.rds" (Design 4 & H0); "simu_4_H1_A.rds"~"simu_4_H1_E.rds" (Design 4 & H1);
			        
			iii. "Table_S2.rds": the summary of calibrated threshold probability

	b) Multiple external control arms: 
			i. Intermediate results (each dataset takes about 8 hours using 15 cores):  
			(***** Notice: The following intermediate results can be downloaded from: "https://doi.org/10.5281/zenodo.13935144" (Version v2).******)
				"simu_5_F.rds"~"simu_5_H.rds" (Design 5, Scenario F-H); 
				"simu_5_I.rds"~"simu_5_K.rds" (Design 5, Scenario I-K)
			
			ii. Final results of operating characteristics used for plotting are all saved for ready use: 
				"simu_5_F_H0.rds"~"simu_5_H_H0.rds" (Design 5, Scenario F-H & H0); "simu_5_I_H0.rds"~"simu_5_K_H0.rds" (Design 5, Scenario I-K & H0)
				"simu_5_F_H1.rds"~"simu_5_H_H1.rds" (Design 5, Scenario F-H & H1); "simu_5_I_H1.rds"~"simu_5_K_H1.rds" (Design 5, Scenario I-K & H1)
			      
			iii. "Table_S3.rds": the summary of calibrated threshold probability

*7  Folder "Application_rst" is used to save the application results in Section 6.
	a) Intermediate results: "rst_case1.rds", "rst_case2.rds"; "rst_varying_ratio.rds"
	b) Final results: "Table_2.rds"
	c) Results for Figure 5: "rst_plot_varying_ratio.rds"

*8  Folder "Figure" is used to save the figures presented in this article and its Supporting information.


****************************************************************************************
****************************************************************************************
### Instructions for the reproduction of figures and tables reported in the article:
****************************************************************************************
****************************************************************************************

*To reproduce the figures in the main text and Supporting information, please direct run the "Plot.R", and the results are saved in the Folder "Figure". 
  Notice that the final operating characteristics are all saved in folders "Sim_rst_for_SingleEC", "Sim_rst_for_MultipleECs", and "Application_rst" for ready use. 

*To reproduce Table 2 and Figure 5 containing the application results, please direct run the "Application.R", and the results are saved in the Folder "Application_rst".

*To reproduce Table S2, S3 containing the calibrated threshold value, there is no need to run the simulation yourself. 
  Please direct run Line 246-270 of "Simulation_study_for_SingleEC.R" and Line 167-185 of "Simulation_study_for_MultipleECs.R", and the results are saved in "Table_S2.rds" in the Folder "Sim_rst_for_SingleEC" and "Table_S3.rds" in the Folder "Sim_rst_for_MultipleECs". 

*To reproduce the simulation results, please run the "Simulation_study_for_SingleEC.R" and "Simulation_study_for_MultipleECs.R". 
  We recommend parallel computation and selecting certain designs or scenarios according to the code annotation to save time.


****************************************************************************************
### Software Information
****************************************************************************************

R version 4.1.3 (2022-03-10)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: OpenCloudOS 8.6

Matrix products: default
BLAS:   /usr/local/lib64/R/lib/libRblas.so
LAPACK: /usr/local/lib64/R/lib/libRlapack.so

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
 [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                  LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] snowfall_1.84-6.2    snow_0.4-4           MASS_7.3-58.3        RBesT_1.6-6          rstantools_2.3.1     rstan_2.21.8         StanHeaders_2.21.0-7
 [8] rlist_0.4.6.2        data.table_1.14.8    lubridate_1.9.2      forcats_1.0.0        stringr_1.5.0        dplyr_1.1.1          purrr_1.0.1         
[15] readr_2.1.4          tidyr_1.3.0          tibble_3.2.1         ggplot2_3.4.1        tidyverse_2.0.0     

loaded via a namespace (and not attached):
 [1] tidyselect_1.2.0   colorspace_2.1-0   vctrs_0.6.1        generics_0.1.3     stats4_4.1.3       loo_2.5.1          utf8_1.2.3         rlang_1.1.0       
 [9] pkgbuild_1.4.0     pillar_1.9.0       glue_1.6.2         withr_2.5.0        matrixStats_0.63.0 lifecycle_1.0.3    munsell_0.5.0      gtable_0.3.3      
[17] mvtnorm_1.1-3      codetools_0.2-18   inline_0.3.19      tzdb_0.3.0         callr_3.7.3        ps_1.7.3           parallel_4.1.3     fansi_1.0.4       
[25] Rcpp_1.0.10        scales_1.2.1       backports_1.4.1    checkmate_2.1.0    RcppParallel_5.1.7 gridExtra_2.3      hms_1.1.3          stringi_1.7.12    
[33] processx_3.8.0     grid_4.1.3         cli_3.6.1          tools_4.1.3        magrittr_2.0.3     Formula_1.2-5      crayon_1.5.2       pkgconfig_2.0.3   
[41] prettyunits_1.1.1  timechange_0.2.0   assertthat_0.2.1   rstudioapi_0.14    R6_2.5.1           compiler_4.1.3    
