-------------------------------------------------------------------------------
 How to produce ROC (Receiver Operatin Characteristic) with cWB
 The ROC used here are plots which reports the efficiency @ hrss50% vs FAR
-------------------------------------------------------------------------------

Consider the following example : 1G vs 2G

Perform the following steps :

1) Produce the standard 1G & 2G background

2) Produce the standard 1G & 2G simulation 

   - NOTE : use the same same MDC set & simulation=1 (user_parameters.C)

   - go into 1G working directory

     - define in config/user_pparameters.C the following parameters :
       pp_eff_vs_thr = "--mode auto --rho_min 4.5 --eff 50";
       pp_rho_max = xx;  // xx must greater than the rho of the loudest background event

     - run 'cwb_report MERGE_LABEL create'  // Ex : MERGE_LABEL = M1 / M1.V_hvetoLHV_cat3LHV
       this command produces the following files under the standard report directory :

       A) eff_50_threshold_*.txt    	    // one for each mdc type : * = MDC_NAME
          is a list of efficiencies vs rho  // start from rho_min , end to pp_rho_max

       B) eff_50_threshold_factors.txt       
          each entry in the list is : MDC_NAME MDC_TYPE FACTOR
          where : MDC_NAME = MDC name
                  MDC_TYPE = identification number (type[0])
                  FACTOR   = factor selected from the list declared in user_parameters.C
                             if --mode auto : factor is automatically computed using the 
                                              post-production thresholds defined in user_pparameters.C
                                              -> it is the nearest factor to the 50% efficiency

          NOTE : the file eff_50_threshold_factors.txt is used as input for the second working directory
                 because the comparison must be done with the the same factors

   - go into 2G working directory

     - define in config/user_pparameters.C the following parameters :
       pp_eff_vs_thr = "--mode eff_50_threshold_factors.txt --eff 50 --quit true";
       pp_rho_max = xx;  // xx must greater than the rho of the loudest background event

     - run 'cwb_report MERGE_LABEL create'  // Ex : MERGE_LABEL = M1 / M1.V_hvetoLHV_cat3LHV
       this command produces the following files under the standard report directory :
       A) eff_50_threshold_*.txt    	    // one for each mdc type : * = MDC_NAME
       B) eff_50_threshold_factors.txt       
     
3) defines the configuration file for MakeROC.C
   The format is described in the header of the MakeROC.C macro    
   
   Ex : tools/cwb/tutorials/MakeROC_Config_1G_vs_2G_Example.C

4) run MakeROC.C

   Ex : root 'MakeROC.C("MakeROC_Config_1G_vs_2G_Example.C",50)'

   The ROC plots and index.html are saved in the output ROC_DIR+"_50"
   The report can be viewed from a WEB browser at the link defined in ROC_WWW
 
