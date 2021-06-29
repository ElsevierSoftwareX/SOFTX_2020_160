{
  #define WATFUN_HH

  bool pp_batch = false;	// if true then the interactive questions are disabled 

  TString pp_label = "";
  TString user_pp_label = "";
  char pp_data_dir[1024] = "data";
  double pp_drho    = (TString(analysis)=="1G") ?  0.1 : 0.2;	// rho step size of bins in FAR histogram
  double pp_rho_min = (TString(analysis)=="1G") ?  2 : 4;
  double pp_rho_max = (TString(analysis)=="1G") ?  10 : 15;
  int    pp_rho_bin = 0; // number of bins for FAR histogram, it is initialized in $CWB_EPPARAMETERS_FILE  
  int pp_max_nloudest_list = 25; // maximum number of loudest event in the report list
  bool pp_rho_log = false;	// if true then rho is plotted in log scale

  // hrss range used in eff_freq plots
  double pp_hrss_min = 1e-23;
  double pp_hrss_max = 1e-20;

  // Create jet_benchmark.png (estimation time plot)
  // (0) -> disabled - (>0) -> max value - (<0) -> auto max
  int pp_jet_benckmark=0;
  // Create mem_benchmark.png (maximum memory)
  // (0) -> disabled - (>0) -> max value - (<0) -> auto max
  int pp_mem_benckmark=0;

  // create job status table : list of jobs not finished (def = false)
  bool pp_job_status = false;

  // default index used for rho = rho[pp_irho] (used by T_cut)
  int pp_irho = (TString(analysis)=="1G") ?  1 : 0;	// 0/1
  // default index used for cc = netcc[pp_inetcc] (used by T_cor)
  int pp_inetcc  = 0;					// 0/1	

  // compute efficiency vs threshold for each mdc type @ specific (custom/auto) factor 
  // used to build ROC plots
  // syntax : "--mode disable/auto/file_name --rho_min xxx"
  // disabled  : not computed
  // auto      : factor is automatically computed -> is the nearest to the 50% efficiency
  // file_name : factor is provided by the input file
  // rho_min   : efficiency is computed for rho>rho_min 
  TString pp_eff_vs_thr = "--mode disabled --quit false";	 

  // options for the PE report 
  // mdc       : used to declare the mdc type[1] to be displayed 
  //             -1 : use the list of mdc declared in the production parameter injectionList
  //              0 : all types are displayed together (default)
  //              k : only type[1]=k is displayed : k = [1:N]
  //      i/j/.../k : display multiple mdc types -> i/j/.../k  : i,j,...,k=[1:N] 
  TString pp_pe = "--mdc 0";

  // pp_sreport is an array used to add auxiliary reports to the current report
  // --link LINK : www link of the auxiliary report
  //               --link <tab> open a sub tab
  //               --link </tab> close sub tab
  // --label LABEL: report label - used in the html tabs
  // --high HIGH: report window height (optional : def=4000)  
  // --name NAME: report html name (optional : def = "header.html, body.html")
  vector<TString> pp_sreport;
  pp_sreport.resize(100);
  for(int i=0;i<100;i++) pp_sreport[i]="";

  // enable/disable draw efficiency fit curves
  // true (default) : show fit curve, eff@hrss10/50/90 are derived from fit
  // false          : show lines connecting points, eff@hrss10/50/90 are not showed 
  bool pp_show_eff_fit_curve = true;

  // enable/disable merge MDC types
  // false (default) : the efficiencies curves are produced individually for each MDC type
  // true            : all MDC types are used to produce only one MDC efficiency curve
  bool pp_merge_types = false;

  // enable/disable draw of rho_vs_rate with no_multiplicity
  // false/true     : disable/enable  (default=false) 
  bool pp_rho_vs_rate_no_multiplicity = false;  

  // enable/disable FAD report (default is disabled)
  // if == "" -> disabled
  // --bkgrep  : directory of the background report (used to read rate_threshold_veto.txt and live.txt)  
  // --hrss    : if it is a number then it is used as normalization constant for all MDC types (def=0)
  //		 if =0 then hrss rescale is not applied
  //             if it is a file then it is the list of hrss used for each MDC type 
  //               (format : for each line -> hrss) 
  // --gfit    : true/false -> enable/disable fit in the output plots (default=false)
  // --rhomin  : minimum rho value selected for plots (default = 5)
  // --nzbins  : (default is 0)
  //             if nzbins=0 the standard FAD statistic is used
  //             if nzbin>0 the FAD statistic is computed until there are nzbins
  //               consecutive bins with zero events inside
  //             if nzbin<0 the FAD statistic is computed with classical FAD and min-hold
  // --units   : K -> Kpc,Kyr : M -> Mpc,Myr  (def=M)
  // --distr   : formula/mdc -> radial distribution is computed from formula or from mdc injections (def=MDC)
  // --nbins   : number of bins used in hist to computed the radial distribution from the mdc injections
  // --header  : if true -> add cwb header to fad html file (def=false)
  // --multi   : if true -> FAD multi plot for each mdc set are created and substituted in the sim report
  //             page to the eff_freq plots (def=false)
  // --title   : title of the html page (default=FAD) : spaces must be filled with *
  // --subtitle: subtitle of the html page (default="") : spaces must be filled with *
  //
  TString pp_fad="";

  // if pp_factor2distance > 0 && simulation==1 the efficiency plots are converted from 'eff vs hrss' to 'eff vs distance'
  // distance = pp_hrss2distance/factor
  // pp_factor2distance is the distance in Kpc when factor=1
  double pp_factor2distance = 0.;

  TString cwb_merge_label;
  if(gSystem->Getenv("CWB_MERGE_LABEL")==NULL) {
    cout << "Error : environment CWB_MERGE_LABEL is not defined!!!" << endl;exit(1);
  } else {
    cwb_merge_label=TString(gSystem->Getenv("CWB_MERGE_LABEL"));
  }

  int cwb_slag_number=-1;  // slag>0
  if(gSystem->Getenv("CWB_SLAG_NUMBER")==NULL) {
    cout << "Error : environment CWB_SLAG_NUMBER is not defined!!!" << endl;exit(1);
  } 
  if(TString(gSystem->Getenv("CWB_SLAG_NUMBER")).IsDigit()) {
    cwb_slag_number=TString(gSystem->Getenv("CWB_SLAG_NUMBER")).Atoi();
  }

  int cwb_lag_number=-1;  // lag>0
  if(gSystem->Getenv("CWB_LAG_NUMBER")==NULL) {
    cout << "Error : environment CWB_LAG_NUMBER is not defined!!!" << endl;exit(1);
  } 
  if(TString(gSystem->Getenv("CWB_LAG_NUMBER")).IsDigit()) {
    cwb_lag_number=TString(gSystem->Getenv("CWB_LAG_NUMBER")).Atoi();
  }

  // ----------------------------------------------------------
  // title & subtitle
  // ----------------------------------------------------------
  char title[1024]=""; 
  char subtitle[1024]=""; 

  // ----------------------------------------------------------
  // thresholds
  // ----------------------------------------------------------
  double T_ifar     = 0.0;        // ifar cut (years) : can be applied only for tree where the ifar is included
  double T_cor      = 0.60;       // cc cut
  double T_cut      = 3.5;        // rho high frequency cut
  double T_scc      = 0.0;        // subcc (netcc[2]) : subnetwork consistency coefficient [2G]
  double T_acor     = 0.0;        // output threshold
  double T_out      = 3.5;        // output threshold
  double T_vED      = 0.0;        // vED threshold
  double T_pen      = 0.;         // penalty threshold
  double T_hrss     = 0.0;        // penalty threshold
  double T_win      = 0.1;        // time window used to skip random coincidence in simplot
  double T_mchirp   = 0.0;        // cut on reconstructed chirp
  double T_sphr     = 0.0;        // cut on sphericity of reconstructed chirp
  double T_efrac    = 0.0;        // cut on Energy fraction
  int    i_hrss1    = 0;
  int    i_hrss2    = 1;
  double T_rms      = 0.;         // penalty threshold
  double hours      = 24;         // bin size in hours for rate vs time plot 

  float ratio_hrss  = 1;          //CORRECTION FREQUENCY
  float pp_freq_offset = 0.;      //CORRECTION FREQUENCY  (SIM)

  char  RunLabel[1024]="";

  // ----------------------------------------------------------
  // Input data files
  // ----------------------------------------------------------

  char sim_file_name[1024];
  char mdc_file_name[1024];
  char net_file_name[1024];
  char liv_file_name[1024];
  char mdc_inj_file[1024];

  if(simulation>0) {	// simulation
    sprintf(sim_file_name,"%s/wave_%s.%s.root",merge_dir,data_label,cwb_merge_label.Data());
    sprintf(mdc_file_name,"%s/mdc_%s.%s.root",merge_dir,data_label,cwb_merge_label.Data());
    cout << sim_file_name << endl;
    cout << mdc_file_name << endl;

    // Declare inj list
    if(TString(injectionList).Sizeof()<=1) {
      cout << "cwb_pparameters.C - Error : injectionList must be declared in simulation mode !!!" << endl;
      exit(1);
    }
    
    TObjArray* token;
    // extract file name
    token = TString(injectionList).Tokenize(TString("/"));
    TString injName = ((TObjString*)token->At(token->GetEntries()-1))->GetString();
    TString fext = "";
    // extract file extention
    if(injName.Contains(".")) { // there is an extention
      token = TString(injName).Tokenize(TString("."));
      fext = TString(".")+((TObjString*)token->At(token->GetEntries()-1))->GetString();
      sprintf(mdc_inj_file,"%s",injName.ReplaceAll(fext,"").Data());
    } else {
      sprintf(mdc_inj_file,"%s",injName.Data());
    }
    // extract the name if mdc injection list
    char _mdc_inj_file[1024];
    sprintf(_mdc_inj_file,"%s/%s.inj",input_dir,mdc_inj_file);
    strcpy(mdc_inj_file,_mdc_inj_file);

    // Check if files exist
    CWB::Toolbox::checkFile(sim_file_name);
    CWB::Toolbox::checkFile(mdc_file_name);
    CWB::Toolbox::checkFile(mdc_inj_file);
  } else if (simulation==0) {	// production
    sprintf(net_file_name,"%s/wave_%s.%s.root",merge_dir,data_label,cwb_merge_label.Data());
    sprintf(liv_file_name,"%s/live_%s.%s.root",merge_dir,data_label,cwb_merge_label.Data());
    cout << net_file_name << endl;
    cout << liv_file_name << endl;

    // Check if files exist
    CWB::Toolbox::checkFile(net_file_name);
    CWB::Toolbox::checkFile(liv_file_name);
  }

  // ----------------------------------------------------------
  // Create output directory
  // ----------------------------------------------------------
  char netdir[1024];   

  // ----------------------------------------------------------
  // templates file for www
  // ----------------------------------------------------------
  char  html_index_template[1024]="";
  char  html_header_template[1024]="";
  char  html_body_prod_template[1024]="";

  if(gSystem->Getenv("CWB_HTML_INDEX")==NULL) {
    cout << "Error : environment CWB_HTML_INDEX is not defined!!!" << endl;exit(1);
  } else {
    strcpy(html_index_template,gSystem->Getenv("CWB_HTML_INDEX"));
  }
  if(gSystem->Getenv("CWB_HTML_HEADER")==NULL) {
    cout << "Error : environment CWB_HTML_HEADER is not defined!!!" << endl;exit(1);
  } else {
    strcpy(html_header_template,gSystem->Getenv("CWB_HTML_HEADER"));
  }
  if(gSystem->Getenv("CWB_HTML_BODY_PROD")==NULL) {
    cout << "Error : environment CWB_HTML_BODY_PROD is not defined!!!" << endl;exit(1);
  } else {
    strcpy(html_body_prod_template,gSystem->Getenv("CWB_HTML_BODY_PROD"));
  }

  CWB::Toolbox::checkFile(html_index_template);
  CWB::Toolbox::checkFile(html_header_template);
  CWB::Toolbox::checkFile(html_body_prod_template);
/*
  // ----------------------------------------------------------
  // ch2 cuts
  // ----------------------------------------------------------

  char  ch2[2048]="";
  if(cwb_lag_number==-1) {
    sprintf(ch2,"netcc[%d]>%f&&lag[%d]>0",pp_inetcc,T_cor,nIFO);
  } else {
    sprintf(ch2,"netcc[%d]>%f&&lag[%d]==%d",pp_inetcc,T_cor,nIFO,cwb_lag_number);
  }
  if (T_vED>0)    sprintf(ch2,"%s&&neted[0]/ecor<%f",ch2,T_vED);
  if (T_pen>0)    sprintf(ch2,"%s&&penalty>%f",ch2,penalty);
  if (T_hrss>0)   sprintf(ch2,"%s&&abs((TMath::Log10(hrss[%i])-TMath::Log10(hrss[%i])))<%f",ch2,i_hrss1,i_hrss2,T_hrss);
*/
  // ----------------------------------------------------------
  // frequency cuts
  // array of frequency cuts [lowFCUT,highFCUT] 
  // ----------------------------------------------------------
  int    nFCUT = 0;           // number of frequency cuts
  double lowFCUT[100];        // low frequency cut
  double highFCUT[100];       // high frequency cut

  // ----------------------------------------------------------
  // VETO cuts
  // dq file list
  // {ifo, dqcat_file, dqcat[0/1/2], shift[sec], inverse[false/true], 4columns[true/false]}
  // ----------------------------------------------------------
  int nVDQF=0;
  dqfile VDQF[100];

  char  veto_vetoed[1024] = "";
  char  veto_not_vetoed[1024]  = "";


  // ----------------------------------------------------------
  // CBC variables
  // ----------------------------------------------------------

  //Fiducial volume times time initialized to -1: in case of redshifted MonteCarlos it has to be set to the 
  //correct value for the injected simulation in the user pp file. A control in cwb_report_cbc will check if it's -1 
  //and, in case, it will exits with an error 
  double VT = -1.0;    

  //Total time covered by the MonteCarlo: in case of redshifted MonteCarlos it has to be set to the 
  //correct value for the injected simulation in the user pp file. A control in cwb_report_cbc will check if it's -1 
  //and, in case, it will exits with an error  
  int TMC = -1.0;   
  
  //By default, the inner shell of the volume is not included. Negligible effects in all relevant cases. 
  bool INCLUDE_INTERNAL_VOLUME = 0;

  //By default, TRIALS is set to 1 
  int TRIALS = 1; 
}
