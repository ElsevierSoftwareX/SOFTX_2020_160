// This macro computes the livetime for each lag using the CWB::Toolbox::getLiveTime method
// NOTE : this macro works only for standard lags !!!
//
// The livetime obtained with this method is not exactly equals to the one computed
// in the CWB pipeline (network::getNetworkPixels). 
//
// In the network::getNetworkPixels the livetime takes into account the time discretization
// used in the TF level with the maximum dt, while in the CWB::Toolbox::getLiveTime method
// the time is not discretized.
// NOTE : in the CWB pipeline only DQ CAT2 are used 
//   
// usage : root -n -l -b $CWB_PARMS_FILES GetLiveTime.C
//
// the DQ lists must be declared in the config/user_parameters.C
// To get the livetime after the CAT1+CAT2            ->  #define DQ_LIVETIME CWB_CAT2
// To get the livetime after the CAT1+CAT2+CAT3       ->  #define DQ_LIVETIME CWB_CAT3
// To get the livetime after the CAT1+CAT2+CAT3+HVETO ->  #define DQ_LIVETIME CWB_HVETO
// livetimes are saved in the OUTPUT_LIVETIME_FILE_NAME 
// NOTE : remember that the DQ files defined in the user_parameters.C have the boolean parameter "inverse"
//        defined with the opposite value respect to the one defined in user_pparameters.C
// 
// Author : Gabriele Vedovato


#define DQ_LIVETIME CWB_CAT2		// livetime after the CAT1+CAT2
//#define DQ_LIVETIME CWB_CAT3		// livetime after the CAT1+CAT2+CAT3
//#define DQ_LIVETIME CWB_HVETO		// livetime after the CAT1+CAT2+CAT3+HVETO

#define OUTPUT_LIVETIME_FILE_NAME	"livetime.txt"

void GetLiveTime() {

  CWB::Toolbox TB;

  // check if parameter files exist
  TB.checkFile(gSystem->Getenv("CWB_ROOTLOGON_FILE"));
  TB.checkFile(gSystem->Getenv("CWB_PARAMETERS_FILE"));
  TB.checkFile(gSystem->Getenv("CWB_UPARAMETERS_FILE"));

  if(simulation) {
    cout << "GetLiveTime : simulation is !=0 -> only zero lag is defined!!!" << endl << endl;
    gSystem->Exit(1);
  }

  // check input user configuration
  CWB::config cfg;
  cfg.Import();
  cfg.Check();

  // define network
  detector* pD[NIFO_MAX];                 	// pointers to detectors
  for(int i=0; i<nIFO; i++) {
    if(strlen(ifo[i])>0) pD[i] = new detector(ifo[i]);        // built in detector
    else                 pD[i] = new detector(detParms[i]);   // user define detector
  }

  network NET;                                  // network
  for(int i=0; i<nIFO; i++) NET.add(pD[i]);     // add deetctors to network object

  wavearray<double> x(segLen*16384);		// dummy array used by setTimeShifts
  pD[0]->TFmap=x;

  // get lags from network class
  int lags = NET.setTimeShifts(lagSize,lagStep,lagOff,lagMax,lagFile,lagMode,lagSite);
  cout<<"lag step: "<<lagStep<<endl;
  cout<<"number of time lags: "<<lags<<endl;

  // compute cat1List 
  // CAT1 list must be defined in user_parameters.C as CWB_CAT1
  vector<waveSegment> cat1List=TB.readSegList(nDQF, DQF, CWB_CAT1);
  // compute cat2List after CAT1+CAT2 
  // CAT2 list must be defined in user_parameters.C as CWB_CAT2
  vector<waveSegment> cat2List=TB.readSegList(nDQF, DQF, CWB_CAT2);
  // compute cat3List after CAT1+CAT2+CAT3 
  // CAT3 list must be defined in user_parameters.C as CWB_CAT3
  vector<waveSegment> cat3List=TB.readSegList(nDQF, DQF, CWB_CAT3);
  // compute havetoList after CAT1+CAT2+CAT3+HVETO 
  // HVETO list must be defined in user_parameters.C as CWB_HVETO
  vector<waveSegment> hvetoList=TB.readSegList(nDQF, DQF, CWB_HVETO);
  // compute jobList
  vector<waveSegment> jobList=TB.getJobList(cat1List, cat2List, segLen, segMLS, segTHR, segEdge);
  cout<<"Number of standard jobs : " << jobList.size() <<endl<<endl;

  TString fName = OUTPUT_LIVETIME_FILE_NAME;
  FILE* fP=NULL; 
  if((fP = fopen(fName.Data(), "w")) == NULL) {
    cout << "cannot open output file " << fName.Data() <<". \n";
    gSystem->Exit(1);
  };
  cout << "Write output file : " << fName.Data() << endl;
  fprintf(fP,"# lagID\t live\t\t\n");

  // compute the approximate value of the livetime obtained with CWB::Toolbox::getLiveTime
  printf("%8s ","lag");
  for(int n=0; n<nIFO; n++) printf("%12.12s%2s","ifo",NET.getifo(n)->Name);
  printf("\n");
  vector<double> shiftList(nIFO);
  for(int m=0;m<lags;m++) {

    printf("%8d ",(int)NET.wc_List[m].shift);
    for(int n=0; n<nIFO; n++) printf("%14.5f",NET.getifo(n)->lagShift.data[m]);
    printf("\n");

    for(int n=0;n<nIFO;n++) shiftList[n]=NET.getifo(n)->lagShift.data[m];
   
    double live = 0;
#if ( DQ_LIVETIME == CWB_CAT2 )
    live = TB.getLiveTime(jobList, cat2List, shiftList);
#endif
#if ( DQ_LIVETIME == CWB_CAT3 )
    live = TB.getLiveTime(jobList, cat3List, shiftList);
#endif
#if ( DQ_LIVETIME == CWB_HVETO )
    live = TB.getLiveTime(jobList, hvetoList, shiftList);
#endif

    cout << endl;
    cout << "-> livetime : " << int(live) << " sec "
         << live/3600. << " h " << live/86400. << " day" << endl;
    cout << endl;

    fprintf(fP,"%8d %10d\n",(int)NET.wc_List[m].shift,int(live));
  }
  if(fP!=NULL) fclose(fP);

  gSystem->Exit(0);
}
