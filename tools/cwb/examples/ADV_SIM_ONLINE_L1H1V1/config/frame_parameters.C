{
  // NET setting
  nIFO = 3;           // number of detectors from the array ifo[]
  strcpy(ifo[0],"L1");
  strcpy(ifo[1],"H1");
  strcpy(ifo[2],"V1");

  slagSize = 1;	     // fix segLen	

  // set the output frame dir
  strcpy(output_dir,"frames");
  
  // simulation parameters
  simulation = 2;     // factor are snr
  nfactor = 1;        // number of strain factors
  double snr = 50.;
  double FACTORS[] ={snr};

  for(int i=0;i<nfactor;i++) factors[i]=FACTORS[i];

  // set the output strain channel names
  strcpy(channelNamesRaw[0],"L1:FAKE_STRAIN");
  strcpy(channelNamesRaw[1],"H1:FAKE_STRAIN");
  strcpy(channelNamesRaw[2],"V1:FAKE_16384Hz");

  // set the output mdc channel names
  strcpy(channelNamesMDC[0],"L1:FAKE-STRAIN_BURST");
  strcpy(channelNamesMDC[1],"H1:FAKE-STRAIN_BURST");
  strcpy(channelNamesMDC[2],"V1:FAKE-STRAIN_BURST");

  // set label of output frame files -> (label-GPS-LEN.gwf)
  strcpy(frFiles[0],"L-L1_llhoft");
  strcpy(frFiles[1],"H-H1_llhoft");
  strcpy(frFiles[2],"V-V1_llhoft");
/*
  // if frFile[i+nIFO]=="ADD_TO_STRAIN" ->  strain=strain+mdc
  strcpy(frFiles[3],"ADD_TO_STRAIN");
  strcpy(frFiles[4],"ADD_TO_STRAIN");
  strcpy(frFiles[5],"ADD_TO_STRAIN");
*/
  // set label of output mdc frame log files -> (Log-X1Y1Z1-label-jobXXX.gwf)
  strcpy(injectionList,"FrameOnline");

  // {ifo, dqcat_file, dqcat[0/1/2], shift[sec], inverse[false/true], 4columns[true/false]}
  nDQF=3;
  dqfile dqf[nDQF]={
                    {"L1" ,"input/frame.in",           CWB_CAT1, 0., false, false},
                    {"H1" ,"input/frame.in",           CWB_CAT1, 0., false, false},
                    {"V1" ,"input/frame.in",           CWB_CAT1, 0., false, false},
                   };
  for(int i=0;i<nDQF;i++) DQF[i]=dqf[i];

  // define frame parameters
  TString frDir       = output_dir;
  TString frLabel     = injectionList;

  network* net = new network;
  detector* pD[NIFO_MAX];      					  //! pointers to detectors
  for(int i=0; i<nIFO; i++) {
    if(strlen(ifo[i])>0) pD[i] = new detector(ifo[i]);        // built in detector
    else                 pD[i] = new detector(detParms[i]);   // user define detector
  }
  for(int i=0; i<nIFO; i++) net.add(pD[i]);

  plugin = TMacro("macro/CWB_Plugin_OnlineFrame.C");        		// Macro source
  configPlugin = TMacro("macro/CWB_Plugin_OnlineFrame_Config.C");  	// Macro config
}
