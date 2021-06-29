{
  // NET setting
  nIFO = 1;           // number of detectors from the array ifo[]
  cfg_search = 'r';       // see description in parameters.C

  strcpy(ifo[0],"H1");

  strcpy(refIFO,"H1"); // reference IFO

  // segments
  segLen     = 600.;  // Segment length [sec]
  segMLS     = 600.;  // Minimum Segment Length after DQ_CAT1 [sec]
  segTHR     = 30.;   // Minimum Segment Length after DQ_CAT2 [sec]
  segEdge    = 32.;   // wavelet boundary offset [sec]

  // simulation parameters
  simulation = 0;     // 1 for simulation, 0 for production
  nfactor = 1;        // number of strain factors
  double FACTORS[] ={1.0}; // array of strain factors
  for(int i=0;i<nfactor;i++) factors[i]=FACTORS[i];

  strcpy(channelNamesRaw[0],"H1:LDAS-STRAIN");
  strcpy(frFiles[0],"input/H1_LDAS_C02_L2_94240000_942500000.frl");

  // dq file list
  // {ifo, dqcat_file, dqcat[0/1/2], shift[sec], inverse[false/true], 4columns[true/false]}
  nDQF=1;
  dqfile dqf[nDQF]={
                    {"H1" ,"input/period.txt",           CWB_CAT0, 0., false, false},
                   };

  for(int i=0;i<nDQF;i++) DQF[i]=dqf[i];

  plugin       = TMacro("macro/CWB_Plugin_S6B_REG_H1_PL_REMOVAL.C");               // Macro source
  configPlugin = TMacro("macro/CWB_Plugin_S6B_REG_H1_PL_REMOVAL_Config.C");        // Macro Config source
}
