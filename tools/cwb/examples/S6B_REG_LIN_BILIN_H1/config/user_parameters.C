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
  simulation = 1;     // 1 for simulation, 0 for production
  nfactor = 12;        // number of strain factors
  double FACTORS[] ={0.0375,
                     0.075,
                     0.15,
                     0.3,
                     0.6,
                     1.2,
                     2.4,
                     4.8,
                     9.6,
                     19.2,
                     38.4,
                     76.8}; // array of strain factors
  for(int i=0;i<nfactor;i++) factors[i]=FACTORS[i];

  strcpy(channelNamesRaw[0],"H1:LDAS-STRAIN");
  strcpy(frFiles[0],"input/S6B_R6_H1_LDAS_C02_L2.frames");

  // dq file list
  // {ifo, dqcat_file, dqcat[0/1/2], shift[sec], inverse[false/true], 4columns[true/false]}
  nDQF=3;
  dqfile dqf[nDQF]={
                    {"H1" ,"input/S6B_OFFLINE_H1SCIENCE.txt",         CWB_CAT0, 0., false, false},
                    {"H1" ,"input/S6B_OFFLINE_H1_DQCAT1SEGMENTS.txt", CWB_CAT1, 0., true, false},
                    {"H1" ,"input/S6B_OFFLINE_H1_DQCAT2SEGMENTS.txt", CWB_CAT2, 0., true, false}
                   };

  for(int i=0;i<nDQF;i++) DQF[i]=dqf[i];

  // Regression Ranking 
//  simulation   = 0;
//  plugin       = TMacro("../../../cwb/plugins/CWB_Plugin_Linear_Bilinear_Regression_Rank.C");               // Macro source
//  configPlugin = TMacro("macro/CWB_Plugin_Linear_Bilinear_Regression_Rank_Config.C");        // Macro Config source

  // Regression Cleaning
  simulation   = 0;
  plugin       = TMacro("../../../cwb/plugins/CWB_Plugin_Linear_Bilinear_Regression.C");               // Macro source
  configPlugin = TMacro("macro/CWB_Plugin_Linear_Bilinear_Regression_Config.C");        // Macro Config source
}
