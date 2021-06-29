{
  strcpy(analysis,"1G");

  nIFO = 1;
  search = 'r';
  strcpy(ifo[0],"");
  strcpy(refIFO,"H2");

  detectorParams detPARAMS[nIFO] = {
                          {"H2", 46.4551, -119.408,  0.0, 0, ( +90-125.998),    0, (+45-125.998 )},   // H2
                         };
  for(int i=0;i<nIFO;i++) detParms[i] = detPARAMS[i];

  lagSize = 1;
  lagOff = 0;
  segTHR = 0;

  strcpy(channelNamesRaw[0],"H2:SIM-STRAIN");

  strcpy(frFiles[0],"input/H2.frames");

  nDQF=1;
  dqfile dqf[nDQF]={
                    {"H2" ,"input/burst.in",           CWB_CAT1, 0., false, false},
                   };
  for(int i=0;i<nDQF;i++) DQF[i]=dqf[i];
  cedRHO = 3.000000;
  fLow=64;

  segLen = 60;
  segMLS = 30;

  gap=5;

  nfactor = 1;
//  simulation = 2;  // snr mode
//  factors[0]=10*sqrt(3);  
  simulation = 1;  // hrss mode
  factors[0]=1;  

  healpix = 7; 

  plugin = TMacro("macro/CWB_Plugin_BRST2.C");               // Macro source
  configPlugin = TMacro("macro/CWB_Plugin_BRST2_Config.C");  // Macro config
}
