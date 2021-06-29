{
  strcpy(analysis,"1G");

  nIFO = 2;
  search = 'r';
  strcpy(ifo[0],"L1");
  strcpy(ifo[1],"H1");
  strcpy(refIFO,"L1");

  strcpy(injectionList,"input/BRST_LF_L1H1V1.inj");

  detectorParams detPARAMS[nIFO] = {
                          {"L1", 30.5629,  -90.7742, 0.0, 0, ( +90-197.716),    0, (    -197.716 )},  // L1 
                          {"H1", 46.4551, -119.408,  0.0, 0, ( +90-125.998),    0, (   -125.998 )},   // H1
                         };
  for(int i=0;i<nIFO;i++) detParms[i] = detPARAMS[i];

  lagSize = 1;
  lagOff = 0;
  segTHR = 0;

  nDQF=2;
  dqfile dqf[nDQF]={
                    {"L1" ,"input/burst.in",           CWB_CAT1, 0., false, false},
                    {"H1" ,"input/burst.in",           CWB_CAT1, 0., false, false},
                   };
  for(int i=0;i<nDQF;i++) DQF[i]=dqf[i];
  cedRHO = 3.000000;
  fLow=64;

  gap=5;

  nfactor = 1;
  simulation = 2;  // snr mode
  factors[0]=10*sqrt(2);  

  // weak regulator 
  delta = 0.0;
  gamma = 0.0;          

  healpix = 7; 

  sprintf(skyMaskCCFile,"input/SkyMaskCC_GWGC_Rev1d8_HPX7_PAD.txt");

  plugin = TMacro("macro/CWB_Plugin_BRST_LF_GWGC.C");        		// Macro source
  configPlugin = TMacro("macro/CWB_Plugin_BRST_LF_GWGC_Config.C");  	// Macro config
}
