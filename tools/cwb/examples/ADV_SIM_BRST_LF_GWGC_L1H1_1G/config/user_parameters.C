{
  strcpy(analysis,"1G");

  nIFO = 2;
  cfg_search = 'r';
  strcpy(ifo[0],"L1");
  strcpy(ifo[1],"H1");
  strcpy(refIFO,"L1");

  strcpy(injectionList,"input/BRST_LF_L1H1V1.inj");

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
  simulation = 2;               // snr mode
  factors[0]=10*sqrt(2);  	// fixed injected network snr ~ 17

  healpix = 7;	   // use healpix sky segmentation
  Psave = true;    // save skymap probability to output root file
  nSky = -99;      // save pixels skymap probability pixels up to cumulative prob < 0.99
                   // NOTE : select a probability threshold not too high to avoid a big number of unusefull pixels
  dump=false;      // disable output to output/*.txt file (only output/*.root file)

  plugin = TMacro("macro/CWB_Plugin_BRST_LF_GWGC.C");        		// Macro source
  configPlugin = TMacro("macro/CWB_Plugin_BRST_LF_GWGC_Config.C");  	// Macro config
}
