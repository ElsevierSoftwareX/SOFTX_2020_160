{
  // NET setting
  nIFO = 3;           // number of detectors from the array ifo[]

  cfg_search = 'I';   // see description in parameters.C (2G)
  optim=true;

  strcpy(ifo[0],"L1");
  strcpy(ifo[1],"H1");
  strcpy(ifo[2],"V1");

  strcpy(refIFO,"L1"); // reference IFO


  l_low    = 4;       // low frequency resolution level
  l_high   = 10;      // high frequency resolution level

  strcpy(wdmXTalk,"wdmXTalk/OverlapCatalog16-1024.bin");

  healpix = 7;

  bpp   = 0.001;      // probability for black pixel selection (netpixel)
  subnet = 0.6;       // [<0.7] subnetwork threshold (supercluster)
  netRHO = 5.0;       // [>4.0] coherent network SNR (supercluster, likelihood)
  netCC = 0.5;        // network correlation (supercluster, likelihood)
  Acore = 1.7;
  Tgap  = 0.0;        // defragmentation time gap between clusters (sec)
  Fgap  = 0.0;        // defragmentation frequency gap between clusters (Hz)
  delta = 0.5;        // [0-1] regularize sky locations with low |fx|^2,
  cfg_gamma = 0.5;    // [0-1] defines threshold on coherent energy: gamma=0 - Ec>0
  LOUD = 300;         // number of pixel per cluster to load TD amplitudes

  segLen = 600.;      // Segment length [sec]

  // simulation parameters
  simulation = 1;          // 1 for simulation, 0 for production
  nfactor = 1;             // number of strain factors
  factors[0]=1.0;

  // channels
  strcpy(channelNamesRaw[0],"L1:LDAS-STRAIN");
  strcpy(channelNamesRaw[1],"H1:LDAS-STRAIN");
  strcpy(channelNamesRaw[2],"V1:h_16384Hz");

  // plugin
  plugin       = TMacro("macro/CWB_Plugin_InjectBigDog.C");         // Macro source
  configPlugin = TMacro("macro/CWB_Plugin_InjectBigDog_Config.C");  // Macro config
}
