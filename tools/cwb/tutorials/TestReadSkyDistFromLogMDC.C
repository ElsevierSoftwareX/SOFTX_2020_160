// Read SkyDistribution From Log MDC file
// Author : Gabriele Vedovato
//

{
  #define N_IFO 3
  #define LOG_FILE "/home/vedovato/S6/coherent/offline/SIMULATION/ADV_SIM_BRST6_S6a_LS2_SGQ3_L1H1V1_run100rws/BurstMDC-BRST6_S6a_LS2_SGQ3-Log.txt"

  TString ifo[N_IFO]={"L1","H1","V1"};
  CWB::mdc MDC(N_IFO,ifo); 

  vector<mdcpar> par;
  MDC.SetSkyDistribution(MDC_LOGFILE,LOG_FILE,par);

  //MDC.DumpLog("TEST_LOG_INJ-Log.root");

  //CWB::mdc MDC2 = MDC; 
  //MDC2.DumpLog("TEST_LOG_INJ-Log2.root");
  


  double theta; double phi; double psi; double rho;
  double xiota; double hrss; int ID; int id;
  for(int n=0;n<10;n++) {
    MDC.GetSourceCoordinates(theta, phi, psi, rho, xiota, hrss, ID, id);
    cout << "Source : " << theta << " " << phi << " " << psi << " " << rho << endl;
  }

  exit(0);
}
