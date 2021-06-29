
// this macro provides an example for the generation of the calibration file for CWB_Plugin_CE.C

#define SPCAL_NPTS	10		// number of frequency points where amp and phase errors are provided
//#define FREQ_MIN	20
#define FREQ_MIN	5		// min frequency 
#define FREQ_MAX	1024		// max frequency

#define nIFO		3		// number of detectors used for calibration file 

#define nENTRY		1		// number of entries in the calibration file

void WriteCalFileForPluginCE(TString clb_file="clb_file_ce.clb")  {

  TString ifo[NIFO] = {"h1","l1","v1"}; 	// WARNING: ifo must be defined with lower-case characters

  // check input name
  ofstream clb;
  if(!clb_file.Contains(".clb")) {
    cout << endl << "Error : calibration file extention must be .clb" << endl;
    exit(1);
  }

  // open file
  clb.open(clb_file.Data(),ios::out);
  if(!clb.good()) {
    cout << endl << "Error Opening Output Calibration File : " << clb_file.Data() << endl;
    exit(1);
  }
  clb.precision(14);

  // write header
  clb << "time\t" << "spcal_active\t" << "spcal_npts\t" << "simulation_id\t";
  for(int n=0;n<nIFO;n++) {
    for(int i=0;i<SPCAL_NPTS;i++) {
      clb << TString::Format("%s_spcal_freq_%d\t",ifo[n].Data(),i);
      clb << TString::Format("%s_spcal_amp_%d\t",ifo[n].Data(),i);
      clb << TString::Format("%s_spcal_phase_%d\t",ifo[n].Data(),i);
    }
  }
  clb << endl;

  // define frequency array in log scale
  double logFactor = (log(FREQ_MAX)-log(FREQ_MIN))/(SPCAL_NPTS-1);
  double logFreq[SPCAL_NPTS];
  logFreq[0] = log(FREQ_MIN);
  for(int i=1;i<SPCAL_NPTS;i++) logFreq[i]=logFreq[i-1]+logFactor;
  for(int i=0;i<SPCAL_NPTS;i++) cout << i << "\tFrequency: " << exp(logFreq[i]) << "\t logFrequency: " << logFreq[i] << endl;

  // init amp/phase calibration errors
  double spcal_amp[nIFO][SPCAL_NPTS];
  double spcal_phase[nIFO][SPCAL_NPTS];
  for(int n=0;n<nIFO;n++) {
    for(int i=0;i<SPCAL_NPTS;i++) {
      if(ifo[n]=="h1") {
        spcal_amp[n][i] = 0.2;				// amplitude error in percentage
        spcal_phase[n][i] = 90.0*TMath::Pi()/180.;	// phase errors in radians
      }
      if(ifo[n]=="l1") {
        spcal_amp[n][i] = 0.3;
        spcal_phase[n][i] = 45.0*TMath::Pi()/180.;
      }
      if(ifo[n]=="v1") {
        spcal_amp[n][i] = 0.0;
        spcal_phase[n][i] = 0.0*TMath::Pi()/180.;
      }
    }
  }

  // write entries
  double time=1267963151.;				// the Plugin_CE uses this GPS time to select the calibration entry
  int spcal_active = 1;					// must be 1
  int simulation_id = 0;				// not used
  for(int k=0;k<nENTRY;k++) {
    clb << time << " \t" << spcal_active << "\t" << SPCAL_NPTS << "\t" << simulation_id << "\t";
    for(int n=0;n<nIFO;n++) {
      for(int i=0;i<SPCAL_NPTS;i++) {
        clb << exp(logFreq[i]) << "\t";
        clb << spcal_amp[n][i] << "\t";
        clb << spcal_phase[n][i] << "\t";
      }
    }
    clb << endl;
  }

  clb.close();

  exit(0);
}
