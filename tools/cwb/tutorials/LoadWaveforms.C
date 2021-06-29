//
// Load waveforms from files
// Author : Gabriele Vedovato

#define WAVEFORM_GA_DIR  "/home/waveburst/WAVEFORMS/GA"
#define WAVEFORM_SG_DIR  "/home/waveburst/WAVEFORMS/SG"
#define WAVEFORM_WNB_DIR "/home/waveburst/WAVEFORMS/WNB"
#define WFRATE	16384.
#define N_IFO    3

//#define READ_GA
#define READ_SG
//#define READ_WNB

#define DRAW_WF	"SG1304Q100"
//#define DRAW_WF	"WNB_1000_1000_0d01"
#define PRINT_WF

waveform GetWaveform(TString fName);
 
CWB::mdc* MDC;

void LoadWaveforms() {

  TString ifo[N_IFO] = {"L1","H1","V1"};

  MDC = new CWB::mdc(N_IFO,ifo);

  MDC->SetInjLength(1.5);

  vector<TString> fileList;

#ifdef READ_GA
  // read GA waveform file names
  fileList = CWB::Toolbox::getFileListFromDir(WAVEFORM_GA_DIR, ".txt", "GA");
  for(int n=0;n<fileList.size();n++) {          // loop over the root output files
    //cout << n << " " << fileList[n].Data() << endl;
    waveform wf = GetWaveform(fileList[n]);	// read waveform
    MDC->AddWaveform(wf); 			// add waveform
  }
#endif

#ifdef READ_SG
  // read SG waveform file names
  fileList = CWB::Toolbox::getFileListFromDir(WAVEFORM_SG_DIR, ".txt", "SG");
  for(int n=0;n<fileList.size();n++) {          // loop over the root output files
    //cout << n << " " << fileList[n].Data() << endl;
    waveform wf = GetWaveform(fileList[n]);	// read waveform
    MDC->AddWaveform(wf); 			// add waveform
  }
#endif

#ifdef READ_WNB
  cout << "read WNB waveform file names ..." << endl;
  fileList = CWB::Toolbox::getFileListFromDir(WAVEFORM_WNB_DIR, ".txt", "WNB");
  for(int n=0;n<fileList.size();n++) {          // loop over the root output files
    //cout << n << " " << fileList[n].Data() << endl;
    waveform wf = GetWaveform(fileList[n]);	// read waveform
    MDC->AddWaveform(wf); 			// add waveform
  }
  cout << "generate WNB with all possible hp,hx permutations ..." << endl;
  for(int i=0;i<MDC->wfList.size();i++) {
    cout << "extract all waveforms belonging to the WNB type " << MDC->wfList[i].name << " ..." << endl;
    int N = MDC->wfList[i].list.size()+1;	// number of waveforms belonging to the same WNB type
    waveform wf,jwf,kwf; 
    vector<waveform> vWF; 
    for(int j=0;j<N;j++) {
      for(int k=0;k<N;k++) {
        if(j!=k) {
  	  jwf = (j==0) ? MDC->wfList[i] : MDC->wfList[i].list[j-1]; 
  	  kwf = (k==0) ? MDC->wfList[i] : MDC->wfList[i].list[k-1]; 
  	  wf.name = jwf.name;
	  wf.type = MDC_USER;
	  wf.hpPath = jwf.hpPath;
	  wf.hxPath = kwf.hpPath;
          wf.hp=jwf.hp;
          wf.hx=kwf.hp;
          wf.par.resize(3);
          wf.par[1].name="hp"; wf.par[1].value=j;
          wf.par[2].name="hx"; wf.par[2].value=k;
          vWF.push_back(wf);
        }
      }
    }
    // add permutations waveform to MDC
    for(int j=0;j<vWF.size();j++) {

      vWF[j].par[0].name="wf"; vWF[j].par[0].value=j+1;	// fill par[0]

      if(j==0) MDC->wfList[i] = vWF[j]; 	// overwrite wf
      if(j>0)  MDC->AddWaveform(vWF[j]); 	// add waveform
    }
  }
#endif

  int nMDC = MDC->wfList.size();			// numer of MDC types

#ifdef PRINT_WF
  MDC->Print(1);
#endif

#ifdef DRAW_WF
  // Draw SG1304Q100 waveform
  for(int i=0;i<nMDC;i++) {
    if(MDC->wfList[i].name == DRAW_WF) {
      cout << i << " " << MDC->wfList[i].name.Data() 
           << " " << MDC->wfList[i].list.size()+1 << endl;

      MDC->Draw(i,0,"hp",MDC_TIME);
      MDC->Draw(i,0,"hx",MDC_TIME,"same",kRed);

      //wavearray<double> hp = MDC->wfList[i].hp;
      //MDC->Draw(hp);
      //MDC->Draw(hp,MDC_FFT);
      //MDC->Draw(hp,MDC_TF);
    }
  }
#endif

}

waveform 
GetWaveform(TString fName) {

  waveform wf;
  wavearray<double> x;
  vector<mdcpar> par(3);

  // extract MDC name 
  TString mdc_name = fName;
  mdc_name.Remove(0,mdc_name.Last('/')+1);
  if(mdc_name.Contains(".")) mdc_name.Remove(mdc_name.Last('.'));
  if(mdc_name.Contains("~")) mdc_name.Remove(mdc_name.Last('~'));
//  cout << "MDC NAME : " << mdc_name << endl;
  // read waveform
  MDC->ReadWaveform(x, fName, WFRATE);
  // build waveform
  wf.name = mdc_name;
  wf.type = MDC_USER;
  wf.hpPath = fName;
  wf.hp = x;
  wf.hxPath = "";
  wf.hx = wf.hp;
  wf.hx = 0;
  wf.par=par;

  return wf;
}
