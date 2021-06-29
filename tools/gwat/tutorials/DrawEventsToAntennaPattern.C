//
// Draw Events & Antenna Pattern
// Author : Gabriele Vedovato


#define RESOLUTION  2 
//#define RESOLUTION  4

//#define COORDINATES "cWB"
#define COORDINATES "Geographic"

#define PROJECTION ""
//#define PROJECTION "hammer"
//#define PROJECTION "sinusoidal"

#define WRITE_PLOT


#define IFILE_NAME "merge/wave_ADV_SIM_BRST_LF_L1H1V1_2GM_run51.M1.root"
#define DRAW_EVENTS

using namespace CWB;

void DrawEventsToAntennaPattern() {

  int nIFO=3;
  TString ifo[3]={"L1","H1","V1"};
  
  gnetwork* gNET = new gnetwork;

  detector* pD[3];
  for(int i=0; i<nIFO; i++) pD[i] = new detector((char*)ifo[i].Data()); // built in detector
  for(int i=0; i<nIFO; i++) gNET->add(pD[i]);

  gskymap* gSM = gNET->GetGskymap();
  gSM->SetOptions(PROJECTION,COORDINATES,RESOLUTION/2);


  TH2D* h2 = (TH2D*)gSM->GetHistogram();
  h2->GetXaxis()->SetTitleSize(0.05);
  h2->GetXaxis()->SetLabelSize(0.05);
  h2->GetYaxis()->SetTitleSize(0.05);
  h2->GetYaxis()->SetLabelSize(0.05);
  h2->GetYaxis()->SetLabelFont(42);
  h2->GetYaxis()->SetLabelFont(42);
  h2->GetXaxis()->SetTitleFont(42);
  h2->GetYaxis()->SetTitleFont(42);

  h2->GetZaxis()->SetRangeUser(0,1.0);

  gNET->DrawAntennaPattern(0);

#ifdef DRAW_EVENTS
  TFile *ifile = TFile::Open(IFILE_NAME);
  if(ifile==NULL) {cout<<"Error opening file : "<<IFILE_NAME<<endl;exit(1);}
  TTree* itree = (TTree *) gROOT->FindObject("waveburst");
  if(itree==NULL) {cout<<"Error opening tree : "<<"waveburst"<<endl;exit(1);}
  itree->Draw("theta[1]:phi[1]","erA[0]>60","goff");
  int isize=itree->GetSelectedRows();
  cout << "isize : " << isize << endl;
  double* itheta = itree->GetV1();
  double* iphi = itree->GetV2();

  for (int i=0;i<isize;i++) {
    if(COORDINATES=="cWB") {
      gSM->DrawMarker(iphi[i], itheta[i], 20, 0.1, kBlack);  // cWB
    }
    if(COORDINATES=="Geographic") {
      double phi,theta;
      CwbToGeographic(iphi[i],itheta[i],phi,theta);
      gSM->DrawMarker(phi,theta, 20, 0.5, kBlack);  // Geographic
    }
  }
#endif

  gSM->GetCanvas()->Update();

#ifdef WRITE_PLOT

  TObjArray* token = TString(IFILE_NAME).Tokenize(TString('/'));
  TObjString* sfile = (TObjString*)token->At(token->GetEntries()-1);
  TString TITLE = sfile->GetString();
  TString ofile = sfile->GetString();
  ofile.ReplaceAll(".root","_EventsVsAntPat.png");

  cout << "Write : " << ofile << endl;
  gSM->Print(ofile);
  exit(0);
#endif
}

