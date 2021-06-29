//
// Draw Time Sky Segmentation
// Author : Gabriele Vedovato


#define ODIR_NAME "plots"
#define OFILE_EXT "png"
//#define OFILE_EXT "root"
//#define WRITE_PLOT

#define SAMPLE_RATE 8192.
//#define SAMPLE_RATE 16384.

#define RESOLUTION  2 
//#define RESOLUTION  4

//#define COORDINATES "cWB"
#define COORDINATES "Geographic"

#define PROJECTION ""
//#define PROJECTION "hammer"
//#define PROJECTION "sinusoidal"

//#define DISPLAY_WORLD_MAP
#define WORLD_MAP_DIR "$WB_GWAT/data/"

// polarization=0 -> |Fx|                     DPF
// polarization=1 -> |F+|                     DPF
// polarization=2 -> |Fx|/|F+|                DPF
// polarization=3 -> sqrt(|F+|^2+|Fx|^2)      DPF
// polarization=4 -> |Fx|^2                   DPF
// polarization=5 -> |F+|^2                   DPF
// polarization=6 -> Fx                       only with 1 detector
// polarization=7 -> F+                       only with 1 detector
// polarization=8 -> F1x/F2x                  only with 2 detectors
// polarization=9 -> F1+/F2+                  only with 2 detectors
// polarization=10 -> sqrt(|F1+|^2+|F1x|^2)/sqrt(|F2+|^2+|F2x|^2)     only with 2 detectors
// polarization=11 -> The same as (10) but averaged over psi          only with 2 detectors

using namespace CWB;
using namespace ROOT::Math;

void DrawTimeSkySegmentation(TString network="L1H1V1", int polarization = 1, int palette = 0, bool btitle = true) {

  if (!gROOT->GetClass("XYZVector")) gSystem->Load("libMathCore");

  int nIFO=0;
  TString ifo[10];
  if(network.Contains("V1")) ifo[nIFO++]="V1";   // VIRGO
  if(network.Contains("H1")) ifo[nIFO++]="H1";   // LHO1
  if(network.Contains("L1")) ifo[nIFO++]="L1";   // LLO
  if(network.Contains("G1")) ifo[nIFO++]="G1";   // GEO
  if(network.Contains("T1")) ifo[nIFO++]="T1";   // TAMA
  if(network.Contains("H2")) ifo[nIFO++]="H2";   // LHO2
  if(network.Contains("A1")) ifo[nIFO++]="A1";   // AIGO
  if(network.Contains("O1")) ifo[nIFO++]="O1";   // AURIGA
  if(network.Contains("N1")) ifo[nIFO++]="N1";   // NAUTILUS
  if(network.Contains("E1")) ifo[nIFO++]="E1";   // EXPLORER
  if(network.Contains("A2")) ifo[nIFO++]="A2";   // AUSTRALIAN 90°
  if(network.Contains("J1")) ifo[nIFO++]="J1";   // JAPANESE

  if(nIFO==0) {cout << "No detectors defined !!! " << endl;exit(1);}

  char ifostr[32]="";
  for(int n=0; n<nIFO; n++) {
    sprintf(ifostr,"%s %s",ifostr,ifo[n].Data());
  }
  cout << "Network : " << ifostr << endl;

  TString title;

  gnetwork* gNET = new gnetwork(nIFO,ifo);
  gskymap* gSM = gNET->GetGskymap();
  gSM->SetOptions(PROJECTION,COORDINATES,RESOLUTION/2);
//  gSM->SetOptions("LVC experiment", 300,40, 1200, 670);

  gNET->setSkyMaps(0.4,0,180,0,360);
  gNET->setAntenna();
  gNET->setDelay(const_cast<char*>(ifo[0].Data()));

#ifdef DISPLAY_WORLD_MAP
  gSM->SetWorldMap();
#endif

  TH2D* h2 = (TH2D*)gSM->GetHistogram();
  h2->GetXaxis()->SetTitleSize(0.05);
  h2->GetXaxis()->SetLabelSize(0.05);
  h2->GetYaxis()->SetTitleSize(0.05);
  h2->GetYaxis()->SetLabelSize(0.05);
// For CHRIS 
  h2->GetYaxis()->SetLabelFont(42);
  h2->GetYaxis()->SetLabelFont(42);
  h2->GetXaxis()->SetTitleFont(42);
  h2->GetYaxis()->SetTitleFont(42);

  if(polarization==3) {
    //if(nIFO>1) h2->GetZaxis()->SetRangeUser(0,1.7); 
    if(nIFO>1) h2->GetZaxis()->SetRangeUser(0,1.0); 
    else       h2->GetZaxis()->SetRangeUser(0,1.0);
  }
  if(polarization==2) h2->GetZaxis()->SetRangeUser(0,1.0);

  gNET->DrawAntennaPattern(polarization,palette,btitle);

  double mTau=gNET->getDelay("MAX");  // maximum time delay
  double ph1,ph2,th1,th2;
  double dt = 1./SAMPLE_RATE;
  int sTau=mTau/dt;
  for (int k=-sTau;k<sTau;k++) {
    for (int h=-sTau;h<sTau;h++) {
      gNET->Delay2Coordinates(0,k*dt,h*dt,ph1,th1,ph2,th2);
//      cout << ph1 << " " << th1 << " " << ph2 << " " << th2 << endl;
      if(COORDINATES=="Geographic") {
        CwbToGeographic(ph1,th1,ph1,th1);
        gSM->DrawMarker(ph1,th1, 1, 0.1, kBlack);  // Geographic
        CwbToGeographic(ph2,th2,ph2,th2);
        gSM->DrawMarker(ph2,th2, 1, 0.1, kBlack);  // Geographic
      }
    }
//exit(0); 
  }

#ifdef DRAW_EQUATORIAL
  double phi,theta;
  CwbToGeographic(EQUATORIAL_PHI,EQUATORIAL_THETA,phi,theta);
cout << phi << " " << theta << endl;
  gSM->DrawMarker(phi,theta, 4, 3.0, kRed);  // Geographic
#endif

#ifdef WRITE_PLOT
  char ofileName[128]=ODIR_NAME;
  sprintf(ofileName,"%s/",ofileName);
  for(int n=0; n<nIFO; n++) {
    sprintf(ofileName,"%s%s",ofileName,ifo[n].Data());
  }
  if(polarization==0) sprintf(ofileName,"%s%s.%s",ofileName,"_Fc",OFILE_EXT);
  if(polarization==1) sprintf(ofileName,"%s%s.%s",ofileName,"_Fp",OFILE_EXT);
  if(polarization==2) sprintf(ofileName,"%s%s.%s",ofileName,"_Fc_over_Fp",OFILE_EXT);
  if(polarization==3) sprintf(ofileName,"%s%s.%s",ofileName,"_Sqrt_Fp2_plus_Fc2",OFILE_EXT);

  cout << "Write : " << ofileName << endl;
  gSM->Print(ofileName);
  exit(0);
#endif
}

