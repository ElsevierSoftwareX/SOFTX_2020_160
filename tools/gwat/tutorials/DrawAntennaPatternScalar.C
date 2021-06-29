//
// Draw Antenna Pattern for Builtin Detectors
// Author : Gabriele Vedovato

#define ODIR_NAME "."
#define OFILE_EXT "png"
//#define OFILE_EXT "root"
//#define WRITE_PLOT

#define RESOLUTION  2 
//#define RESOLUTION  4

//#define COORDINATES "cWB"
#define COORDINATES "Geographic"

#define PROJECTION ""
//#define PROJECTION "hammer"
//#define PROJECTION "sinusoidal"

#define DISPLAY_WORLD_MAP
#define WORLD_MAP_DIR "$CWB_GWAT/data/"

//#define HEALPIX

POLARIZATION GW=SCALAR;

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

TCanvas* 
DrawAntennaPatternScalar(TString network="L1H1V1", int polarization = 3, int palette = 0, bool btitle = true) {

  int nIFO=0;
  TString ifo[10];
  if(network.Contains("H1")) ifo[nIFO++]="H1";   // LHO1
  if(network.Contains("L1")) ifo[nIFO++]="L1";   // LLO
  if(network.Contains("G1")) ifo[nIFO++]="G1";   // GEO
  if(network.Contains("V1")) ifo[nIFO++]="V1";   // VIRGO
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

  gnetwork* gNET = new gnetwork;

  detector* pD[NIFO_MAX];
  for(int i=0; i<nIFO; i++) pD[i] = new detector((char*)ifo[i].Data()); // built in detector
  for(int i=0; i<nIFO; i++) gNET->add(pD[i]);

  gskymap* gSM = gNET->GetGskymap();
  gSM->SetOptions(PROJECTION,COORDINATES,RESOLUTION/2);
//  gSM->SetOptions("LVC experiment", 300,40, 1200, 670);

#ifdef DISPLAY_WORLD_MAP
  TString world_map = gSystem->ExpandPathName(WORLD_MAP_DIR);
  gSM->SetWorldMapPath(world_map.Data());
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
    else h2->GetZaxis()->SetRangeUser(0,1.0);
  }
  if(polarization==2) h2->GetZaxis()->SetRangeUser(0,1.0);

  if(GW==SCALAR) {
    for(int n=0;n<(int)gNET->ifoListSize();n++) {
      detector *d = gNET->getifo(n);
      d->setPolarization(SCALAR);
    }
#ifdef HEALPIX
    int healpix_order = 2;  // set the resolution of the HEALPix skygrid
    gNET->setSkyMaps((int)healpix_order);
    gNET->setAntenna();
#endif
  }
  //detector* det = net->getifo(0);
  //det->rotate(-70); 
  gNET->DrawAntennaPattern(polarization,palette,btitle);
  //gNET->DrawDelay("H1","L1");
  //cout << gNET->GetDelay("H1","L1",0,50) << endl;
  //cout << gNET->GetDelay("H1","L1",0,120) << endl;
  //gNET->DrawCircles(100,60,kWhite);
  //gNET->ClearCircles();
  //gNET->DrawSites(kBlue,1.0);
  //gNET->DrawSitesLabel(kBlue,0.05);
  //gNET->DrawSites(kBlack,2.0);
  //gNET->DrawSites(kBlack,2.5);
  gNET->DrawSitesShortLabel(kBlack);
  //gNET->DrawSitesLabel(kWhite,0.05);
  gNET->DrawSites(kBlack,2.0);
  gNET->DrawSitesArms(1000000,kWhite,3.0);
  if(GW==SCALAR) {
    if(polarization==1) {
      TString pTitle=TString("|F_{o}|");
      h2->SetTitle(TString("Network = ")+network+"          Antenna Pattern =  "+pTitle);
    }
    if(polarization==5 || polarization==3 || polarization==4) {
      TString pTitle=TString("|F_{o}|^{2}");
      h2->SetTitle(TString("Network = ")+network+"          Antenna Pattern =  "+pTitle);
    }
  }

#ifdef WRITE_PLOT
  char ofileName[128]=ODIR_NAME;
  sprintf(ofileName,"%s/",ofileName);
  for(int n=0; n<nIFO; n++) {
    sprintf(ofileName,"%s%s",ofileName,ifo[n].Data());
  }
  if(GW==SCALAR) { 
    if(polarization==1) sprintf(ofileName,"%s%s.%s",ofileName,"_Fo",OFILE_EXT);
    if(polarization==3) sprintf(ofileName,"%s%s.%s",ofileName,"_Fo",OFILE_EXT);
    if(polarization==4) sprintf(ofileName,"%s%s.%s",ofileName,"_Fo2",OFILE_EXT);
    if(polarization==5) sprintf(ofileName,"%s%s.%s",ofileName,"_Fo2",OFILE_EXT);
  } else {
    if(polarization==0) sprintf(ofileName,"%s%s.%s",ofileName,"_Fc",OFILE_EXT);
    if(polarization==1) sprintf(ofileName,"%s%s.%s",ofileName,"_Fp",OFILE_EXT);
    if(polarization==2) sprintf(ofileName,"%s%s.%s",ofileName,"_Fc_over_Fp",OFILE_EXT);
    if(polarization==3) sprintf(ofileName,"%s%s.%s",ofileName,"_Sqrt_Fp2_plus_Fc2",OFILE_EXT);
    if(polarization==4) sprintf(ofileName,"%s%s.%s",ofileName,"_Fc2",OFILE_EXT);
    if(polarization==5) sprintf(ofileName,"%s%s.%s",ofileName,"_Fp2",OFILE_EXT);
  }
  cout << "Write : " << ofileName << endl;
  gSM->Print(ofileName);
  exit(0);
#endif

  return gSM->GetCanvas();
}

