//
// Draw Delay Time skymap for two detectors
// Author : Gabriele Vedovato

#define ODIR_NAME "plots"
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

//#define DISPLAY_WORLD_MAP
#define WORLD_MAP_DIR "/aufs/data9/users/vedovato/S5-VSR1-Y2/lib/skyplot/"

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


#define SKYRES 0.4

#define SAMPLE_RATE 16384.

#define WAT_LEVEL 3
#define WAT_FILTER "up2"
#define WAT_FILTER_LENGTH 32

void DrawDelayIndex(TString network="H1V1") {

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

  gnetwork* gNET = new gnetwork;

  detector* pD[3];
  for(int i=0; i<nIFO; i++) pD[i] = new detector((char*)ifo[i].Data()); // built in detector
  for(int i=0; i<nIFO; i++) gNET->add(pD[i]);

  gNET->setSkyMaps(SKYRES,0,180,0,360);
  gNET->setAntenna();
  gNET->setDelay(const_cast<char*>(ifo[0].Data()));

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
  h2->GetYaxis()->SetLabelFont(42);
  h2->GetYaxis()->SetLabelFont(42);
  h2->GetXaxis()->SetTitleFont(42);
  h2->GetYaxis()->SetTitleFont(42);

  //h2->GetZaxis()->SetRangeUser(0,1.0);

  //gNET->DrawAntennaPattern(polarization,palette,btitle);
  //gNET->DrawDelay("H1","L1");
  //cout << gNET->GetDelay("H1","L1",0,50) << endl;
  //cout << gNET->GetDelay("H1","L1",0,120) << endl;
  //gNET->DrawCircles(100,60,kWhite);
  //gNET->ClearCircles();
  //gNET->DrawSites(kBlue,1.0);
  //gNET->DrawSitesLabel(kBlue,0.05);
  //gNET->DrawSites(kBlack,2.0);
  //gNET->DrawSites(kBlack,2.5);
//  gNET->DrawSitesShortLabel(kBlack);
  //gNET->DrawSitesLabel(kWhite,0.05);
//  gNET->DrawSites(kBlack,2.0);

  skymap sm(SKYRES,0,180,0,360);
  int L = sm.size();

  double mTau=gNET->getDelay("MAX");
  cout<<"maximum time delay between detectors: "<<mTau<<endl;
  char tdf00[256];
  sprintf(tdf00,"%s/Meyer1024wat482_00%s_L%1d.dat",gSystem->ExpandPathName("$HOME_WAT_FILTERS"),WAT_FILTER,WAT_LEVEL);
  gNET->setDelayFilters(tdf00);
  gNET->setDelayIndex();
  gNET->setIndexMode(0);
  for(int n=0;n<nIFO;n++) pD[n]=gNET->getifo(n);

//  cout << pD[0]->index.data[0] << endl;

  for(int l=0;l<L;l++) sm.set(l,pD[1]->index.data[l]);

  gskymap* smd = new gskymap(SKYRES,0,180,0,360);
  for(int l=0;l<L;l++) {
    double th = sm.getTheta(l);
    double ph = sm.getPhi(l);
    double delay=pD[1]->getTau(th,ph)-pD[0]->getTau(th,ph);
    //double sample_rate =  SAMPLE_RATE/(1<<WAT_LEVEL);
    //double sample_rate =  SAMPLE_RATE/WAT_LEVEL;
    double sample_rate =  SAMPLE_RATE;
    double idelay=pD[1]->index.data[l]-WAT_FILTER_LENGTH;
    idelay/=sample_rate;
    //smd->set(l,-delay);
    //smd->set(l,idelay);
    //smd->set(l,idelay-mTau);
    smd->set(l,idelay-mTau+delay);
  }
  smd->Draw();

//  gNET->DrawDelay("V1","H1");

  gNET->DrawSitesShortLabel(kBlack);
  //gNET->DrawSitesLabel(kWhite,0.05);
  gNET->DrawSites(kBlack,2.0);

}

