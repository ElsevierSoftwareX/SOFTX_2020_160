//
// Draw Antenna Pattern for Builtin Detectors
// Author : Gabriele Vedovato

{
  #define WORLD_MAP_DIR "$CWB_GWAT/data/"

  #define RESOLUTION  2
  #define COORDINATES "Geographic"
  #define PROJECTION "hammer"

  int nIFO=3;
  TString ifo[3]={"L1","H1","J1"};

  int polarization=3;
  bool btitle=true;
  int palette=0;

  gnetwork gNET;

  detector* pD[3];
  for(int i=0; i<nIFO; i++) pD[i] = new detector((char*)ifo[i].Data()); // built in detector
  for(int i=0; i<nIFO; i++) gNET.add(pD[i]);
//  gNET.setSkyMaps(int(7));
//  gNET.setAntenna();
//  gNET.setDelay(const_cast<char*>(ifo[0].Data()));
  
  gskymap* gSM = gNET.GetGskymap();
  gSM->SetOptions(PROJECTION,COORDINATES,RESOLUTION/2);
//  gSM->SetOptions("LVC experiment", 300,40, 1200, 670);

  TString world_map = gSystem->ExpandPathName(WORLD_MAP_DIR);
  gSM->SetWorldMapPath(world_map.Data());
  gSM->SetWorldMap();
/*
  TH2D* h2 = (TH2D*)gSM->GetHistogram();
  h2->GetXaxis()->SetTitleSize(0.05);
  h2->GetXaxis()->SetLabelSize(0.05);
  h2->GetYaxis()->SetTitleSize(0.05);
  h2->GetYaxis()->SetLabelSize(0.05);
  h2->GetYaxis()->SetLabelFont(42);
  h2->GetYaxis()->SetLabelFont(42);
  h2->GetXaxis()->SetTitleFont(42);
  h2->GetYaxis()->SetTitleFont(42);

  if(polarization==3) h2->GetZaxis()->SetRangeUser(0,1.0); 
  if(polarization==2) h2->GetZaxis()->SetRangeUser(0,1.0);
*/
  gNET.DrawAntennaPattern(polarization,palette,btitle);
  gNET.DrawSitesShortLabel(kBlack);
  gNET.DrawSites(kBlack,2.0);
  gNET.DrawSitesArms(1000000,kWhite,3.0);

}

