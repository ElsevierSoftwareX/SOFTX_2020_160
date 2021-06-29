//
// Draw Antenna Pattern for Builtin Detectors dump/load snetwork object
// Author : Gabriele Vedovato

{
  #define RESOLUTION  2
  #define COORDINATES "Geographic"
  #define PROJECTION "hammer"

  #define ROOT_FILE "DrawGnetwork3.root"

  //#define DUMP_OBJECT
  #define LOAD_OBJECT

  #ifdef DUMP_OBJECT

  int nIFO=3;
  TString ifo[3]={"L1","H1","J1"};

  int polarization=3;
  bool btitle=true;
  int palette=0;

  gnetwork gNET(3,ifo);
  
  gskymap* gSM = gNET.GetGskymap();
  gSM->SetOptions(PROJECTION,COORDINATES,RESOLUTION/2);
//  gSM->SetOptions("LVC experiment", 300,40, 1200, 670);

  gSM->SetWorldMap();
  gNET.DrawAntennaPattern(polarization,palette,btitle);
  gNET.DrawSitesShortLabel(kBlack);
  gNET.DrawSites(kBlack,2.0);
  gNET.DrawSitesArms(1000000,kWhite,3.0);

  gNET.DumpObject(const_cast<char*>(ROOT_FILE));
#endif

#ifdef LOAD_OBJECT
  gnetwork iNET;
  iNET.LoadObject(const_cast<char*>(ROOT_FILE));
  gskymap* iSM = iNET.GetGskymap();
  cout << "iSM size " << iSM->size() << endl;
  iSM->SetWorldMap();
  iSM->SetOptions(PROJECTION,COORDINATES,RESOLUTION/2);
  iSM->Draw();
  iNET.DrawSitesShortLabel(kBlack);
  iNET.DrawSites(kBlack,2.0);
  iNET.DrawSitesArms(1000000,kWhite,2.0);
#endif
}

