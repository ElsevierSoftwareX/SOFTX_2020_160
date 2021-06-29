//
// Draw Antenna Pattern Componets for a single detector
// Author : Gabriele Vedovato

{
  #define PROJECTION ""
  //#define PROJECTION "hammer"
  #define RESOLUTION  2
  #define COORDINATES "Geographic"

  #define DISPLAY_WORLD_MAP
  #define WORLD_MAP_DIR "$CWB_GWAT/data/"

  //#define IFO	"L1"
  //#define IFO	"H1"
  #define IFO	"V1"

  #define FPLUS

  //#define OFILE_NAME	"DetectorAntennaPattern.png" 

  gskymap* gSM = new gskymap(int(6));
  gSM->SetOptions(PROJECTION,COORDINATES,RESOLUTION/2);

  detector D((char*)IFO);

  int L = gSM->size();
  for(int l=0;l<L;l++) {
    double th = gSM->getTheta(l);
    double ph = gSM->getPhi(l);
    wavecomplex F = D.antenna(th,ph,0);
#ifdef FPLUS
    gSM->set(l,F.real());
#else
    gSM->set(l,F.imag());
#endif
  }

#ifdef DISPLAY_WORLD_MAP
  TString world_map = gSystem->ExpandPathName(WORLD_MAP_DIR);
  gSM->SetWorldMapPath(world_map.Data());
  gSM->SetWorldMap();
#endif

#ifdef FPLUS
  TString title = TString("F_{+}");
#else
  TString title = TString("F_{x}");
#endif

  gSM->SetTitle(title);
  gSM->SetGridxColor(kWhite);
  gSM->SetGridyColor(kWhite);
  gSM->Draw(0);

#ifdef OFILE_NAME
  cout << "Write : " << OFILE_NAME << endl;
  gSM->Print(OFILE_NAME);
  exit(0);
#endif

}
