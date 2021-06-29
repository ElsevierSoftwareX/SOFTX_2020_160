//
// This example show how to display skymap probability from fits file
// Author : Gabriele Vedovato


//#define COORDINATES "cWB"
//#define COORDINATES "Geographic"
#define COORDINATES "Celestial"

#define PROJECTION ""
//#define PROJECTION "hammer"
//#define PROJECTION "sinusoidal"

//#define PRINT_SKYMAP

void DrawFits(TString fitsName) {

  gskymap* gSM = new gskymap(fitsName.Data());
  gSM->SetOptions(PROJECTION,COORDINATES);
  gSM->Draw();
#ifdef PRINT_SKYMAP
  gSM->Print("probability_skymap.png");
#endif
}
