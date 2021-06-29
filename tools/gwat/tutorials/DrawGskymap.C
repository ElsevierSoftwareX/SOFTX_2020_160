// This example shows how to produce a skymap plot

//#define SAVE_PLOT

#define RESOLUTION  1
//#define RESOLUTION  2

//#define COORDINATES "cwb"
#define COORDINATES "geographic"


TCanvas* DrawGskymap(TString projection="cartesian") {

  // create gskymap with HEALPix order=7
  gskymap* gSM = new gskymap((int)7);

  // set gskymap options
  gSM->SetOptions(projection,COORDINATES,RESOLUTION);

  // set title
  TString title;
  if(projection=="cartesian")  title="Projection : cartesian"; 
  if(projection=="hammer")     title="Projection : hammer"; 
  if(projection=="parabolic")  title="Projection : parabolic"; 
  if(projection=="sinusoidal") title="Projection : sinusoidal"; 

  if(TString(COORDINATES)=="cwb")        title=title+"  -  Coordinates : cwb"; 
  if(TString(COORDINATES)=="geographic") title=title+"  -  Coordinates : geographic"; 
  if(TString(COORDINATES)=="celestial")  title=title+"  -  Coordinates : celestial"; 

  gSM->SetTitle(title);

  // set world map
  gSM->SetWorldMap();

  // draw skymap (exclude TPaletteAxis)
  gSM->Draw(0,"col");

#ifdef SAVE_PLOT
  TString fname;
  if(projection=="cartesian")  fname="gskymap_cartesian_plot.png"; 
  if(projection=="hammer")     fname="gskymap_hammer_plot.png"; 
  if(projection=="parabolic")  fname="gskymap_parabolic_plot.png"; 
  if(projection=="sinusoidal") fname="gskymap_sinusoidal_plot.png"; 

  if(TString(COORDINATES)=="cwb")        fname.ReplaceAll("gskymap_","gskymap_cwb_"); 
  if(TString(COORDINATES)=="geographic") fname.ReplaceAll("gskymap_","gskymap_geographic_"); 
  if(TString(COORDINATES)=="celestial")  fname.ReplaceAll("gskymap_","gskymap_celestial_"); 

  gSM->Print(fname);
  cout << "Write : " << fname << endl;
  exit(0);
#endif

  return gSM->GetCanvas();  // used by THtml
}
