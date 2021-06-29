//
// Draw Antenna Pattern for user defined detectors
// Author : Gabriele Vedovato

#define ODIR_NAME "."
#define OFILE_EXT "png"
//#define OFILE_EXT "root"
//#define WRITE_PLOT

#define RESOLUTION  2 
//#define RESOLUTION  4

//#define COORDINATES "cWB"
#define COORDINATES "Geographic"

//#define PROJECTION ""
#define PROJECTION "hammer"
//#define PROJECTION "sinusoidal"

#define DISPLAY_WORLD_MAP
#define WORLD_MAP_DIR "$CWB_GWAT/data/"


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

void DrawAntennaPatternUdD(int polarization = 1, int palette = 0, bool btitle = true) {

  /*
  // ------------------------------------------------
  // 3 ET at cartesian vertexes
  // ------------------------------------------------
  int nIFO=9;

  detectorParams dP[9] = {{"X1", 0.0,  0.0, 0.0, 0,   90, 0,   30},
                          {"X2", 0.0,  0.0, 0.0, 0,  150, 0, -150},
                          {"X3", 0.0,  0.0, 0.0, 0,  -90, 0,  -30},

                          {"Y1", 0.0, 90.0, 0.0, 0,   90, 0,   30},
                          {"Y2", 0.0, 90.0, 0.0, 0,  150, 0, -150},
                          {"Y3", 0.0, 90.0, 0.0, 0,  -90, 0,  -30},

                          {"Z1", 90.0, 0.0, 0.0, 0,   90, 0,   30},
                          {"Z2", 90.0, 0.0, 0.0, 0,  150, 0, -150},
                          {"Z3", 90.0, 0.0, 0.0, 0,  -90, 0,  -30},
                         };
  */
  /*
  // ------------------------------------------------
  // 3 ET located on big circle (triangulation simmetry)
  // ------------------------------------------------
  int nIFO=5;
  
  detectorParams dP[9] = {{"X1", 0.0,  0.0, 0.0, 0,   90, 0,   30},
                          {"X2", 0.0,  0.0, 0.0, 0,  150, 0, -150},
                          {"X3", 0.0,  0.0, 0.0, 0,  -90, 0,  -30},
 
                          {"Y1", 0.0, 120.0, 0.0, 0,   90, 0,   30},
                          {"Y2", 0.0, 120.0, 0.0, 0,  150, 0, -150},
                          {"Y3", 0.0, 120.0, 0.0, 0,  -90, 0,  -30},
  
                          {"Z1", 0.0,-120.0, 0.0, 0,   90, 0,   30},
                          {"Z2", 0.0,-120.0, 0.0, 0,  150, 0, -150},
                          {"Z3", 0.0,-120.0, 0.0, 0,  -90, 0,  -30},
                         };  
  */
  /* 
  // ------------------------------------------------
  // 3 ET at L1,H1,V1 sites
  // ------------------------------------------------
  int nIFO=5;
  
                         // V1
  detectorParams dP[9] = {{"X1", 43.6314,  10.5045, 0.0, 0,   90, 0,   30},
                          {"X2", 43.6314,  10.5045, 0.0, 0,  150, 0, -150},
                          {"X3", 43.6314,  10.5045, 0.0, 0,  -90, 0,  -30},
  
                         // L1
                          {"Y1", 30.5629, -90.7742, 0.0, 0,   90, 0,   30},
                          {"Y2", 30.5629, -90.7742, 0.0, 0,  150, 0, -150},
                          {"Y3", 30.5629, -90.7742, 0.0, 0,  -90, 0,  -30},
  
                         // H1
                          {"Z1", 46.4551, -119.408, 0.0, 0,   90, 0,   30},
                          {"Z2", 46.4551, -119.408, 0.0, 0,  150, 0, -150},
                          {"Z3", 46.4551, -119.408, 0.0, 0,  -90, 0,  -30},
                         };
  
  */ 
  /*
  // ------------------------------------------------
  // 3 LCT at L1,H1,V1 sites
  // ------------------------------------------------
  int nIFO=9;
  
                         // V1
  detectorParams dP[9] = {{"X1", 43.6314,  10.5045, 0.0, 0, ( +90-70.5675), 0, (    -70.5675)},
                          {"X2", 43.6314,  10.5045, 0.0, 0, (+180-70.5675), 0, (+135-70.5675)},
                          {"X3", 43.6314,  10.5045, 0.0, 0, ( -45-70.5675), 0, ( -90-70.5675)},
  
                         // L1
                          {"Y1", 30.5629, -90.7742, 0.0, 0, ( +90-197.716), 0, (    -197.716)},
                          {"Y2", 30.5629, -90.7742, 0.0, 0, (+180-197.716), 0, (+135-197.716)},
                          {"Y3", 30.5629, -90.7742, 0.0, 0, ( -45-197.716), 0, ( -90-197.716)},
  
                         // H1
                          {"Z1", 46.4551, -119.408, 0.0, 0, ( +90-125.999), 0, (    -125.999)},
                          {"Z2", 46.4551, -119.408, 0.0, 0, (+180-125.999), 0, (+135-125.999)},
                          {"Z3", 46.4551, -119.408, 0.0, 0, ( -45-125.999), 0, ( -90-125.999)},
                         };
  */
  // ------------------------------------------------
  // 3 LCI at L1,H1,V1 sites
  // ------------------------------------------------
  /*
  int nIFO=5;
  
                         // V1
  detectorParams dP[9] = {{"X1", 43.6314,  10.5045, 0.0, 0, ( +90-70.5675), 0, (    -70.5675)},
                          {"X2", 43.6314,  10.5045, 0.0, 0, ( +90-70.5675), 0, ( +45-70.5675)},
                          {"X3", 43.6314,  10.5045, 0.0, 0, ( +45-70.5675), 0, (    -70.5675)},
  
                         // L1
                          {"Y1", 30.5629, -90.7742, 0.0, 0, ( +90-197.716), 0, (    -197.716)},
                          {"Y2", 30.5629, -90.7742, 0.0, 0, ( +90-197.716), 0, ( +45-197.716)},
                          {"Y3", 30.5629, -90.7742, 0.0, 0, ( +45-197.716), 0, (    -197.716)},
  
                         // H1
                          {"Z1", 46.4551, -119.408, 0.0, 0, ( +90-125.999), 0, (    -125.999)},
                          {"Z2", 46.4551, -119.408, 0.0, 0, ( +90-125.999), 0, ( +45-125.999)},
                          {"Z3", 46.4551, -119.408, 0.0, 0, ( +45-125.999), 0, (    -125.999)},
                         };

  */
  /*
  // ------------------------------------------------
  // 3 LCI at L1,H1,V1 sites
  // ------------------------------------------------
  int nIFO=5;

                         // V1
  detectorParams dP[9] = {{"X1", 43.6314,  10.5045, 0.0, 0, ( +90-70.5675), 0, (    -70.5675)},
                          {"X2", 43.6314,  10.5045, 0.0, 0, ( +90-70.5675), 0, ( +45-70.5675)},
                          {"X3", 43.6314,  10.5045, 0.0, 0, ( +45-70.5675), 0, (    -70.5675)},

                         // L1
                          {"Y1", 30.5629, -90.7742, 0.0, 0, ( +90-197.716), 0, (    -197.716)},
                          {"Y2", 30.5629, -90.7742, 0.0, 0, ( +90-197.716), 0, ( +45-197.716)},
                          {"Y3", 30.5629, -90.7742, 0.0, 0, ( +45-197.716), 0, (    -197.716)},

                         // H1
                          {"Z1", 46.4551, -119.408, 0.0, 0, ( +90-125.999), 0, (    -125.999)},
                          {"Z2", 46.4551, -119.408, 0.0, 0, ( +90-125.999), 0, ( +45-125.999)},
                          {"Z3", 46.4551, -119.408, 0.0, 0, ( +45-125.999), 0, (    -125.999)},
                         };
  */
  /*
  int nIFO=4;
  detectorParams dP[4] = {
                          {"L2", 30.5629,  -90.7742, 0.0, 0, ( +90-197.716),    0, ( +45-197.716 )},  // L2 LCI
                          {"L3", 30.5629,  -90.7742, 0.0, 0, ( +45-197.716),    0, (    -197.716 )},  // L3 LCI

                          {"H1", 46.4551,  -119.408, 0.0, 0, ( +90-125.999),    0, (     -125.999)},  // H1 LCI

                          {"V1", 43.6314,   10.5045, 0.0, 0, ( +90-70.5675),    0, (     -70.5675)},  // V1 LCI
                         };
  */
  /*
  int nIFO=4;
  detectorParams dP[4] = {
                          {"L1", 30.5629,  -90.7742, 0.0, 0, ( +90-197.716),    0, (    -197.716 )},  // L1

                          {"H2", 46.4551,  -119.408, 0.0, 0, ( +90-125.999),    0, ( +45-125.999)},   // H2 LCI
                          {"H3", 46.4551,  -119.408, 0.0, 0, ( +45-125.999),    0, (    -125.999)},   // H3 LCI

                          {"V1", 43.6314,   10.5045, 0.0, 0, ( +90-70.5675),    0, (    -70.5675)},   // V1
                         };
  */
  /*
  int nIFO=3;
  detectorParams dP[3] = {
                          {"L1", 30.5629,  -90.7742, 0.0, 0, ( +90-197.716),    0, (    -197.716 )},  // L1

                          {"H1", 46.4551,  -119.408, 0.0, 0, ( +90-125.999),    0, (    -125.999)},   // H1

                          {"V1", 43.6314,   10.5045, 0.0, 0, ( +90-70.5675),    0, (    -70.5675)},   // V1
                         };
  */
  /*
  int nIFO=2;
  detectorParams dP[2] = {
                          {"L2", 30.5629,  -90.7742, 0.0, 0, ( +90-197.716),    0, ( +45-197.716 )},  // L2 LCI
                          {"L3", 30.5629,  -90.7742, 0.0, 0, ( +45-197.716),    0, (    -197.716 )},  // L3 LCI
                         };
  */
  /*
  int nIFO=2;
  detectorParams dP[2] = {
                          {"L1", 30.5629,  -90.7742, 0.0, 0, ( +90-197.716),    0, (    -197.716 )},  // L1 LCI
                          {"L3", 30.5629,  -90.7742, 0.0, 0, ( 90+45-197.716),    0, (   45 -197.716 )},  // L3 LCI
                         };
  */
   /*  
  int nIFO=6;
  detectorParams dP[6] = {

                          {"H1", 46.4551,  -119.408, 0.0, 0, ( +90-125.999),    0, (    -125.999)},   // H2 LCI
                          {"H4", 46.4551,  -119.408, 0.0, 0, ( +90- 80.999),    0, ( +45-125.999)},   // H3 LCI

                          {"V1", 43.6314,   10.5045, 0.0, 0, ( +90-70.5675),    0, (    -70.5675)},   // V1
                          {"V4", 43.6314,   10.5045, 0.0, 0, ( +90-25.5675),    0, (    -25.5675)},   // V1

                          {"I1", 14.4,      76.4,    0.0, 0, ( +90+0.0    ),    0, (    +0.0    )},   // I1
                          {"I4", 14.4,      76.4,    0.0, 0, ( +90+45.0    ),    0, (  +45.0    )},   // I1
                         };
   */  
   /*
  int nIFO=4;
  detectorParams dP[4] = {

                          {"L1", 30.5629,  -90.7742, 0.0, 0, ( +90-197.716),    0, (    -197.716 )},  // L1 LCI
                          {"H1", 46.4551,  -119.408, 0.0, 0, ( +90-125.999),    0, (    -125.999)},   // H2 LCI
			  {"H4", 46.4551,  -119.408, 0.0, 0, ( +90- 80.999),    0, ( +45-125.999)},   // H3 LCI
                          {"V1", 43.6314,   10.5045, 0.0, 0, ( +90-70.5675),    0, (    -70.5675)},   // V1
                          {"V4", 43.6314,   10.5045, 0.0, 0, ( +90-25.5675),    0, (    -25.5675)},   // V1
                          {"I1", 14.4,      76.4,    0.0, 0, ( +90+0.0    ),    0, (    +0.0    )},   // I1
                          {"I4", 14.4,      76.4,    0.0, 0, ( +90+55.0    ),    0, (  +55.0    )},   // I1
                         };  
   */ 
   
  int nIFO=4;
  detectorParams dP[4] = {
     {"L1", 30.5629,  -90.7742, 0.0, 0, ( +90-197.716),    0, (    -197.716 )},  // L1	

     {"S1", 55.4,      82.6,    0.0, 0, ( +90-10.0    ),    0, (   -10.0    )},   // S1 
			 
     {"S4", 55.4,      82.6,    0.0, 0, ( +90+35.0    ),   0, (  +35.0     )},   // S1
			  
     {"H1", 46.4551,  -119.408,  0.0, 0, ( +90-125.999),   0, (     -125.999)},  // H1
     /*
     {"V1", 43.6314,   10.5045, 0.0, 0, ( +90-70.5675),    0, (    -70.5675)},   // V1
     */
  };
   
     //			  {"S4", 55.4,      82.6,    0.0, 0, ( +90+45.0    ),   0, (  +45.0     )},   // S1
  //			 {"H4", 46.4551,  -119.408,  0.0, 0, ( +90- 80.999),    0, ( +45-125.999)}.   // H3 LCI

 
  detector* pD[NIFO_MAX];
  for(int n=0;n<nIFO;n++) pD[n] = new detector(dP[n]);
  for(int n=0;n<nIFO;n++) cout << n << " " << pD[n]->Name << endl;


  gnetwork* gNET = new gnetwork(nIFO,NULL,dP);

  gskymap* gSM = gNET->GetGskymap();
  gSM->SetOptions(PROJECTION,COORDINATES,RESOLUTION/2);
//  gSM->SetOptions("LVC experiment", 300,40, 1200, 670);

#ifdef DISPLAY_WORLD_MAP
  TString world_map = gSystem->ExpandPathName(WORLD_MAP_DIR);
  gSM->SetWorldMapPath(world_map.Data());
  gSM->SetWorldMap();
#endif

  TString title;

  TH2D* h2 = (TH2D*)gSM->GetHistogram();
  h2->GetXaxis()->SetTitleSize(0.05);
  h2->GetXaxis()->SetLabelSize(0.05);
  h2->GetYaxis()->SetTitleSize(0.05);
  h2->GetYaxis()->SetLabelSize(0.05);
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

  //detector* det = gNET->getifo(0);
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

  //gNET->GetSite("X1");

  cout << "SkyMap Mean : " << gSM->mean() << endl;

#ifdef WRITE_PLOT
  char ofileName[128]=ODIR_NAME;
  sprintf(ofileName,"%s/",ofileName);
  for(int n=0; n<nIFO; n++) {
    sprintf(ofileName,"%s%s",ofileName,pD[n]->Name);
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

