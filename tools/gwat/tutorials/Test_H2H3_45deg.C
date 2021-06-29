//
// Test L1H2H3V1 with H2 H3 @ 45deg
// Author : Gabriele Vedovato

#define WORLD_MAP_DIR "$CWB_GWAT/data/"

#define RESOLUTION  2
//#define RESOLUTION  4

//#define COORDINATES "cWB"
#define COORDINATES "Geographic"

//#define PROJECTION ""
#define PROJECTION "hammer"
//#define PROJECTION "sinusoidal"

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

#define H2_ELEVATION       1000000	// H2 elevation in meters

#define N_IFO 4 

void Test_H2H3_45deg(int polarization=3) {

  // define network

  char ifo[NIFO_MAX][8];

  strcpy(ifo[0],"L1");
  strcpy(ifo[1],"");
  strcpy(ifo[2],"");
  strcpy(ifo[3],"V1");

  detectorParams detParms[N_IFO] = {
              {"L1", 30.5629,  -90.7742, 0.0,          0, ( +90-197.716),    0, (    -197.716 )},  // L1
              {"H2", 46.4551, -119.408,  H2_ELEVATION, 0, ( +90-125.998),    0, (+45-125.998 )},   // H2
              {"H3", 46.4551, -119.408,  0.0,          0, ( +45-125.998),    0, (   -125.998 )},   // H3
              {"V1", 43.6314,   10.5045, 0.0,          0, ( +90-70.5675),    0, (    -70.5675)},   // V1
             };

  detector* pD[NIFO_MAX];                                     //! pointers to detectors
  for(int i=0; i<N_IFO; i++) {
    if(strlen(ifo[i])>0) pD[i] = new detector(ifo[i]);        // built in detector
    else                 pD[i] = new detector(detParms[i]);   // user define detector
  }

  gnetwork* gNET = new gnetwork;
  for(int i=0; i<N_IFO; i++) gNET->add(pD[i]);
  for(int n=0;n<N_IFO;n++) gNET->GetSite(detParms[n].name);

  // display network antenna pattern

  gskymap* gSM = gNET->GetGskymap();
  gSM->SetOptions(PROJECTION,COORDINATES,RESOLUTION/2);

  TString world_map = gSystem->ExpandPathName(WORLD_MAP_DIR);
  gSM->SetWorldMapPath(world_map.Data());
  gSM->SetWorldMap();

  gNET->DrawAntennaPattern(polarization,0,true);
  gNET->DrawSitesShortLabel(kBlack);
  gNET->DrawSites(kBlack,2.0);
  gNET->DrawSitesArms(1000000,kWhite,3.0);

  // print delays (theta,phi is the H2,H3 zenit direction)

  double theta=90-46.4551;      // theta
  double phi=360-119.408;       // phi

  for(int i=0;i<N_IFO;i++) {
    for(int j=i+1;j<N_IFO;j++) {
      cout << detParms[i].name << " " << detParms[j].name << " -> "
           << gNET->GetDelay(detParms[i].name,detParms[j].name,phi,theta) << " sec " << endl;
    }
  }
}

