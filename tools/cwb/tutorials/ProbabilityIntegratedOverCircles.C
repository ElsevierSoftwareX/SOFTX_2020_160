// Integrate likelihood skymap over the big circles
// Author : Gabriele Vedovato
//

{
  //#define SKYMAP_FILE "data/ced_968654036_60_sg554q8d9_obj_20_job1/L1H1V1_968654042.750_968654042.750_968654042.750/probability.root"
  #define SKYMAP_FILE "/home/vedovato/waveburst/coherence/offline/SIMULATION/ADV_SIM_SGQ9_L1H1V1_run5/data/ced_931158378_44_ADV_SIM_SGQ9_L1H1V1_run5_30_job1/L1H1V1_931158395.102_931158395.102_931158395.102/probability.root"

  #define SpeedOfLightInVacuo 2.99792458000000000e+08

  #define ALPHA_FACTOR 4
  #define BETA_FACTOR 4

  #define H1L1_AXIS

  using namespace ROOT::Math;

  gSystem->Load("libMathCore");

  gskymap* iSM = new gskymap(TString(SKYMAP_FILE));


  int nIFO=3;
  TString ifo[3]={"L1","H1","V1"};
  gnetwork* gNET = new gnetwork(nIFO,ifo);
  gskymap* gSM = gNET->GetGskymap();
  gSM->SetOptions("","Geographic");
  *gSM = *iSM;

  skymap* osm[3];
  for(int n=0;n<nIFO;n++) {
    osm[n] = new skymap(*gSM);
    int L = osm[n]->size();
    for(int l=0;l<L;l++) osm[n]->set(l,0);
  }

  detector *pD[3];
  for(int n=0; n<nIFO; n++) pD[n]=gNET->getifo(n);

  XYZVector L1(pD[0]->Rv[0],pD[0]->Rv[1],pD[0]->Rv[2]);
  XYZVector H1(pD[1]->Rv[0],pD[1]->Rv[1],pD[1]->Rv[2]);
  XYZVector V1(pD[2]->Rv[0],pD[2]->Rv[1],pD[2]->Rv[2]);

  L1=L1/SpeedOfLightInVacuo;
  H1=H1/SpeedOfLightInVacuo;
  V1=V1/SpeedOfLightInVacuo;

  XYZVector D1 = H1-L1;
  XYZVector D3 = D1.Cross(V1-L1);
  XYZVector D2 = D1.Cross(D3);

  XYZVector eD1 = D1.Unit();
  XYZVector eD2 = D2.Unit();
  XYZVector eD3 = D3.Unit();

  Rotation3D R(eD1,eD2,eD3);

  XYZVector H1L1 = H1-L1;
  H1L1.SetRho(1);
  XYZVector V1L1 = V1-L1;
  V1L1.SetRho(1);
  XYZVector V1H1 = V1-H1;
  V1H1.SetRho(1);

  TVector3 vH1L1(H1L1.X(),H1L1.Y(),H1L1.Z());
  TVector3 vV1L1(V1L1.X(),V1L1.Y(),V1L1.Z());
  TVector3 vV1H1(V1H1.X(),V1H1.Y(),V1H1.Z());
  TVector3 veD3(eD3.X(),eD3.Y(),eD3.Z());

  double r2d = 180./TMath::Pi();

  double _PI=TMath::Pi(); 

  int N=360*BETA_FACTOR;
  wavearray<double> th(N);
  wavearray<double> ph(N);
  double dalpha = _PI/180./(double)ALPHA_FACTOR;
  for(int n=0;n<nIFO;n++) {
    for(int k=0;k<ALPHA_FACTOR*180;k++) {

      double alpha = k*dalpha;
      TRotation oR;
      oR.Rotate(alpha,veD3.Unit());
      TVector3 roSD;
      if(n==0) roSD = oR*vH1L1;
      if(n==1) roSD = oR*vV1L1;
      if(n==2) roSD = oR*vV1H1;

      double like=0;
      for (int i=0;i<N;i++) {

        double angle = TMath::Pi()*i/180./BETA_FACTOR;
        TRotation vR;
        if(n==0) vR.Rotate(angle,vH1L1.Unit());
        if(n==1) vR.Rotate(angle,vV1L1.Unit());
        if(n==2) vR.Rotate(angle,vV1H1.Unit());

        TVector3 rvSD = vR*roSD;
        double theta = rvSD.Theta()*r2d;
        double phi = rvSD.Phi()*r2d;
        if (phi<0) phi+=360; if (phi>360) phi-=360;
        if (theta<0) theta+=180; if (theta>180) theta-=180;

        th[i]=theta; 
        ph[i]=phi; 
      }

      for (int i=0;i<N;i++) {
        int l=gSM->getSkyIndex(th[i],ph[i]);
        like+=gSM->get(l);
      }
      like/=N;
//like*=alpha;
//if(alpha>_PI/2) like*=(_PI-alpha);
      for (int i=0;i<N;i++) {
        int l=osm[n]->getSkyIndex(th[i],ph[i]);
        osm[n]->set(l,like);
      }
    }
  }

  int L = gSM->size();
  for(int l=0;l<L;l++) {
    double like=0;
    for (int n=0;n<nIFO;n++) like+=osm[n]->get(l);
    //like=osm[0]->get(l);
    gSM->set(l,like);
  }

  gSM->Draw(0);
}
