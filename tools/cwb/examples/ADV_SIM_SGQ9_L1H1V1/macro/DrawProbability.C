//
// Draw Probability SkyMap
// Author : Gabriele Vedovato


#define RESOLUTION  2 
//#define RESOLUTION  4

//#define COORDINATES "cWB"
#define COORDINATES "Geographic"

#define PROJECTION ""
//#define PROJECTION "hammer"
//#define PROJECTION "sinusoidal"

//#define WRITE_PLOT

#define SKYMAP_DIR  "data/ced_931158100_600_ADV_SIM_SGQ9_L1H1V1_30_job1"
#define SKYMAP_FILE "L1H1V1_931158128.695_931158128.695_931158128.695/probability.root"

using namespace CWB;

void DrawProbability() {

  int nIFO=3;
  TString ifo[3]={"L1","H1","V1"};

  char ifostr[32]="";
  for(int n=0; n<nIFO; n++) {
    sprintf(ifostr,"%s %s",ifostr,ifo[n].Data());
  }

  SkyPlot* skyplot = new SkyPlot("skymap",PROJECTION,COORDINATES,RESOLUTION);
  TH2D* h2 = (TH2D*)skyplot->GetHistogram();
//  h2->GetZaxis()->SetRangeUser(0.78,1.0);
  skyplot->SetTitle(ifostr);

  TString skymap_path = TString(SKYMAP_DIR)+"/"+TString(SKYMAP_FILE);
  skymap sm(skymap_path);
  int L = sm.size();
  double F=1e16*double(L);
//cout << F << " " << log10(F) << endl;
  double norm=0.;
  for(int l=0;l<L;l++) {
    double val = log10(F)*sm.get(l);
    val = pow(10,val)/F;
    sm.set(l,val);
    norm+=val;
  }
  sm *= 1./norm;
  skyplot->FillData(sm);
  skyplot->Draw(0);

#ifdef WRITE_PLOT
  cout << "Write : " << ofileName << endl;
  skyplot->Print(ofileName);
#endif
}

