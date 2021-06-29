//
// This example show how to display skymap probability from fits file
// Author : Gabriele Vedovato
//
// Ex : root 'DrawFitsGraceDB.C("skyprobcc.fits","_skymap",true,false)'

//#define COORDINATES "cWB"
//#define COORDINATES "Geographic"
#define COORDINATES "Celestial"

//#define PROJECTION ""
#define PROJECTION "hammer"
//#define PROJECTION "sinusoidal"

#define GRACEDB_LABEL	"GXXXXXX"
#define EM_LABEL	"EM_FOLLOWUP"
#define EM_PHI		360.-65.5
#define EM_THETA	-82.06

void DrawFitsGraceDB(TString fitsName, TString label="", bool cumulative=false, bool save=false) {

  gskymap* gSM = new gskymap((char*)fitsName.Data());

  int L     = gSM->size();                // number of pixels in the sphere
  double pi = TMath::Pi();
  double S  = 4*pi*pow(180/pi,2);         // solid angle of a sphere
  double dS = S/L;                        // solid angle of a pixel

  if(cumulative) {

    double cumul=0;
    int L=gSM->size();
    int* index = new int[L];              // sorted index
    double* prob = new double[L];         // probability array
    for(int l=0;l<L;l++) {                // loop over the sky grid
      prob[l] = gSM->get(l);              // get skymap prob
    }

    TMath::Sort(L,prob,index,false);      // sort prob
    cumul=0;
    for(int l=0;l<L;l++) {                // loop over the sky grid
      int m=index[l];                     // sorted index
      double dec = gSM->getTheta(m);      // get dec  
      double ra  = gSM->getPhi(m);        // get ra  
      cumul+=prob[m];
      gSM->set(m,cumul);                  // set !=0 to force dark blue background
/*
      if(prob[index[l]]>0) {
        cout << m << "\tdec : " << dec << "\tra : " << ra << "\tprob : "
             << prob[m] << "\tcumul prob : " << cumul << endl;
      }
*/
    }
    if(label!="") {
      TString tag=label;
      tag.ReplaceAll("_"," "); 
      gSM->SetTitle(TString::Format("%s : SkyMap Cumulative Probability Distribution (%s)",GRACEDB_LABEL,tag.Data()));
    } else {
      gSM->SetTitle(TString::Format("%s : SkyMap Cumulative Probability Distribution",GRACEDB_LABEL));
    }
  } else {
    if(label!="") {
      TString tag=label;
      tag.ReplaceAll("_"," "); 
      gSM->SetTitle(TString::Format("%s : SkyMap Probability Distribution (%s)",GRACEDB_LABEL,tag.Data()));
    } else {
      gSM->SetTitle(TString::Format("%s : SkyMap Probability Distribution",GRACEDB_LABEL));
    }
    gSM->SetZaxisTitle("prob. per deg^{2}");
  }

  for(int l=0;l<L;l++) {                  // loop over the sky grid
    double P = gSM->get(l);               // get probability per pixel
    P/=dS;                                // normalize probability to prob per deg^2
    if(P==0) gSM->set(l,1e-40);           // set !=0 to force dark blue background
  }


  gSM->SetOptions(PROJECTION,COORDINATES,2);
  gSM->Draw(57);

  //gSM->DrawMarker(EM_PHI,EM_THETA,29,2.5,kWhite);
  //gSM->DrawMarker(EM_PHI,EM_THETA,30,2.5,kBlack);

  //gSM->DrawMarker(EM_PHI,EM_THETA,30,2.0,kWhite);
  //gSM->DrawText(EM_PHI-10,EM_THETA+10, EM_LABEL, 0.04, kWhite);

  //gSM->DrawMarker(EM_PHI,EM_THETA,30,2.0,kBlack);
  //gSM->DrawText(EM_PHI-10,EM_THETA+10, EM_LABEL, 0.04, kWhite);

  if(save) {
    if(cumulative) {
      gSM->Print(TString::Format("%s_cumulative_probability_skymap%s.png",GRACEDB_LABEL,label.Data()));
    } else {
      gSM->Print(TString::Format("%s_probability_skymap%s.png",GRACEDB_LABEL,label.Data()));
    }
    exit(0);
  }
}
