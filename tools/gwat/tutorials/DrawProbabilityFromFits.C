//
// This example show how to display skymapcc probability from fits file
// Author : Gabriele Vedovato


//#define COORDINATES "cWB"
//#define COORDINATES "Geographic"
#define COORDINATES "Celestial"

#define PROJECTION ""
//#define PROJECTION "hammer"
//#define PROJECTION "sinusoidal"

#define RESOLUTION 2

#define DRAW_SKYMAP
//#define DUMP_SKYMAP
//#define SAVE_SKYMAP

void DrawProbabilityFromFits(TString fitsName) {

  gskymap* gSM = new gskymap(fitsName.Data());
  gSM->SetOptions(PROJECTION,COORDINATES,RESOLUTION);

  int order=gSM->getOrder();		// Get HEALPix order
  int L=gSM->size();			// number of pixels in the sphere

  // get pixel area
  double pi = TMath::Pi();
  double S = 4*pi*pow(180/pi,2);	// solid angle of a sphere 
  double dS = S/L;			// solid angle of a pixel	
  cout << "solid angle of a pixel : " << dS << endl;  

  // normalize probability
  for(int l=0;l<L;l++) {		// loop over the sky grid
    double prob = gSM->get(l);		// get probability per pixel
    prob/=dS;				// normalize probability to prob per deg^2
    if(prob==0) gSM->set(l,1e-40);      // set !=0 to force dark blue background
  }

#ifdef DRAW_SKYMAP
  TGaxis::SetMaxDigits(3);
  gSM->SetZaxisTitle("prob. per deg^{2}");

  gSM->Draw(1);

  TH2D* h2 = gSM->GetHistogram();
  h2->GetZaxis()->SetTitleOffset(0.85);
  h2->GetZaxis()->SetTitleSize(0.03);
#endif

#ifdef DUMP_SKYMAP
  double cumulative=0;
  int L=gSM->size();
  int* index = new int[L];		// sorted index
  double* prob = new double[L];		// probability array
  for(int l=0;l<L;l++) {		// loop over the sky grid
    prob[l] = gSM->get(l);		// get skymap prob
  }

  TMath::Sort(L,prob,index,true);	// sort prob
  cumulative=0;
  for(int l=0;l<L;l++) {		// loop over the sky grid
    int m=index[l];			// sorted index
    double dec = gSM->getTheta(m);	// get dec  
    double ra  = gSM->getPhi(m);	// get ra  
    cumulative+=prob[m];
    if(prob[index[l]]>0) {
      cout << m << "\tdec : " << dec << "\tra : " << ra << "\tprob : " 
           << value[m] << "\tcumulative prob : " << cumulative << endl; 
    }
  }
#endif

#ifdef SAVE_SKYMAP
  gSM->Print("skyprobcc.png");
#endif

}
