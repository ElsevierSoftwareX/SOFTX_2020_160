//
// Read event probability skymap from the online cWB probability and draw
// Author : Gabriele Vedovato
//
// root -l 'DrawOnlineProbability.C("error_region_1058523803.4934.txt",true)'
//
// online probability file format 
//
// sid skyID  theta   DEC     step   phi     R.A    step  probability    cumulative
// 

#define PROJECTION ""
//#define PROJECTION "hammer"
#define RESOLUTION  2
#define COORDINATES "Geographic"

#define HEALPix	7

void DrawOnlineProbability(TString fName, bool save=false) {

  gskymap* gSM = new gskymap(int(HEALPix));
  gSM->SetOptions(PROJECTION,COORDINATES,RESOLUTION/2);

  int    sid, skyid;
  double theta, DEC, step_DEC, phi, RA, step_RA, probability, cumulative;

  // Read probability file
  ifstream in;
  in.open(fName.Data(),ios::in);
  if (!in.good()) {cout << "Error Opening File : " << fName.Data() << endl;exit(1);}

  char str[1024];
  int fpos=0;

  while(true) {
    fpos=in.tellg();
    in.getline(str,1024);
    if(str[0] == '#') continue;
    in.seekg(fpos, ios::beg);
    in >> sid >> skyid >> theta >> DEC >> step_DEC >> phi >> RA >> step_RA >> probability >> cumulative;
    if(!in.good()) break;
    gSM->set(skyid, probability);	// fill skymap
    fpos=in.tellg();
    in.seekg(fpos+1, ios::beg);
  }
  in.close();

  TString title = TString("Probability SkyMap");
  gSM->SetTitle(title);
  gSM->Draw(0);

  if(save) {				// save skymap
    TString ofName=fName;
    ofName.ReplaceAll(".txt",".png"); 
    cout << "Write : " << ofName.Data() << endl;
    gSM->Print(ofName.Data());
    exit(0);
  }
}
