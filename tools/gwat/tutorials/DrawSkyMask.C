//
// Draw SkyMask
// Author : Gabriele Vedovato

#define IFILE_NAME "SkyMaskCC_BRST50MPC_S6_Rev1d7_R0d4_NSOURCES.txt" 
//#define IFILE_NAME "SkyMaskCC_BRST50MPC_S6_Rev1d7_R0d4_PROB.txt" 
//#define IFILE_NAME "SkyMaskCC_BRST50MPC_S6_Rev1d7_R0d4.txt" 
//#define IFILE_NAME "SkyMaskCC_BRST50MPC_S6_Rev1d7_R0d4_PHPAD.txt" 

void DrawSkyMask() {

  gskymap* gSM = new gskymap(0.4,0,180,0,360);
  //gSM->SetOptions("sinusoidal","Geographic");
  gSM->SetOptions("","Geographic");
  //gSM->SetOptions("","cWB");
  gSM->SetTitle(IFILE_NAME);
  int L = gSM->size();

  ifstream in;
  in.open(IFILE_NAME,ios::in);
  if (!in.good()) {cout << "Error Opening File : " << IFILE_NAME << endl;exit(1);}

  int nPixels=0;
  int l;
  float mask;
  while (1) {
    in >> l >> mask;
    if (!in.good()) break;
    gSM->set(l,mask);
    if(mask) nPixels++;
  }
  in.close();
  cout << "% : " << 100.*(double)nPixels/(double)L << endl;

  gSM->Draw(0);

}

