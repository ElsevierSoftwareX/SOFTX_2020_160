//
// Read Frame & Produce Bicoherence
// Author : Gabriele Vedovato
//
// to run macro do :
// root BicoConfig.C MakeBicoherence.C
// BicoConfig.C defines the setup
//

#include <vector>

void
MakeBicoherence() {

  // set & get frame file list
  CWB::frame fr_1(FRLIST_NAME_1);
  int nfrFiles_1=fr_1.getNfiles();
  cout << "nfrFiles_1 : " << nfrFiles_1 << endl;
  CWB::frame fr_2(FRLIST_NAME_2);
  int nfrFiles_2=fr_2.getNfiles();
  cout << "nfrFiles_2 : " << nfrFiles_2 << endl;

  int bsize = BSIZE;
  double srate = SRATE;

  CWB::Bicoherence* blc = new CWB::Bicoherence(CHNAME_1, CHNAME_2, srate, bsize,
                                               BCM_X_NUM_SLICES, BCM_Y_NUM_SLICES,
                                               BCM_X_SLICE_INDEX, BCM_Y_SLICE_INDEX,
                                               BCM_ORDER);

#ifdef SEGLIST
  int sizeList = blc->readSegList(SEGLIST);
  cout << "sizeList : " << sizeList << endl;
#endif

  double pi = TMath::Pi();
  double dt = 1./srate;
  TRandom normal;

  double start = START;

  wavearray<double> x(bsize);
  wavearray<double> y; y.resize(0);
  x.rate(srate);
//  for (int cnt=0;cnt<NLOOP;cnt++) {
  int cnt_rejected=0;
  int cnt=0;
  while (cnt<NLOOP) {
    double stop = start+bsize/srate;
    if(blc->segListCheck(start,stop)) {
      cout << "READ Data Channel " << CHNAME_1 << endl;
      frfile FRF_1 = fr_1.getFrList(start, stop, 0);
      fr_1.readFrames(FRF_1,(char*)CHNAME_1,x);  
      cout << "READ Data Channel " << CHNAME_2 << endl;
      frfile FRF_2 = fr_2.getFrList(start, stop, 0);
      fr_2.readFrames(FRF_2,(char*)CHNAME_2,y);  
      cout << "X Channel rate : " << x.rate() << endl;
      cout << "Y Channel rate : " << y.rate() << endl;
      cout << "Loop : " << blc->GetAverages() << endl;
      if(blc->MakeBicoherence(x,y)) cnt++; else cout << "Warning - buffer rejected : " << cnt_rejected++ << endl;
    } else cout << "Warning - buffer rejected : " << cnt_rejected++ << endl;
    start=stop; 
  }

  cout << "Warning - buffer rejected : " << cnt_rejected << endl;

/*
  vector<bico> listBico = blc->GetBicoherence(BLC_THRESHOLD,REBIN);
  cout << "list size : " << listBico.size() << endl;

  char ofile[256];
  sprintf(ofile,"%s/bico_%s_%s.txt",ODIR_NAME,CHNAME_1,chname_2.Data());
  ofstream out;
  out.open(ofile,ios::out);

  for(int i=0;i<listBico.size();i++) {
    printf("%d -> %3.2f - %3.2f : %1.2f\n",i,listBico[i].x,listBico[i].y,listBico[i].c);
    char ostr[256];
    sprintf(ostr,"%d -> %3.2f - %3.2f : %1.2f\n",i,listBico[i].x,listBico[i].y,listBico[i].c);
    out << ostr <<endl;
  }
  out.close();
*/

  int fDeltaX  = srate/BCM_X_NUM_SLICES/2; 
  int fOffsetX = fDeltaX*BCM_X_SLICE_INDEX; 
  int fDeltaY  = srate/BCM_Y_NUM_SLICES/2; 
  int fOffsetY = fDeltaY*BCM_Y_SLICE_INDEX; 

#ifdef BATCH
  bool batch=true;
#else
  bool batch=false;
#endif

  int graphId=0;
#ifdef GRAPHID
  graphId=GRAPHID;
#endif

  char ofname[256]="";
#ifdef SAVE
  if(graphId==0 || graphId==2) {
    sprintf(ofname,"%s/bico_O%d_%dS_%dL_%dR_%s_%dHz_%dHz_%s_%dHz_%dHz.png",
            ODIR_NAME,BCM_ORDER,START,NLOOP,REBIN,CHNAME_1,fOffsetX,fDeltaX,CHNAME_2,fOffsetY,fDeltaY);
  } else {
    sprintf(ofname,"%s/bico_O%d_%dS_%dL_%dR_%s_%dHz_%dHz_%s_%dHz_%dHz.png",
            ODIR_NAME,BCM_ORDER,START,NLOOP,REBIN,CHNAME_2,fOffsetY,fDeltaY,CHNAME_1,fOffsetX,fDeltaX);
  } 
  TString tmp(ofname);tmp.ReplaceAll(":","_");
  strcpy(ofname,tmp.Data());
#endif

  blc->DrawBicoherence(REBIN,graphId,ofname,batch);
//  blc->Reset();

}
