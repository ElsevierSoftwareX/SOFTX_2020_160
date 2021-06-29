//
// This macro is used to test the command cwb_combine_cbc
//
// Author : Gabriele Vedovato
//
// Combined IFAR -> the IFAR of any event reconstructed by the BBH and IMBHB searches is defined in the table below
//
// -----------------------------------------------------------------------------------------------------------------
// Case 1 -> IMBHB band (IMBHB BBH), BBH band (         ) -> IFAR = IFAR(IMBHB)                 TRIALS=1
// Case 2 -> IMBHB band (         ), BBH band (IMBHB BBH) -> IFAR = IFAR(BBH  )                 TRIALS=1
// Case 3 -> IMBHB band (IMBHB    ), BBH band (      BBH) -> IFAR = max IFAR(IMBHB,BBH)/2       TRIALS=2
// Case 4 -> IMBHB band (      BBH), BBH band (IMBHB    ) -> IFAR = max IFAR(IMBHB,BBH)/2       TRIALS=2
// Case 5 -> IMBHB band (IMBHB    ), BBH band (         ) -> IFAR = max IFAR(IMBHB,BBH)         TRIALS=1
// Case 6 -> IMBHB band (      BBH), BBH band (         ) -> IFAR = max IFAR(IMBHB,BBH)/2       TRIALS=2 
// Case 7 -> IMBHB band (         ), BBH band (IMBHB    ) -> IFAR = max IFAR(IMBHB,BBH)/2       TRIALS=2
// Case 8 -> IMBHB band (         ), BBH band (      BBH) -> IFAR = max IFAR(IMBHB,BBH)         TRIALS=1
// -----------------------------------------------------------------------------------------------------------------
//
// Ex: root -b -l 'MakeWaveTestCombineCBC.C("BBH","test_combine",8)'
//

// must be used a wave output root file, after the pp-cuts, the unique selection and the ifar set
#define IWAVE_FILE	"merge/wave_O3_K02_C00_LH_BBH_SIM_NR-MIX_tst1.M1.C_U.S_bin1_cut_nodq_sk8.root"

#define MAX_TREE_SIZE   100000000000LL

#define nIFO		2

#define YEAR		(24.*3600.*365.)

void MakeWaveTestCombineCBC(TString osearch, TString otag, int Case) {

  if(Case<1 || Case>8) {
    cout << "Error - The allowed cases are 1...8" << endl << endl;
    gSystem->Exit(1);
  }

  // check if cwb_combine_osearch is defined
  if(osearch!="IMBHB" && osearch!="BBH") {
    cout << "Error - osearch not defined, valid values are: IMBHB,BBH" << endl << endl;
    gSystem->Exit(1);
  }

  // check it ifwave exist
  TString ifwave = IWAVE_FILE;
  CWB::Toolbox::checkFile(ifwave);

  // check if output files exist
  TString ofwave = IWAVE_FILE;
  ofwave.ReplaceAll(".root","_"+osearch+"_"+otag+".root");
  bool overwrite = CWB::Toolbox::checkFile(ofwave,true);
  if(!overwrite) gSystem->Exit(1);

  TFile *iwwave = TFile::Open(ifwave);
  TTree* iwtree = (TTree*)iwwave->Get("waveburst");

  int isize = iwtree->GetEntries();
  cout << "isize " << isize << endl;

  bool check_ifar=false;
  TBranch* branch;
  TIter next(iwtree->GetListOfBranches());
  while ((branch=(TBranch*)next())) {
    if(TString(branch->GetName())=="ifar") {
      check_ifar=true;
      TLeaf* leaf = branch->FindLeaf("ifar");
      //cout << branch->GetName() << "\t" << leaf->GetTypeName() << endl;
    }
  }
  next.Reset();
  if(check_ifar==false) {
    cout << "Error: ifar not present in : " << IWAVE_FILE << endl;
    exit(1);
  }

  float   iIFAR;
  double* iTIME   = new double[2*nIFO];
  float*  iFREQ   = new float[nIFO];
  iwtree->SetBranchAddress("ifar",&iIFAR);
  iwtree->SetBranchAddress("time",iTIME);
  iwtree->SetBranchAddress("frequency",iFREQ);

  // create output wave root file
  TFile* owroot = new TFile(ofwave,"RECREATE");
  TTree* owtree = (TTree*)iwtree->CloneTree(0);
  owtree->SetMaxTreeSize(MAX_TREE_SIZE);

  float   oIFAR;
  double* oTIME   = new double[2*nIFO];
  float*  oFREQ   = new float[nIFO];
  owtree->SetBranchAddress("ifar",&oIFAR);
  owtree->SetBranchAddress("time",oTIME);
  owtree->SetBranchAddress("frequency",oFREQ);

  for(int i=0;i<isize;i++) {
    iwtree->GetEntry(i);

    oIFAR = iIFAR;
    for(int n=0;n<2*nIFO;n++) oTIME[n] = iTIME[n];
    for(int n=0;n<nIFO;n++)   oFREQ[n] = iFREQ[n];

// Case 1 -> IMBHB band (IMBHB BBH), BBH band (         ) -> IFAR = IFAR(IMBHB)                 TRIALS=1
    if(Case==1) {
        if(osearch=="IMBHB") {
        oIFAR = 100.;
        for(int n=0;n<nIFO;n++) oFREQ[n] = 70;
      }
      if(osearch=="BBH") {
        oIFAR = 200.;
        for(int n=0;n<nIFO;n++) oFREQ[n] = 70;
      }
    }

// Case 2 -> IMBHB band (         ), BBH band (IMBHB BBH) -> IFAR = IFAR(BBH  )                 TRIALS=1
    if(Case==2) {
      if(osearch=="IMBHB") {
        oIFAR = 100.;
        for(int n=0;n<nIFO;n++) oFREQ[n] = 90;
      }
      if(osearch=="BBH") {
        oIFAR = 200.;
        for(int n=0;n<nIFO;n++) oFREQ[n] = 90;
      }
    }

// Case 3 -> IMBHB band (IMBHB    ), BBH band (      BBH) -> IFAR = max IFAR(IMBHB,BBH)/2       TRIALS=2
    if(Case==3) {
      if(osearch=="IMBHB") {
        oIFAR = 100.;
        for(int n=0;n<nIFO;n++) oFREQ[n] = 70;
      }
      if(osearch=="BBH") {
        oIFAR = 200.;
        for(int n=0;n<nIFO;n++) oFREQ[n] = 90;
      }
    }

// Case 4 -> IMBHB band (      BBH), BBH band (IMBHB    ) -> IFAR = max IFAR(IMBHB,BBH)/2       TRIALS=2
    if(Case==4) {
      if(osearch=="IMBHB") {
        oIFAR = 100.;
        for(int n=0;n<nIFO;n++) oFREQ[n] = 90;
      }
      if(osearch=="BBH") {
        oIFAR = 200.;
        for(int n=0;n<nIFO;n++) oFREQ[n] = 70;
      }
    }

// Case 5 -> IMBHB band (IMBHB    ), BBH band (         ) -> IFAR = max IFAR(IMBHB,BBH)         TRIALS=1
    if(Case==5) {
      if(osearch=="IMBHB") {
        oIFAR = 100.;
        for(int n=0;n<nIFO;n++) oFREQ[n] = 70;
      }
      if(osearch=="BBH") {	// not detected
        oIFAR = 0.;
        for(int n=0;n<nIFO;n++) oFREQ[n] = 70;
      }
    }

// Case 6 -> IMBHB band (      BBH), BBH band (         ) -> IFAR = max IFAR(IMBHB,BBH)/2       TRIALS=2 
    if(Case==6) {
      if(osearch=="IMBHB") {	// not detected
        oIFAR = 0.;
        for(int n=0;n<nIFO;n++) oFREQ[n] = 70;
      }
      if(osearch=="BBH") {
        oIFAR = 200.;
        for(int n=0;n<nIFO;n++) oFREQ[n] = 70;
      }
    }

// Case 7 -> IMBHB band (         ), BBH band (IMBHB    ) -> IFAR = max IFAR(IMBHB,BBH)/2       TRIALS=2
    if(Case==7) {
      if(osearch=="IMBHB") {
        oIFAR = 100.;
        for(int n=0;n<nIFO;n++) oFREQ[n] = 90;
      }
      if(osearch=="BBH") {	// not detected
        oIFAR = 0.;
        for(int n=0;n<nIFO;n++) oFREQ[n] = 90;
      }
    }

// Case 8 -> IMBHB band (         ), BBH band (      BBH) -> IFAR = max IFAR(IMBHB,BBH)         TRIALS=1
    if(Case==8) {
      if(osearch=="IMBHB") {	// not detected
        oIFAR = 0.;
        for(int n=0;n<nIFO;n++) oFREQ[n] = 90;
      }
      if(osearch=="BBH") {
        oIFAR = 200.;
        for(int n=0;n<nIFO;n++) oFREQ[n] = 90;
      }
    }

    owtree->Fill(); 
  }

  owroot->cd();
  owtree->Write();
  owroot->Close();

  // create output mdc file -> symbolic link 
  TString ofmdc;
  TString odir = gSystem->DirName(ofwave);
  // create symbolic link to live root file
  TString ifmdc = ifwave;
  ifmdc.ReplaceAll("wave_","mdc_");
  ifmdc.Remove(0,ifmdc.Last('/')+1);    // strip path
  ofmdc = ofwave;
  ofmdc.ReplaceAll("wave_","mdc_");
  ofmdc.Remove(0,ofmdc.Last('/')+1);    // strip path
  cout << odir << endl;
  //cout << ofmdc << endl;
  Long_t id,size,flags,mt;
  int estat = gSystem->GetPathInfo(odir+"/"+ifmdc,&id,&size,&flags,&mt);
  if(estat==0) {
    char cmd[1024];
    sprintf(cmd,"cd %s;ln -sf %s %s",odir.Data(),ifmdc.Data(),ofmdc.Data());
    //cout << cmd << endl;
    gSystem->Exec(cmd);
  }

  // create output merge list combine file name
  TString ofmerge = ofwave;
  ofmerge.ReplaceAll("wave_","merge_");
  ofmerge.ReplaceAll(".root",".lst");
  char cmd[1024]; sprintf(cmd,"touch %s",ofmerge.Data());
  gSystem->Exec(cmd);

  // write history intooutput  wave file
  CWB::History* history = (CWB::History*)iwwave->Get("history");
  if(history!=NULL) {
    TFile owfile(ofwave,"UPDATE");
    history->Write("history");
    owfile.Close();
  }
  delete history;

  // output files 
  cout << endl;
  cout << "Output Combined files : " << endl;
  cout << ofwave << endl;
  cout << odir+"/"+ofmdc << endl;
  cout << ofmerge << endl;
  cout << endl;

  exit(0); 
}
