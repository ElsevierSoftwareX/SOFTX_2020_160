//
// Author : Gabriele Vedovato
// this macro is used to redefine the rho statistic

void ChangeRhoStatistic(TString iwfile, bool simulation=false) {
// -------------------------------------------------------------------------------------------------------------

   #define TREE_NAME 	"waveburst"

   char cmd[1024];

   cout << endl; 
   TString owfile = iwfile;
   owfile.ReplaceAll(".root",".C_newrho.root");
   cout << "owave : " << owfile << endl;
   cout << "iwave : " << iwfile << endl;

   cout << endl; 
   TString omfile = owfile;
   omfile.ReplaceAll(".root",".lst");
   omfile.ReplaceAll("wave_","merge_");
   cout << "omerge : " << omfile << endl;

   TString imfile = omfile;
   imfile.ReplaceAll(".C_newrho.lst",".lst");
   cout << "imerge : " << imfile << endl;

   sprintf(cmd,"cd merge;ln -sf ../%s ../%s",imfile.Data(),omfile.Data());
   cout << cmd << endl;
   gSystem->Exec(cmd);

   cout << endl; 
   TString olfile = owfile;
   olfile.ReplaceAll("wave_","live_");
   cout << "olive : " << olfile << endl;

   TString ilfile = olfile;
   ilfile.ReplaceAll(".C_newrho.root",".root");
   cout << "ilive : " << ilfile << endl;

   sprintf(cmd,"cd merge;ln -sf ../%s ../%s",ilfile.Data(),olfile.Data());
   cout << cmd << endl;
   if(!simulation) gSystem->Exec(cmd);

   cout << endl; 
   TString oMfile = owfile;
   oMfile.ReplaceAll("wave_","mdc_");
   cout << "omdc : " << oMfile << endl;

   TString iMfile = oMfile;
   iMfile.ReplaceAll(".C_newrho.root",".root");
   cout << "imdc : " << iMfile << endl;
   cout << endl; 

   sprintf(cmd,"ln -sf %s %s",iMfile.Data(),oMfile.Data());
   sprintf(cmd,"cd merge;ln -sf ../%s ../%s",iMfile.Data(),oMfile.Data());
   cout << cmd << endl;
   if(simulation) gSystem->Exec(cmd);

   CWB::Toolbox::checkFile(iwfile);
   TFile f1(iwfile);
   TTree *itree = (TTree*)f1.Get(TREE_NAME);
   itree->Draw("rho[1]:penalty","","goff");
   Int_t nentries = (Int_t)itree->GetSelectedRows();
 
   float* irho = new float[2];
   itree->SetBranchAddress("rho", irho);
   float* inetcc = new float[4];
   itree->SetBranchAddress("netcc", inetcc);
   float ipenalty;
   itree->SetBranchAddress("penalty", &ipenalty);

   TFile f(owfile,"recreate");
   TTree *otree = (TTree*)itree->CloneTree(0);
 
   float* orho = new float[2];
   otree->SetBranchAddress("rho", orho);

   for (Int_t i=0;i<nentries;i++) {
      if (i%100000==0 && i>0) cout << "write entry : " << i << "/" << nentries << endl;
      itree->GetEntry(i);
      //cout << irho[0] << " " << irho[1] << " " << ipenalty << endl;
      //orho[0] = irho[0];
      orho[0] = irho[0]*sqrt(ipenalty);
      orho[1] = irho[1]*sqrt(ipenalty);
      otree->Fill();
   } 
   otree->Write();

   gSystem->Exit(0);
}

