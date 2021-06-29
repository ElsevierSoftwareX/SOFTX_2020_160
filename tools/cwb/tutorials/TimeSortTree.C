//
// Sort waveburst root file IFILE_NAME according time[0] parameters and to the ISELECTION rules
// write the sorted root file to OFILE_NAME
// write the sorted selected parameters OEXPRESSION to file LOG_NAME
// Author : Gabriele Vedovato

{

   #define IFILE_NAME 	"merge/wave_S6A_R4_SIM_SGQ9_L1H1V1_1G_run2.M1.root"
   #define OFILE_NAME 	"wave_S6A_R4_SIM_SGQ9_L1H1V1_1G_run2.M1.sorted_factor_0d3.root"
   #define TREE_NAME 	"waveburst"

   #define ISELECTION	"abs(time[0]-time[3])<0.1 && abs(factor-0.3)<0.1"
   #define OEXPRESSION     "run:time[0]:rho[1]:netcc[0]"

   #define LOG_NAME	"log.txt"

   CWB::Toolbox::checkFile(IFILE_NAME);
   TFile f(IFILE_NAME);
   TTree *tree = (TTree*)f.Get(TREE_NAME);
   //Drawing variable t with no graphics option.
   //variable pz stored in array fV1 (see TTree::Draw)
   tree->Draw("time[0]:Entry$",ISELECTION,"goff");
   Int_t nentries = (Int_t)tree->GetSelectedRows();
   Int_t *_index = new Int_t[nentries];
   double* entry = tree->GetV2(); 
   //sort array containing pz in decreasing order
   //The array index contains the entry numbers in decreasing order of t
   TMath::Sort(nentries,tree->GetV1(),_index,false);

   //open new file to store the sorted Tree
   TFile f2(OFILE_NAME,"recreate");
   //Create an empty clone of the original tree
   TTree *tsorted = (TTree*)tree->CloneTree(0);
   for (Int_t i=0;i<nentries;i++) {
      if (i%1000==0 && i>0) cout << "write entry : " << i << endl;
      tree->GetEntry(entry[_index[i]]);
      tsorted->Fill();
   }
   tsorted->Write();
   delete [] _index;

#ifdef LOG_NAME
   tsorted->SetScanField(0);
   ((TTreePlayer*)(tsorted->GetPlayer()))->SetScanRedirect(true); 
   ((TTreePlayer*)(tsorted->GetPlayer()))->SetScanFileName(LOG_NAME);
   tsorted->Scan(OEXPRESSION,"","");
#endif

   gSystem->Exit(0);
}

