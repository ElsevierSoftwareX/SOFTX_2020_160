// this example show how to extract the first 2000 jobs from the merge root file

{
  #define DATA_LABEL "S6D_R11_SIM_BRST_L1H1V1_2G_MP_run5.M1"


  char ifile_wave_name[256];
  sprintf(ifile_wave_name,"merge/wave_%s.root",DATA_LABEL);
  char ofile_wave_name[256];
  sprintf(ofile_wave_name,"merge/wave_%s.run2000.root",DATA_LABEL);

  TFile *ifile_wave = TFile::Open(ifile_wave_name);
  if(ifile_wave==NULL) {
    cout << "Error : File " << ifile_wave_name << " not exist !!!" << endl;
    gSystem->Exit(1);
  }
  TTree* itree_wave = (TTree *) gROOT->FindObject("waveburst");

  TFile ofile_wave(ofile_wave_name,"RECREATE");
  TTree* otree_wave = itree_wave->CopyTree("run<=2000");
  otree_wave->SetDirectory(&ofile_wave);
  otree_wave->Write();
  ofile_wave.Close();

  char ifile_mdc_name[256];
  sprintf(ifile_mdc_name,"merge/mdc_%s.root",DATA_LABEL);
  char ofile_mdc_name[256];
  sprintf(ofile_mdc_name,"merge/mdc_%s.run2000.root",DATA_LABEL);

  TFile *ifile_mdc = TFile::Open(ifile_mdc_name);
  TTree* itree_mdc = (TTree *) gROOT->FindObject("mdc");

  TFile ofile_mdc(ofile_mdc_name,"RECREATE");
  TTree* otree_mdc = itree_mdc->CopyTree("run<=2000");
  otree_mdc->SetDirectory(&ofile_mdc);
  otree_mdc->Write();
  ofile_mdc.Close();


  exit(0);
}
