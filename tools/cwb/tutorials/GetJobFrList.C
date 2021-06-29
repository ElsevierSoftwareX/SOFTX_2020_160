
void GetJobFrList(TString fName) {
  //
  // read the frame list from the job file
  // Author : Gabriele Vedovato
  //
  // Example : root 'GetJobFrList("jobFileName.root")'

  // read cwb object from job file
  TFile* ifile = new TFile(fName);
  if(ifile==NULL||!ifile->IsOpen())
    {cout << "GetJobFrList - Error opening root file : " << fName.Data() << endl;gSystem->Exit(1);}
  ifile->cd();
  cwb* iCWB = (cwb*)ifile->Get("cwb");
  if(iCWB==NULL) {
    cout << "GetJobFrList - Error : cwb is not contained in root file " << fName.Data() << endl;
    gSystem->Exit(1);
  }

  //vector<frfile> frlist = iCWB->GetFrList("L1");	// get L1 frlist
  //vector<frfile> frlist = iCWB->GetFrList(1);		// get second detector frlist
  vector<frfile> frlist = iCWB->GetFrList();		// get full frame list

  for(int i=0;i<frlist.size();i++) {
     cout << "start : " << frlist[i].start << "(sec) stop : " << frlist[i].stop 
          << "(sec) length : " << frlist[i].length << "(sec)" << endl;
     vector<TString> flist = frlist[i].file;
     for(int j=0;j<flist.size();j++) cout << j << " " << flist[j] << endl;
  }

  gSystem->Exit(0);
}
