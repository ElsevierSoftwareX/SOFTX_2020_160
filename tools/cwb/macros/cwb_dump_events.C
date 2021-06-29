/*
# Copyright (C) 2019 Gabriele Vedovato
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/


//  dump/view the event parameters : used by the cwb_dump events

{
  CWB::Toolbox TB;

  TB.checkFile(gSystem->Getenv("CWB_ROOTLOGON_FILE"));

  TString cwb_dump_evt_file_name;
  if(gSystem->Getenv("CWB_DUMP_EVT_FILE_NAME")==NULL) {
    cout << "Error : environment CWB_DUMP_EVT_FILE_NAME is not defined!!!" << endl;exit(1);
  } else {
    cwb_dump_evt_file_name=TString(gSystem->Getenv("CWB_DUMP_EVT_FILE_NAME"));
  }
  if(cwb_dump_evt_file_name.Contains(".root")==0) {
    cout << "Error : " << cwb_dump_evt_file_name.Data() << " is not a root file!!!" << endl;exit(1);
  }

  TString cwb_dump_evt_mode="view";
  if(gSystem->Getenv("CWB_DUMP_EVT_MODE")!=NULL) {
    cwb_dump_evt_mode=TString(gSystem->Getenv("CWB_DUMP_EVT_MODE"));
  }

  TFile *ifile = TFile::Open(cwb_dump_evt_file_name);
  if(ifile==NULL) {cout << "Failed to open " << cwb_dump_evt_file_name.Data() << endl;exit(-1);}

  // get waveburst tree
  TTree* itree = (TTree *) ifile->Get("waveburst");
  if(itree==NULL) {cout << "Error : tree waveburst not present !!!" << endl;exit(1);}

  // get detector list
  TList* list = itree->GetUserInfo();
  int nifo=list->GetSize();
  if (nifo==0) {cout << "Error : no ifo present in the tree" << endl;exit(1);}
  detector** pD = new detector*[nifo];
  for (int n=0;n<nifo;n++) {
    pD[n] = (detector*)list->At(n);
    detectorParams dParams = pD[n]->getDetectorParams();
    //pD[n]->print();
  }

  netevent EVT(itree,nifo);
  for (int n=0;n<nifo;n++) EVT.ifoList.push_back(pD[n]);

  // build dump file name
  char evtFile[512];
  TObjArray* token = TString(cwb_dump_evt_file_name).Tokenize(TString("/"));
  sprintf(evtFile,"%s/%s",dump_dir,
          TString(((TObjString*)token->At(token->GetEntries()-1))->GetString()).ReplaceAll(".root",".evt").Data());

  // dump header
  if(cwb_dump_evt_mode.CompareTo("view")!=0) {
    cout << "Write : " << evtFile << endl;
    EVT.dopen(evtFile,const_cast<char*>("w"),true); 
    EVT.dclose();
  }
  // dump events
  int nEVT = EVT.GetEntries();
  for(int i=0;i<nEVT;i++) {

    EVT.GetEntry(i);

    if(cwb_dump_evt_mode.CompareTo("view")==0) EVT.fP = stdout;
    else EVT.dopen(evtFile,const_cast<char*>("a"),false);
    EVT.Dump();
    if(cwb_dump_evt_mode.CompareTo("view")!=0) EVT.dclose();
  }

  exit(0);
}
