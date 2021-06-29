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


//  dump/view the history infos : used by the cwb_dump command

{
  CWB::Toolbox TB;

  TB.checkFile(gSystem->Getenv("CWB_ROOTLOGON_FILE"));

  TString cwb_dump_hist_file_name;
  if(gSystem->Getenv("CWB_DUMP_HIST_FILE_NAME")==NULL) {
    cout << "Error : environment CWB_DUMP_HIST_FILE_NAME is not defined!!!" << endl;exit(1);
  } else {
    cwb_dump_hist_file_name=TString(gSystem->Getenv("CWB_DUMP_HIST_FILE_NAME"));
  }
  if(cwb_dump_hist_file_name.Contains(".root")==0) {
    cout << "Error : " << cwb_dump_hist_file_name.Data() << " is not a root file!!!" << endl;exit(1);
  }

  TString cwb_dump_hist_mode="view";
  if(gSystem->Getenv("CWB_DUMP_HIST_MODE")!=NULL) {
    cwb_dump_hist_mode=TString(gSystem->Getenv("CWB_DUMP_HIST_MODE"));
  }

  TFile *ifile = TFile::Open(cwb_dump_hist_file_name);
  if(ifile==NULL) {cout << "Failed to open " << cwb_dump_hist_file_name.Data() << endl;exit(-1);}

  CWB::History* ihistory = (CWB::History*)ifile->Get("history");
  if(ihistory==NULL) {
    cout << "Error : history is not present!!!" << endl;exit(1);
  }

  if(cwb_dump_hist_mode.CompareTo("view")==0) {
    ihistory->Print();
  } else {
    char historyFile[512];
    TObjArray* token = TString(cwb_dump_hist_file_name).Tokenize(TString("/"));
    sprintf(historyFile,"%s/%s",dump_dir,
            TString(((TObjString*)token->At(token->GetEntries()-1))->GetString()).ReplaceAll(".root",".history").Data());

    ihistory->DumpToTextFile(historyFile);
    cout << "Write : " << historyFile << endl;
  }

  exit(0);
}
