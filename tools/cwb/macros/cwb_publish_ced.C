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


// obsolete

{
  CWB::Toolbox TB;
  TB.checkFile(gSystem->Getenv("CWB_ROOTLOGON_FILE"));
  TB.checkFile(gSystem->Getenv("CWB_PARAMETERS_FILE"));
  TB.checkFile(gSystem->Getenv("CWB_UPARAMETERS_FILE"));


  // data label
  char www_label[512];
  TObjArray* token = TString(work_dir).Tokenize(TString("/"));
  sprintf(www_label,((TObjString*)token->At(token->GetEntries()-1))->GetString().Data());

  char cmd[256];
  sprintf(cmd,"mkdir %s/%s",www_dir,www_label);
  cout << endl;
  cout << cmd << endl;
  gSystem->Exec(cmd);

  sprintf(cmd,"ln -s %s/%s %s/%s",work_dir,ced_dir,www_dir,www_label);
  cout << cmd << endl;
  cout << endl;
  gSystem->Exec(cmd);

  exit(0);
}
