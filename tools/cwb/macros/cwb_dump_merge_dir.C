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


// dump the list of labels of the merged files in the merge dir
// used by the commands : cwb_setveto, cwb_setmulti, cwb_setcuts, cwb_report, cwb_merge

{
  CWB::Toolbox TB;

  vector<TString> dirList = TB.getFileListFromDir(merge_dir, ".root","","wave_");

  cout << endl;
  cout << " --------------------" << endl; 
  cout << " List of merge labels" << endl; 
  cout << " --------------------" << endl; 
  cout << endl;
  for(int i=0;i<dirList.size();i++) {
    if(dirList[i].Contains("wave")) {
      dirList[i].ReplaceAll(data_label,"");
      dirList[i].ReplaceAll(merge_dir,"");
      dirList[i].ReplaceAll("/","");
      dirList[i].ReplaceAll("wave_.","");
      dirList[i].ReplaceAll(".root","");
      cout << " - " << dirList[i].Data() << endl;
    }
  }
  cout << endl;

  exit(0);
}
