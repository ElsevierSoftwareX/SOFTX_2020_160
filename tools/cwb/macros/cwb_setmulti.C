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


{
  TString cwb_merge_label=TString(gSystem->Getenv("CWB_MERGE_LABEL"));
  
  sprintf(net_file_name,"wave_%s.%s.root",data_label,cwb_merge_label.Data());
  cout << net_file_name << endl;

  CWB::Toolbox TB;
  TB.setMultiplicity(net_file_name,merge_dir,merge_dir,nIFO,Tgap);

  exit(0);
}
