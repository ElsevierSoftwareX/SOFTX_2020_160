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


// apply veto to merged root file : used by the cwb_setveto command

{
  if(nIFO==2 && TString(ifo[1])==ifo[0]) nIFO=1;	// set single detector mode when nIFO=2 and ifo[1]=ifo[0] 

  vector<TString> ifos(nIFO);
  for(int i=0;i<nIFO;i++) ifos[i]=ifo[i];

  cwb_merge_label=TString(gSystem->Getenv("CWB_MERGE_LABEL"));
  
  sprintf(net_file_name,"wave_%s.%s.root",data_label,cwb_merge_label.Data());
  cout << net_file_name << endl;

  CWB::Toolbox::setVeto(net_file_name,merge_dir,merge_dir,ifos,nVDQF,VDQF,nDQF,DQF,segLen,segMLS,segEdge);

  exit(0);
}
