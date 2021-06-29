/*
# Copyright (C) 2019 Francesco Salemi
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


void AddChip(TString filein,TString treename) { 
	CWB::CBCTool cbcTool;	
	TFile *f = new TFile(filein,"update"); 
	TTree *T = (TTree*)f->Get(treename);
	float mass[2];
  	float spin[6]; 
	float chip; 
        TBranch *bchip = T->Branch("chip",&chip,"chip/F");
	T->SetBranchAddress("mass", mass);
  	T->SetBranchAddress("spin", spin); 
	Long64_t nentries = T->GetEntries(); 
	for (Long64_t i=0;i<nentries;i++) { 
		T->GetEntry(i); 
		chip = cbcTool.chip(mass[0], mass[1], spin[0], spin[1], spin[2], spin[3], spin[4], spin[5]);
		//cout <<  mass[0] << " " << mass[1] << " " << spin[0] << " " << spin[1] << " " <<  spin[2] << " " << spin[3] << " " << spin[4] << " " <<spin[5] << endl;
		bchip->Fill(); 
	} 
	T->Print(); 
	T->Write(); 
	delete f; 
}
