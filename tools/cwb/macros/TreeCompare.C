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


/*#include "ObjComparer.h"
#include "TFile.h"
#include "TTree.h"
#include "TSystem.h"
#include <iostream>
#include "/work/salemi/git/cWB/library/wat/netevent.hh"
*/
/*void printCDBEntry(const char * fileName) {
  std::cout << "printing CDB file " << fileName << "\n";
  TFile *file = TFile::Open(fileName);
  if (!file) return;
  TObject *obj = (TObject*)file->Get("AliCDBEntry");

  ROOTObjComparer cmp;
  cmp.Print(obj);

  // for debugging purposes, can dump XML here via ROOT serialization
  // dump them to XML just to have a text comparison also
  TString outxml(fileName); outxml+=".xml";
  TFile *file1xml = TFile::Open(outxml.Data(), "RECREATE");
  if (file1xml){
    std::cout << "Writing XML for comparision\n";
    file1xml->WriteTObject(obj);
    file1xml->Close();
  }

  file->Close();
}
*/

void TreeCompare(const char * filename,  const char *postfix1, const char *postfix2, TString treename="waveburst", float tollerance = 1e-4) {
  TString fileName1 = filename;
  TString fileName2 = fileName1;
  fileName2.ReplaceAll(postfix1,postfix2);
  std::cout << "processing cWB file1 " << fileName1.Data() << "\n";
  std::cout << "processing cWB file2 " << fileName2.Data() << "\n";
  cout.precision(17);
  TFile *file1 = TFile::Open(fileName1);
  if (!file1) return;
  TTree *obj1 = (TTree*)file1->Get(treename);
 // netevent EVT1(obj1,2);  

  TFile *file2 = TFile::Open(fileName2);
  if (!file2) return;
  TTree *obj2 = (TTree*)file2->Get(treename);
 // netevent EVT2(obj2,2);
/*
  ROOTObjComparer cmp;
  cmp.SetVerbose(false);
  cmp.SetDebug(false);
  cmp.Diff(obj1, obj2);
*/
  // dump events
  int nEVT1 = obj1->GetEntries(); 
  int nEVT2 = obj2->GetEntries(); 
  
  //check if the number of events is matching between the 2 TTrees
  if(nEVT1 != nEVT2){
    std::cout << "Unequal number of events in the 2 TTrees!!! nEVT1=" << nEVT1 <<" nEVT2=" << nEVT2 << " Exiting now..." << "\n";	
    exit(1);
  }  
  //check if the TTrees are empty
  if(nEVT1 == 0) {std::cout << "Empty ntuple!!!" << endl; exit(0);}
  int index;
  //sort the TTrees
//  if(treename!=liveTime) {obj1->BuildIndex("time[0]","time[1]");}
 // else {
  obj1->BuildIndex("start[0]","start[1]");  //}

  TTreeIndex *I1=(TTreeIndex*)obj1->GetTreeIndex(); // get the tree index
  Long64_t* index1=I1->GetIndex(); //create an array of entries in sorted order
  //obj2->BuildIndex("time[0]","time[1]");
  obj2->BuildIndex("start[0]","start[1]");
  TTreeIndex *I2=(TTreeIndex*)obj2->GetTreeIndex(); // get the tree index
  Long64_t* index2=I2->GetIndex(); //create an array of entries in sorted order
 
  TObjArray *branchList; 
  branchList  = obj1->GetListOfBranches();
  int nBranch     = obj1->GetNbranches();
  TString varnames[63];
  TLeaf *var1[63];
  TLeaf *var2[63];
  double value_1;
  double value_2;
 for (Long64_t jentry=0; jentry<nEVT1;jentry++) { 
    obj1->GetEntry(jentry);
    obj2->GetEntry(jentry);
    //obj1->Show(jentry);
    for(int i=0;i<nBranch;i++){
	varnames[i] = branchList->At(i)->GetName();
	var1[i] = obj1->GetBranch(varnames[i])->GetLeaf(varnames[i]);
	var2[i] = obj2->GetBranch(varnames[i])->GetLeaf(varnames[i]);
	for(int k=0;k<obj1->GetBranch(varnames[i])->GetLeaf(varnames[i])->GetLen();k++){
	    value_1 = var1[i]->GetValue(k);
	    value_2 = var2[i]->GetValue(k);
	    if (TMath::Abs(value_1-value_2) > tollerance){
		cout <<jentry <<" " << varnames[i].Data() << "["<<k <<"] " << value_1 << " " << value_2 << " " << TMath::Abs(value_1-value_2) <<endl;
	    }
	}
    }
 }
  file1->Close();
  file2->Close();
}

