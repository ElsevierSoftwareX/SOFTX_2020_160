/*
# Copyright (C) 2019 Marco Drago
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

   char ifile[1024];
   char ofile[1024];
   int start,stop;
   cin >> ifile >> ofile >> start >> stop;

   TFile *f_ifile = TFile::Open(ifile);

   TTree* livetree = (TTree *) f_ifile->Get("liveTime");
   TTree* otree = (TTree*) livetree->CloneTree(0);
   double live;
   livetree->SetBranchAddress("live",&live);
   otree->SetBranchAddress("live",&live);

   TFile* f_ofile = TFile::Open(ofile,"RECREATE");
   otree->SetDirectory(f_ofile);
   for (int i=0; i<livetree->GetEntries(); i++)
     {
       livetree->GetEntry(i);
       live=stop-start;
       otree->Fill();
     }
   otree->Write();
   f_ofile->Close();
   f_ifile->Close();







}
