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


void AdjustSlag()
{

   char ifile[1024];
   char ofile[1024];
   char ss_lag[1024];
   char sshift[1024];
   cin >> ifile >> ofile >> ss_lag>>sshift;
   TObjArray* token = TString(ss_lag).Tokenize(TString(","));
   TObjArray* token2 = TString(sshift).Tokenize(TString(","));
   int const nifo=token2->GetEntries();
   int shift[nifo];
   float slag[nifo+1];
   for (int i=0; i<nifo; i++)
     { 
        TString s=((TObjString*)token->At(i))->GetString();
        slag[i]=s.Atoi();
        TString s2=((TObjString*)token2->At(i))->GetString();
        shift[i]=s2.Atoi();
     }
   TString s=((TObjString*)token->At(nifo))->GetString();
   slag[nifo]=s.Atoi();
   double time[2*nifo];

   TFile *f_ifile = TFile::Open(ifile);

   TTree* witree = (TTree *) f_ifile->Get("waveburst");
   witree->SetBranchAddress("time",time);
   TTree* wotree = witree->CloneTree(0);
   TFile f_ofile(ofile,"RECREATE");
   wotree->SetDirectory(&f_ofile);
   wotree->SetBranchAddress("slag",slag);
   wotree->SetBranchAddress("time",time);
 
   for (int i=0; i<witree->GetEntries(); i++)
    {
     witree->GetEvent(i);
     for (int j=0; j<nifo; j++) time[j]-=shift[j];
     wotree->Fill();
    }
   cout << wotree->GetEntries() << endl;

   TTree* itree = (TTree *) f_ifile->Get("liveTime");
   TTree* otree = itree->CloneTree(0);
   otree->SetDirectory(&f_ofile);
   otree->SetBranchAddress("slag",slag);
 
   for (int i=0; i<itree->GetEntries(); i++)
    {
     itree->GetEvent(i);
     otree->Fill();
    }
   cout << otree->GetEntries() << endl;

   CWB::History* hist = (CWB::History *) f_ifile->Get("history");

   hist->Write("history");
   wotree->Write();
   otree->Write();
   f_ofile.Close();

   f_ifile->Close();
}
