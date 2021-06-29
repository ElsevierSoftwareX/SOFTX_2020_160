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


void SelectTriggers()
{
  char file[1024];
  char name[1024];
  char start[1024];
  char stop[1024];
  char cut[2048];

  cin >> file >> name >> start >> stop;
  printf("%s %s %s %s \n",file, name, start, stop);
  sprintf(cut,"time[0]>%s && time[0]<%s",start,stop);
  cout << cut << endl;
  TFile *ifile = TFile::Open(file);
  TTree* itree = (TTree *) ifile->Get("waveburst");
  sprintf(file,"trigger_%s.root",name);
  TFile ofile(file,"RECREATE");
  TTree* otree = itree->CopyTree(cut);
  otree->SetDirectory(&ofile);
  ifile->Close();
  otree->Write();
  ofile.Close();
  
}
