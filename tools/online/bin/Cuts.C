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


void Cuts()
{

TString IFILE;
TString CUTFILE;
TString CUT;
TString TIME;
cin >> IFILE >> CUTFILE >> CUT>> TIME;

gROOT->LoadMacro(CUTFILE.Data());

TFile* f=TFile::Open(IFILE.Data());
TTree* t=(TTree*)f->Get("waveburst");

TGlobal *global=(TGlobal*)gROOT->GetListOfGlobals()->FindObject(CUT.Data());
TCut *tcut = (TCut*)global->GetAddress();
cout << tcut->GetTitle() << endl;

//t->Draw("run","","goff");
//cout << t->GetSelectedRows() << endl;
char sel[2048];
sprintf(sel,"%s && abs(time[0]-%s)<.1",tcut->GetTitle(),TIME.Data());
cout << sel << endl;
t->Draw("run",sel,"goff");
cout << t->GetSelectedRows() << endl;


}
