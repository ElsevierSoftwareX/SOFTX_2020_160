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


void addHistory()
{
   char ifile[1024];
   char ofile[2014];
   cin >> ofile >> ifile;
   TFile *f_ifile = TFile::Open(ifile);

   CWB::History* hist = (CWB::History *) f_ifile->Get("history");

   TFile f_ofile(ofile,"UPDATE");
   hist->Write("history");
   f_ofile.Close();

   f_ifile->Close();
}
