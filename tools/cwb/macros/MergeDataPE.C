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


int ReadConfig(TString ifconfig, vector<TString> &gwname, vector<TString> &wpath, vector<TString> &rdir);

void MergeDataPE(TString side, TString odir, TString ifconfig) {

  cout<<"MergeDataPE.C init ..."<<endl;

//  CWB::Toolbox::checkFile(odir);	// check if odir exist

  gSystem->Exec("date");

  vector<TString> gwname; 
  vector<TString> wpath;
  vector<TString> rdir;
  int gwsize = ReadConfig(ifconfig, gwname, wpath, rdir);

  vector<TString> rpath(gwsize);

//  for(int i=0;i<gwsize;i++) rpath[i] = wpath[i]+"/report/dump/"+rdir[i]+"/"+side;
//  for(int i=0;i<gwsize;i++) cout << rpath[i] << endl;

  for(int i=0;i<gwsize;i++) {
    rpath[i] = wpath[i]+"/report/dump/"+rdir[i]+"/"+side+"/distributions";
    char cmd1[1024];
    char cmd2[1024];
    char cmd3[1024];
    if(i==0) {
      sprintf(cmd1,"cat %s/*_ResidualEnergy.txt > %s/ResidualEnergy.txt",rpath[i].Data(),odir.Data());
      sprintf(cmd2,"cat %s/*_FittingFactor.txt  > %s/FittingFactor.txt",rpath[i].Data(),odir.Data());
      sprintf(cmd3,"cat %s/*_OverlapFactor.txt  > %s/OverlapFactor.txt",rpath[i].Data(),odir.Data());
    } else {
      sprintf(cmd1,"cat %s/*_ResidualEnergy.txt >> %s/ResidualEnergy.txt",rpath[i].Data(),odir.Data());
      sprintf(cmd2,"cat %s/*_FittingFactor.txt  >> %s/FittingFactor.txt",rpath[i].Data(),odir.Data());
      sprintf(cmd3,"cat %s/*_OverlapFactor.txt  >> %s/OverlapFactor.txt",rpath[i].Data(),odir.Data());
    }
    cout << cmd1 << endl;
    gSystem->Exec(cmd1);
    cout << cmd2 << endl;
    gSystem->Exec(cmd2);
    cout << cmd3 << endl;
    gSystem->Exec(cmd3);
  }

  for(int i=0;i<gwsize;i++) {
    rpath[i] = wpath[i]+"/report/dump/"+rdir[i]+"/"+side+"/distributions";
    char cmd1[1024];
    char cmd2[1024];
    char cmd3[1024];
    if(i==0) {
      sprintf(cmd1,"cat %s/*_ResidualEnergy_sensitivity.txt > %s/ResidualEnergy_sensitivity.txt",rpath[i].Data(),odir.Data());
      sprintf(cmd2,"cat %s/*_FittingFactor_sensitivity.txt  > %s/FittingFactor_sensitivity.txt",rpath[i].Data(),odir.Data());
      sprintf(cmd3,"cat %s/*_OverlapFactor_sensitivity.txt  > %s/OverlapFactor_sensitivity.txt",rpath[i].Data(),odir.Data());
    } else {
      sprintf(cmd1,"cat %s/*_ResidualEnergy_sensitivity.txt >> %s/ResidualEnergy_sensitivity.txt",rpath[i].Data(),odir.Data());
      sprintf(cmd2,"cat %s/*_FittingFactor_sensitivity.txt  >> %s/FittingFactor_sensitivity.txt",rpath[i].Data(),odir.Data());
      sprintf(cmd3,"cat %s/*_OverlapFactor_sensitivity.txt  >> %s/OverlapFactor_sensitivity.txt",rpath[i].Data(),odir.Data());
    }
    cout << cmd1 << endl;
    gSystem->Exec(cmd1);
    cout << cmd2 << endl;
    gSystem->Exec(cmd2);
    cout << cmd3 << endl;
    gSystem->Exec(cmd3);
  }

  for(int i=0;i<gwsize;i++) {
    rpath[i] = wpath[i]+"/report/dump/"+rdir[i]+"/"+side+"/distributions";
    char cmd1[1024];
    char cmd2[1024];
    char cmd3[1024];
    if(i==0) {
      sprintf(cmd1,"cat %s/*_ResidualEnergy_onsource_prob.txt > %s/ResidualEnergy_onsource_prob.txt",rpath[i].Data(),odir.Data());
      sprintf(cmd2,"cat %s/*_FittingFactor_onsource_prob.txt  > %s/FittingFactor_onsource_prob.txt",rpath[i].Data(),odir.Data());
      sprintf(cmd3,"cat %s/*_OverlapFactor_onsource_prob.txt  > %s/OverlapFactor_onsource_prob.txt",rpath[i].Data(),odir.Data());
    } else {
      sprintf(cmd1,"cat %s/*_ResidualEnergy_onsource_prob.txt >> %s/ResidualEnergy_onsource_prob.txt",rpath[i].Data(),odir.Data());
      sprintf(cmd2,"cat %s/*_FittingFactor_onsource_prob.txt  >> %s/FittingFactor_onsource_prob.txt",rpath[i].Data(),odir.Data());
      sprintf(cmd3,"cat %s/*_OverlapFactor_onsource_prob.txt  >> %s/OverlapFactor_onsource_prob.txt",rpath[i].Data(),odir.Data());
    }
    cout << cmd1 << endl;
    gSystem->Exec(cmd1);
    cout << cmd2 << endl;
    gSystem->Exec(cmd2);
    cout << cmd3 << endl;
    gSystem->Exec(cmd3);
  }

  exit(0);
}

int ReadConfig(TString ifconfig, vector<TString> &gwname, vector<TString> &wpath, vector<TString> &rdir) {

  ifstream in;
  in.open(ifconfig.Data(),ios::in);
  if (!in.good()) {cout << "MergeDataPE::ReadConfig Error Opening File : " << ifconfig.Data() << endl;exit(1);}

  bool start=false;
  char log[1024];
  while(true) {
    in.getline(log,1024);
    if (!in.good()) break;
    if(log[0]=='#') continue;
    TString line = log;
    if(line.Contains("ifeq") && line.Contains("$(GW_NAME)")) {
      start=true;
      TString gw = line(line.Index(',')+1, line.Last(')')-line.Index(',')-1);
      gwname.push_back(gw);
      //cout << gw << endl;
    }
    if(line.Contains("endif") && start==true) start=false;
    if(start) {
      if(line.Contains("PE_OFFPATH")) {
        TString path = line(line.Index('/'), line.Sizeof()-line.Index(',')-1);
        wpath.push_back(path);
        //cout << path << endl;
      }
      if(line.Contains("PE_RDIR")) {
        TString dir = line(line.Index('=')+1, line.Sizeof()-line.Index(',')-2);
	dir.ReplaceAll(" ","");
        rdir.push_back(dir);
        //cout << dir << endl;
      }
    }
  }

  in.close();

  if((wpath.size() != gwname.size())) {
    cout << "MergeDataPE::ReadConfig - missing PE_OFFPATH config declarations: MISSED = "
         << gwname.size()-wpath.size() << " in file " << ifconfig << endl;exit(1);
  }
  if((rdir.size() != gwname.size())) {
    cout << "MergeDataPE::ReadConfig - missing PE_RDIR config declarations: MISSED = " 
         << gwname.size()-rdir.size() << " in file " << ifconfig << endl;exit(1);
  }

  for(int i=0;i<gwname.size();i++) cout << i << "\t" << gwname[i] << endl;
  for(int i=0;i<gwname.size();i++) cout << i << "\t" << wpath[i] << endl;
  for(int i=0;i<gwname.size();i++) cout << i << "\t" << rdir[i] << endl;

  return gwname.size();
}

