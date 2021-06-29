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


// dump the list of injections : used by the cwb_dump command

{
  #include <vector>

  CWB::Toolbox TB;

  TB.checkFile(gSystem->Getenv("CWB_ROOTLOGON_FILE"));
  TB.checkFile(gSystem->Getenv("CWB_PARAMETERS_FILE"));
  TB.checkFile(gSystem->Getenv("CWB_UPARAMETERS_FILE"));

  network NET;         // network

  if(!simulation) {     
    cout << "cwb_dump_inj.C : dump injection type list works only in simulation mode !!!" << endl;
    exit(1);
  }

  if(TString(frFiles[nIFO])!="") {   // MDC from frame files : read MDC list from MDC frames log file

    cout << "cwb_dump_inj.C : Opening MDC frame file ... " << frFiles[nIFO] << endl;
    TB.checkFile(frFiles[nIFO]);
    CWB::frame ifr(frFiles[nIFO]);
    int nfrFiles=ifr.getNfiles();
    cout << "MDC " << " -> nfrFiles : " << nfrFiles << endl;
    waveSegment mdc_range  = ifr.getFrRange();
    cout << "mdc_range : " << mdc_range.start << " " << mdc_range.stop << endl;
    double GPS = (mdc_range.stop+mdc_range.start)/2;
    int nINJ = NET.readMDClog(injectionList,GPS); 
    printf("GPS: %16.6f saved,  injections: %d\n",GPS,nINJ);
    //for(int i=0;i<NET.mdcType.size();i++) cout << i << " " << NET.mdcType[i] << endl; 
    //for(int i=0;i<NET.mdcTime.size();i++) cout << i << " " << NET.mdcTime[i] << endl; 

  } else { 	// MDC from cWB engine : read MDC list from configPlugin

    cout << "cwb_dump_inj.C : read MDC list from configPlugin ... " << configPlugin.GetTitle() << endl;
    TB.checkFile(configPlugin.GetTitle());

    // export to CINT net,cfg
    char cmd[128];
    sprintf(cmd,"network* net = new network;");
    gROOT->ProcessLine(cmd);
    CWB::config* cfg = new CWB::config;
    cfg->Import();
    sprintf(cmd,"CWB::config* cfg = (CWB::config*)%p;",cfg);
    gROOT->ProcessLine(cmd);
    sprintf(cmd,"int gIFACTOR=1;");
    gROOT->ProcessLine(cmd);
    // exec config plugin
    configPlugin.Exec();
//    MDC.Print();

    // fill mdcType list
    NET.mdcType.clear();
    if(MDC.GetInspiral()!="") {		// inspiral MDC type
      NET.mdcType.push_back(MDC.GetInspName.Data());
    } else {                            // burst MDC type
      int size = MDC.wfList.size();
      for(int i=0;i<size;i++) {
        bool save=true;
        for(int j=0; j<(int)NET.mdcType.size(); j++){
          if(MDC.wfList[i].name.CompareTo(NET.mdcType[j])==0) {save = false; break;}
        }
        if(save) {
          NET.mdcType.push_back(MDC.wfList[i].name.Data());
        }
      }
    }  
    //cout.precision(14);
    //for(int k=0;k<(int)NET.mdcType.size();k++) cout << k << " mdcType " << NET.mdcType[k] << endl;
  }

  if(NET.mdcType.size()==0) {
    cout << endl << "cwb_dump_inj.C - Error : injection types not found !!! " << endl << endl;
    exit(1);
  }

  // dump the inj list to file

  char mdc_inj_file[1024];
  sprintf(mdc_inj_file,"%s/%s.inj",dump_dir,data_label);
  bool overwrite=TB.checkFile(mdc_inj_file,true);
  if(!overwrite) gSystem->Exit(1);

  FILE* fP;
  if((fP = fopen(mdc_inj_file, "w")) == NULL) {
    cout << "cwb_dump_inj.C : cannot open output file " << mdc_inj_file <<". \n";
    exit(1);
  };
  cout << "Write output file : " << mdc_inj_file << endl << endl;
  for(int i=0;i<NET.mdcType.size();i++) {
    char ostr[256];
    sprintf(ostr,"MDC_SET\t%4d\t%20s\t%6.0f\t%6.0f",i,NET.mdcType[i].data(),(fLow+fHigh)/2.,(fHigh-fLow)/2.);  
    cout << ostr << endl;
    fprintf(fP,"%s\n",ostr);
  }
  fclose(fP);

  char cmd[256];
  sprintf(cmd,"sort -k 3 -d %s > %s",mdc_inj_file,TString(mdc_inj_file).ReplaceAll(".inj",".inj.tmp").Data());
  //cout << cmd << endl;
  gSystem->Exec(cmd);
  sprintf(cmd,"mv %s %s",TString(mdc_inj_file).ReplaceAll(".inj",".inj.tmp").Data(),mdc_inj_file);
  //cout << cmd << endl;
  gSystem->Exec(cmd);
  cout << endl;
  cout << "Output inj MDC file : " << mdc_inj_file << endl << endl;

  exit(0);
}
