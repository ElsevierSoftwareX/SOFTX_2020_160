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


// merge plugins into one plugin
{

  #define nPLUGIN 6

  CWB::Toolbox TB;

  TB.checkFile(gSystem->Getenv("CWB_ROOTLOGON_FILE"));
  TB.checkFile(gSystem->Getenv("CWB_PARAMETERS_FILE"));
  TB.checkFile(gSystem->Getenv("CWB_NETC_FILE"));

  // get used defined plugins up to nPLUGIN plugins
  TString cwb_mplugin[nPLUGIN];
  TMacro  cwb_mplugin_macro[nPLUGIN];
  for(int i=0;i<nPLUGIN;i++) {
    char env_label[32];sprintf(env_label,"CWB_MPLUGIN_%d",i+1);
    cwb_mplugin[i]=TString(gSystem->Getenv(env_label));
    if(i==0) { 		// this is the output plugin
      if(cwb_mplugin[i].CompareTo("")!=0) {
        TString DirName  = gSystem->DirName(cwb_mplugin[0]);
        TString BaseName = gSystem->BaseName(cwb_mplugin[0]);
        TB.checkFile(DirName);
        bool overwrite=TB.checkFile(cwb_mplugin[i],true);
        if(!overwrite) gSystem->Exit(1);
      } else {
        cout << "cwb_mplugin.C error : output file name not declared !!!" << endl;
        gSystem->Exit(1);
      }
    } else {
      if(cwb_mplugin[i].CompareTo("")!=0) {
        TB.checkFile(cwb_mplugin[i]);
        cwb_mplugin_macro[i] = TMacro(cwb_mplugin[i]); 
      }
    }
  }

  // print plugin list
  cout << "Output Plugin Name : " << endl << endl; 
  cout << 0 << " -> : " << cwb_mplugin[0] << endl;
  cout << endl << "Input Plugin Names : " << endl << endl; 
  for(int i=1;i<nPLUGIN;i++) {
    if(cwb_mplugin[i].CompareTo("")!=0) {
      cout << i << " -> " << cwb_mplugin[i] << endl;
    }
  }
  cout << endl;
  char answer[256];
  strcpy(answer,"");
  do {
    cout << "Do you want to continue it ? (y/n) ";
    cin >> answer;
    cout << endl << endl;
  } while ((strcmp(answer,"y")!=0)&&(strcmp(answer,"n")!=0));
  if (strcmp(answer,"n")==0) gSystem->Exit(0);


  // write plugin, if plugin is defined then we merge cwb_inet_plugin & plugin
  ofstream out;
  out.open(cwb_mplugin[0].Data(),ios::out);
  if(!out.good()) {cout << "cwb_mplugin.C - Error : Opening File : " 
                        << cwb_mplugin[0] << endl;gSystem->Exit(1);}

  // write warning message
  out << "// -------------------------------------------------------------------------" << endl;
  out << "// WARNING, NOT EDIT : This is a multi plugin generated with cwb_mplugin !!!" << endl;
  out << "// NOTE              : The main is listed on the bottom                      " << endl;
  out << "// " << endl;
  out << "// INPUT PLUGINS     :                                                       " << endl;
  out << "// " << endl;
  for(int i=1;i<nPLUGIN;i++) {
    if(cwb_mplugin[i].CompareTo("")!=0) {
      out << "// " << i << " - " << gSystem->BaseName(cwb_mplugin[i]) << endl; 
    }
  }
  out << "// " << endl;
  out << "// -------------------------------------------------------------------------" << endl;

  // write user user plugin code
  for(int i=1;i<nPLUGIN;i++) {
    if(cwb_mplugin[i].CompareTo("")!=0) {
      bool include=true;
      char user_plugin_name[32];sprintf(user_plugin_name,"CWB_UserPlugin_%d(",i);
      char user_plugin_namespace[32];sprintf(user_plugin_namespace,"CWB_UserPluginNamespace_%d",i);
      TList* fLines = cwb_mplugin_macro[i].GetListOfLines();
      TObjString *obj;
      TIter next(fLines);
      out << endl << endl;
      out << "// -------------------------------------------------------------------------" << endl;
      out << "// --> BEGIN CWB_USER PLUGIN CODE " << i << endl;
      out << "// " << endl; 
      out << "// " << gSystem->BaseName(cwb_mplugin[i]) << endl; 
      out << "// -------------------------------------------------------------------------" << endl;
      out << endl;
      while ((obj = (TObjString*) next())) {
        TString line = obj->GetName();
        TString sline = line;
        sline.ReplaceAll(" ","");
        if(sline.BeginsWith("#include")||sline.BeginsWith("#define")||sline.BeginsWith("#pragma")) {
          if(!include) out << "}" << endl;
          include=true; 
        } else {
          if(sline!="") {
            if(include) out << "namespace " << user_plugin_namespace << " {" << endl << endl;
            include=false;
          }
        }
        out <<  line.Data() << endl;
      }
      out << "}" << endl << endl;
      out << "// -------------------------------------------------------------------------" << endl;
      out << "// --> END CWB_USER PLUGIN CODE "   << i << endl;
      out << "// " << endl; 
      out << "// " << gSystem->BaseName(cwb_mplugin[i]) << endl; 
      out << "// -------------------------------------------------------------------------" << endl;
    }
  }

  // write user plugin declarations
  out << endl << endl;
  out << "// -------------------------------------------------------------------------" << endl;
  out << "// --> MAIN CWB_USER PLUGIN CODE " << endl;
  out << "// " << endl;
  out << "// INPUT PLUGINS     :                                                       " << endl;
  out << "// " << endl;
  for(int i=1;i<nPLUGIN;i++) {
    if(cwb_mplugin[i].CompareTo("")!=0) {
      out << "// " << i << " - " << gSystem->BaseName(cwb_mplugin[i]) << endl; 
    }
  }
  out << "// " << endl;
  out << "// -------------------------------------------------------------------------" << endl;
  out << endl;
  out << "#define XIFO 4" << endl;
  out << "#pragma GCC system_header" << endl;
  out << "#include \"cwb.hh\"" << endl;

  // write mplugin main
  out << endl << endl;
  out << "void" << endl;
  out << "CWB_Plugin(TFile* jfile, CWB::config* cfg, network* net, WSeries<double>* x, TString ifo, int type)  {" << endl;
  out << endl;
  for(int i=1;i<nPLUGIN;i++) {
    if(cwb_mplugin[i].CompareTo("")!=0) {
      char user_plugin_namespace[32];sprintf(user_plugin_namespace,"CWB_UserPluginNamespace_%d",i);
      out << "  " << user_plugin_namespace;
      out << "::CWB_Plugin(jfile, cfg, net, x, ifo, type); // CALL USER PLUGIN CODE " << i << endl;
    }
  }
  out << "}" << endl;

  out.close();

  cout << "Output Plugin Name : " << cwb_mplugin[0] << endl << endl;

  gSystem->Exit(0);
}
