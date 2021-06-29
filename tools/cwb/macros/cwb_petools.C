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


// this macro uses the following json package to read json file 
// https://github.com/nlohmann/json
// git clone https://github.com/nlohmann/json.git

R__ADD_INCLUDE_PATH(/home/waveburst/git/json/include)
#include <nlohmann/json.hpp>
using json = nlohmann::json;

#include "MACROS/ReadChunkList.C"

#define GW_MAX_SIZE		1000
#define CHUNK_FILE_LIST         "Chunk_List.txt"

#define PE_GIT_DIR		"/home/waveburst/git/pe"

struct evtparameters {
  TString name;
  double  time;
  int     chunk;
  int     start;
  int     stop;
  TString search;
  TString calib;
  TString run;
  TString tag;
};

struct peparameters {
  evtparameters evt;
  TString jsonFile;
  TString run;
  TString approx;
  float   fref;
  float   flow;
  TString webdir;
  TString pspath;
};

vector<evtparameters> ReadInputEventsFile(TString ifile);
vector<evtparameters> GetEventsParameters(TString ifile, TString run );
vector<TString> GetPreferredPE(vector<evtparameters> evtparms, TString pedir, TString run);
vector<peparameters> GetPreferredParmsPE(vector<evtparameters> evtparms, TString pedir, TString run, bool posteriors=false);

void DumpUpdateGIT(vector<evtparameters> evtparms, TString ofile="");
void DumpMkDirPE(TString ofile, vector<peparameters> vpeparms, TString run);
void DumpConfigPE(TString ofile, vector<peparameters> vpeparms, TString run);
void DumpXML(vector<peparameters> vpeparms, TString xmldir, bool execute=false);
void DumpMacroPosterior2XML(TString sample, peparameters peparms, TString odir);

void cwb_petools(TString ifile, TString option="", TString opath="", TString pe_git_dir=PE_GIT_DIR, TString run="O3") {
//
// this macro provides a set of methods used for the PE analysis
//

  // cwb_petools options
  if(ifile=="") {
    cout << endl;
    cout << "cwb_petools parameters" << endl;
    cout << endl;
    cout << " COMMAND: cwb_petools ifile option opath" << endl;
    cout << endl;
    cout << " IFILE: " << endl;
    cout << endl;
    cout << "  event list file" << endl;
    cout << endl;
    cout << "   file format: GW-NAME         GW-TIME         SEARCH  CALIB   PE-RUN  TAG  " << endl;
    cout << "                S190408an       1238782700.29   BBH     C00     NONE    dev1 " << endl;
    cout << endl;
    cout << " OPTIONS: " << endl;
    cout << endl;
    cout << "  1 - dump_update_git_pe       : Generates a script used to update the local git repository (remote: https://git.ligo.org/pe/O3)" << endl;
    cout << "  2 - dump_peparms             : Dump to file the PE parameters of the preferred posterior samples" << endl;
    cout << "  3 - dump_pe_xml_files        : Generates XML files for map, maxl & offsource samples" << endl;
    cout << "  4 - dump_pe_mkdirs           : Dump to file the instructions used to create the working directories" << endl;
    cout << "  5 - dump_cwb_pereport_config : Write the config used as imput for Makefile.pe" << endl;
    cout << endl;
    cout << " OPATH: " << endl;
    cout << endl;
    cout << "  output file/dir" << endl;
    cout << endl;
    exit(0);
  }

  vector<evtparameters> evtparms = GetEventsParameters(ifile, run); 

  vector<TString> jsonFile = GetPreferredPE(evtparms, pe_git_dir, run);
  for(int n=0;n<jsonFile.size();n++) cout << n << "\t" << jsonFile[n] << endl; cout << endl;

  if(option!="dump_pe_xml_files" && option!="dump_peparms" && option!="dump_pe_mkdirs" && 
     option!="dump_cwb_pereport_config" && option!="dump_update_git_pe") exit(0);

  TSystemFile odir(opath.Data(),"");
  if(option=="dump_pe_xml_files") {
    if(!odir.IsDirectory()) {cout << "cwb_petools - Error: \"" << opath << "\" is not a directory" << endl << endl; exit(1);}
  } else {
    if(opath!="") {
      bool overwrite = CWB::Toolbox::checkFile(opath,true); 
      if(!overwrite) exit(1);
    } else {
      cout << "cwb_petools - Error : output file not declared" << endl << endl;exit(1);
    }
  }

  if(option=="dump_update_git_pe") {
    DumpUpdateGIT(evtparms, opath);
    cout << endl << "Output file: " << opath << endl << endl;
    exit(0);
  }

  cout << endl;
  vector<peparameters> vpeparms = GetPreferredParmsPE(evtparms, pe_git_dir, run, false);
  if(option=="dump_peparms") {
    ofstream outpe;
    outpe.open(opath.Data());
    if (!outpe.good()) {cout << "Error Opening File : " << opath << endl;exit(1);}
    for(int n=0;n<vpeparms.size();n++) {
      outpe << endl;
      outpe << n+1 << "/" << vpeparms.size() << "\t" << vpeparms[n].evt.name << endl;
      outpe << endl;
      outpe << "jsonFile = " << vpeparms[n].jsonFile << endl;
      outpe << "run      = " << vpeparms[n].run << endl;
      outpe << "approx   = " << vpeparms[n].approx << endl;
      outpe << "fref     = " << vpeparms[n].fref << endl;
      outpe << "flow     = " << vpeparms[n].flow << endl;
      outpe << "webdir   = " << vpeparms[n].webdir << endl;
      outpe << "pspath   = " << vpeparms[n].pspath << endl;
      outpe << endl;
    }
    cout << endl << "Output file: " << opath << endl << endl;
  }

  if(option=="dump_cwb_pereport_config") {
    DumpConfigPE(opath, vpeparms, run);
    cout << endl << "Output file: " << opath << endl << endl;
  }

  if(option=="dump_pe_mkdirs") {
    DumpMkDirPE(opath, vpeparms, run);
    cout << endl << "Output file: " << opath << endl << endl;
  }

  if(option=="dump_pe_xml_files") {
    DumpXML(vpeparms, odir.GetName(),true);
  }

  exit(0);
}

vector<evtparameters> 
GetEventsParameters(TString ifile, TString run) {
//
// merge the user input list of event parameters with chunk parameters
//
// Input: 
//        ifile	- user input file: list of event parameters
//          run - O3 - Implemented only for O3
// 
// Returns the vector of event parameters + chunk parameters
//

  // get CWB_CONFIG
  char cwb_config_env[1024] = "";
  if(gSystem->Getenv("CWB_CONFIG")!=NULL) {
    strcpy(cwb_config_env,TString(gSystem->Getenv("CWB_CONFIG")).Data());
  }

  TString search="";
  char chunk_file_list[1024];
  sprintf(chunk_file_list,"%s/%s/CHUNKS/%s/%s",cwb_config_env,run.Data(),search.Data(),CHUNK_FILE_LIST);
  cout << chunk_file_list << endl;

  int    chunk[CHUNK_MAX_SIZE];
  double start[CHUNK_MAX_SIZE];
  double stop[CHUNK_MAX_SIZE];
  int nChunks = ReadChunkList(chunk_file_list,chunk,start,stop);

  vector<evtparameters> evtparms = ReadInputEventsFile(ifile);

  for(int n=0;n<evtparms.size();n++) {
    bool found=false;
    for(int i=0;i<nChunks;i++) {
      if(evtparms[n].time>start[i] && evtparms[n].time<=stop[i]) {
        printf("\t%d\t%s\t%.2f\tK%02d\n",n+1,evtparms[n].name.Data(),evtparms[n].time,chunk[i]);
        evtparms[n].chunk=chunk[i];
        evtparms[n].start=start[i];
        evtparms[n].stop=stop[i];
	found=true;
      }
    }
    if(!found) {
      evtparms[n].chunk=-1;
      evtparms[n].start=-1;
      evtparms[n].stop=-1;
    }
  }
  cout << endl;

  return evtparms;
}

void 
DumpUpdateGIT(vector<evtparameters> evtparms, TString ofile) {
//
// Generates a script used to update the local git repository (remote: https://git.ligo.org/pe/O3)
//
// Input: 
//     evtparms	- the vector of event parameters
//        ofile - output file name 
//
// output file format :
//
//  if [ ! -d S190408an ]; then 
//    git clone git@git.ligo.org:pe/O3/S190408an.git
//  fi
//  git lfs install --skip-smudge
//  cd S190408an/Preferred/PESummary_metafile
//  git lfs pull -I  posterior_samples.json
//  cd -
//

  ofstream outpe;
  outpe.open(ofile.Data());
  if (!outpe.good()) {cout << "Error Opening File : " << ofile << endl;exit(1);}
  outpe << "if [ ! -d O3 ]; then" << endl; 
  outpe << " mkdir O3" << endl;
  outpe << "fi" << endl;
  outpe << "cd O3" << endl;
  outpe << "" << endl;
  for(int k=0;k<evtparms.size();k++) {
    outpe << "if [ ! -d " << evtparms[k].name << " ]; then" << endl; 
    outpe << " git clone git@git.ligo.org:pe/O3/"<< evtparms[k].name << ".git" << endl;
    outpe << "fi" << endl;
    outpe << "git lfs install --skip-smudge" << endl;
    outpe << "cd " << evtparms[k].name << "/Preferred/PESummary_metafile" << endl;
    outpe << "git lfs pull -I  posterior_samples.json" << endl;
    outpe << "cd -" << endl;
    outpe << "" << endl;
  }
  outpe.close();
}

vector<peparameters> 
GetPreferredParmsPE(vector<evtparameters> vevtparms, TString pedir, TString run, bool posteriors) {
// 
// extract the PE parameters (fref, flow, pe-run, approximant) from json files
//
// Input:
//       
//     evtparms	- the vector of event parameters
//        pedir - the local PE git directory (/home/waveburst/git/pe)
//          run - O3 - Implemented only for O3
//   posteriors - true/false, if true it retrive the path of posterior_samples.dat file
//
// Return the vector of pe parameters
//

  vector<peparameters> vpeparms(vevtparms.size());

  for(int k=0;k<vevtparms.size();k++) {

    cout << "\t" << k+1 << "/" << vevtparms.size() << "\t" << "GetPreferredParmsPE\t" << vevtparms[k].name << endl;

    vpeparms[k].evt=vevtparms[k]; 

    // Get json file name
    TString jsonFile = TString::Format("%s/%s/%s/Preferred/PESummary_metafile/posterior_samples.json",
                                       pedir.Data(),run.Data(),vevtparms[k].name.Data());

    Long_t id,size=0,flags,mt;
    int estat = gSystem->GetPathInfo(jsonFile.Data(),&id,&size,&flags,&mt);
    if (estat!=0) continue;

    vpeparms[k].jsonFile=jsonFile; 

    // fix json file -> remove NaN, Infinity 
    unsigned int Pid = gSystem->GetPid();  // used to tag in a unique way the temporary files
    TString jsonFileTmp = TString::Format("/tmp/%d_%s",Pid,gSystem->BaseName(jsonFile));

    char cmd[1024];

    // copy json file to temporary dir
    sprintf(cmd,"cp %s %s",jsonFile.Data(),jsonFileTmp.Data());
    //cout << cmd << endl;
    gSystem->Exec(cmd);

    // replace NaN with 0, (avoid parser crash) 
    sprintf(cmd,"sed -i 's/NaN/0/g' %s",jsonFileTmp.Data());
    //cout << cmd << endl;
    gSystem->Exec(cmd);

    // replace Infinity with 0, (avoid parser crash) 
    sprintf(cmd,"sed -i 's/Infinity/0/g' %s",jsonFileTmp.Data());
    //cout << cmd << endl;
    gSystem->Exec(cmd);

    // open temp json file 
    std::ifstream ifile(jsonFileTmp.Data());
    json jroot;
    ifile >> jroot;

    json japprox = jroot["approximant"];

    // extract run (Ex: EXP2)
    TString run="";
    for (auto& el : japprox.items()) {
      //cout << el.key() << " : " << el.value() << "\n";
      string key = el.key();
      string value = el.value();
      if(TString(value.c_str()).Contains("IMRPhenomPv2")) run = key.c_str();
    }

    // substitute non alphanumeric characters with 'x'
    vpeparms[k].run = run;
    for(int i=0;i<run.Sizeof()-1;i++) if(!TString(run[i]).IsAlnum()) vpeparms[k].run[i]='x';
    // overwrite with user defined run (if != NONE)
    if(vevtparms[k].run!="NONE") vpeparms[k].run=vevtparms[k].run;

    // extract approximant
    try {
      string capprox = jroot["config_file"][run.Data()]["engine"]["approx"];
      TString sapprox = capprox.c_str();
      vpeparms[k].approx=sapprox;
    } catch (...) {
      cout << "extract approximant exception ...." << endl;
      vpeparms[k].approx="";
    }

    // extract fref
    try {
      string cfref = jroot["config_file"][run.Data()]["engine"]["fref"];
      vpeparms[k].fref = TString(cfref.c_str()).Atof();
    } catch (...) {
      cout << "extract fref exception ...." << endl;
      vpeparms[k].fref=-1;
    }

    // extract flow
    try {
      string cflow = jroot["config_file"][run.Data()]["lalinference"]["flow"];
      TString sflow = cflow.c_str();
      sflow = sflow(sflow.Last(':')+1,sflow.Last('}')-sflow.Last(':')-1);
      vpeparms[k].flow = sflow.Atof();
    } catch (...) {
      cout << "extract flow exception ...." << endl;
      vpeparms[k].flow=-1;
    }

    // extract webdir
    try {
      string cwebdir = jroot["config_file"][run.Data()]["paths"]["webdir"];
      vpeparms[k].webdir = cwebdir.c_str();
    } catch (...) {
      cout << "extract webdir exception ...." << endl;
      vpeparms[k].webdir="";
    }

    // search posterior_samples.dat file in webdir
    if(posteriors && vpeparms[k].webdir!="") {
      sprintf(cmd,"find %s -name posterior_samples.dat",vpeparms[k].webdir.Data());
      TString ocmd = gSystem->GetFromPipe(cmd);
      TObjArray* token = ocmd.Tokenize(TString('\n'));
      TString pspath="";
      for(int j=0;j<token->GetEntries();j++){
        TObjString* tok = (TObjString*)token->At(j);
        TString stok = tok->GetString();
        if(stok.Contains("/L1/")) continue;
        if(stok.Contains("/H1/")) continue;
        if(stok.Contains("/V1/")) continue;
        if(stok.Contains("/K1/")) continue;
        pspath=stok;
        //cout << j << " " << stok << endl;
      }
      if(pspath=="") pspath=vpeparms[k].webdir;
      vpeparms[k].pspath=pspath;
    } else {	// use converted json file
      //vpeparms[k].pspath=TString::Format("pesummary_%s.dat",vpeparms[k].run.Data());
      // save the non alpha numeric run 
      vpeparms[k].pspath=TString::Format("pesummary_%s.dat",run.Data());
    }

    // remove temporary json file 
    sprintf(cmd,"rm %s",jsonFileTmp.Data());
    //cout << cmd << endl;
    gSystem->Exec(cmd);
  }

  return vpeparms;
}

void 
DumpMkDirPE(TString ofile, vector<peparameters> vpeparms, TString run) {
// 
// Dump to file the instructions used to create the working directories
//
// Input:
//       
//        ofile - output file name 
//      peparms	- the vector of pe parameters
//          run - O3 - Implemented only for O3
//

  ofstream out;
  out.open(ofile.Data());
  if (!out.good()) {cout << "DumpMkDirPE: Error Opening File : " << ofile << endl;exit(1);}

  for(int k=0;k<vpeparms.size();k++) {

    cout << "\t" << k+1 << "/" << vpeparms.size() << "\t" << vpeparms[k].evt.name << endl << endl;
     out << "\t" << k+1 << "/" << vpeparms.size() << "\t" << vpeparms[k].evt.name << endl << endl;

    cout  << endl;
     out  << endl;

    TString cmd1 = TString::Format("cwb_mkchunk --run %s --chunk %02d --cal %s --net LH --search %s --type SIM/ONSPE --tag %s_%s_maxl_%s",
                                   run.Data(), vpeparms[k].evt.chunk, vpeparms[k].evt.calib.Data(), vpeparms[k].evt.search.Data(), 
				   vpeparms[k].evt.name.Data(), vpeparms[k].run.Data(), vpeparms[k].evt.tag.Data());
    cout  << cmd1 << endl;
     out  << cmd1 << endl;

    TString cmd2 = TString::Format("cwb_inet @ %.2f ced",vpeparms[k].evt.time);
    cout.precision(12);
    cout  << cmd2 << endl << endl;
     out.precision(12);
     out  << cmd2 << endl << endl;

    TString cmd3 = TString::Format("cwb_mkchunk --run %s --chunk %02d --cal %s --net LH --search %s --type SIM/OFSPE --tag %s_%s_%s",
                                   run.Data(), vpeparms[k].evt.chunk, vpeparms[k].evt.calib.Data(), vpeparms[k].evt.search.Data(), 
				   vpeparms[k].evt.name.Data(), vpeparms[k].run.Data(), vpeparms[k].evt.tag.Data());
    cout  << cmd3 << endl;
     out  << cmd3 << endl;

    TString cmd4 = TString::Format("cwb_ppchunk --run %s --chunk %02d --cal %s --net LH --search %s --type SIM/OFSPE --tag %s_%s_%s --opt=all",
                                   run.Data(), vpeparms[k].evt.chunk, vpeparms[k].evt.calib.Data(), vpeparms[k].evt.search.Data(), 
				   vpeparms[k].evt.name.Data(), vpeparms[k].run.Data(), vpeparms[k].evt.tag.Data());
    cout  << cmd4 << endl << endl;
     out  << cmd4 << endl << endl;

    TString cmd5 = TString::Format("cwb_mkchunk --run %s --chunk %02d --cal %s --net LH --search %s --type SIM/ONSPE --tag %s_%s_pe_%s",
                                   run.Data(), vpeparms[k].evt.chunk, vpeparms[k].evt.calib.Data(), vpeparms[k].evt.search.Data(), 
				   vpeparms[k].evt.name.Data(), vpeparms[k].run.Data(), vpeparms[k].evt.tag.Data());
    cout  << cmd5 << endl;
     out  << cmd5 << endl;

    TString cmd6 = TString::Format("cwb_ppchunk --run %s --chunk %02d --cal %s --net LH --search %s --type SIM/ONSPE --tag %s_%s_pe_%s --opt=all",
                                   run.Data(), vpeparms[k].evt.chunk, vpeparms[k].evt.calib.Data(), vpeparms[k].evt.search.Data(), 
				   vpeparms[k].evt.name.Data(), vpeparms[k].run.Data(), vpeparms[k].evt.tag.Data());
    cout  << cmd6 << endl << endl;
     out  << cmd6 << endl << endl;
  }

  out.close();
}

void 
DumpConfigPE(TString ofile, vector<peparameters> vpeparms, TString run) {
// 
// Write the config used as imput for Makefile.AnalizePE
//
// Input:
//       
//        ofile - output file name 
//      peparms	- the vector of pe parameters
//          run - O3 - Implemented only for O3
//

  ofstream out;
  out.open(ofile.Data());
  if (!out.good()) {cout << "DumpConfigPE: Error Opening File : " << ofile << endl;exit(1);}

  out << endl << "STATUS=FALSE" << endl << endl;

  for(int k=0;k<vpeparms.size();k++) {

    TString onpath  = TString::Format("%s/%s_K%02d_%s_LH_%s_SIM_ONSPE_%s_%s_maxl_%s",
                      vpeparms[k].evt.name.Data(), run.Data(), vpeparms[k].evt.chunk, vpeparms[k].evt.calib.Data(), 
                      vpeparms[k].evt.search.Data(), vpeparms[k].evt.name.Data(), vpeparms[k].run.Data(), vpeparms[k].evt.tag.Data());

    TString offpath = TString::Format("%s/%s_K%02d_%s_LH_%s_SIM_OFSPE_%s_%s_%s",
                      vpeparms[k].evt.name.Data(), run.Data(), vpeparms[k].evt.chunk, vpeparms[k].evt.calib.Data(), 
                      vpeparms[k].evt.search.Data(), vpeparms[k].evt.name.Data(), vpeparms[k].run.Data(), vpeparms[k].evt.tag.Data());

    out << "ifeq ($(GW_NAME)," << vpeparms[k].evt.name << ")" << endl;
    out << "GW_GPS        = " << vpeparms[k].evt.name << endl;
    out << "GW_TSTEP      = 150." << endl;
    out << "GW_ONPATH     = " << onpath << endl;
    out << "GW_OFFPATH    = " << offpath << endl;
    out << "GW_RDIR       = pestat_maxl_" << vpeparms[k].evt.tag << endl;
    out << "STATUS        = TRUE" << endl;
    out << "endif" << endl;
    out << endl;
  }

  out << "GW_WDIR = $(shell basename $(GW_OFFPATH))" << endl;
  out << "GW_MTAG = M1.C_U.S_bin1_cut" << endl;
  out << endl;
  out << "GW_DTAG = $(GW_WDIR).$(GW_MTAG)" << endl;
  out << endl;
  out << "WWW_DIR = " << endl;
  out << endl;

  out.close();
}

void 
DumpXML(vector<peparameters> vpeparms, TString xmldir, bool execute) {
// 
// Generates XML files for map, maxl & offsource samples
//
// Input:
//       
//      peparms	- the vector of pe parameters
//       xmldir - output xml directory 
//      execute - false -> creates directories & Posteriors2XML macros
//                 true -> creates & execute Posteriors2XML macros   
//

  // write Posterior2XML macros
  for(int k=0;k<vpeparms.size();k++) {
    if(vpeparms[k].jsonFile=="" || vpeparms[k].run=="") continue;
    TString odir = xmldir+"/"+vpeparms[k].evt.name+"/PE/"+vpeparms[k].run;
    char cmd[1024];

    // create output directory
    sprintf(cmd,"mkdir -p %s",odir.Data());
    cout << cmd << endl;
    gSystem->Exec(cmd);

    // copy json file
    sprintf(cmd,"cp %s %s",vpeparms[k].jsonFile.Data(),odir.Data());
    cout << cmd << endl;
    gSystem->Exec(cmd);

    //convert json file to dat 
    TString json2dat = TString("python ")+gSystem->ExpandPathName("$CWB_SCRIPTS/cwb_json2dat.py");
    TString setup = "source /home/vedovato/virtualenv/pycwb/bin/activate";
    sprintf(cmd,"cd %s;%s;%s %s",odir.Data(),setup.Data(),json2dat.Data(),vpeparms[k].jsonFile.Data());
    cout << cmd << endl;
    gSystem->Exec(cmd);

    // rename non alphanumeric pspath to alphanumeric
    TString pspath = TString::Format("pesummary_%s.dat",vpeparms[k].run.Data());
    sprintf(cmd,"cd %s;mv \"%s\" %s",odir.Data(),vpeparms[k].pspath.Data(),pspath.Data());
    cout << cmd << endl;
    gSystem->Exec(cmd);
    vpeparms[k].pspath=pspath;

    // create Posterior2XML macros
    DumpMacroPosterior2XML("map", vpeparms[k], odir);
    DumpMacroPosterior2XML("maxl", vpeparms[k], odir);
    DumpMacroPosterior2XML("pe", vpeparms[k], odir);
    DumpMacroPosterior2XML("offsource", vpeparms[k], odir);

    // excecute Posterior2XML macros
    odir = odir+"/xml";
    sprintf(cmd,"mkdir -p %s",odir.Data());
    gSystem->Exec(cmd);
    if(execute) {	// generate XML files
      sprintf(cmd,"cd %s; root -l -b ../Posterior2XML_MAP.C",odir.Data()); gSystem->Exec(cmd);
      sprintf(cmd,"cd %s; root -l -b ../Posterior2XML_MAXL.C",odir.Data()); gSystem->Exec(cmd);
      //sprintf(cmd,"cd %s; root -l -b ../Posterior2XML_ONSOURCE.C",odir.Data()); gSystem->Exec(cmd);
      sprintf(cmd,"cd %s; root -l -b ../Posterior2XML_OFFSOURCE.C",odir.Data()); gSystem->Exec(cmd);
    }
    cout << endl << "output xml directory: " << odir << endl << endl;
  }
}

vector<TString> 
GetPreferredPE(vector<evtparameters> evtparms, TString pedir, TString run) {
// 
// extract the json files names
//
// Input:
//       
//     evtparms	- the vector of event parameters
//        pedir - the local PE git directory (/home/waveburst/git/pe)
//          run - O3 - Implemented only for O3
//
// Return the vector of json files
//

  vector<TString> jsonFile;
  for(int k=0;k<evtparms.size();k++) {
    TString json_file_path = TString::Format("%s/%s/%s/Preferred/PESummary_metafile/posterior_samples.json",pedir.Data(),run.Data(),evtparms[k].name.Data());
    jsonFile.push_back(json_file_path);
  }
  return jsonFile;
}

void 
DumpMacroPosterior2XML(TString sample, peparameters peparms, TString odir) {
// 
// creates Posteriors2XML macros
//
// Input:
//       sample - map/maxl/onsource/offsource
//      peparms	- the vector of pe parameters
//         odir - output directory 
//
// Output:
//		Posterior2XML_MAP.C
//		Posterior2XML_MAXL.C					
//		Posterior2XML_ONSOURCE.C
//		Posterior2XML_OFFSOURCE.C
//

  TString ofile;
  TString xmlFile;
  TString clbFile;

  if(sample=="map") {
    xmlFile = TString::Format("posterior_samples_%s_%s_map.xml",peparms.evt.name.Data(),peparms.run.Data());
    ofile = odir+"/Posterior2XML_MAP.C";
  }
  if(sample=="maxl") {
    xmlFile = TString::Format("posterior_samples_%s_%s_maxl.xml",peparms.evt.name.Data(),peparms.run.Data());
    ofile = odir+"/Posterior2XML_MAXL.C";
  }
  if(sample=="pe") {
    xmlFile = TString::Format("posterior_samples_%s_%s_pe.xml",peparms.evt.name.Data(),peparms.run.Data());
    ofile = odir+"/Posterior2XML_ONSOURCE.C";
  }
  if(sample=="offsource") {
    xmlFile = TString::Format("posterior_samples_%s_%s_offsource.xml",peparms.evt.name.Data(),peparms.run.Data());
    ofile = odir+"/Posterior2XML_OFFSOURCE.C";
  }
  clbFile = xmlFile; clbFile.ReplaceAll(".xml",".clb");

  ofstream out;
  out.open(ofile.Data());
  if (!out.good()) {cout << "DumpMacroPosterior2XML- Error Opening File : " << ofile << endl;exit(1);}

  out << "{" << endl;
  out << "  TString options = \"\";" << endl;
  out << endl;
  out << "  TString posteriorFile = \"" << "../" << peparms.pspath << "\";" << endl;
  out << endl;
  if(sample=="map") {
    out << "  options += \"--sample map \";" << endl;
    out << "  options += \"--ninjections 1 \";" << endl;
  }
  if(sample=="maxl") {
    out << "  options += \"--sample maxl \";" << endl;
    out << "  options += \"--ninjections 1 \";" << endl;
  }
  if(sample=="offsource") {
    out << "  options += \"--gps_start_time " << TString::Format("%d",peparms.evt.start) << " \";" << endl;
    out << "  options += \"--gps_stop_time  " << TString::Format("%d",peparms.evt.stop)  << " \";" << endl;
    out << endl;
    out << "  options += \"--time_step 150 \";" << endl;
    out << "  options += \"--seed 1 \";" << endl;
  }
  if(sample=="pe") {
    out << "  // extract from posteriorFile the number of entries		" << endl;
    out << "  char cmd[1024];							" << endl;
    out << "  sprintf(cmd,\"wc %s\",posteriorFile.Data());			" << endl;
    out << "  TString ocmd = gSystem->GetFromPipe(cmd);				" << endl;
    out << "  TObjArray* token = ocmd.Tokenize(TString(' '));			" << endl;
    out << "  TObjString* tok = (TObjString*)token->At(0);			" << endl;
    out << "  TString stok = tok->GetString();					" << endl;
    out << endl;
    out << "  cout << \"ocmd \" << ocmd << endl;				" << endl;
    out << "  int max_sampleID = stok.Atoi()-1;					" << endl;
    out << "  cout << \"max_sampleID \" << max_sampleID << endl;		" << endl;
    out << endl;
    out << "  gRandom->SetSeed(150914);						" << endl;
    out << endl;
    out << "  CWB::Toolbox::mkDir(\"pe\",true);				" << endl;
    out << endl;
    out << "  for(int i=0;i<1000;i++) {						" << endl;
    out << endl;
    out << "    int sampleID = gRandom->Uniform(0,max_sampleID);		" << endl;
    out << endl;
    out << "    cout << endl;							" << endl;
    out << "    cout << i << \" -> sampleID = \" << sampleID << endl << endl;	" << endl;
    out << endl;
    out << "    options = \"\";							" << endl;
    out << endl;
    out << "    options += TString::Format(\"--sample %d \",sampleID);		" << endl;
    out << "    options += \"--ninjections 1 \";				" << endl;
    out << endl;
    out << "    options += \"--source "   << peparms.evt.name 	<< " \";" << endl;
    out << "    options += \"--waveform " << peparms.approx   	<< " \";" << endl;
    out << "    options += \"--f_ref "    << peparms.fref     	<< " \";" << endl;
    out << "    options += \"--f_lower "  << peparms.flow     	<< " \";" << endl;
    out << endl;
    out << "    TString xmlFile = \"" << "pe/" << xmlFile 	<< "\";" << endl;
    out << "    TString clbFile = \"" << "pe/" << clbFile 	<< "\";" << endl;
    out << "    TString xml = TString::Format(\"_%d.xml\",i+1);			" << endl;
    out << "    xmlFile.ReplaceAll(\".xml\",xml);				" << endl;
    out << "    TString clb = TString::Format(\"_%d.clb\",i+1);			" << endl;
    out << "    clbFile.ReplaceAll(\".clb\",clb);				" << endl;
    out << "    options += \"--clb_file \"+clbFile+\" \";			" << endl;
    out << endl;
    out << "    CWB::mdc::Posterior2XML(posteriorFile, xmlFile, options);	" << endl;
    out << "  }									" << endl;
    out << endl;
  } else {
    out << endl;
    out << "  options += \"--source "   << peparms.evt.name << " \";" << endl;
    out << "  options += \"--waveform " << peparms.approx   << " \";" << endl;
    out << "  options += \"--f_ref "    << peparms.fref     << " \";" << endl;
    out << "  options += \"--f_lower "  << peparms.flow     << " \";" << endl;
    out << endl;
    out << "  TString xmlFile       = \"" << xmlFile << "\";" << endl;
    out << "  options += \"--clb_file   " << clbFile << " \";" << endl;
    out << endl;
    out << "  CWB::mdc::Posterior2XML(posteriorFile, xmlFile, options);" << endl;
    out << endl;
  }
  out << "  exit(0);" << endl;
  out << endl;
  out << "}" << endl;

  out.close();
}

vector<evtparameters> 
ReadInputEventsFile(TString ifile) {
//
// read the user input file: list of event parameters
//
// file format: GW-NAME         GW-TIME         SEARCH  CALIB   PE-RUN  TAG
//              S190408an       1238782700.29   BBH     C00     NONE    dev1
//
// Input: 
//        ifile	- user input file: list of event parameters
//          run - O3 - Implemented only for O3
// 
// Returns the vector of event parameters
//

  CWB::Toolbox::checkFile(ifile);

  ifstream in;
  in.open(ifile.Data(),ios::in);
  if (!in.good()) {cout << "cwb_petools: Error Opening File : " << ifile << endl;exit(1);}

  int isize=0;
  char str[1024];
  int fpos=0;
  while(true) {
    in.getline(str,1024);
    if (!in.good()) break;
    if(str[0] != '#') isize++;
  }
  in.clear(ios::goodbit);
  in.seekg(0, ios::beg);
  if(isize==0)  {cout << "cwb_petools: Error : File " << ifile << " is empty" << endl;exit(1);}
  if(isize>GW_MAX_SIZE) {cout << "cwb_petools: Error : File " << ifile << " > " << GW_MAX_SIZE << endl;exit(1);}

  char   name[256];
  double time;
  char   search[256];
  char   calib[256];
  char   run[256];
  char   tag[256];

  cout.precision(10); 

  cout << endl;
  int k=0;
  vector<evtparameters> evtparms(isize);
  while (1) {
    fpos=in.tellg();
    in.getline(str,1024);
    if(str[0] == '#') continue;
    in.seekg(fpos, ios::beg);

    in >> name >> time >> search >> calib >> run >> tag;
    if (!in.good()) break;
    fpos=in.tellg();
    in.seekg(fpos+1, ios::beg);

    evtparms[k].name=name;
    evtparms[k].time=time;
    evtparms[k].search=search;
    evtparms[k].calib=calib;
    evtparms[k].run=run;
    evtparms[k].tag=tag;
    k++;
  }
  evtparms.resize(k);
  in.close();

  return evtparms;
}

