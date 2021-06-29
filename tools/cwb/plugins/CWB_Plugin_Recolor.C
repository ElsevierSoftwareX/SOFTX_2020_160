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


#define XIFO 4

#pragma GCC system_header

#include "cwb.hh"
#include "config.hh"
#include "network.hh"
#include "wavearray.hh"
#include "TString.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "TRandom.h"
#include "Toolbox.hh"

void ReadUserOptions(CWB::config* cfg);
void PrintUserOptions(CWB::config* cfg, int ifoId);

// ---------------------------------------------------------------------------------
// plugin parameters
// ---------------------------------------------------------------------------------

TString frDir[NIFO_MAX]; 	// output directory for frame files
TString frName[NIFO_MAX];	// name to tag the frame 
TString chName[NIFO_MAX]; 	// name to tag the frame channel data
TString frLabel[NIFO_MAX];	// label used to tag the output frame file name
                                // Ex: frDir/ifo_frLabel-GPS-LENGTH.gwf 
TString psdFile[NIFO_MAX]; 	// file used to recolor the data

// NOTE : cwb whitened data is 0 for freq<16
// for freq<padFreq the amplitude in frequency domain is padded 
// with a value given by -> (famplitude at padFreq) * padFactor

int padFreq[NIFO_MAX];		
int padFactor[NIFO_MAX];

bool gausNoise[NIFO_MAX];	// true/false -> gaussian/real data colored noise

/* Example : parPlugin definition in user_parameters.C file for L1,H1 network
   Note    : definitions must be repeated for each detector          

  TString optrc = "";           // NOTE : add space at the end of each line
  optrc += "frDir=/home/vedovato/WP/RECOLOR/O1_C01_RECOLOR_tst1/recolored/L1 ";
  optrc += "frDir=/home/vedovato/WP/RECOLOR/O1_C01_RECOLOR_tst1/recolored/H1 ";
  optrc += "frLabel=HOFT_C01 ";
  optrc += "frLabel=HOFT_C01 ";
  optrc += "chName=L1:RECOLOR-CALIB_STRAIN_C01 ";
  optrc += "chName=H1:RECOLOR-CALIB_STRAIN_C01 ";
  optrc += "frName=LLO_4k ";
  optrc += "frName=LHO_4k ";
  optrc += "psdFile=$HOME_CWB/plugins/strains/advLIGO_NSNS_Opt_8khz_one_side.txt ";
  optrc += "psdFile=$HOME_CWB/plugins/strains/advLIGO_NSNS_Opt_8khz_one_side.txt ";
  optrc += "padFreq=20 ";
  optrc += "padFreq=20 ";
  optrc += "padFactor=1000 ";
  optrc += "padFactor=1000 ";
  optrc += "gausNoise=false ";
  optrc += "gausNoise=false ";

  strcpy(parPlugin,optrc.Data());      // set plugin parameters
*/
/* Note : define DQF array with CAT0 DQ files 
  
  // dq file list
  // {ifo, dqcat_file, dqcat[0/1/2], shift[sec], inverse[false/true], 4columns[true/false]}
  nDQF=2;
  dqfile dqf[nDQF]={
    {"L1" ,"/home/vedovato/O1/DQvetos/O1_12Sep19Jan_C01/L1Cat0.txt",   CWB_CAT0, 0., false, false},
    {"H1" ,"/home/vedovato/O1/DQvetos/O1_12Sep19Jan_C01/H1Cat0.txt",   CWB_CAT0, 0., false, false}
  };
  for(int i=0;i<nDQF;i++) DQF[i]=dqf[i];
*/

void 
CWB_Plugin(TFile* jfile, CWB::config* cfg, network* net, WSeries<double>* x, TString ifo, int type)  {
//!MISCELLANEA
// Plugin to produce recolored frames of whitened data
// First data are whitened then are recolored and save to file frames 
//

  cout << endl;
  cout << "-----> CWB_Plugin_Recolor.C" << endl;
  cout << "ifo " << ifo.Data() << endl;
  cout << "type " << type << endl;
  cout << endl;

  CWB::Toolbox TB;

  if(type==CWB_PLUGIN_CONFIG) {
    cfg->levelR=0; 		// no resampling, data are whitened with yjr original sample rate
    cfg->simulation = 0;	// force non simulation mode
    cfg->nfactor=1;

    ReadUserOptions(cfg);	// user config options : read from parPlugin
  }

  if(type==CWB_PLUGIN_INIT_JOB) {
    // get ifo id
    int id=-1;
    for(int n=0;n<cfg->nIFO;n++) if(ifo==net->ifoName[n]) {id=n;break;}
    if(id<0) {cout << "Plugin : Error - bad ifo id" << endl; gSystem->Exit(1);}

    PrintUserOptions(cfg, id);  // print config options
  }

  if(type==CWB_PLUGIN_ODATA_CONDITIONING) {

    // get ifo id
    int id=-1;
    for(int n=0;n<cfg->nIFO;n++) if(ifo==net->ifoName[n]) {id=n;break;}
    if(id<0) {cout << "Plugin : Error - bad ifo id" << endl; gSystem->Exit(1);}

    TB.mkDir(frDir[id],false,false);	// create output frame dir

    // data recolored
    int size=x->size();
    double start=x->start();
    if(gausNoise[id]) {
      int seed = 1000+id;
      TB.getSimNoise(*x, psdFile[id], seed, net->nRun);			// gaussian colored noise
    } else {
      TB.getSimNoise(*x, psdFile[id], -padFreq[id], padFactor[id]);	// real data colored noise
    }
 
    x->resize(size);
    x->start(start);

    // save cleaned data into frame
    wavearray<double> X = *x;
    // remove scratch
    int os = cfg->segEdge*x->rate();
    X.start(x->start()+cfg->segEdge);
    X.stop(x->stop()-cfg->segEdge);
    for(int i=0;i<X.size()-2*os;i++) X[i]=X[i+os];
    X.resize(x->size()-2*os);

    // create frame file
    char frFile[1024];
    sprintf(frFile,"%s/%s_%s-%d-%d.gwf",
            frDir[id].Data(),ifo.Data(),frLabel[id].Data(),
            int(X.start()),(int)(X.stop()-X.start()));

    // open frame file
    CWB::frame fr(frFile,chName[id],"WRITE");

    // write frame to file
    fr.writeFrame(X, frName[id], chName[id]);
    cout << "CWB_Plugin_Recolor.C : write " << frFile << endl;

    fr.close();

    int nIFO=net->ifoListSize();
    if(TString(ifo).CompareTo(net->ifoName[nIFO-1])==0) gSystem->Exit(0);  // last ifo
  }

  return;
}

void ReadUserOptions(CWB::config* cfg) {

  TString options = cfg->parPlugin;

  int n_frDir=0;
  int n_frName=0;
  int n_frLabel=0;
  int n_chName=0;
  int n_psdFile=0;
  int n_padFreq=0;
  int n_padFactor=0;
  int n_gausNoise=0;
  for(int i=0;i<cfg->nIFO;i++) {
    frDir[i]="";
    frName[i]="";
    frLabel[i]="";
    chName[i]="";
    psdFile[i]="";
    padFreq[i]=2;
    padFactor[i]=0;
    gausNoise[i]=false;
  }
  if(options.CompareTo("")!=0) {
    cout << options << endl;
    if(!options.Contains("--")) {  // parameters are used only by cwb_inet

      TObjArray* token = TString(options).Tokenize(TString(' '));
        for(int j=0;j<token->GetEntries();j++){

        TObjString* tok = (TObjString*)token->At(j);
        TString stok = tok->GetString();

        if(stok.Contains("frDir=")) {
          frDir[n_frDir]=stok;
          frDir[n_frDir].Remove(0,frDir[n_frDir].Last('=')+1);
          frDir[n_frDir]=gSystem->ExpandPathName(frDir[n_frDir].Data());
          if(n_frDir<(NIFO_MAX-1)) n_frDir++;
        }

        if(stok.Contains("frName=")) {
          frName[n_frName]=stok;
          frName[n_frName].Remove(0,frName[n_frName].Last('=')+1);
          if(n_frName<(NIFO_MAX-1)) n_frName++;
        }

        if(stok.Contains("frLabel=")) {
          frLabel[n_frLabel]=stok;
          frLabel[n_frLabel].Remove(0,frLabel[n_frLabel].Last('=')+1);
          if(n_frLabel<(NIFO_MAX-1)) n_frLabel++;
        }

        if(stok.Contains("chName=")) {
          chName[n_chName]=stok;
          chName[n_chName].Remove(0,chName[n_chName].Last('=')+1);
          if(n_chName<(NIFO_MAX-1)) n_chName++;
        }

        if(stok.Contains("padFreq=")) {
          TString opt=stok;
          opt.Remove(0,opt.Last('=')+1);
          if(opt.IsDigit()) padFreq[n_padFreq]=opt.Atoi();
          if(padFreq[n_padFreq]<2) padFreq[n_padFreq]=2;
          if(n_padFreq<(NIFO_MAX-1)) n_padFreq++;
        }

        if(stok.Contains("padFactor=")) {
          TString opt=stok;
          opt.Remove(0,opt.Last('=')+1);
          if(opt.IsDigit()) padFactor[n_padFactor]=opt.Atoi();
          if(padFactor[n_padFactor]<0) padFactor[n_padFactor]=0;
          if(n_padFactor<(NIFO_MAX-1)) n_padFactor++;
        }

        if(stok.Contains("psdFile=")) {
          psdFile[n_psdFile]=stok;
          psdFile[n_psdFile].Remove(0,psdFile[n_psdFile].Last('=')+1);
          psdFile[n_psdFile]=gSystem->ExpandPathName(psdFile[n_psdFile].Data());
          if(n_psdFile<(NIFO_MAX-1)) n_psdFile++;
        }

        if(stok.Contains("gausNoise=")) {
          TString opt=stok;
          opt.Remove(0,opt.Last('=')+1);
          if(opt=="false")  gausNoise[n_gausNoise]=false;
          if(opt=="true")   gausNoise[n_gausNoise]=true;
          if(n_gausNoise<(NIFO_MAX-1)) n_gausNoise++;
        }
      }
    }
  }

  for(int i=0;i<cfg->nIFO;i++) {
    if(frDir[i]=="")     {cout << "CWB_Plugin_Recolor : Error - frDir["<<cfg->ifo[i]<<"]     not defined" << endl; gSystem->Exit(1);}
    if(frName[i]=="")    {cout << "CWB_Plugin_Recolor : Error - frName["<<cfg->ifo[i]<<"]    not defined" << endl; gSystem->Exit(1);}
    if(chName[i]=="")    {cout << "CWB_Plugin_Recolor : Error - chName["<<cfg->ifo[i]<<"]    not defined" << endl; gSystem->Exit(1);}
    if(frLabel[i]=="")   {cout << "CWB_Plugin_Recolor : Error - frLabel["<<cfg->ifo[i]<<"]   not defined" << endl; gSystem->Exit(1);}
    if(psdFile[i]=="")   {cout << "CWB_Plugin_Recolor : Error - psdFile["<<cfg->ifo[i]<<"]   not defined" << endl; gSystem->Exit(1);}
    if(i>=n_gausNoise)   {cout << "CWB_Plugin_Recolor : Error - gausNoise["<<cfg->ifo[i]<<"] not defined" << endl; gSystem->Exit(1);}
  }
}

void PrintUserOptions(CWB::config* cfg, int ifoId) {

    cout << "-----------------------------------------"     << endl;
    cout << "Recolor config options : ifo = " << cfg->ifo[ifoId] << endl;
    cout << "-----------------------------------------"     << endl << endl;
    cout << endl;
    cout << "frDir        : " << frDir[ifoId]     << endl;
    cout << "chName       : " << chName[ifoId]    << endl;
    cout << "frName       : " << frName[ifoId]    << endl;
    cout << "frLabel      : " << frLabel[ifoId]   << endl;
    cout << "psdFile      : " << psdFile[ifoId]   << endl;
    cout << "padFreq      : " << padFreq[ifoId]   << endl;
    cout << "padFactor    : " << padFactor[ifoId] << endl;
    cout << "gausNoise    : " << gausNoise[ifoId] << endl;
    cout << endl;
}
