/*
# Copyright (C) 2020 Gabriele Vedovato
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
#include "gwavearray.hh"
#include "TString.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "TRandom.h"
#include "mdc.hh"
#include "gwavearray.hh"
#include <vector>

// ---------------------------------------------------------------------------------
// WHAT IS?
// this plugin can be used in simulation mode to apply calibration errors to the injected signals 
// HOW TO CONFIGURE THE CE PLUGIN
// the following is an example : must be included in the config/user_parameters.C file
// see the DumpUserOptions function for the full description of parameters
// ---------------------------------------------------------------------------------
/*
  plugin = TMacro(gSystem->ExpandPathName("$HOME_CWB/plugins/CWB_Plugin_CE.C"));        // Macro source

  TString optce = "";                  		// NOTE : add space at the end of each line
  optce += "ce_type=uniform ";                  // calibration error type                                             -> det 1 
  optce += "ce_eamp=0.1 ";			// max percentage -> amplitude miscalibration  (uniform in -0.9,+1.1) -> det 1
  optce += "ce_ephs=10 ";			// max phase (degrees) -> phase miscalibration (uniform in -10,+10)   -> det 1
  optce += "ce_etim=0.005 ";			// max time (sec) -> time miscalibration (uniform in -0.005,+0.005)   -> det 1
  optce += "ce_type=fixed ";                    // ...                                                                -> det 2
  optce += "ce_eamp=0.1 ";			// ...                                                                -> det 2
  optce += "ce_ephs=10 ";			// ...                                                                -> det 2     
  optce += "ce_etim=0 ";			// ...                                                                -> det 2     
  optce += "ce_seed=1234 ";			// seed used by CE for ramdom generation

  strcpy(parPlugin,optce.Data());      		// set CE plugin parameters
  strcpy(comment,"ce configuration example");

*/

// ---------------------------------------------------------------------------------
// DEFINES
// ---------------------------------------------------------------------------------

#define CE_eAMP 		1		// max percentage of amplitude error  (1 -> no amplitude error)
#define CE_ePHS 		0		// max phase (degrees) phase error (0 -> no phase error)
#define CE_eTIM 		0		// max time (sec) error (0 -> no time error)
#define CE_SEED 		150914		// seed used by CE for ramdom generation
#define CE_TYPE			"fixed"		// CE type: fixed(default)/uniform/gussian
						//    fixed:    the ce_(eamp/ephs/etim) is used as calibration error for all IFOs 
						//    uniform:  the ce_(eamp/ephs/etim) is used to randomly extract the calibration error 
                                                //              from a uniform distribution in the range [-ce_(eamp/ephs/etim),+ce_(eamp/ephs/etim)] 
						//    guassian: the ce_(eamp/ephs/etim) is used to randomly extract the calibration error 
                                                //              from a gaussian distribution with median (1/0) and sigma=ce_(eamp/ephs/etim)
						//    file:     use the ce_file
						//              the file format is the same used by the posterior samples output PE analysis
#define CE_CFILE		""		//		is the file that contains the list of eamp/epha vs frequency
						//              the file format is the same used by the posterior samples output PE analysis
#define CE_GDUMP                false           // if true ->   dump in/out injections to root file in the working dir 'report/dump', used only for debugging
#define CE_VERBOSE              false           // if true ->   output debug infos

// ---------------------------------------------------------------------------------
// FUNCTIONS
// ---------------------------------------------------------------------------------

void SetCalibrationErrors(wavearray<double>& x, int ifoid);
void SetCalibrationErrors(wavearray<double>& x, TString ifo, double time);
double GetNearestCalibrationEntry(double time, int& simulation_id);

// ---------------------------------------------------------------------------------
// USER CONFIG OPTIONS
// ---------------------------------------------------------------------------------

struct uoptions {
  float   eamp[NIFO_MAX];
  float   ephs[NIFO_MAX];
  float   etim[NIFO_MAX];
  TString type[NIFO_MAX];
  TString cfile;
  int     seed;
  bool    gdump;
  bool    verbose;
};

void ResetUserOptions();
void ReadUserOptions(TString options);
void PrintUserOptions(CWB::config* cfg);
void DumpUserOptions(TString odir, CWB::config* cfg);

// ---------------------------------------------------------------------------------
// Global Variables
// ---------------------------------------------------------------------------------

uoptions	gOPT;                        	// global User Options
TString		gFDUMP;				// file name used to dump inj comparison into temporary root file
CWB::mdc*       gMDC;
  
void 
CWB_Plugin(TFile* jfile, CWB::config* cfg, network* NET, WSeries<double>* x, TString ifo, int type)  {
//!MISCELLANEA
// Plugin used for Calibration Errors

  if(type==CWB_PLUGIN_CONFIG) {  

    ResetUserOptions(); 				// set default config options
    ReadUserOptions(cfg->parPlugin); 			// user config options : read from parPlugin

    gFDUMP="";

    cout << endl;
    for(int n=0;n<cfg->nIFO;n++) {

      if(gOPT.type[n]=="file") {
        // Check if calibration file exist
        cout << "CWB_Plugin_CE : check if calibration file defined in plugin parameter ce_cfile=" << gOPT.cfile << " exist" << endl;
        CWB::Toolbox::checkFile(gOPT.cfile);
      }
    }
    cout << endl;
  }

  if(type==CWB_PLUGIN_NETWORK) {
    PrintUserOptions(cfg);      			// print config options

    gMDC = new CWB::mdc(NET);				// init gMDC
    for(int n=0;n<cfg->nIFO;n++) {
      if(gOPT.type[n]=="file") {
        gMDC->SetInspiralCLB(gOPT.cfile);		// setup calibration file
      }
    }
  }

  if(type==CWB_PLUGIN_MDC) {				// stage after reading injections

    // import ifactor
    int gIFACTOR=-1; IMPORT(int,gIFACTOR)
    int ifactor = (cfg->simulation==4) ? gIFACTOR : 0;

    // get ifo id
    int ifoid=-1;
    for(int n=0;n<cfg->nIFO;n++) if(ifo==NET->ifoName[n]) {ifoid=n;break;}
    if(ifoid<0) {cout << "CWB_Plugin_CE : Error - bad ifo ifoid" << endl; gSystem->Exit(1);}

    gwavearray<double> ginj_in;				// injection before calibration error (in)
    gwavearray<double> ginj_out;			// injection after calibration error (out)

    // init the output root file name used to store plots in/out injections
    if((gOPT.gdump) && (gFDUMP=="")) {
      TString jname = jfile->GetPath();
      jname.ReplaceAll(":/","");
      jname.ReplaceAll(".root","_ce_gdump.root");
      gFDUMP=jname;
      TFile *froot = new TFile(gFDUMP, "RECREATE");
      froot->Close();
      cout << endl << "CWB_Plugin_CE : Created file: " << gFDUMP << endl << endl; 
    }

    // loop over the injected events and apply amp/phs/tim errors to the injections
    for(int nmdc=0; nmdc<(int)NET->mdcListSize(); nmdc++) {	

      int seed = gOPT.seed+NET->nRun+ifoid+ifactor+nmdc; 
      gRandom->SetSeed(seed);				// for each injection/run/ifo/factor set a unique seed (used CE random generation)

      if(gOPT.verbose) cout << "CWB_Plugin_CE : injection_id=" << nmdc << "\trun=" << NET->nRun
                            << "\tifoid=" << ifoid << "\tfactor=" << ifactor << "\tseed=" << seed << endl;

      TString mdcstring(NET->getmdcList(nmdc));
      TObjArray* token = mdcstring.Tokenize(' ');
      TObjString* iname = (TObjString*)token->At(11);
      TString wavename = iname->GetString();
      TObjString* itime = (TObjString*)token->At(10);
      TString wavetime = itime->GetString();
      double mdctime = wavetime.Atof();
      if (mdctime<x->start()||mdctime>x->start()+x->size()/x->rate()) continue;
      //cout << mdcstring.Data() << endl;
      //printf(" Time : %s %f %f %f\t", wavetime.Data(), mdctime, x->start(), x->start()+x->size()/x->rate());
      //cout << "String: " << wavename.Data() << " time : " << wavetime.Data() << " " <<mdctime << endl;
      int istart = (mdctime - x->start()-cfg->iwindow/2.)*x->rate();
      int istop  = (mdctime - x->start()+cfg->iwindow/2.)*x->rate();
      if(istart<0) istart=0;
      if(istop>(int)x->size()) istop=x->size();
      int isize = istop-istart;	

      // save in injection into ginj_in 
      ginj_in.resize(isize);
      ginj_in.start(0);
      ginj_in.rate(x->rate());
      for(int i=0; i<isize; i++) ginj_in[i] = x->data[i+istart];

      // copy jnjection into a temporary wavearray X
      wavearray<double> X = ginj_in;

      if(gOPT.type[ifoid]=="file") {

        // Apply Amplitude/Phase Calibration Errors from file gOPT.cfile
        SetCalibrationErrors(X, ifo, mdctime);

      } else {

        // Apply Amplitude/Phase Calibration Errors
        SetCalibrationErrors(X, ifoid); 
      }

      // copy calibrated data into the original data wavearray x
      for(int i=0; i<isize; i++) x->data[i+istart] = X[i];

      // save in injection into ginj_out
      ginj_out.resize(isize);
      ginj_out.start(0);
      ginj_out.rate(x->rate());
      for(int i=0; i<isize; i++) ginj_out[i] = x->data[i+istart];

      // store in/out injection into root file
      if(gOPT.gdump) {

        TFile *froot = new TFile(gFDUMP, "UPDATE");

        watplot* plot;
        TString name;  

        ginj_in.Draw(GWAT_TIME);
        ginj_in.Draw(&ginj_out,GWAT_TIME,"SAME",kRed);
        plot = ginj_in.GetWATPLOT();              // get pointer to watplot object  
        name = TString::Format("%s_INJ_Time_InOut_%d",ifo.Data(),TMath::Nint(mdctime));
        plot->canvas->Write(name.Data());

        ginj_in.Draw(GWAT_FFT);
        ginj_in.Draw(&ginj_out,GWAT_FFT,"SAME",kRed);
        plot = ginj_in.GetWATPLOT();              // get pointer to watplot object  
        name = TString::Format("%s_INJ_FFT_InOut_%d",ifo.Data(),TMath::Nint(mdctime));
        plot->canvas->Write(name.Data());

        froot->Close();

        // move root file from tmp to report/dump directory
        if(ifoid==cfg->nIFO-1) {
          TString ofdump = gFDUMP;
          ofdump.ReplaceAll(TString::Format("%s/",cfg->tmp_dir).Data(),TString::Format("%s/",cfg->dump_dir).Data());
          char cmd[1024];
          sprintf(cmd,"mv %s %s",gFDUMP.Data(),ofdump.Data());
          cout << endl << cmd << endl << endl;
          gSystem->Exec(cmd);
        }
      }
    }
  }

  return;
}

void PrintUserOptions(CWB::config* cfg) {

  cout << "-----------------------------------------"     << endl;
  cout << "CE config options                        "     << endl;
  cout << "-----------------------------------------"     << endl << endl;
  for(int n=0;n<cfg->nIFO;n++) {
    cout << "CE_TYPE              " << cfg->ifo[n] << "   " << gOPT.type[n] << endl;
    cout << "CE_eAMP              " << cfg->ifo[n] << "   " << gOPT.eamp[n] << endl;
    cout << "CE_ePHS              " << cfg->ifo[n] << "   " << gOPT.ephs[n] << endl;
    cout << "CE_eTIM              " << cfg->ifo[n] << "   " << gOPT.etim[n] << endl << endl;
  }
  cout << "CE_CFILE             " << gOPT.cfile           << endl;
  cout << "CE_SEED              " << gOPT.seed            << endl;
  cout << "CE_GDUMP             " << gOPT.gdump           << endl;
  cout << "CE_VERBOSE           " << gOPT.verbose         << endl;

  cout << endl;
}

void DumpUserOptions(TString odir, CWB::config* cfg) {

  TString ofName = odir+"/ce_config.txt";

  ofstream out;
  out.open(ofName,ios::out);
  if(!out.good()) {cout << "DumpUserOptions : Error Opening File : " << ofName << endl;exit(1);}

  TString info="";
  TString tabs="\t\t\t\t";

  char version[128];
  sprintf(version,"WAT Version : %s - SVN Revision : %s - Tag/Branch : %s",watversion('f'),watversion('r'),watversion('b'));

  out << endl;
  out << "--------------------------------"     << endl;
  out << "CE config options               "     << endl;
  out << "--------------------------------"     << endl;
  out << endl;

  for(int n=0;n<cfg->nIFO;n++) {
    out << "ce_eamp  " << cfg->ifo[n] << "   " << gOPT.eamp[n] << endl;
  }
  info = "// max percentage of amplitude miscalibration : def(0) -> disabled";
  out << tabs << info << endl;

  for(int n=0;n<cfg->nIFO;n++) {
    out << "ce_ephs  " << cfg->ifo[n] << "   " << gOPT.ephs[n] << endl;
  }
  info = "// max phase (degrees) miscalibration : def(0) -> disabled";
  out << tabs << info << endl;

  for(int n=0;n<cfg->nIFO;n++) {
    out << "ce_etim  " << cfg->ifo[n] << "   " << gOPT.etim[n] << endl;
  }
  info = "// max time (sec) miscalibration : def(0) -> disabled";
  out << tabs << info << endl;

  for(int n=0;n<cfg->nIFO;n++) {
    out << "ce_type         " << cfg->ifo[n] << "   " << gOPT.type[n] << endl;
  }

  out << "ce_seed              " << gOPT.seed << endl;
  info = "// seed used by CE for random generation - 0(def) -> random seed";
  out << tabs << info << endl;

  out << "ce_gdump         " << gOPT.gdump << endl;
  out << tabs << info << endl;

  out << "ce_verbose       " << gOPT.verbose << endl;
  out << tabs << info << endl;

  out << "ce_cfile         " << gOPT.cfile << endl;
  out << tabs << info << endl;

  out.close();
}

void ResetUserOptions() {

  for(int n=0;n<NIFO_MAX;n++) {
    gOPT.eamp[n] = CE_eAMP;
    gOPT.ephs[n] = CE_ePHS;
    gOPT.etim[n] = CE_eTIM;
    gOPT.type[n] = CE_TYPE;
  }
  gOPT.cfile   = CE_CFILE;
  gOPT.seed    = CE_SEED;
  gOPT.gdump   = CE_GDUMP;
  gOPT.verbose = CE_VERBOSE;
}

void ReadUserOptions(TString options) {

  int n_type=0;
  int n_eamp=0;
  int n_ephs=0;
  int n_etim=0;
  if(options.CompareTo("")!=0) {
    cout << options << endl;
    if(!options.Contains("--")) {  // parameters are used only by cwb_inet

      TObjArray* token = TString(options).Tokenize(TString(' '));
        for(int j=0;j<token->GetEntries();j++){

        TObjString* tok = (TObjString*)token->At(j);
        TString stok = tok->GetString();

        if(stok.Contains("ce_eamp=")) {
          TString ce_eamp=stok;
          ce_eamp.Remove(0,ce_eamp.Last('=')+1);
          if(ce_eamp.IsFloat()) gOPT.eamp[n_eamp]=ce_eamp.Atof();
          if(n_eamp<(NIFO_MAX-1)) n_eamp++;
        }

        if(stok.Contains("ce_ephs=")) {
          TString ce_ephs=stok;
          ce_ephs.Remove(0,ce_ephs.Last('=')+1);
          if(ce_ephs.IsFloat()) gOPT.ephs[n_ephs]=ce_ephs.Atof();
          if(n_ephs<(NIFO_MAX-1)) n_ephs++;
        }

        if(stok.Contains("ce_etim=")) {
          TString ce_etim=stok;
          ce_etim.Remove(0,ce_etim.Last('=')+1);
          if(ce_etim.IsFloat()) gOPT.etim[n_etim]=ce_etim.Atof();
          if(n_etim<(NIFO_MAX-1)) n_etim++;
        }

        if(stok.Contains("ce_type=")) {
          TString type=stok;
          type.Remove(0,type.Last('=')+1);
          gOPT.type[n_type]=type;
          if(n_type<(NIFO_MAX-1)) n_type++;
        }

        if(stok.Contains("ce_cfile=")) {
          TString cfile=stok;
          cfile.Remove(0,cfile.Last('=')+1);
          gOPT.cfile=cfile;
        }

        if(stok.Contains("ce_seed=")) {
          TString ce_seed=stok;
          ce_seed.Remove(0,ce_seed.Last('=')+1);
          if(ce_seed.IsDigit()) gOPT.seed=ce_seed.Atoi();
        }

        if(stok.Contains("ce_gdump=")) {
          TString gdump=stok;
          gdump.Remove(0,gdump.Last('=')+1);
          gOPT.gdump=gdump;
        }

        if(stok.Contains("ce_verbose=")) {
          TString verbose=stok;
          verbose.Remove(0,verbose.Last('=')+1);
          gOPT.verbose=verbose;
        }
      }
    }
  }
}

//______________________________________________________________________________
void
SetCalibrationErrors(wavearray<double>& x, TString ifo, double time) {
//
// apply amplitude/phase calibration errors
//
//
// Input: x       - wavearray which contains the waveform data
//        ifo     - eg: L1/H1/V1
//        time    - event gps time 
//

  // get from the calibration file the nearest time entry 
  int simulation_id=0;
  time = GetNearestCalibrationEntry(time, simulation_id);

  // search begin,end of non zero data
  int ibeg=0; int iend=0;
  for(int i=0;i<(int)x.size();i++) {
    if(x[i]!=0 && ibeg==0) ibeg=i;
    if(x[i]!=0) iend=i;
  }
  int ilen=iend-ibeg+1;
  // create temporary array for FFTW & add scratch buffer + tShift
  int isize = 2*ilen;
  isize = isize + (isize%4 ? 4 - isize%4 : 0); // force to be multiple of 4
  wavearray<double> w(isize);
  w.rate(x.rate()); w=0;
  // copy x data !=0 in the middle of w array & set x=0
  for(int i=0;i<ilen;i++) {w[i+isize/4]=x[ibeg+i];x[ibeg+i]=0;}

  gMDC->CalibrateInspiral(w, ifo, time, simulation_id); // apply calibration errors

  // copy shifted data to input x array
  for(int i=0;i<(int)w.size();i++) {
    int j=ibeg-isize/4+i;
    if((j>=0)&&(j<(int)x.size())) x[j]=w[i];
  }

  return;
}

double GetNearestCalibrationEntry(double time, int& simulation_id) {

  ifstream in;
  in.open(gOPT.cfile.Data(),ios::in);
  if(!in.good()) {
    cout << endl << "CWB_Plugin_CE - "
         << "Error Opening Calibration File : " << gOPT.cfile.Data() << endl;
    exit(1);
  }

  std::map<TString, int> hsample;               // header sample

  // get header line
  char hline[4*1024];
  in.getline(hline,4*1024);

  TObjArray* token = TString(hline).Tokenize('\t');
  TString first_sheader="";
  for(int j=0;j<token->GetEntries();j++) {
     TString sheader = ((TObjString*)token->At(j))->GetString();
     sheader.ReplaceAll(" " ,"");
     // the first entry is used to check is an item is present. hsample["whatever itime not in the list"] returns 0
     // we move the first item to the end of the list
     if(j==0) {
       hsample[""] = j;
       first_sheader=sheader;
     } else {
       hsample[sheader] = j;
     }
  }
  hsample[first_sheader] = token->GetEntries();

  // extract parameters at GPS=time with simulation_id value
  bool spcal = false;
  int hsize=token->GetEntries();
  std::vector<double> sample(hsample.size());
  double min_diff_time=1e20;
  double nearest_time=time;
  while (1) {
    for(int i=0;i<hsize;i++) in >> sample[i];
    sample[hsize]=sample[0];    // we move the first item to the end of the list (see comment above)
    if(!in.good()) break;
    if(fabs(time-sample[hsample["time"]])<min_diff_time) {
      min_diff_time=fabs(time-sample[hsample["time"]]);
      nearest_time=sample[hsample["time"]];
      simulation_id=sample[hsample["simulation_id"]];
    } 
  }
  in.close();

  return nearest_time;
}

//______________________________________________________________________________
void
SetCalibrationErrors(wavearray<double>& x, int ifoid) {
//
// apply amplitude/phase calibration errors
//
//
// Input: x       - wavearray which contains the waveform data
//        ifod    - eg: L1/H1/V1 -> 0/1/2
//

  if(gOPT.verbose) cout << endl;

  // select amplitude calibration error eamp
  double eamp=1.;
  if(gOPT.eamp[ifoid]!=0) {
    if(gOPT.type[ifoid]=="fixed")    eamp = 1.+gOPT.eamp[ifoid];
    if(gOPT.type[ifoid]=="uniform")  eamp = gRandom->Uniform(1-fabs(gOPT.eamp[ifoid]),1+fabs(gOPT.eamp[ifoid]));
    if(gOPT.type[ifoid]=="gaussian") eamp = gRandom->Gaus(1,fabs(gOPT.eamp[ifoid]));
    cout.precision(3);
    if(gOPT.verbose) cout << "CWB_Plugin_CE : selected amplitude error (" << gOPT.eamp[ifoid] <<"%) -> eamp : " << eamp << endl;
  }

  // select phase calibration error ephs
  double ephs=0.;
  if(gOPT.ephs[ifoid]!=0) {
    if(gOPT.type[ifoid]=="fixed")    ephs = gOPT.ephs[ifoid];
    if(gOPT.type[ifoid]=="uniform")  ephs = gRandom->Uniform(-fabs(gOPT.ephs[ifoid]),fabs(gOPT.ephs[ifoid]));
    if(gOPT.type[ifoid]=="gaussian") ephs = gRandom->Gaus(0,fabs(gOPT.ephs[ifoid]));
    cout.precision(3);
    if(gOPT.verbose) cout << "CWB_Plugin_CE : selected phase error (" << gOPT.ephs[ifoid] <<" deg) -> ephs : " << ephs << " deg" << endl;
  }

  // select time calibration error etim
  double etim=0.;
  if(gOPT.etim[ifoid]!=0) {
    if(gOPT.type[ifoid]=="fixed")    etim = gOPT.etim[ifoid];
    if(gOPT.type[ifoid]=="uniform")  etim = gRandom->Uniform(-fabs(gOPT.etim[ifoid]),fabs(gOPT.etim[ifoid]));
    if(gOPT.type[ifoid]=="gaussian") etim = gRandom->Gaus(0,fabs(gOPT.etim[ifoid]));
    cout.precision(5);
    if(gOPT.verbose) cout << "CWB_Plugin_CE : selected time error (" << gOPT.etim[ifoid] <<" sec) -> etim : " << etim << " sec" << endl;
  }

  if(gOPT.verbose) cout << endl;

  // Apply Amplitude Calibration Error
  if(eamp!=1) x *= eamp;

  // Apply Phase Calibration Error
  if(ephs!=0) CWB::mdc::PhaseShift(x,ephs);

  // Apply Time Calibration Error
  if(etim!=0) CWB::mdc::TimeShift(x,etim);

  return;
}
