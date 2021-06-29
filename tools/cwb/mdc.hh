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


/**********************************************************
 * Package:      mdc Class Library
 * File name:    mdc.hh
 * Author:       Gabriele Vedovato (vedovato@lnl.infn.it)
 **********************************************************/


#ifndef MDC_HH
#define MDC_HH

#include "TObjString.h"
#include "TObjArray.h"
#include "TString.h"
#include "TROOT.h"
#include "TMath.h"
#include "TSystem.h"
#include "TLatex.h"
#include "TTree.h"
#include "TFile.h"
#include "TNamed.h"
#include "TGlobal.h"
#include "TRandom.h"
#include "TComplex.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TColor.h"
#include "TF1.h"
#include "TF3.h"
#include "TMacro.h"
#include "Math/VectorUtil.h"

#include <string>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <vector>

#include "wavecomplex.hh"
#include "wavearray.hh"
#include "watplot.hh"
#include "skymap.hh"
#include "network.hh"
#include "detector.hh"
#include "wat.hh"

#include "History.hh"
#include "Toolbox.hh"
#include "config.hh"
#include "time.hh"
#include "gskymap.hh"
#include "STFT.hh"

#include "FrameL.h"

#ifndef __CINT__
#ifdef _USE_LAL

#define LAL_USE_OLD_COMPLEX_STRUCTS
#define restrict  // this is to getrid of restrict which is defined only in c99

#ifdef _USE_ROOT6
// this instruction is necessary in ROOT6 because it seems to be defined
// by ROOT or other packages and it forbid to LAL to include Random.h 
#undef _RANDOM_H        
#endif

#include <lal/LIGOLwXMLInspiralRead.h>
#include <lal/LALConfig.h>
#include <lal/LALStdio.h>
#include <lal/LALStdlib.h>
#include <lal/LALError.h>
#include <lal/LALDatatypes.h>
#include <lal/LIGOMetadataUtils.h>
#include <lal/LIGOMetadataTables.h>
#include <lal/AVFactories.h>
#include <lal/NRWaveIO.h>
#include <lal/NRWaveInject.h>
//#include <lal/Inject.h>	// obsolete
#include <lal/FileIO.h>
#include <lal/Units.h>
#include <lal/FrequencySeries.h>
#include <lal/TimeSeries.h>
#include <lal/TimeFreqFFT.h>
#include <lal/VectorOps.h>
#include <lal/LALDetectors.h>
#include <lal/FindChirp.h>
#include <lal/Random.h>
#include <lal/LALNoiseModels.h>
#include <lal/Date.h>
#include <lal/LALStatusMacros.h>
#include <lal/LALSimulation.h>
#include <lal/LALInspiral.h>
#include <lal/LALSimInspiralWaveformFlags.h>
#include <lal/LALSimIMR.h>
#if LAL_VERSION_MAJOR >   6 || (LAL_VERSION_MAJOR ==  6 && \
   (LAL_VERSION_MINOR >  14 || (LAL_VERSION_MINOR == 14 && \
    LAL_VERSION_MICRO >=  0                             ))) 	// LAL_VERSION >= 6.14.0
#include <lal/SnglBurstUtils.h>
#else
#include <lal/LIGOMetadataBurstUtils.h>
#endif
#include <lal/LIGOLwXMLBurstRead.h>
#include <lal/GenerateBurst.h>
#include <lal/LIGOLwXML.h>
#if LAL_VERSION_MAJOR >   6 || (LAL_VERSION_MAJOR ==  6 && \
   (LAL_VERSION_MINOR >  19 || (LAL_VERSION_MINOR == 19 && \
    LAL_VERSION_MICRO >=  2                             ))) 	// LAL_VERSION >= 6.19.2
#include <lal/LIGOLwXMLlegacy.h>
#endif

#endif
#endif

#ifdef _USE_EBBH
#include "eBBH.hh"
#endif

#define MDC_INJ_RATE       0.01   // sec^-1
#define MDC_INJ_JITTER     10.    // sec
#define MDC_INJ_LENGTH     1.     // sec
#define MDC_INJ_HRSS       2.5e-21

#define MDC_SAMPLE_RATE    16384. // sample/sec

#define EPZOOM		   0.9999999

enum MDC_TYPE {
  MDC_SG            = 0,
  MDC_SGL           = 0,
  MDC_SGC           = 1,
  MDC_SGE           = 2,
  MDC_RD            = 3,
  MDC_RDL           = 3,
  MDC_RDC           = 4,
  MDC_RDE           = 5,
  MDC_WNB           = 6,
  MDC_GA            = 7,
  MDC_EBBH          = 8,
  MDC_GA_LAL        = 9,
  MDC_SGE_LAL       = 10,
  MDC_WNB_LAL       = 11,
  MDC_CG            = 12,
  MDC_CGL           = 12,
  MDC_CGC           = 13,
  MDC_CGE           = 14,
  MDC_SC_LAL        = 15,
  MDC_USER          = 16	  // must be declared as the last one
};

enum MDC_COORDINATES {
  MDC_EARTH         = 0,
  MDC_CELESTIAL     = 1
};

enum MDC_DISTRIBUTION {
  MDC_RANDOM        = -1000,
  MDC_EARTH_FIX     = 0,
  MDC_CELESTIAL_FIX = 1,
  MDC_MNGD          = 2,
  MDC_GWGC          = 3,
  MDC_LOGFILE       = 4,
  MDC_CUSTOM        = 5,
  MDC_XMLFILE       = 6
};

inline const char* DistributionToString(MDC_DISTRIBUTION n) {
  switch (n) {
    case MDC_RANDOM:   		return "MDC_RANDOM";
    case MDC_EARTH_FIX:   	return "MDC_EARTH_FIX";
    case MDC_CELESTIAL_FIX:   	return "MDC_CELESTIAL_FIX";
    case MDC_MNGD:   		return "MDC_MNGD";
    case MDC_GWGC:   		return "MDC_GWGC";
    case MDC_LOGFILE:   	return "MDC_LOGFILE";
    case MDC_CUSTOM:   		return "MDC_CUSTOM";
    case MDC_XMLFILE:   	return "MDC_XMLFILE";
    default:      		return "MDC_RANDOM";
  }
}

enum MDC_DRAW {
  MDC_TIME          = 0,
  MDC_FFT           = 1,
  MDC_TF            = 2
};

struct mdcid {
  TString  name;     // name of waveform
  int      ID;       // major id of waveform list
  int      id;       // minor id of waveform list (Ex : the list WNB with different random waveforms)
};

struct mdcpar {
  TString  name;
  double   value;
  TString  svalue;
};

struct waveform {             // waveform structure
  MDC_TYPE          type;     // mdc type
  bool              status;   // status flag used for internal mdc processing
  TString           name;     // name of waveform
  TString           hpPath;   // path of user define hp file
  TString           hxPath;   // path of user define hx file
  vector<mdcpar>    par;      // waveform parameters 
  wavearray<double> hp;       // waveform hp component
  wavearray<double> hx;       // waveform hx component;
  vector<waveform>  list;     // list of waveforms belonging to the same waveform class
};

struct source {
  double   gps;               // waveform source structure 
  double   theta;             // latitude     (degrees)
  double   phi;               // longitude    (degrees)
  double   psi;               // polarization (degrees)
  double   rho;               // distance
  double   iota;              // elliptical inclination angle (degrees)
  double   hrss;              // sqrt(hp^2+hx^2)
  int      ID;                // major id of waveform list
  int      id;                // minor id of waveform list (Ex : the list WNB with different random waveforms)
  waveform wf;                // waveform
};

static vector<double> DEFAULT_VECtOR_DOUBLE;

namespace CWB {

class mdc : public TNamed {

public:

  mdc(); 
  mdc(int nIFO, TString* ifo);
  mdc(int nIFO, detector** pD);
  mdc(network* net);
  mdc(const CWB::mdc& value);
  ~mdc();  

  // operators
  mdc& operator = (const mdc&);

  TString Get(wavearray<double>& x, TString ifo);

  mdcid AddWaveform(MDC_TYPE mdc_type, vector<mdcpar> par, TString uname="");    // built-in waveform
  void AddWaveform(TString mdc_name, TString hp_fName, TString hx_fName, 
                   double srate, vector<mdcpar> par=vector<mdcpar>());
  void AddWaveform(TString mdc_name, TString hp_fName, 
                   double srate, vector<mdcpar> par=vector<mdcpar>());
  void AddWaveform(TString mdc_name, TString hp_fName, TString hx_fName);
  void AddWaveform(TString mdc_name, TString hp_fName);
  mdcid AddWaveform(waveform wf);
#ifdef _USE_LAL
#ifndef __CINT__
  // LAL SimBurst waveforms
  mdcid AddWaveform(MDC_TYPE mdc_type, SimBurst* sim_burst, vector<mdcpar> par, TString uname=""); 
#endif
#endif

  // read ascii waveform with fixed sample rate (1 column [h])
  void ReadWaveform(wavearray<double>& x, TString fName, double srate);
  // read ascii waveform with variable sample rate (2 columns [t,h])
  void ReadWaveform(wavearray<double>& x, TString fName);

  // get size of wfList
  inline size_t wfListSize() { return wfList.size(); }

  void SetCoordinatesSystem(MDC_COORDINATES mdc_coordinates=MDC_EARTH) {this->mdc_coordinates=mdc_coordinates;}
  MDC_COORDINATES GetCoordinatesSystem() {return mdc_coordinates;}

  void SetSkyDistribution(MDC_DISTRIBUTION sky_distribution, vector<mdcpar> par, int seed=0, bool add=false);
  void SetSkyDistribution(MDC_DISTRIBUTION sky_distribution, TString fName, vector<mdcpar> par, int seed=0, bool add=false);
  void SetSkyDistribution(MDC_DISTRIBUTION sky_distribution, TString fName, int seed, bool add=false);
  MDC_DISTRIBUTION GetSkyDistribution() {return sky_distribution;}
  TString GetSkyFile() {return sky_file;} 
  vector<mdcpar> GetSkyParms() {return sky_parms;}
  void DrawSkyDistribution(TString name = "skymap", TString projection = "", 
                           TString coordinate = "Geographic", double resolution = 2, bool background=true);

  void GetSourceCoordinates(double& theta, double& phi, double& psi, double& rho, 
                            double& iota, double& hrss, int& ID, int& id);
  void GetSourceCoordinates(double gps, double& theta, double& phi, double& psi, double& rho, 
                            double& iota, double& hrss, int& ID, int& id);

  void   SetInjLength(double inj_length=MDC_INJ_LENGTH) {
         this->inj_length = inj_length>0 ? inj_length : MDC_INJ_LENGTH;}
  double GetInjLength() {return inj_length;}
  void   SetInjHrss(double inj_hrss=MDC_INJ_HRSS) { 
         // if inj_hrss=0 the hrss @10Kpc is the one which is defined by hp,hc waveforms
         this->inj_hrss = inj_hrss>=0 ? inj_hrss : MDC_INJ_HRSS;}
  double GetInjHrss() {return inj_hrss;}
  void   SetInjRate(double inj_rate=MDC_INJ_RATE) {
         this->inj_rate = inj_rate>0 ? inj_rate : MDC_INJ_RATE;}
  double GetInjRate() {return inj_rate;}
  void   SetInjOffset(double inj_offset=0) {
         this->inj_offset = inj_offset>0 ? inj_offset : 0;}
  double GetInjOffset() {return inj_offset;}
  void   SetInjJitter(double inj_jitter=MDC_INJ_JITTER) {
         this->inj_jitter = inj_jitter>=0 ? inj_jitter : MDC_INJ_JITTER;}
  double GetInjJitter() {return inj_jitter;}
  void   SetSourceListSeed(unsigned int srcList_seed) {
         this->srcList_seed = srcList_seed;}
  double GetSourceListSeed() {return srcList_seed;}
  double GetSampleRate() {return MDC_SAMPLE_RATE;}

#ifdef _USE_LAL
  void SetInspiral(TString inspName, TString inspOptions="");
  TString GetInspiral() {return inspOptions;}
  wavearray<double> GetInspiral(TString pol, int gps_start_time, int gps_end_time);
  TString GetInspiralOption(TString inspOption);	
#endif

  void Dump(TString fname, int ID, int id, TString polarization);
  void Dump(TString fname, TString name, int id, TString polarization);
  void Dump(TString fname, wavearray<double>& x);
  void Dump(TString dname);    // dump all waveforms in dname dir

  watplot* Draw(TString name, int id=0, TString polarization="hp", MDC_DRAW type=MDC_TIME, 
                TString options = "ALP", Color_t color=kBlack);
  watplot* Draw(int ID, int id=0, TString polarization="hp", MDC_DRAW type=MDC_TIME, 
                TString options = "ALP", Color_t color=kBlack);
  watplot* Draw(wavearray<double>& x, MDC_DRAW type=MDC_TIME, 
                TString options = "ALP", Color_t color=kBlack);
  watplot* Draw(TString ifo, double gpsStart, double gpsEnd, int id,
                MDC_DRAW type=MDC_TIME, TString options = "ALP", Color_t color=kBlack);

  void DumpLog(TString fName, TString label="", bool append=false);
  void DumpLogHeader(TString fName, TString label="", int size=0);

  waveform GetWaveform(int ID, int id=0);
  waveform GetWaveform(TString name, int id=0);
  int      GetWaveformID(TString name);
  void     GetWaveform(waveform& wf); 

  static double GetCentralTime(wavearray<double> x);
  static double GetCentralTime(waveform wf);
  static double GetCentralFrequency(wavearray<double> x);
  static double GetCentralFrequency(waveform wf);
  static double GetTimeRange(wavearray<double> x, double& tMin, double& tMax, double efraction = EPZOOM);
  static void   TimeShift(wavearray<double>& x, double tShift=0.);
  static void   PhaseShift(wavearray<double>& x, double pShift=0.);

  double GetAntennaPattern(TString ifo, double phi, double theta, double psi=0., TString polarization="hp");
  double GetDelay(TString ifo1, TString ifo2, double phi, double theta);

  void Print(int level=0);  // *MENU*
  virtual void Browse(TBrowser *b) {Print();}

  wavearray<double> GetWNB(double frequency, double bandwidth, double duration, int seed=0, bool mode=0);
  wavearray<double>  GetRD(double frequency, double tau, double iota, bool polarization=0);
  wavearray<double> GetSGQ(double frequency, double Q);
  wavearray<double> GetCGQ(double frequency, double Q);
  wavearray<double>  GetGA(double duration);

  static void AddGauss(wavearray<double> &td, double v, double u=0.);
  static void AddExp(wavearray<double> &td, double v, int M);
  static void AddSGBurst(wavearray<double> &td, double a, double f, double s, double d=0.);
  static void AddCGBurst(wavearray<double> &td, double a, double f, double s, double d=0.);
  static void AddWGNoise(wavearray<double> &td, double a, double s);

  TString WriteFrameFile(TString frDir, TString frLabel, size_t gps, 
                         size_t length=1000, bool log=false, vector<TString> chName=vector<TString>());

  network* GetNetwork() {return net;}

  vector<waveform> GetWaveformList() {return wfList;}

  vector<waveform> wfList;

  std::vector<std::string>  mdcList;   // list of injections
  std::vector<std::string>  mdcType;   // list of injection types
  std::vector<std::string>  xmlType;   // list of xml injection types (used with MDC_XMLFILE)
  std::vector<double>       mdcTime;   // gps time of selected injections
  std::vector<std::string>  mdcName;   // list of injection names used with MDC_LOG
  std::vector<source>       srcList;   // list of source MDC parameters

  double  GetPar(TString name, vector<mdcpar> par, bool& error); 
  TString GetParString(TString name, vector<mdcpar> par, bool& error); 

  waveform GetSourceWaveform(int& ID, int& id);
  vector<source> GetSourceList(double start, double stop);
  TString GetBurstLog(source src, double FrameGPS, double SimHpHp, double SimHcHc, double SimHpHc);

  watplot* DrawTime(wavearray<double>& x, TString options = "ALP", Color_t color=kBlack);
  watplot* DrawFFT(wavearray<double>& x, TString options = "ALP", Color_t color=kBlack);
  void     DrawTF(wavearray<double>& x, TString options = "");

  TString GetBurst(wavearray<double>& x, TString ifo);

#ifdef _USE_LAL
  TString GetInspiral(wavearray<double>& x, TString ifo);
  TString GenInspiralXML(int gps_start_time, int gps_end_time, bool rmFile=false);
  TString GetInspName() {return inspName;}
  TString GetWaveName() {return waveName;}
#ifndef __CINT__
  TString GetInspiralLog(TString inspName, double FrameGPS, SimInspiralTable *thisInj);
  void FindChirpInjectSignals(wavearray<double>& w, SimInspiralTable *injections, TString ifo);
#endif
  void CreateBurstXML(TString fName, vector<mdcpar> xml_parms);
  void FillBurstXML(bool precision=false);
  void CloseBurstXML();
  static TString GetBurstNameLAL(TString options);
  int XLALInspiralTDWaveformFromSimInspiral(REAL8TimeSeries **hplus, REAL8TimeSeries **hcross, SimInspiralTable *thisRow, REAL8 deltaT);
  static double SimIMRSEOBNRv4ROMTimeOfFrequency(double freq, double m1, double m2, double chi1, double chi2);
  static double SimIMRSEOBNRv4ROMFrequencyOfTime(double time, double m1, double m2, double chi1, double chi2);
#endif
  static double e2cosi(double e);
  static double cosi2e(double cosi);

  TString GetTemporaryFileName(TString tag="mdc", TString ext="txt", TString dir="/tmp", bool mkdir=false);

  void SetZoom(double epzoom=EPZOOM) {if(epzoom<0||epzoom>1) epzoom=EPZOOM; this->epzoom=epzoom;}
  double GetZoom() {return epzoom;}

  gskymap*   GetGSkyMap() {return psp;}
  CWB::STFT* GetSTFT()    {return stft;}
  watplot*   GetWatPlot() {return pts;}

#ifdef _USE_LAL
  // posterioris sample functions
  static void Posterior2XML(TString sampleFile, TString xmlFile, TString options="");
  static std::map<TString, double> GetPsample(std::map<TString, int> hsample, vector<double> sample, 
                                              float if_lower=-1, float if_ref=-1);
  // calibration functions
  void CalibrateInspiral(wavearray<double>& w, TString ifo, double time, int simulation_id);
  void SetInspiralCLB(TString inspCLB) {
    TString extension = inspCLB(inspCLB.Last('.'),inspCLB.Sizeof()-inspCLB.Last('.')-1);
    if(extension!=".clb") {
      cout << "CWB::mdc::SetInspiralCLB - Error : Calibration file extension must be .clb" << endl;
      exit(1);
    }
    this->inspCLB=inspCLB;
  }
  TString GetInspiralCLB() {return this->inspCLB;}
#endif

  // wavearray methods
  static double 
	 PhaseSync(wavearray<double>& w1, wavearray<double>& w2, double& sync_phase);
  static double 
	 TimeSync(wavearray<double>& w1, wavearray<double>& w2, double& sync_time);
  static double 
	 TimePhaseSync(wavearray<double>& w1, wavearray<double>& w2, double& sync_time, double& sync_phase);
  static void 
	 PhaseSync(wavearray<double>& w, double sync_phase);
  static void 
	 TimeSync(wavearray<double>& w, double sync_time);
  static void 
	 TimePhaseSync(wavearray<double>& w, double sync_time, double sync_phase);
  static wavearray<double> 
	 GetAligned(wavearray<double>* w1, wavearray<double>* w2);
  static wavearray<double>
	 GetAdd(wavearray<double>* w1, wavearray<double>* w2);
  static wavearray<double>
	 GetDiff(wavearray<double>* w1, wavearray<double>* w2);
  static int 
	 Align(wavearray<double>& w1, wavearray<double>& w2);
  static double 
	 GetTimeBoundaries(wavearray<double> x, double P, double& bT, double& eT, double T=-1., double Q=-1.);
  static double
         GetFrequencyBoundaries(wavearray<double> x, double P, double& bF, double& eF);
  static wavearray<double> 
	 GetXCorr(wavearray<double>& w1, wavearray<double>& w2);
  static double
         GetMatchFactor(TString match, vector<wavearray<double> >& w1, vector<wavearray<double> >& w2,
                        vector<double> tstart=DEFAULT_VECtOR_DOUBLE, vector<double> tstop=DEFAULT_VECtOR_DOUBLE);
  static wavearray<double>
         GetBandpass(wavearray<double> x, double bF, double eF);
  static wavearray<double>
         GetEnvelope(wavearray<double>* x);
  static wavearray<double>
         GetSpectrum(wavearray<double>* x, bool oneside=false);


private:

  void Init(int seed=0);
  void exit(int err) {gSystem->Exit(err);}	// overwrite exit system function

  network* net;     //!
  injection* inj;   //!
  TTree* inj_tree;  //! 

  MDC_COORDINATES  mdc_coordinates;
  MDC_DISTRIBUTION sky_distribution;
  TString sky_file;
  TString xml_filename;	
  vector<mdcpar> sky_parms;

  std::vector<std::string>  nameList;
  std::vector<float>        thList;
  std::vector<float>        phList;
  std::vector<float>        psiList;
  std::vector<float>        rhoList;
  std::vector<float>        iotaList;
  std::vector<double>       hrssList;
  std::vector<double>       gpsList;
  std::vector<int>          IDList;
  std::vector<int>          idList;

  double inj_rate;
  double inj_offset;
  double inj_jitter;
  double inj_hrss;
  double inj_length;  		//length sec
  unsigned int srcList_seed;	// value added to the seed used in GetSourceList (default is 0)

  // inspiral parameters
  TString inspCLB;
  TString inspXML;
  TString inspDIR;
  TString inspName;
  TString inspOptions;
  TString waveName;     //- PN approximant to be used in computing the waveform
                        //  (see enum Approximant in LALSimInspiral.h)

  skymap sm;

  double epzoom;	// used to zoom plots : energy percentage displayed

  gskymap*      psp;  //!
  CWB::STFT*    stft; //!
  watplot*      pts;  //!

#ifdef _USE_LAL
#ifndef __CINT__
  // used by the CreateBurstXML,FillBurstXML,CloseBurstXML functions
  ProcessTable 		*xml_process_table_head;		//!
  ProcessTable 		*xml_process;				//!
  ProcessParamsTable 	*xml_process_params_table_head;		//!
  SearchSummaryTable 	*xml_search_summary_table_head;    	//!
  SearchSummaryTable 	*xml_search_summary;			//!
  TimeSlide 		*xml_time_slide_table_head;		//!
  SimBurst 		*xml_sim_burst_table_head;		//!
  SimBurst 		**xml_sim_burst;			//!
#endif
#endif

  ClassDef(mdc,9)
};  

} // end namespace

#endif
