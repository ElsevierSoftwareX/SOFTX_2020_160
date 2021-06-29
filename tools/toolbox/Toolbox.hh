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
 * Package:      Toolbox Class Library
 * File name:    Toolbox.hh
 * Author:       Gabriele Vedovato (vedovato@lnl.infn.it)
 **********************************************************/


#ifndef TOOLBOX_HH
#define TOOLBOX_HH

#include "TROOT.h"
#include "TMath.h"
#include "TSystem.h"
#include "TString.h"
#include "TTree.h"
#include "TFile.h"
#include "TRandom.h"
#include "TRandom3.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "TSystemDirectory.h"
#include "TChain.h"
#include "TApplication.h"
#include "TH1I.h"
#include "TF1.h"
#include "TH1D.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TGraphSmooth.h"
#include "TMacro.h"

#include <string>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <vector>
#include <string.h>

#include "wavecomplex.hh"
#include "wavearray.hh"
#include "network.hh"
#include "Meyer.hh"
#include "netevent.hh"

#include "FrameL.h"

#include "History.hh"
#include "Toolfun.hh"

#define LST_TREE_NAME "frl"

enum CWB_CAT {
  CWB_CAT0  = 0,
  CWB_CAT1  = 1,
  CWB_CAT2  = 2,
  CWB_CAT3  = 3,
  CWB_HVETO = 10,
  CWB_PEM   = 11,
  CWB_EXC   = 12,
  CWB_USER  = 13
};

struct dqfile {
  char ifo[32];
  char file[1024];
  CWB_CAT cat;
  double shift;
  bool invert;
  bool c4;
};

struct frfile {
  int start;
  int stop;
  int length;
  vector<TString> file;
};

struct slag {
  int jobId;		// job id : sequential progressive number 
  vector<int> slagId;	// slag id vector : [0]=jobId - [1]=1/0 1=header slag - [2,..,nIFO+1] ifo slag 
  vector<int> segId;	// seg id vector : [0,..,nIFO-1] ifo segment number
};

struct mdcshift {
  double startMDC;
  double stopMDC;
  double offset;
};

struct ifoparms {
  TString name;
  double latitude;
  double longitude;
  double elevation;
  double AltX;       // elevation of the x arm
  double AzX;        // azimut of the x arm  (angle-deg from nord)
  double AltY;       // elevation of the y arm
  double AzY;        // azimut of the y arm (angle-deg from nord)
};

static vector<double> DEFAULT_DOUBLE_VECTOR;

using namespace std;

namespace CWB {

class Toolbox {

public:
  
  /* ************************ */
  /* * segment methods      * */
  /* ************************ */

  static vector<waveSegment> 
	 readSegments(TString ifile); 
  static vector<waveSegment> 
	 unionSegments(vector<waveSegment>& ilist); 
  static vector<waveSegment> 
	 sortSegments(vector<waveSegment>& ilist); 
  static vector<waveSegment> 
	 invertSegments(vector<waveSegment>& ilist); 
  static vector<waveSegment> 
         mergeSegLists(vector<waveSegment>& ilist1, vector<waveSegment>& ilist2);
  static waveSegment 
         getMaxSeg(vector<waveSegment> list);
  static double 
         getTimeSegList(vector<waveSegment> list);
  static int 
         dumpSegList(vector<waveSegment> list, TString fName, bool c4=false);
  static vector<waveSegment> 
         getSlagJobList(vector<waveSegment> ilist, int seglen=600);

  /* ************************ */
  /* * dq methods           * */
  /* ************************ */

  static vector<waveSegment> 
         readSegList(dqfile DQF);
  static vector<waveSegment> 
         readSegList(int nDQF, dqfile* DQF, CWB_CAT dqcat);
  static void 
         setSlagShifts(slag SLAG, vector<TString> ifos, double segLen, int nDQF, dqfile* DQF);

  /* ************************ */
  /* * job methods          * */
  /* ************************ */

  static vector<waveSegment> 
         getJobList(vector<waveSegment> ilist, double segLen=600., 
         double segMLS=300., double segEdge=8.);
  static vector<waveSegment> 
         getJobList(vector<waveSegment> cat1List, vector<waveSegment> cat2List,
         double segLen=600., double segMLS=300., double segTHR=30., double segEdge=8.);
  static void 
         dumpJobList(vector<waveSegment> ilist, TString fName, double segLen=600., 
         double segMLS=300., double segEdge=8.);
  static vector<waveSegment> 
         getSegList(int jobId, int nIFO, double segLen, double segMLS, double segEdge,
         vector<waveSegment> dqList); 
  static double 
         getLiveTime(vector<waveSegment>& jobList, vector<waveSegment>& dqList);
  static double 
         getLiveTime(vector<waveSegment>& jobList, vector<waveSegment>& dqList,
         vector<double> shiftList);

  /* ************************ */
  /* * slag methods         * */
  /* ************************ */

  static vector<slag> 
         getSlagList(size_t  nIFO,   size_t slagSize, int slagSegs, 
         int slagOff, size_t slagMin, size_t slagMax, 
         size_t* slagSite, char* slagFile=NULL);
  static void 
         dumpSlagList(vector<slag> slagList, TString slagFile, bool slagOnly=false);
  static slag 
         getSlag(vector<slag> slagList, int jobid); 
  static vector<slag> 
         getSlagList(vector<slag> slagList, vector<TString> ifos,  double segLen, 
         double segMin, double segEdge, int nDQF, dqfile* iDQF, CWB_CAT dqcat);
  static void 
         getSlagList(vector<slag>& slagList, int slagSize, 
         int slagRank, int nifo, vector<int>& id);
  static vector<waveSegment> 
         getSegList(slag SLAG, vector<waveSegment> jobList, double segLen, 
         double segMLS, double segEdge, vector<waveSegment> dqList); 

  /* ************************ */
  /* * MDC methods          * */
  /* ************************ */

  static double 
         getMDCShift(mdcshift mshift, double time);
  static int 
         shiftBurstMDCLog(std::vector<std::string>& mdcList, 
         vector<TString> ifos, double mdc_shift);
  static TString 
         SetMDCLog(TString log, int pos, TString val);
  static TString 
         SetMDCLog(TString log, int pos, double val);
  static TString 
         GetMDCLog(TString log, int pos); 
  static int 
         GetMDCLogSize(TString log); 

  /* ************************ */
  /* * condor methods       * */
  /* ************************ */

  static int 
         createDagFile(vector<waveSegment> jobList, TString condor_dir, TString label, 
         int jobmin=0, int jobmax=0); 
  static int 
         createDagFile(vector<slag> slagList, TString condor_dir, TString label, 
         int jobmin=0, int jobmax=0);  
  static int 
         createDagFile(vector<int> jobList, TString condor_dir, TString label, 
         vector<TString> jobFiles, TString stage="CWB_STAGE_FULL", int jobmin=0, int jobmax=0); 
  static int 
         createDagFile(vector<waveSegment> jobList, TString condor_dir, TString label, 
         vector<TString> jobFiles, TString stage="CWB_STAGE_FULL", int jobmin=0, int jobmax=0); 
  static int 
         createDagFile(vector<slag> slagList, TString condor_dir, TString label, 
         vector<TString> jobFiles, TString stage="CWB_STAGE_FULL", int jobmin=0, int jobmax=0);  
  static int 
         createSubFile(TString label, TString condor_dir, TString out_dir, TString err_dir, 
         TString log_dir, TString ext="", TString condor_tag=""); 
  static vector<int> 
         getCondorJobList(TString condor_dir, TString label); 
  static vector<float> 
         getJobBenchmark(TString ifName, int stageID, TString bench); 
  static vector<float> 
         getJobBenchmark(TString ifName, int stageID, int resID, int factorID, TString bench); 
  static TString
         DAG2LSF(char* dagFile, char* data_label, char* nodedir, char* data_dir,
         char* condor_dir, char* log_dir, char* output_dir, char* work_dir);

  /* ************************ */
  /* * merge methods        * */
  /* ************************ */

  static vector<TString> 
         mergeCWBTrees(TString dir_name, bool simulation, TString odir, 
         TString label, bool brms=false, bool bvar=false, bool bpsm=false);
  static void 
         mergeCWBTrees(vector<TString> fileList, bool simulation, TString odir, 
         TString label, bool brms=false, bool bvar=false, bool bpsm=false);
  static vector<TString> 
         mergeCWBTrees(int nthreads, TString dir_name, bool simulation, TString odir, 
         TString label, bool brms=false, bool bvar=false, bool bpsm=false);
  static void 
         mergeCWBTrees(int nthreads, vector<TString> fileList, bool simulation, TString odir, 
         TString label, bool brms=false, bool bvar=false, bool bpsm=false);
  static void 
         mergeTrees(vector<TString> fileList, TString treeName, TString odir, TString ofName, bool bhistory);
  static vector<TString> 
         readFileList(TString ifName);
  static void 
         dumpFileList(vector<TString> fileList, TString ofName);
  static vector<int> 
         getMergeJobList(TString merge_dir, TString label, int version); 
  static vector<int> 
         getMergeJobList(TString ifname, vector<TString>& jobFileList); 
  static vector<int> 
         getMergeJobList(TString ifname); 

  /* ************************ */
  /* * history methods      * */
  /* ************************ */

  static char* 
         readFile(TString ifName);
  static char* 
         getEnvCWB();

  /* ************************ */
  /* * spectra methods      * */
  /* ************************ */

  static void 
         makeSpectrum(wavearray<double>& psd, wavearray<double> x, double chuncklen=8, 
         double scratchlen=0, bool oneside=true);
  static void 
         makeSpectrum(TString ofname, wavearray<double> x, double chuncklen=8, 
         double scratchlen=0, bool oneside=true);
  static void 
         getSimNoise(wavearray<double>& u, TString fName, int seed, int run);
  static void 
         convertSampleRate(wavearray<double> iw, wavearray<double> ow);
  static void 
         convertSampleRate(wavearray<double> ix, wavearray<double> iy, 
         wavearray<double> ox, wavearray<double>& oy);
  static wavearray<double> 
         GetDetectorPSD(TString fName, double fWidth=8192., double dFreq=1.);

  /* ************************ */
  /* * system methods       * */
  /* ************************ */

  static vector<TString> 
         getFileListFromDir(TString dir_name, TString endString="", 
         TString beginString="", TString containString="",bool fast=false);
  static std::map<int,TString> 
         getJobFileMapFromDir(TString dir_name, TString endString="", 
         TString beginString="", TString containString="",bool fast=false);
  static bool 
         isFileExisting(TString fName);
  static bool 
         checkFile(TString fName, bool question=false, TString message="");
  static void 
         mkDir(TString dir, bool question=false, bool remove=true);
  static bool 
         rmDir(TString dir, bool question=true);
  static void 
         PrintProcInfo(TString str="");
  static TString 
         getTemporaryFileName(TString label="temporary", TString extension="tmp");
  static int 
         mksTemp(char *fTemplate) {return mkstemp(fTemplate);}  // wrapper to system funcion (used in CINT)
  static vector<TString> 
         sortStrings(vector<TString> ilist);

  /* ******************************** */
  /* * post processing methods      * */
  /* ******************************** */

  static int 
         setVeto(TString ifName, TString idir, TString odir, vector<TString> ifos, 
         int nVDQF, dqfile* VDQF, int nDQF, dqfile* DQF, double segLen, 
         double segMLS, double segEdge);
  static int 
         setCuts(TString ifName, TString idir, TString odir, 
                 TString trname, TString cuts, TString olabel);
  static int 
         setIFAR(TString ifName, TString idir, TString odir, TString trname, 
                 TString sels, TString farFile, int irho, TString olabel, bool inclusive=true);
  static int 
         setChunk(TString ifName, TString idir, TString odir, TString trname, int chunk);
  static void 
         setUniqueEvents(TString ifwave, TString ofwave, int nIFO, int pp_irho);
  static int 
         CombineCBC(vector<TString> ifwave, vector<TString> ifmdc, TString ofwave, TString ofmdc, 
                    int nIFO, float fthr, TString msearch, TString infos="", int lag=0, int slag=0, float ifarthr=0.);
  static double 
         getLiveTime(int nIFO, TChain& liv, wavearray<double>& Trun, 
         wavearray<double>* Wlag, wavearray<double>* Wslag, wavearray<double>& Tlag, 
         wavearray<double>& Tdlag, int lag_number=-1, int slag_number=-1, int dummy=0);
  static double 
         getZeroLiveTime(int nIFO, TChain& liv, int dummy=0);
  static void 
         doPoissonPlot(int nIFO, wavearray<double>* Wlag, 
         wavearray<double> Tlag, wavearray<double> Rlag, TString odir);
  static bool 
         isLeafInTree(TTree* itree, TString leaf);
  static int 
         setMultiplicity(TString ifName, TString idir, TString odir, int nIFO, double dTime);
  static wavecomplex 
         getRate(double rho, double Tgap, int nIFO, TChain& wav, wavearray<int>& Wsel,
         wavearray<double>* Wlag, wavearray<double>* Wslag, wavearray<double> Tlag);
  static int GetStepFunction(TString fName, vector<double>&  x, vector<double>&  y,
         vector<double>& ex = DEFAULT_DOUBLE_VECTOR, vector<double>& ey = DEFAULT_DOUBLE_VECTOR);
  static double GetStepFunction(TString option, TString fName, double V=0,
         vector<double>&  x = DEFAULT_DOUBLE_VECTOR, vector<double>&  y = DEFAULT_DOUBLE_VECTOR,
         vector<double>& ex = DEFAULT_DOUBLE_VECTOR, vector<double>& ey = DEFAULT_DOUBLE_VECTOR);
  static wavearray<double> GetPhase(wavearray<double> hi, wavearray<double> hr, wavearray<double>& fi, 
                                    wavearray<double>& fr, wavearray<double>& s, wavearray<double>& tt);

  /* ************************ */
  /* * utility methods      * */
  /* ************************ */

  static void 
         resampleToPowerOfTwo(wavearray<double>& w);
  static int  
         getJobId(TString file, TString fext="root");    
  static TString 
         getParameter(TString options, TString param="");   
  static TString 
         getFileName(FILE* fp);      
  static TString 
         getFileName(char* symlink); 
  static void	
         getUniqueFileList(TString ifile, TString ofile); 
  static TString 
         WriteFrameFile(wavearray<double> x, TString chName, TString frName, 
         TString frLabel, TString frDir="");
  static bool 
         question(TString question);
  static TString
         addCWBFlags(TMacro macro, TString ofname="");

  /* ************************************************** */
  /* * Hilbert & Wigner-Ville Transforms methods      * */
  /* ************************************************** */

  static wavearray<double> 
         getHilbertTransform(wavearray<double> x);		
  static wavearray<double> 
         getHilbertEnvelope(wavearray<double> x);		
  static wavearray<double> 
         getHilbertIFrequency(wavearray<double> x);		
  static wavearray<double>  
         getWignerVilleTransform(wavearray<double> x);	
  static void 
         unWrapPhase(wavearray<float>& p);			
  static void 
         getSineFittingParams(double a, double b, double c, double rate,
         double& amplitude, double& omega, double& phase);

  /* ************************************************** */
  /* * Color blind set methods    			      * */
  /* ************************************************** */

  static Int_t
  	 getTableau10BlindColor( Int_t index );
  static Int_t
  	 getTableau10BlindColor( TString name );
  static void
  	 getTableau10BlindColorPalette( const int nSteps, Int_t *palette );


private:

  static void 
         blackmanharris (double* window, int n);

  ClassDef(Toolbox,1)
};  

} // end namespace

#endif
