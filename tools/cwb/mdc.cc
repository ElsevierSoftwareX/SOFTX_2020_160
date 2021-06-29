/*
# Copyright (C) 2019 Gabriele Vedovato, Francesco Salemi
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

//#define NEW_CODE_FOR_SEOBNRv4P_AND_SEOBNRv4PHM	// uncomment this line to enable the SEOBNRv4P_AND_SEOBNRv4PHM
                                                // corrently available only in the following branch:
						// https://git.ligo.org/serguei.ossokine/lalsuite/tree/SEOBNRv4PHM-review-rebase-Sep2019

#include "mdc.hh"
#include "constants.hh"
#include "TSpline.h"

#ifdef _USE_LAL

#include <lal/Interpolate.h>

// LAL Stuff
#define FAILMSG( stat, func, file, line, id )                                  \
  do {                                                                         \
    if ( lalDebugLevel & LALERROR )                                            \
    {                                                                          \
      LALPrintError( "Error[0]: file %s, line %d, %s\n"                        \
          "\tLAL_CALL: Function call `%s' failed.\n", file, line, id, func );  \
    }                                                                          \
  } while( 0 )

#define LAL_CALL( function, statusptr ) \
  ((function),lal_errhandler(statusptr,#function,__FILE__,__LINE__,"$Id$"))

LALStatus blank_status;

typedef int ( *lal_errhandler_t )(
    LALStatus  *,
    const char *func,
    const char *file,
    const int   line,
    volatile const char *id
    );

int LAL_ERR_ABRT(
    LALStatus  *stat,
    const char *func,
    const char *file,
    const int   line,
    volatile const char *id
    )
{
  if ( stat->statusCode )
  {
    FAILMSG( stat, func, file, line, id );
    abort();
  }
  return 0;
}

int LAL_ERR_EXIT(
    LALStatus  *stat,
    const char *func,
    const char *file,
    const int   line,
    volatile const char *id
    )
{
  if ( stat->statusCode )
  {
    FAILMSG( stat, func, file, line, id );
    exit( 1 );
  }
  return stat->statusCode;
}

//#define LAL_ERR_DFLT LAL_ERR_ABRT
lal_errhandler_t lal_errhandler = LAL_ERR_ABRT;

// see LAL binj.c
static 
ProcessParamsTable **add_process_param(ProcessParamsTable **proc_param, const ProcessTable *process, 
                                       const char *type, const char *param, const char *value) {
  *proc_param = XLALCreateProcessParamsTableRow(process);
  snprintf((*proc_param)->program, sizeof((*proc_param)->program), "%s", "cWB");
  snprintf((*proc_param)->type, sizeof((*proc_param)->type), "%s", type);
  snprintf((*proc_param)->param, sizeof((*proc_param)->param), "--%s", param);
  snprintf((*proc_param)->value, sizeof((*proc_param)->value), "%s", value);
  return &(*proc_param)->next;
}

#define ADD_PROCESS_PARAM(proc_param, process, type, param, value)				\
  do { proc_param = add_process_param(proc_param, process, type, param, value); } while(0)

#endif

////////////////////////////////////////////////////////////////////////////////
/* BEGIN_HTML
<p>The mdc class is designed to easily produce mdc data for simulation

Overview:
<ol style="list-style-type: upper-roman;">
  <li><a href="#init">Init</a></li>
  <li><a href="#conf">Configuration</a>
  <ol>
     <li><a href="#conf:setparms">Set Injection Params</a></li>
     <li><a href="#conf:waveform">Add Waveforms</a></li>
     <li><a href="#conf:skydistribution">Set SkyDistribution</a></li>
  </ol>
  </li>
  <li><a href="#getdata">Get Data</a></li>
  <li><a href="#writeframe">Write Frame</a></li>
</ol>

<h3><a name="init">I. Init</a></h3>
Object mdc must be instantiate setting the network detectors:
<br><br>
<pre>
    root[] char wf_name[256];
    root[] waveform wf;
    root[] vector<mdcpar> par;
    root[] int nIFO = 3;               // number of detectors
    root[] TString ifo[3] = {"L1","H1","V1"}; 
    root[] mdc MDC(nIFO,ifo);          // create a mdc object
</pre>
<h3><a name="conf">II. Configuration</a></h3>
In the configuration phase must be define injection parameters, waveforms and sky distribution:
<br>
<h3><a name="conf:setparms">II.1 Set Injection Params</a></h3>
<br>
<pre>
    root[] // --------------------------------------------------------
    root[] // define injection parameters
    root[] // --------------------------------------------------------
    root[] MDC.SetInjHrss(2.5e-21);
    root[] MDC.SetInjRate(0.0333333);
    root[] MDC.SetInjJitter(10.0);
</pre>
<h3><a name="conf:setparms">II.2 Add Waveforms</a></h3>
<br>
<pre>
    root[] // --------------------------------------------------------
    root[] // Add SinGaussian waveform 
    root[] // --------------------------------------------------------
    root[] par.resize(2);
    root[] par[0].name="frequency"; par[0].value=235.;
    root[] par[1].name="Q"; par[1].value=8.9;
    root[] MDC.AddWaveform(MDC_SG, par);
    root[] MDC.Print();   // list defined waveforms
</pre>
<h3><a name="conf:setparms">II.3 Set SkyDistribution</a></h3>
<br>
<pre>
    root[] // --------------------------------------------------------
    root[] // define sky distribution
    root[] // --------------------------------------------------------
    root[] par.resize(3);
    root[] par[0].name="entries";par[0].value=100000;    // pool of events
    root[] par[1].name="rho_min";par[1].value=1;         // min rho // Kpc
    root[] par[2].name="rho_max";par[2].value=100;       // max rho // Kpc
    root[] MDC.SetSkyDistribution(MDC_RANDOM,par,seed);
</pre>
<h3><a name="getdata">III. Get Data</a></h3>
<br>
<pre>
    root[] // --------------------------------------------------------
    root[] // get data 
    root[] // --------------------------------------------------------
    root[] waveform x;
    root[] x.start(0);x.rate(16384);x.size(600*x.rate());
    root[] MDC.Get(x,ifo[0]);
</pre>
<h3><a name="writeframe">IV. Write Frame</a></h3>
<br>
<pre>
    root[] // --------------------------------------------------------
    root[] // write frame 
    root[] // --------------------------------------------------------
    root[] TString frDir   = "frames";
    root[] TString frLabel = "TEST";
    root[] size_t gps      = 123456789;
    root[] size_t length   = 1000;  // sec
    root[] bool log        = true;  
    root[] TString ofName = MDC.WriteFrameFile(frDir, frLabel, gps, length, log);
</pre>

</p>

END_HTML */
////////////////////////////////////////////////////////////////////////////////

ClassImp(CWB::mdc)

//______________________________________________________________________________
CWB::mdc::mdc() {
//
// mdc default constructor 
//

  Init();
  net=NULL;
}

//______________________________________________________________________________
CWB::mdc::mdc(int nIFO, TString* ifo) {
//
// mdc constructor 
//
// Input: nIFO        - number of detectors
//        ifo         - array of detector names       
//

  Init();
  net = new network;

  for(int n=0;n<nIFO;n++) net->add(new detector(const_cast<char*>(ifo[n].Data())));
}

//______________________________________________________________________________
CWB::mdc::mdc(int nIFO, detector** pD) {
//
// mdc constructor 
//
// Input: nIFO        - number of detectors
//        pD          - array of detector objects       
//

  Init();
  net = new network;

  for(int n=0; n<nIFO; n++) {
    detector* _pD = NULL;
    if(pD[n]->isBuiltin()) _pD = new detector(pD[n]->Name);
    else                   _pD = new detector(pD[n]->getDetectorParams());
    this->net->add(_pD);
  }
}

//______________________________________________________________________________
CWB::mdc::mdc(network* net) {
//
// mdc constructor 
//
// Input: net         - network object 
//

  Init();

  if(net!=NULL) { 
    
    this->net = new network;

    int nIFO=net->ifoListSize();
    for(int n=0; n<nIFO; n++) {
      detector* pD = NULL;
      if(net->getifo(n)->isBuiltin()) pD = new detector(net->getifo(n)->Name);
      else                            pD = new detector(net->getifo(n)->getDetectorParams());
      POLARIZATION polarization = net->getifo(n)->getPolarization(); 
      pD->setPolarization(polarization);
      this->net->add(pD);
    }
  } else {
    this->net = NULL;
  }
}

//______________________________________________________________________________
CWB::mdc::mdc(const CWB::mdc& value) {
//
// mdc copy constructor 
//

   *this = value;
}

//______________________________________________________________________________
CWB::mdc::~mdc() {
//
// mdc destructor 
//

  if(net!=NULL)  {
    for(int n=0;n<(int)net->ifoListSize();n++) delete net->getifo(n);
    delete net;
  }
  if(inj!=NULL)  delete inj;
  if(sky_distribution!=MDC_LOGFILE) 
    if(inj_tree!=NULL)  delete inj_tree;

  if(psp!=NULL)  delete psp;
  if(stft!=NULL) delete stft;
  if(pts!=NULL)  delete pts;
}

//______________________________________________________________________________
CWB::mdc& 
CWB::mdc::operator=(const CWB::mdc& value) {
//
// mdc copy constructor 
//

  if(value.net!=NULL) {
     if(net==NULL) net = new network;
     else {for(int n=0;n<(int)net->ifoListSize();n++) delete net->getifo(n);}
     *(net)=*(value.net);
     net->ifoList.resize(0);
     int nIFO=value.net->ifoListSize();
     for(int n=0; n<nIFO; n++) {
       detector* pD = value.net->getifo(n);
       if(value.net->getifo(n)->isBuiltin()) net->add(new detector(pD->Name));
       else                                  net->add(new detector(pD->getDetectorParams()));
       POLARIZATION polarization = value.net->getifo(n)->getPolarization(); 
       net->getifo(n)->setPolarization(polarization);
     }
  }

  if(value.inj!=NULL) {
     int nIFO=net->ifoListSize();
     if(inj!=NULL) delete inj;
     inj = new injection(nIFO);
     inj_tree = inj->setTree();
     if(value.inj_tree!=NULL) {
       inj->output(inj_tree,net,1,false);
       delete inj;
       inj = new injection(inj_tree,nIFO);
     }
  }

//  if(psp!=NULL)  delete psp;  psp=NULL; 
//  if(stft!=NULL) delete stft; stft=NULL;
//  if(pts!=NULL)  delete pts;  pts=NULL;

  mdc_coordinates  = value.mdc_coordinates;
  sky_distribution = value.sky_distribution;

  inspCLB         = value.inspCLB;
  inspXML         = value.inspXML;
  inspDIR         = value.inspDIR;
  inspName        = value.inspName;
  inspOptions     = value.inspOptions;
  waveName        = value.waveName;

  inj_rate        = value.inj_rate;
  inj_offset      = value.inj_offset;
  inj_jitter      = value.inj_jitter;
  inj_hrss        = value.inj_hrss;
  srcList_seed    = value.srcList_seed;

   wfList.clear(); this->wfList  = value.wfList;
  mdcList.clear(); this->mdcList = value.mdcList;
  mdcType.clear(); this->mdcType = value.mdcType;
  xmlType.clear(); this->xmlType = value.xmlType;
  mdcTime.clear(); this->mdcTime = value.mdcTime;
  mdcName.clear(); this->mdcName = value.mdcName;
  srcList.clear(); this->srcList = value.srcList;
 nameList.clear(); this->nameList  = value.nameList;
   thList.clear(); this->thList  = value.thList;
   phList.clear(); this->phList  = value.phList;
  psiList.clear(); this->psiList = value.psiList;
  rhoList.clear(); this->rhoList = value.rhoList;
 iotaList.clear(); this->iotaList = value.iotaList;
 hrssList.clear(); this->hrssList = value.hrssList;
  gpsList.clear(); this->gpsList = value.gpsList;
   IDList.clear(); this->IDList = value.IDList;
   idList.clear(); this->idList = value.idList;

  return *this;
}

//______________________________________________________________________________
void
CWB::mdc::Init(int seed) {
//
// mdc init 
//
// Input:  seed - seed for random generation
//

  inj=NULL;
  inj_tree=NULL;
  psp=NULL;
  stft=NULL;
  pts=NULL;
  inspCLB="";
  inspXML="";
  inspDIR="";
  inspName="";
  inspOptions="";
  waveName="";
  xml_filename="";

  epzoom = EPZOOM;

  sky_distribution = MDC_RANDOM;
  inj_jitter   = MDC_INJ_JITTER;
  inj_rate     = MDC_INJ_RATE;
  inj_offset   = 0.;
  inj_hrss     = MDC_INJ_HRSS;
  inj_length   = MDC_INJ_LENGTH;
  srcList_seed = 0;
  gRandom->SetSeed(seed);

   wfList.clear();
  mdcList.clear();
  mdcType.clear();
  xmlType.clear();
  mdcTime.clear();
  mdcName.clear();
  srcList.clear();
 nameList.clear();
   thList.clear();
   phList.clear();
  psiList.clear();
  rhoList.clear();
 iotaList.clear();
 hrssList.clear();
  gpsList.clear();
   IDList.clear();
   idList.clear();
}

//______________________________________________________________________________
double 
CWB::mdc::GetPar(TString name, vector<mdcpar> par, bool& error) {
//
// Return Waveform parameter
//
//
// Input: name     - parameter name
//        par      - vector of parameters
//
// Return: parameter value 
//

  for(int i=0;i<(int)par.size();i++) {
    if(par[i].name.CompareTo(name)==0) return par[i].value;
  }
  error=true;  
  return 0.;
}

//______________________________________________________________________________
TString 
CWB::mdc::GetParString(TString name, vector<mdcpar> par, bool& error) {
//
// Return Waveform string parameter
//
//
// Input: name     - parameter name
//        par      - vector of parameters
//
// Return: parameter string 
//

  for(int i=0;i<(int)par.size();i++) {
    if(par[i].name.CompareTo(name)==0) return par[i].svalue;
  }
  error=true;  
  return "";
}

//______________________________________________________________________________
mdcid 
CWB::mdc::AddWaveform(MDC_TYPE mdc_type, vector<mdcpar> par, TString uname) {
//
// Add Waveform to this object
//
//
// Input: mdc_type - name of mdc  
//        par      - input mdc parameters
//
//                   mdc_type & par can be one of the following choices :  
//
//            MDC_SG/MDC_SGL = Linear SinGaussian
//                   MDC_SGC = Circular SinGaussian        
//                   MDC_SGE = Elliptical SinGaussian   // created as SGC (ellipticity must be defined by iota) 
//                             vector<mdcpar> par(2);
//  `                          par[0].name="frequency"; par[0].value=XXX; // Hz
//                             par[1].name="Q";         par[1].value=YYY;
//
//            MDC_CG/MDC_CGL = Linear CosGaussian
//                   MDC_CGC = Circular CosGaussian        
//                   MDC_CGE = Elliptical CosGaussian   // created as CGC (ellipticity must be defined by iota) 
//                             vector<mdcpar> par(2);
//  `                          par[0].name="frequency"; par[0].value=XXX; // Hz
//                             par[1].name="Q";         par[1].value=YYY;
// 
//            MDC_RD/MDC_RDL = Linear RingDown          
//                   MDC_RDC = Circular RingDown          
//                   MDC_RDE = Elliptical RingDown      // created as RGC (ellipticity must be defined by iota) 
//                             vector<mdcpar> par(2);
//  `                          par[0].name="frequency"; par[0].value=XXX; // Hz
//                             par[1].name="tau";       par[1].value=YYY;
// 
//               MDC_SGE_LAL = LAL SinGaussian (see definition in LALSimBurst.c)  
//                             Difference from MDC_SGE : LAL use CG instead of SG !!!    
//                             force eccentricity=0 --> circularly polarized 
//                             force polarization=0  
//                             NOTE : created as SGC (ellipticity must be defined by iota)
//                                    polarization=0 (polarization is defined by psi)
//                             vector<mdcpar> par(2);
//  `                          par[0].name="frequency"; par[0].value=XXX; // Hz
//                             par[1].name="Q";         par[1].value=YYY;
// 
//                   MDC_WNB = White Noise Burst          
//                             vector<mdcpar> par(x);
//                             par[0].name="frequency"; par[0].value=XXX; // Hz : initial frequency
//                             par[1].name="bandwidth"; par[1].value=YYY; // Hz 
//                             par[2].name="duration";  par[2].value=ZZZ; // sec
//                             par[3].name="pseed";     par[3].value=ppp; // (opt) seed fo hp component (def=0)
//                             par[4].name="xseed";     par[4].value=xxx; // (opt) seed fo hx component (def=0)
//                             par[5].name="mode";      par[5].value=mmm; // (opt) 0/1 a/symmetric respect central freq
//                                                                        // 0 = (default)
// 
//               MDC_WNB_LAL = LAL White Noise Burst (see definition in LALSimBurst.c)      
//                             Difference from MDC_WNB : LAL has a gaussian envelope also in frequency    
//                             NOtE : LAL seed = waveform_number = hseed*10000000000+lseed
//                             vector<mdcpar> par(5);
//                             par[0].name="frequency"; par[0].value=XXX; // Hz : central frequency 
//                             par[1].name="bandwidth"; par[1].value=YYY; // Hz
//                             par[2].name="duration";  par[2].value=ZZZ; // sec
//                             par[3].name="lseed";     par[3].value=ppp; //  low seed (max 10 int digits)
//                             par[4].name="hseed";     par[4].value=xxx; // high seed (max 10 int digits)
// 
//                   MDC_GA  = Gaussian Burst          
//                             vector<mdcpar> par(1);
//                             par[0].name="duration";  par[0].value=XXX; // sec
// 
//                MDC_GA_LAL = LAL Gaussian (see definition in LALSimBurst.c)  
//                             vector<mdcpar> par(1);
//                             par[0].name="duration";  par[0].value=XXX; // sec
// 
//                MDC_SC_LAL = LAL Generates cosmic string cusp waveforms (see definition in XLALGenerateStringCusp, LALSimBurst.c)  
//                             vector<mdcpar> par(1);
//                             par[0].name="frequency"; par[0].value=XXX; // High frequency cut-off (Hertz)
//                             par[1].name="amplitude"; par[1].value=YYY; // amplitude (strain)
// 
//                  MDC_EBBH = Eccentric Binary Black Holes // created as Circular (ellipticity must be defined by iota)
//                                                          // are generated with hrss @ distance GM/c^2 meters
//                                                          // custom hrss is disabled
//                             vector<mdcpar> par(1);
//                             par[0].name="file_name";     // ebbh list file name [.lst/.root] 
//                                                          // .lst  format : one event for each line 
//                                                          // -> m1 m2 rp0 e0 (hp,hx are generated 'On The Fly')
//                                                          // .root format : one event for each entry 
//                                                          // -> id m1 m2 rp0 e0 hp hx 
//                             par[1].name="tree cuts";     // optional : used .root file to select tree entries 
//
//                             or
//
//                             vector<mdcpar> par(1);
//                             par[0].name="ebbh parameters"; // ebbh parameters : id m1 m2 rp0 e0 dist(Kpc) redshift 
//
//        uname    - is the name of waveform, if uname="" then uses the builtin name (default) 
//                   not enable for MDC_EBBH
//
// Return: the name of waveform  
//
//
// Note : For the elliptically polarized MDC's, they were produced as follows:
//
//        h(t) = F_+ h_+(t) + F_x h_x(t)
//
//        where
// 
//        h_+ = 0.5*(1+cos(iota) * A(t) cos(Phi(t))
//        h_x = cos(iota) * A(t) sin(Phi(t))
//
//        See Patrick's document here: (https://dcc.ligo.org/LIGO-P1000041) 
//

  // check if mdc_type is in the MDC_TYPE list
  if(mdc_type<0 || mdc_type>=MDC_USER) {
     cout << "CWB::mdc::AddWaveform - Error : mdc type " << mdc_type << " not allowed  !!!" << endl;
     exit(1);
  }

  char wf_name[128];
  waveform wf;
  mdcid waveid = {"",0,0};

  // get number of decimals to be used to format the waveform name
  bool error=false;
  int decimals = (int)GetPar("decimals",par,error);
  if(error) decimals = -1;		// default format (like BurstMDC)

  if((mdc_type==MDC_RDE) ||
     (mdc_type==MDC_RDC) ||
     (mdc_type==MDC_RD)) {

    error=false;
    double frequency = GetPar("frequency",par,error);
    double tau       = GetPar("tau",par,error);
    double iota=0.;
    if(mdc_type==MDC_RDE) iota = 0.;	// iota can be customized with SetSkyDistribution
    if(mdc_type==MDC_RDC) iota = 0.;
    if(mdc_type==MDC_RD)  iota = 90.;

    if(error) {
      cout << "CWB::mdc::AddWaveform - "
           << "Error : num par must be at least 2 [frequency,tau]" << endl;
      exit(1);
    }

    int d1 = decimals==-1 ? 0 : decimals;
    int d2 = decimals==-1 ? 1 : decimals;

    if(mdc_type==MDC_RD)  sprintf(wf_name, "RD_%.*f_%.*f",d1,frequency,d2,tau);
    if(mdc_type==MDC_RDC) sprintf(wf_name,"RDC_%.*f_%.*f",d1,frequency,d2,tau);
    if(mdc_type==MDC_RDE) sprintf(wf_name,"RDE_%.*f_%.*f",d1,frequency,d2,tau);

    wf.type = mdc_type;
    wf.name = wf_name;
    wf.name.ReplaceAll(".","d");
    if(uname!="") wf.name = uname;
    wf.hp = GetRD(frequency,tau,iota,0);
    wf.hx = GetRD(frequency,tau,iota,1);
    wf.hpPath = wf.name;
    wf.hxPath = wf.name;
    wf.par=par;
    return AddWaveform(wf);

  } else
  if((mdc_type==MDC_SGC) ||
     (mdc_type==MDC_SGE) ||
     (mdc_type==MDC_SG)) {

    error=false;
    double frequency = GetPar("frequency",par,error);
    double Q         = GetPar("Q" ,par,error);

    if(error) {
      cout << "CWB::mdc::AddWaveform - Error : num par must be 2 [frequency,Q]" << endl;
      exit(1);
    }

    int d1 = decimals==-1 ? 0 : decimals;
    int d2 = decimals==-1 ? 1 : decimals;

    if(mdc_type==MDC_SG)  sprintf(wf_name, "SG%.*fQ%.*f",d1,frequency,d2,Q);
    if(mdc_type==MDC_SGC) sprintf(wf_name,"SGC%.*fQ%.*f",d1,frequency,d2,Q);
    if(mdc_type==MDC_SGE) sprintf(wf_name,"SGE%.*fQ%.*f",d1,frequency,d2,Q);

    wf.type = mdc_type;
    wf.name = wf_name;
    wf.name.ReplaceAll(".","d");
    if(decimals==-1) wf.name.ReplaceAll("d0","");
    if(uname!="") wf.name = uname;
    wf.hp = GetSGQ(frequency,Q);
    if(mdc_type==MDC_SG) {   // iota = 90
      wf.hx = wf.hp;
      wf.hx = 0;
      wf.hpPath = wf.name;
      wf.hxPath = "";
    } else {	             // iota = 0; for MDC_SGE iota can be customized with SetSkyDistribution
      wf.hx = GetCGQ(frequency,Q);
      wf.hpPath = wf.name;
      wf.hxPath = wf.name;
    }
    wf.par=par;
    return AddWaveform(wf);

  } else
  if((mdc_type==MDC_CGC) ||
     (mdc_type==MDC_CGE) ||
     (mdc_type==MDC_CG)) {

    error=false;
    double frequency = GetPar("frequency",par,error);
    double Q         = GetPar("Q" ,par,error);

    if(error) {
      cout << "CWB::mdc::AddWaveform - Error : num par must be 2 [frequency,Q]" << endl;
      exit(1);
    }

    int d1 = decimals==-1 ? 0 : decimals;
    int d2 = decimals==-1 ? 1 : decimals;

    if(mdc_type==MDC_CG)  sprintf(wf_name, "CG%.*fQ%.*f",d1,frequency,d2,Q);
    if(mdc_type==MDC_CGC) sprintf(wf_name,"CGC%.*fQ%.*f",d1,frequency,d2,Q);
    if(mdc_type==MDC_CGE) sprintf(wf_name,"CGE%.*fQ%.*f",d1,frequency,d2,Q);

    wf.type = mdc_type;
    wf.name = wf_name;
    wf.name.ReplaceAll(".","d");
    if(decimals==-1) wf.name.ReplaceAll("d0","");
    if(uname!="") wf.name = uname;
    wf.hp = GetCGQ(frequency,Q);
    if(mdc_type==MDC_CG) {   // iota = 90
      wf.hx = wf.hp;
      wf.hx = 0;
      wf.hpPath = wf.name;
      wf.hxPath = "";
    } else {	             // iota = 0; for MDC_CGE iota can be customized with SetSkyDistribution
      wf.hx = GetSGQ(frequency,Q);
      wf.hpPath = wf.name;
      wf.hxPath = wf.name;
    }
    wf.par=par;
    return AddWaveform(wf);

  } else
  if(mdc_type==MDC_WNB) {

    error=false;
    double frequency = GetPar("frequency",par,error);
    double bandwidth = GetPar("bandwidth",par,error);
    double duration  = GetPar("duration" ,par,error);

    if(error) {
      cout << "CWB::mdc::AddWaveform - "
           << "Error : WNB par must at least  "
           << "[frequency,bandwidth,duration]" << endl;
      exit(1);
    }

    error=false; int pseed = (int)GetPar("pseed",par,error); if(error) pseed=0;
    error=false; int xseed = (int)GetPar("xseed",par,error); if(error) xseed=0;
    error=false; bool mode = (bool)GetPar("mode",par,error); if(error) mode=0;

    int d1 = decimals==-1 ? 0 : decimals;
    int d2 = decimals==-1 ? 0 : decimals;
    int d3 = decimals==-1 ? 3 : decimals;

    sprintf(wf_name,"WNB%.*f_%.*f_%.*f",d1,frequency,d2,bandwidth,d3,duration);

    wf.type = mdc_type;
    wf.name = wf_name;
    wf.name.ReplaceAll(".","d");
    if(uname!="") wf.name = uname;
    wf.hp = GetWNB(frequency,bandwidth,duration,pseed,mode);
    wf.hpPath = wf.name;
    wf.hx = GetWNB(frequency,bandwidth,duration,xseed,mode);
    wf.hxPath = wf.name;
    wf.par=par;
    return AddWaveform(wf);
  
  } else
  if(mdc_type==MDC_GA) {

    error=false;
    double duration  = GetPar("duration" ,par,error);

    if(error) {
      cout << "CWB::mdc::AddWaveform - "
           << "Error : num par must at least 1 "
           << "[duration]" << endl;
      exit(1);
    }

    int d1 = decimals==-1 ? 1 : decimals;

    sprintf(wf_name,"GA%.*f",d1,1000.*duration);

    wf.type = mdc_type;
    wf.name = wf_name;
    wf.name.ReplaceAll(".","d");
    if(uname!="") wf.name = uname;
    wf.hp = GetGA(duration);
    wf.hpPath = wf.name;
    wf.hx = wf.hp;
    wf.hx = 0;
    wf.hxPath = "";
    wf.par=par;
    return AddWaveform(wf);

  } else
  if(mdc_type==MDC_GA_LAL) {	// use LAL GA waveform 

#if LAL_VERSION_MAJOR >   6 || (LAL_VERSION_MAJOR ==  6 && \
   (LAL_VERSION_MINOR >  15 || (LAL_VERSION_MINOR == 15 && \
    LAL_VERSION_MICRO >=  0                             )))     // LAL_VERSION >= 6.15.0
#ifdef _USE_LAL
    error=false;
    double duration = GetPar("duration",par,error);

    if(error) {
      cout << "CWB::mdc::AddWaveform - "
           << "Error : wrong input LAL GA parameters : "
           << "[duration,decimals(opt)]" << endl;
      exit(1);
    }

    // create & populate SimBurst structure
    SimBurst *sim_burst = XLALCreateSimBurst();
    strcpy(sim_burst->waveform, "Gaussian");
    sim_burst->duration          = duration;
    mdcpar mpar={"normalization",1};par.push_back(mpar);	// force normalization
    waveid = AddWaveform(MDC_GA_LAL, sim_burst, par, uname);
    XLALDestroySimBurst(sim_burst); 
    return waveid;
#else
    cout << "CWB::mdc::AddWaveform - "
         << "Error : MDC_GA_LAL is enabled only with LAL" << endl;
    exit(1);
#endif
#else
    cout << "CWB::mdc::AddWaveform - MDC_GA_LAL can not be used with LAL ver < 6.15.0" << endl; 
    exit(1);
#endif

  } else
  if(mdc_type==MDC_SGE_LAL) {	// use LAL SG waveform 

#ifdef _USE_LAL
    error=false;
    double frequency = GetPar("frequency",par,error);
    double Q         = GetPar("Q",par,error);

    if(error) {
      cout << "CWB::mdc::AddWaveform - "
           << "Error : wrong input LAL SG parameters : "
           << "[frequency,Q,decimals(opt)]" << endl;
      exit(1);
    }

    // create & populate SimBurst structure
    SimBurst *sim_burst = XLALCreateSimBurst();
    strcpy(sim_burst->waveform, "SineGaussian");
    sim_burst->frequency         = frequency;
    sim_burst->q                 = Q;
    sim_burst->hrss              = 1;
    // force eccentricity=0 --> circularly polarized (it is a SG parameter)
    // the ellipticity is applied elsewhere using iota
    sim_burst->pol_ellipse_e=0.;
    // force polarization=0 (it is a SG parameter)
    // the polarization is applied elsewhere using psi
    sim_burst->pol_ellipse_angle=0.;
    mdcpar mpar={"normalization",1};par.push_back(mpar);	// force normalization
    waveid = AddWaveform(MDC_SGE_LAL, sim_burst, par, uname);
    XLALDestroySimBurst(sim_burst); 
    return waveid;
#else
    cout << "CWB::mdc::AddWaveform - "
         << "Error : MDC_SGE_LAL is enabled only with LAL" << endl;
    exit(1);
#endif

  } else
  if(mdc_type==MDC_WNB_LAL) {	// use LAL WNB waveform != from cWB WNB

#ifdef _USE_LAL
    error=false;
    double frequency =      GetPar("frequency",par,error);
    double bandwidth =      GetPar("bandwidth",par,error);
    double duration  =      GetPar("duration" ,par,error);
    // LAL waveform_number = hseed*10000000000+lseed
    int    lseed     = (int)GetPar("lseed"    ,par,error);
    int    hseed     = (int)GetPar("hseed"    ,par,error);

    if(error) {
      cout << "CWB::mdc::AddWaveform - "
           << "Error : wrong input LAL WNB parameters : "
           << "[frequency,bandwidth,duration,lseed,hseed,decimals(opt)]" << endl;
      exit(1);
    }

    // create & populate SimBurst structure
    SimBurst *sim_burst = XLALCreateSimBurst();
    strcpy(sim_burst->waveform, "BTLWNB");
    sim_burst->frequency         = frequency;
    sim_burst->duration          = duration;
    sim_burst->bandwidth         = bandwidth;
    sim_burst->egw_over_rsquared = 1;
    sim_burst->waveform_number   = hseed*10000000000+lseed;
    mdcpar mpar={"normalization",1};par.push_back(mpar);	// force normalization
    waveid = AddWaveform(MDC_WNB_LAL, sim_burst, par, uname);
    XLALDestroySimBurst(sim_burst); 
    return waveid;
#else
    cout << "CWB::mdc::AddWaveform - "
         << "Error : MDC_WNB_LAL is enabled only with LAL" << endl;
    exit(1);
#endif
  
  } else
  if(mdc_type==MDC_SC_LAL) {	// use LAL cosmic string cusp waveform

#if LAL_VERSION_MAJOR >   6 || (LAL_VERSION_MAJOR ==  6 && \
   (LAL_VERSION_MINOR >  15 || (LAL_VERSION_MINOR == 15 && \
    LAL_VERSION_MICRO >=  0                             )))     // LAL_VERSION >= 6.15.0
#ifdef _USE_LAL
    error=false;
    double frequency = GetPar("frequency",par,error);
    double amplitude = GetPar("amplitude",par,error);

    if(error) {
      cout << "CWB::mdc::AddWaveform - "
           << "Error : wrong input LAL SC parameters : "
           << "[frequency,amplitude,decimals(opt)]" << endl;
      exit(1);
    }

    // create & populate SimBurst structure
    SimBurst *sim_burst = XLALCreateSimBurst();
    strcpy(sim_burst->waveform, "StringCusp");
    sim_burst->frequency = frequency;
    if(amplitude>0) {
      sim_burst->amplitude = amplitude;
    } else {
      sim_burst->amplitude = 1;
      mdcpar mpar={"normalization",1};par.push_back(mpar);	// force normalization
    }
    waveid = AddWaveform(MDC_SC_LAL, sim_burst, par, uname);
    XLALDestroySimBurst(sim_burst); 
    return waveid;
#else
    cout << "CWB::mdc::AddWaveform - "
         << "Error : MDC_SC_LAL is enabled only with LAL" << endl;
    exit(1);
#endif
#else
    cout << "CWB::mdc::AddWaveform - MDC_SC_LAL can not be used with LAL ver < 6.15.0" << endl; 
    exit(1);
#endif

  } else
  if(mdc_type==MDC_EBBH) {

#ifdef _USE_EBBH

    if(par.size()<1) {
      cout << "CWB::mdc::AddWaveform - "
           << "Error : num par must at least 1 "
           << "[file name (.lst/.root) : list of eBBH mdc]" << endl; 
      exit(1);
    }

    TString fName = par[0].name;	// get list file name of eBBH mdc

    if(fName.EndsWith(".lst")) {	// read lst file

      ifstream in;
      in.open(fName,ios::in);
      if (!in.good()) {
        cout << "CWB::mdc::AddWaveform - Error Opening File : " << fName << endl;
        exit(1);
      }

      // get number of entries
      int entries=0;
      char str[1024];
      while(true) {
        in.getline(str,1024);
        if (!in.good()) break;
        if(str[0] != '#') entries++;
      }
      cout << "entries " << entries << endl;
      in.clear(ios::goodbit);
      in.seekg(0, ios::beg);
  
      double m1,m2,rp0,e0,dist,redshift;
      int id;
      int fpos=0;
      while (1) {
        fpos=in.tellg();
        in.getline(str,1024);
        if(str[0] == '#') continue;
        if (!in.good()) break;

        dist=0; 
        redshift=0;
        std::stringstream linestream(str);
        if(!(linestream >> id >> m1 >> m2 >> rp0 >> e0 >> dist >> redshift)) {
          linestream.str(str); 
          linestream.clear();		// clear stringstream error status 
          if(!(linestream >> id >> m1 >> m2 >> rp0 >> e0 >> dist)) {
            linestream.str(str); 
            linestream.clear();		// clear stringstream error status 
            if(!(linestream >> id >> m1 >> m2 >> rp0 >> e0)) {
              cout << "CWB::mdc::AddWaveform - Wrong Format for File : " << fName << endl;
              cout << "input line : " << endl;
              cout << str << endl;
              cout << "must be : " << endl;
              cout << "event# " << "id " << " m1 " << " m2 " << " rp0 " << " e0 " << endl;
              cout << "or : " << endl;
              cout << "event# " << "id " << " m1 " << " m2 " << " rp0 " << " e0 " << " dist " << endl;
              cout << "or : " << endl;
              cout << "event# " << "id " << " m1 " << " m2 " << " rp0 " << " e0 " << " dist " << " redshift " << endl;
              exit(1);
            }
          }
        }

        wf.type = mdc_type;
        wf.name = "eBBH";
        wf.hpPath = "eBBH";
        wf.hxPath = "eBBH";
        wf.par.resize(1);
        wf.par[0].name=str;
        wf.par[0].value=MDC_EBBH;
        waveid = AddWaveform(wf);
      }

      in.close();
      return waveid;

    } else
    if(fName.EndsWith(".root")) {	// read root file

      TFile* efile = new TFile(fName);
      if(efile==NULL) {
        cout << "CWB::mdc::AddWaveform - Error opening root file : " << fName.Data() << endl;
        exit(1);
      }

      int id;
      double m1,m2,rp0,e0,dist,redshift;
      wavearray<double>* hp = new wavearray<double>;
      wavearray<double>* hx = new wavearray<double>;

      TTree* etree = (TTree *) efile->Get("ebbh");
      if(etree==NULL) {
        cout << "CWB::mdc::AddWaveform - file : " << fName.Data()  
             << " not contains tree ebbh" << endl;
        exit(1);
      }

      etree->SetBranchAddress("id",&id);
      etree->SetBranchAddress("m1",&m1);
      etree->SetBranchAddress("m2",&m2);
      etree->SetBranchAddress("rp0",&rp0);
      etree->SetBranchAddress("e0",&e0);
      etree->SetBranchAddress("hp",&hp);
      etree->SetBranchAddress("hx",&hx);
      int dstatus = (etree->GetBranch("dist")==NULL)     ? 0 : etree->SetBranchAddress("dist",&dist);
      int rstatus = (etree->GetBranch("redshift")==NULL) ? 0 : etree->SetBranchAddress("redshift",&redshift);

      TString ecut = "";
      if(par.size()==2) ecut = par[1].name;	// get tree selection cut of eBBH mdc
      etree->Draw("Entry$",ecut,"goff");
      double* entry = etree->GetV1();
      int esize = etree->GetSelectedRows();
      for(int i=0;i<esize;i++) {
        etree->GetEntry(entry[i]); 
        char str[256];
        if(rstatus)      sprintf(str,"%d %f %f %f %f %f %f",id,m1,m2,rp0,e0,dist,redshift);
        else if(dstatus) sprintf(str,"%d %f %f %f %f %f",id,m1,m2,rp0,e0,dist);
        else             sprintf(str,"%d %f %f %f %f",id,m1,m2,rp0,e0);
        
        //cout << id << " " << m1 << " " << m2 << " " << rp0 << " " << e0 << endl;

        wf.type = mdc_type;
        wf.name = "eBBH";
        wf.hpPath = "eBBH";
        wf.hxPath = "eBBH";
        wf.par.resize(2);
        wf.par[0].name=str;
        wf.par[0].value=MDC_EBBH;
        wf.par[1].name=fName;
        wf.par[1].value=entry[i];
        waveid = AddWaveform(wf);
      }
      delete hp;
      delete hx;
      delete efile;
      return waveid;

    } else {	// eBBH values are defined in the par[0].name

      char str[1024];
      sprintf(str,par[0].name.Data());
      double m1,m2,rp0,e0,dist,redshift;
      int id;
      dist=0; 
      redshift=0;
      std::stringstream linestream(str);
      if(!(linestream >> id >> m1 >> m2 >> rp0 >> e0 >> dist >> redshift)) {
        linestream.str(str); 
        linestream.clear();		// clear stringstream error status 
        if(!(linestream >> id >> m1 >> m2 >> rp0 >> e0 >> dist)) {
          linestream.str(str); 
          linestream.clear();		// clear stringstream error status 
          if(!(linestream >> id >> m1 >> m2 >> rp0 >> e0)) {
            cout << "CWB::mdc::AddWaveform - Wrong Input Parameter Format : " << str << endl;
            cout << "input line : " << endl;
            cout << str << endl;
            cout << "must be : " << endl;
            cout << "event# " << "id " << " m1 " << " m2 " << " rp0 " << " e0 " << endl;
            cout << "or : " << endl;
            cout << "event# " << "id " << " m1 " << " m2 " << " rp0 " << " e0 " << " dist " << endl;
            cout << "or : " << endl;
            cout << "event# " << "id " << " m1 " << " m2 " << " rp0 " << " e0 " << " dist " << " redshift " << endl;
            exit(1);
          }
        }
      }

      wf.type = mdc_type;
      wf.name = "eBBH";
      wf.hpPath = "eBBH";
      wf.hxPath = "eBBH";
      wf.par.resize(1);
      wf.par[0].name=str;
      wf.par[0].value=MDC_EBBH;
      waveid = AddWaveform(wf);
    }

#else
    cout << "CWB::mdc::AddWaveform - Error : MDC_EBBH not enabled !!!" << endl;
    exit(1);
#endif

    return waveid;
  }

  cout << "CWB::mdc::AddWaveform - Warning : waveform not added !!!" << endl;
  return waveid;
}

#ifdef _USE_LAL
//______________________________________________________________________________
mdcid 
CWB::mdc::AddWaveform(MDC_TYPE mdc_type, SimBurst* sim_burst, vector<mdcpar> par, TString uname) {
//
// Add LAL SimBurst Waveform to this object
//
//
// Input: mdc_type  - name of mdc  
//
//               MDC_SGE_LAL = LAL SinGaussian   
//                             Difference from MDC_SGE : LAL use CG instead of SG !!!
//                                                       Uses eccentricity instead of ellipticity !!!
// 
//               MDC_WNB_LAL = LAL White Noise Burst          
//                             Difference from MDC_WNB : LAL has a gaussian envelope also in frequency
// 
//               MDC_GA_LAL  = LAL Gaussian Burst          
//
//               MDC_SC_LAL  = LAL Cosmic String Cusp          
//
//      sim_burst  - The LAL SimBurst structure describes a burst injection 
//                   (see LIGOMetadataTables.h) 
//
//           par   - number of decimals used to format the waveform name
//                   if not defined ; used default format
//                   vector<mdcpar> par(1);
//                   par[0].name="decimals"; par[0].value=XXX; 
//                   par[1].name="normalization"; par[1].value=X;   // X= 0/1 = disabled/enabled
//
//        uname    - is the name of waveform, if uname="" then uses the builtin name (default) 
// 
// Return: the name of the waveform  
//

  // check if mdc_type is in the MDC_TYPE list
  if(mdc_type<0 || mdc_type>=MDC_USER) {
     cout << "CWB::mdc::AddWaveform - Error : mdc type " << mdc_type << " not allowed  !!!" << endl;
     exit(1);
  }
  if(sim_burst==NULL) {
    cout << endl << "CWB::mdc::AddWaveform - "
         << "Error in LAL AddWaveform : SimBurst is NULL" << endl;
    exit(1);
  }

  char wf_name[128];
  waveform wf;

  // get number of decimals to be used to format the waveform name
  bool error=false;
  int decimals = (int)GetPar("decimals",par,error);
  if(error) decimals = -1;		// default format (like BurstMDC)
  error=false;
  int normalization = (int)GetPar("normalization",par,error);
  if(error) normalization = 0;		// default : use normalization defined in sim_burst

  REAL8TimeSeries *hp=NULL;
  REAL8TimeSeries *hx=NULL;

  REAL8 deltaT = 1./MDC_SAMPLE_RATE;

  /* Get waveform from LAL : see GenerateBurst.c */
  int ret = XLALGenerateSimBurst(&hp, &hx, sim_burst, deltaT);
  if( ret==XLAL_FAILURE ) {
    cout << endl << "CWB::mdc::AddWaveform - "
         << "Error in LAL AddWaveform : check sim burst parameters" << endl;
    exit(1);
  }
  if(hp->data->length!=hx->data->length) {
     cout << "CWB::mdc::AddWaveform - Error : LAL hp,hx size not equal !!!" << endl;
     exit(1);
  }

  // create waveform arrays hp,hx
  wf.hp.resize(hp->data->length);
  wf.hp.rate(MDC_SAMPLE_RATE);
  wf.hx=wf.hp;

  double sum;
  // fill hp
  wf.hp=0; sum=0; 
  for(int i=0;i<wf.hp.size();i++) {wf.hp[i]=hp->data->data[i]; sum+=pow(wf.hp[i],2);}
  if(normalization&&sum>0) wf.hp *= sqrt(MDC_SAMPLE_RATE/sum/2.);	// normalization -> 1/sqrt(2)
  // fill hx
  wf.hx=0; sum=0; 
  for(int i=0;i<wf.hx.size();i++) {wf.hx[i]=hx->data->data[i]; sum+=pow(wf.hx[i],2);}
  if(normalization&&sum>0) wf.hx *= sqrt(MDC_SAMPLE_RATE/sum/2.);	// normalization -> 1/sqrt(2)

  XLALDestroyREAL8TimeSeries(hp);
  XLALDestroyREAL8TimeSeries(hx);

  if(mdc_type==MDC_SGE_LAL) {

    int d1 = decimals==-1 ? 0 : decimals;
    int d2 = decimals==-1 ? 1 : decimals;

    wf.par.resize(4);
    wf.par[0].name="frequency";    wf.par[0].value=sim_burst->frequency;
    wf.par[1].name="Q";            wf.par[1].value=sim_burst->q;
    wf.par[2].name="eccentricity"; wf.par[2].value=sim_burst->pol_ellipse_e;
    wf.par[3].name="phase"; 	   wf.par[3].value=sim_burst->pol_ellipse_angle;
  
    sprintf(wf_name, "LAL_SGE%.*fQ%.*f",d1,wf.par[0].value,d2,wf.par[1].value);

    wf.type = mdc_type;
    wf.name = wf_name;
    wf.name.ReplaceAll(".","d");
    if(decimals==-1) wf.name.ReplaceAll("d0","");
    if(uname!="") wf.name = uname;
    wf.hpPath = wf.name;
    wf.hxPath = wf.name;
    mdcid waveid = AddWaveform(wf);
    return waveid;

  } else
  if(mdc_type==MDC_WNB_LAL) {

    int d1 = decimals==-1 ? 0 : decimals;
    int d2 = decimals==-1 ? 0 : decimals;
    int d3 = decimals==-1 ? 3 : decimals;

    wf.par.resize(3);
    wf.par[0].name="frequency"; wf.par[0].value=sim_burst->frequency;
    wf.par[1].name="bandwidth"; wf.par[1].value=sim_burst->bandwidth;
    wf.par[2].name="duration";  wf.par[2].value=sim_burst->duration;

    sprintf(wf_name,"LAL_WNB%.*f_%.*f_%.*f",d1,wf.par[0].value,d2,wf.par[1].value,d3,wf.par[2].value);

    wf.type = mdc_type;
    wf.name = wf_name;
    wf.name.ReplaceAll(".","d");
    if(uname!="") wf.name = uname;
    wf.hpPath = wf.name;
    wf.hxPath = wf.name;
    mdcid waveid = AddWaveform(wf);
    return waveid;
  
  } else
  if(mdc_type==MDC_GA_LAL) {

    int d1 = decimals==-1 ? 1 : decimals;

    wf.par.resize(1);
    wf.par[0].name="duration"; wf.par[0].value=sim_burst->duration;

    sprintf(wf_name,"LAL_GA%.*f",d1,1000.*wf.par[0].value);

    wf.type = mdc_type;
    wf.name = wf_name;
    wf.name.ReplaceAll(".","d");
    if(uname!="") wf.name = uname;
    wf.hpPath = wf.name;
    wf.hxPath = wf.name;
    mdcid waveid = AddWaveform(wf);
    return waveid;
  } else
  if(mdc_type==MDC_SC_LAL) {

    int d1 = decimals==-1 ? 1 : decimals;

    wf.par.resize(2);
    wf.par[0].name="frequency"; wf.par[0].value=sim_burst->frequency;
    wf.par[1].name="amplitude"; wf.par[1].value=sim_burst->amplitude;

    sprintf(wf_name,"LAL_SC%.*f",d1,wf.par[0].value);

    wf.hx=0;	// string cup waveforms are linearly polarized
    wf.type = mdc_type;
    wf.name = wf_name;
    wf.name.ReplaceAll(".","d");
    if(uname!="") wf.name = uname;
    wf.hpPath = wf.name;
    wf.hxPath = wf.name;
    mdcid waveid = AddWaveform(wf);
    return waveid;
  }

  cout << "CWB::mdc::AddWaveform - Warning : LAL waveform not added !!!" << endl;
  mdcid waveid = {"",0,0};
  return waveid;
}
#endif

//______________________________________________________________________________
void
CWB::mdc::AddWaveform(TString mdc_name, TString hp_fName) {
//
// Add Waveform to this object
//
//
// Input: mdc_name - name of mdc  
//        hp_fName - name of the input text file which contains hp component
//
// NOTE: the input text file is composed by two columns of ascii values (time hp/hx) 
//

  AddWaveform(mdc_name, hp_fName, "");
  return;
}

//______________________________________________________________________________
void
CWB::mdc::AddWaveform(TString mdc_name, TString hp_fName, TString hx_fName) {
//
// Add Waveform to this object
//
//
// Input: mdc_name - name of mdc  
//        hp_fName - name of the input text file which contains hp component
//        hx_fName - name of the input text file which contains hx component
//
// NOTE: the input text file is composed by two columns of ascii values (time hp/hx) 
//

  waveform wf;
  wf.type=MDC_USER;
  wf.name=mdc_name;
  ReadWaveform(wf.hp, hp_fName);
  wf.hpPath=hp_fName;
  if(hx_fName.Sizeof()>1) {
    ReadWaveform(wf.hx, hx_fName);
    wf.hxPath=hx_fName;
  } else {
    wf.hx=wf.hp; wf.hx=0;
    wf.hxPath="";
  }

  if(wf.hp.size()!=wf.hx.size()) {
     cout << "CWB::mdc::AddWaveform - Error : hp,hx size not equal !!!" << endl;
     exit(1);
  }
  if(wf.hp.rate()!=wf.hx.rate()) {
     cout << "CWB::mdc::AddWaveform - Error : hp,hx rate not equal  !!!" << endl;
     exit(1);
  }

  // check if waveform is already declared in the wfList
  int ID=-1;
  for(int i=0;i<(int)wfList.size();i++) if(wfList[i].name.CompareTo(mdc_name)==0) {ID=i;break;}
   
  if(ID==-1) {  // waveform is not in the list
    wfList.push_back(wf);
  } else {
    wfList[ID].list.push_back(wf);
  }

  return;
}

//______________________________________________________________________________
void
CWB::mdc::AddWaveform(TString mdc_name, TString hp_fName, double srate, vector<mdcpar> par) {
//
// Add Waveform to this object
//
//
// Input: mdc_name - name of mdc  
//        hp_fName - name of the input text file which contains hp component
//        srate    - sample rate of the input waveform (Hz)
//        par      - these parameters are only for infos and are stored in the waveform structure
//
// NOTE: the input text file s composed by a column of ascii values (hp/hx) 
//       at constant sample rate (srate)
//

  AddWaveform(mdc_name, hp_fName, "", srate, par);
  return;
}

//______________________________________________________________________________
void
CWB::mdc::AddWaveform(TString mdc_name, TString hp_fName, TString hx_fName, 
                      double srate, vector<mdcpar> par) {
//
// Add Waveform to this object
//
//
// Input: mdc_name - name of mdc  
//        hp_fName - name of the input text file which contains hp component
//        hx_fName - name of the input text file which contains hx component
//        srate    - sample rate of the input waveform (Hz)
//        par      - these parameters are only for infos and are stored in the waveform structure
//                   only hrss is used to modify the input waveform :
//                   - hrss<0  : wf-hrss is not modify
//                   - hrss=0  : wf-hrss is normalized to 1	(default)
//                   - hrss>0  : wf-hrss is normalized to hrss
//
// NOTE:  the input text file s composed by a column of ascii values (hp/hx) 
//        at constant sample rate (srate)
//

  waveform wf;
  wf.type=MDC_USER;
  wf.name=mdc_name;
  ReadWaveform(wf.hp, hp_fName, srate);
  wf.hpPath=hp_fName;
  if(hx_fName.Sizeof()>1) {
    ReadWaveform(wf.hx, hx_fName, srate);
    wf.hxPath=hx_fName;
  } else {
    wf.hx=wf.hp; wf.hx=0;
    wf.hxPath="";
  }

  if(wf.hp.size()!=wf.hx.size()) {
     cout << "CWB::mdc::AddWaveform - Error : hp,hx size not equal !!!" << endl;
     exit(1);
  }
  if(wf.hp.rate()!=wf.hx.rate()) {
     cout << "CWB::mdc::AddWaveform - Error : hp,hx rate not equal  !!!" << endl;
     exit(1);
  }

  // extract hrss from parameters
  bool error=false;
  double hrss_par = GetPar("hrss",par,error);
  if(error) hrss_par=0.;

  // normalization 
  double hrssp=0; for (int i=0;i<(int)wf.hp.size();i++) hrssp+=wf.hp[i]*wf.hp[i];
  double hrssc=0; for (int i=0;i<(int)wf.hx.size();i++) hrssc+=wf.hx[i]*wf.hx[i];
  hrssp=sqrt(hrssp/wf.hp.rate());
  hrssc=sqrt(hrssc/wf.hx.rate());
  double hrss=sqrt(hrssp*hrssp+hrssc*hrssc);

  // if(hrss_par=0) hrss is normalized to 1
  if(hrss_par==0) {
    for(int i=0;i<(int)wf.hp.size();i++) {wf.hp[i]/=hrss;wf.hx[i]/=hrss;}
    hrssp /= hrss;
    hrssc /= hrss;
    hrss = 1.;
  }
  // if(hrss_par>0) hrss is normalized to hrss_par
  if(hrss_par>0) {
    for(int i=0;i<(int)wf.hp.size();i++) {wf.hp[i]/=hrss/hrss_par;wf.hx[i]/=hrss/hrss_par;}
    hrssp /= hrss/hrss_par;
    hrssc /= hrss/hrss_par;
    hrss = hrss_par;
  }

  // set htss,hrssp,hrssc and add par to waveform
  bool bhrss=false;
  vector<mdcpar> wfpar;
  for(int i=0;i<par.size();i++) {
    if(par[i].name=="hrss") {par[i].value=hrss;bhrss=true;}
    wfpar.push_back(par[i]);
  }
  mdcpar upar = {"hrss",1.,""}; if(!bhrss) wfpar.push_back(upar);
  mdcpar ppar = {"hrssp",hrssp,""}; wfpar.push_back(ppar);
  mdcpar cpar = {"hrssc",hrssc,""}; wfpar.push_back(cpar);
  wf.par=wfpar;
 
  // check if waveform is already declared in the wfList
  int ID=-1;
  for(int i=0;i<(int)wfList.size();i++) if(wfList[i].name.CompareTo(mdc_name)==0) {ID=i;break;}
   
  if(ID==-1) {  // waveform is not in the list
    wfList.push_back(wf);
  } else {
    wfList[ID].list.push_back(wf);
  }

  return;
}

//______________________________________________________________________________
mdcid
CWB::mdc::AddWaveform(waveform wf) {
//
// Add Waveform to this object
//
//
// Input: wf       - waveform structure

  if(wf.hp.size()!=wf.hx.size()) {
     cout << "CWB::mdc::AddWaveform - Error : hp,hx size not equal !!!" << endl;
     exit(1);
  }
  if(wf.hp.rate()!=wf.hx.rate()) {
     cout << "CWB::mdc::AddWaveform - Error : hp,hx rate not equal  !!!" << endl;
     exit(1);
  }

  // check if wf mdc type is in the MDC_TYPE list
  if(wf.type<0 || wf.type>MDC_USER) {
     cout << "CWB::mdc::AddWaveform - Error : mdc type not allowed  !!!" << endl;
     exit(1);
  }

  // check if waveform is already declared in the wfList
  int ID=-1;
  for(int i=0;i<(int)wfList.size();i++) if(wfList[i].name.CompareTo(wf.name)==0) {ID=i;break;}

  if(ID==-1) {  // waveform is not in the list
    wfList.push_back(wf);
  } else {
    wfList[ID].list.push_back(wf);
  }

  mdcid waveid;
  waveid.name = wf.name;
  waveid.ID = ID==-1 ? wfList.size()-1 : ID;
  waveid.id = ID==-1 ? 0 : wfList[ID].list.size();
  return waveid;
}

//______________________________________________________________________________
TString 
CWB::mdc::Get(wavearray<double>& x, TString ifo) {
//
// Get burst/inspiral mdc data of the detector = ifo
//
//
// Input:  x       - the input start/stop gps time are obtained from the wavearray values
//                   start = x.start();  stop = x.start()+x.size()/x.rate()
//         ifo     - name of the detector defined in the network 
//
// Output: x       - x.data contains the mdc data
//
// Return: log     - ascii string with mdc parameters 
//

  TString listLog;
#ifdef _USE_LAL
  if(inspOptions!="") 	listLog=GetInspiral(x, ifo);   // new inspiral mdc
  else 			listLog=GetBurst(x, ifo);      // built-in mdc
#else
  listLog=GetBurst(x, ifo);      // built-in mdc
#endif

  return listLog;
}

//______________________________________________________________________________
TString 
CWB::mdc::GetBurst(wavearray<double>& x, TString ifo) {
//
// Get burst mdc data of the detector = ifo
//
//
// Input:  x       - the input start/stop gps time are obtained from the wavearray values
//                   start = x.start();  stop = x.start()+x.size()/x.rate()
//         ifo     - name of the detector defined in the network 
//
// Output: x       - x.data contains the mdc data
//
// Return: log     - ascii string with mdc parameters 
//

  if(net==NULL) {
    cout << "CWB::mdc::GetBurst - Error : Dummy method : network is not initialized " << endl; 
    exit(1);
  }
  if(x.rate()!=MDC_SAMPLE_RATE) {
    cout << "CWB::mdc::GetBurst - Error : x.rate() != " << MDC_SAMPLE_RATE << endl;
    exit(1);
  }

  TString listLog="";
  double deg2rad = TMath::Pi()/180.;
  double rad2deg = 180./TMath::Pi();

  double dt = 1./x.rate();

  double start = x.start();
  double stop  = x.start()+x.size()*dt;

  x=0.;

  mdcList.clear();
  mdcType.clear();
  mdcTime.clear();
  srcList.clear();

  srcList = GetSourceList(start, stop);

  // fill mdcType list
  if(sky_distribution==MDC_XMLFILE) {
    mdcType=xmlType;
  } else {
    for(int i=0;i<(int)wfList.size();i++) {
      bool save=true;
      for(int j=0; j<(int)mdcType.size(); j++){
        if(wfList[i].name.CompareTo(mdcType[j])==0) {save = false; break;}
      }
      if(save) {
        mdcType.push_back(wfList[i].name.Data());
      }
    }
  }

  for(int k=0;k<(int)srcList.size();k++) {

    double gps   = srcList[k].gps;
    double theta = srcList[k].theta;
    double phi   = srcList[k].phi;
    double psi   = srcList[k].psi;
    double rho   = srcList[k].rho;
    double iota  = srcList[k].iota;
    double hrss  = srcList[k].hrss;

    //cout.precision(14);
    //cout << "wf : " << srcList[k].wf.name.Data() << " type :" << srcList[k].wf.type 
    //     << " gps : " << gps << " theta : " << theta << " phi : " << phi 
    //     << " psi : " << psi << " rho : " << rho << " iota : " << iota << endl;

    double fPlus  = GetAntennaPattern(ifo, phi, theta, psi, "hp");
    double fCross = GetAntennaPattern(ifo, phi, theta, psi, "hx");
    double tShift = 0.;
    if(sky_distribution==MDC_XMLFILE) {
      // Time Delay is computed respect to geocenter
      // this is required to be LAL compliant
      tShift = GetDelay(ifo,"",phi,theta);
    } else {
      // Time Delay is computed respect to the first detector in the network
      tShift = GetDelay(ifo,net->ifoName[0],phi,theta);
    }
 
    // if allowed then set ellipticity
    bool ellipticity=false;
    // if MDC_XMLFILE than hp,hx already contains the eccentricity rescaling
    if((srcList[k].wf.type==MDC_SGE_LAL)&&(sky_distribution!=MDC_XMLFILE))  ellipticity=true;
    if(srcList[k].wf.type==MDC_SGE)  ellipticity=true;
    if(srcList[k].wf.type==MDC_CGE)  ellipticity=true;
    if(srcList[k].wf.type==MDC_RDE)  ellipticity=true;
    if(srcList[k].wf.type==MDC_EBBH) ellipticity=true;
    if(srcList[k].wf.type==MDC_USER) ellipticity=true;
    double ePlus  = ellipticity ? (1+cos(iota*deg2rad)*cos(iota*deg2rad))/2 : 1.;
    double eCross = ellipticity ? cos(iota*deg2rad) : 1.;
    // if ellipticity=false we force iota to be 90
    if(!ellipticity) srcList[k].iota=90;
    // for circular waves we force iota to be 0
    if(srcList[k].wf.type==MDC_RDC) srcList[k].iota=0;
    if(srcList[k].wf.type==MDC_SGC) srcList[k].iota=0;
    if(srcList[k].wf.type==MDC_CGC) srcList[k].iota=0;

    waveform wf = srcList[k].wf;

    // build waveform vector 
    int iShift = fabs(tShift)*wf.hp.rate();
    wavearray<double> w(wf.hp.size()+iShift);	// add iShift to take into account time shift
    w.rate(wf.hp.rate());w=0;
    double SimHpHp=0;
    double SimHcHc=0;
    double SimHpHc=0;
    iShift = tShift<0 ? iShift : 0;
    for(int i=0;i<(int)wf.hp.size();i++) {
      w[i+iShift] = ePlus*fPlus*wf.hp[i]+eCross*fCross*wf.hx[i];
      SimHpHp+=wf.hp[i]*wf.hp[i];
      SimHcHc+=wf.hx[i]*wf.hx[i];
      SimHpHc+=wf.hp[i]*wf.hx[i];
    }
    SimHpHp*=dt;
    SimHcHc*=dt;
    SimHpHc*=dt;
    double SrcHrss=sqrt(SimHpHp+SimHcHc);
    std::stringstream linestream;
    int id; double m1,m2,rp0,e0,dist=0.;
    switch(srcList[k].wf.type) {
    case MDC_EBBH :
      // check if distance is already defined by the user
      linestream.str(wf.par[0].name.Data());
      if(!(linestream >> id >> m1 >> m2 >> rp0 >> e0 >> dist)) {
        // distance is not defined
        // for eBBH hp,hx are already scaled to 10 Kpc
        // add distance in Kpc to par string  
        char pars[256];
        sprintf(pars,"%s %g",srcList[k].wf.par[0].name.Data(),rho); 
        srcList[k].wf.par[0].name = pars;
      }
      break;
    case MDC_WNB_LAL :
      if(sky_distribution==MDC_XMLFILE) {
        // for MDC_WNB_LAL the hrss is not present in the sim_table
        hrss=SrcHrss;
        srcList[k].hrss=hrss;
      }
      break;
    case MDC_SGE_LAL :
      if(sky_distribution==MDC_XMLFILE) {
        // if MDC_XMLFILE than iota is computed from the eccentricity 
        // (see SetSkyDistribution MDC_XMLFILE)
        double eccentricity = iota>1e-10 ? iota : 1e-10;
        // in GetBurstLog the iota is reverted to eccentricity 
        srcList[k].iota = (1.-sqrt(1-eccentricity*eccentricity))/eccentricity;
        srcList[k].iota = acos(srcList[k].iota)*rad2deg;
        // if MDC_XMLFILE than hp,hx contain the eccentricity rescaling
        // SimHpHp,SimHcHc,SimHpHc must be rescaled according to the eccentricity
/*
        double ePlus  = (1+cos(iota*deg2rad)*cos(iota*deg2rad))/2;
        double eCross = cos(iota*deg2rad);
        SimHpHp *= 1./(ePlus*ePlus);
        SimHcHc *= fabs(eCross)>1e-5 ? 1./(eCross*eCross) : 0.;
        SimHpHc *= fabs(eCross)>1e-5 ? 1./(ePlus*eCross)  : 0.; 
*/
        //double eccentricity = cosi2e(cos(iota*deg2rad));
        double E = eccentricity;
        double A = 1./sqrt(2-E*E);
        double B = A*sqrt(1-E*E);
        SimHpHp *= 1./(A*A)/2.;
        SimHcHc *= fabs(B)>1e-5 ? 1./(B*B)/2. : 0.;
        SimHpHc *= fabs(B)>1e-5 ? 1./(A*B)/2. : 0.; 
        if(SimHcHc==0) SimHcHc=SimHpHp;

        //cout << "A " << A << " B " << B << " E " << E << endl;
        //cout << "SimHpHp : " << SimHpHp << " " << " SimHcHc " << SimHcHc << " SimHpHc " << SimHpHc << endl;
        //cout << "hrss : " << hrss << " " << " sqrt(SimHpHp+SimHcHc) " << sqrt(SimHpHp+SimHcHc) << endl;
      }
      break;
    default :
      if(inj_hrss>0) {
        // scaled source to hrss (variable for each source) or inj_hrss (fixed for all source) : @ 10Kpc
        SimHpHp *= hrss>0 ? pow(hrss/SrcHrss,2) : pow(inj_hrss/SrcHrss,2);
        SimHcHc *= hrss>0 ? pow(hrss/SrcHrss,2) : pow(inj_hrss/SrcHrss,2);
        SimHpHc *= hrss>0 ? pow(hrss/SrcHrss,2) : pow(inj_hrss/SrcHrss,2); 
        w *= hrss>0 ? hrss/SrcHrss : inj_hrss/SrcHrss;
      } else {
        // the hrss @ 10Kpc is the one which is defined by hp,hc waveforms 
        srcList[k].hrss = SrcHrss;
      }
    }
    // scale amplitude with the inverse of distance (standard candle @ 10Kpc)
    if(rho>0) {
      SimHpHp*=100./pow(rho,2);
      SimHcHc*=100./pow(rho,2);
      SimHpHc*=100./pow(rho,2);
      w*=10./rho;    
    }

    TimeShift(w, tShift);

    // compute the number of samples (sT) of the mdc central time T
    double T = GetCentralTime(wf);
    int   sT = TMath::Nint(T*w.rate());

    // offset (oS) of inj time respect to beginning of array
    // iShift must be subtracted !!!
    int oS = (srcList[k].gps-start)*x.rate()-iShift;

    // add waveform to mdc vector : T @ gps time
    for(int i=0;i<(int)w.size();i++) {
      int j=i+oS-sT;
      if (j>=(int)x.size()) break;
      if (j>=0) x[j]=w[i];
    }

    // add waveform log string to mdc log string
    TString log = GetBurstLog(srcList[k], start, SimHpHp, SimHcHc, SimHpHc);
    listLog = listLog+log;

    // add infos to lists
    mdcList.push_back(log.Data());
//    mdcType.push_back(srcList[k].wf.name.Data());
    mdcTime.push_back(srcList[k].gps);
/*
    bool save=true;
    for(int j=0; j<(int)mdcType.size(); j++){
      if(srcList[k].wf.name.CompareTo(mdcType[j])==0) {save = false; break;}
    }
    if(save) {
      mdcType.push_back(srcList[k].wf.name.Data());
    }
*/
  }

  return listLog;
}

//______________________________________________________________________________
void 
CWB::mdc::ReadWaveform(wavearray<double>& x, TString fName) {
//
// Read Waveform from input ascii file
//
//
// Input:  fName   - name of the input text file which contains hp component
// 
// Output: x       - wavearray x.data contains waveform data
//
// NOTE: the input text file is composed by two columns of ascii values  
//

  Long_t id,fsize,flags,mt;
  int estat = gSystem->GetPathInfo(fName.Data(),&id,&fsize,&flags,&mt);
  if (estat!=0) {
    cout << "CWB::mdc::ReadWaveform - File : " << fName.Data() << " Not Exist" << endl;
    exit(1);
  }

  // read Waveform
  ifstream in;
  in.open(fName.Data(),ios::in);
  if (!in.good()) {cout << "CWB::mdc::ReadWaveform - Error Opening File : " << fName.Data() << endl;exit(1);}

  int size=0;
  char* str = new char[fsize+1];
  TObjArray* tok;
  while(true) {
    in.getline(str,fsize+1);
    if (!in.good()) break;
    if(str[0] != '#') size++;
    else continue;
    tok = TString(str).Tokenize(TString(' '));
    if(tok->GetEntries()!=2) {
      cout << "CWB::mdc::ReadWaveform - Input file with bad format : must be 2 columns " << endl;
      exit(1);
    }
    delete tok;
  }
  in.clear(ios::goodbit);
  in.seekg(0, ios::beg);

  wavearray<double> w;
  w.resize(size);
  wavearray<double> t=w;

  int cnt=0;
  double tMin=1.e30;
  double tMax=0.;
  while (1) {
    int fpos = in.tellg();
    in.getline(str,1024);
    if (!in.good()) break;
    if(str[0] != '#') {
      in.seekg(fpos, ios::beg);
      in >> t.data[cnt] >> w.data[cnt];
      if (!in.good()) break;
      if(t[cnt]<tMin) tMin=t[cnt];
      if(t[cnt]>tMax) tMax=t[cnt];
      cnt++;
      fpos=in.tellg();
      in.seekg(fpos+1, ios::beg);
    }
  }
  in.close();
  delete [] str;

  // convert to rate MDC_SAMPLE_RATE

  x.rate(MDC_SAMPLE_RATE);
  x.start(0.);
  int offset = 0.05*x.rate();
  x.resize(offset+int(x.rate()*(tMax-tMin)));
  x=0.;

  double dt=1./x.rate();
  wavearray<double> r=x;
  for (int i=0;i<(int)r.size();i++) r[i] = tMin+(i-offset)*dt;

  CWB::Toolbox::convertSampleRate(t, w, r, x);

  return;
}

//______________________________________________________________________________
void 
CWB::mdc::ReadWaveform(wavearray<double>& x, TString fName, double srate) {
//
// Read Waveform from input ascii file
//
//
// Input: fName    - name of the input text file which contains hp component
//        srate    - sample rate of the input waveform (Hz)
// 
// Output: x       - wavearray x.data contains waveform data
//
// NOTE: the input text file s composed by a column of ascii values  
//       at constant sample rate (srate)
//

  Long_t id,fsize,flags,mt;
  int estat = gSystem->GetPathInfo(fName.Data(),&id,&fsize,&flags,&mt);
  if (estat!=0) {
    cout << "CWB::mdc::ReadWaveform - File : " << fName.Data() << " Not Exist" << endl;
    exit(1);
  }

  // read Waveform
  ifstream in;
  in.open(fName.Data(),ios::in);
  if (!in.good()) {cout << "CWB::mdc::ReadWaveform - Error Opening File : " << fName.Data() << endl;exit(1);}

  int size=0;
  char* str = new char[fsize+1];
  TObjArray* tok;
  while(true) {
    in.getline(str,fsize+1);
    if (!in.good()) break;
    if(str[0] != '#') size++;
    tok = TString(str).Tokenize(TString(' '));
    if((tok->GetEntries()!=1)&&(size>1)) {
      cout << "CWB::mdc::ReadWaveform - Input file with bad format : must be 1 column " << endl;
      exit(1);
    }
    if((tok->GetEntries()>1)&&(size==1)) {	// all data are in one line
      x.resize(tok->GetEntries());
      x.rate(srate);
      for(int i=0;i<x.size();i++) {		// get data from line
        TString stok  = ((TObjString*)tok->At(i))->GetString();
        x[i]=stok.Atof();			// fill array
      }
      size=0;break;
    }
    delete tok;
  }
  in.clear(ios::goodbit);
  in.seekg(0, ios::beg);

  if(size>0) {					// data are written with 1 column format
    x.resize(size);
    x.rate(srate);

    int cnt=0;
    while (1) {
      in.getline(str,1024);
      if (!in.good()) break;
      if(str[0] != '#') {
        x[cnt]=TString(str).Atof();
        cnt++;
      }
    }
  }
  in.close();

  // resample data if is not equals to the default
  if(srate!=MDC_SAMPLE_RATE) {
    wavearray<double> y=x;
    y.resample(x,MDC_SAMPLE_RATE);
    y*=x.rms()/y.rms();	// set y energy = x energy
    x=y;
  }

  delete [] str;

  return;
}

//______________________________________________________________________________
void
CWB::mdc::GetSourceCoordinates(double gps, double& theta, double& phi, double& psi, double& rho, 
                               double& iota, double& hrss, int& ID, int& id) {
//
// Get source coordinates (earth system) from the user defined sky distribution
//
//
// Output: gps     - time         (sec)
//         theta   - latitude     (degrees)
//         phi     - longitude    (degrees)
//         psi     - polarization (degrees)
//         rho     - distance     (KPc)
//         iota    - elliptical inclination angle (degrees)
//         hrss    - sqrt(hp^2+hx^2)
//

  GetSourceCoordinates(theta, phi, psi, rho, iota, hrss, ID, id);

  if((sky_distribution==MDC_GWGC)||
     (sky_distribution==MDC_MNGD)||
     (sky_distribution==MDC_CELESTIAL_FIX)) {

    if(gps>0) phi=sm.RA2phi(phi,gps);  // celestial 2 earth coordinates
  }

  return;
}

//______________________________________________________________________________
void
CWB::mdc::GetSourceCoordinates(double& theta, double& phi, double& psi, double& rho, 
                               double& iota, double& hrss, int& ID, int& id) {
//
// Get source coordinates (earth system) from the user defined sky distribution
//
//
// Output: theta   - latitude     (degrees)
//         phi     - longitude    (degrees)
//         psi     - polarization (degrees)
//         rho     - distance     (KPc)
//         iota    - elliptical inclination angle (degrees)
//         hrss    - sqrt(hp^2+hx^2)
//

  if((sky_distribution==MDC_GWGC)||
     (sky_distribution==MDC_MNGD)||
     (sky_distribution==MDC_RANDOM)||
     (sky_distribution==MDC_CELESTIAL_FIX)||
     (sky_distribution==MDC_EARTH_FIX)) {

     if(thList.size()==0) {
       cout << "CWB::mdc::GetSourceCoordinates - Error : injections not defined  !!!" << endl;
       exit(1);
     }

    int id = gRandom->Uniform(0,thList.size());

    theta = thList[id];
    phi   = phList[id];
    psi   = psiList[id];
    rho   = rhoList[id];
    iota  = iotaList[id];
    hrss  = hrssList[id];
      ID  = IDList[id];
      id  = idList[id];

    GeographicToCwb(phi,theta,phi,theta);

  } else {

    if(inj==NULL) {
      cout << "CWB::mdc::GetSourceCoordinates - Error : distribution not defined !!!" << endl;
      exit(1);
    }
  }

  return;
}

//______________________________________________________________________________
waveform 
CWB::mdc::GetSourceWaveform(int& ID, int& id) {
//
// Return waveform randomly selected from the input list
//

  ID = (int)gRandom->Uniform(0,wfList.size());
  id = (int)gRandom->Uniform(0,1+wfList[ID].list.size());

  waveform wf = id==0 ? wfList[ID] : wfList[ID].list[id-1]; 
  GetWaveform(wf);   // if hp,hx are empty -> fill waveforms
  return wf;
}

//______________________________________________________________________________
vector<source>
CWB::mdc::GetSourceList(double start, double stop) {
//
// Get a list of waveforms in the interval [start, stop]
//
//
// Input: start    - start time
//        stop     - stop  time
//
// Return the source list
//

  if((stop-start)<=2*inj_length) {
    cout << "CWB::mdc::GetSourceList - Warning : buffer too small (stop-start)<=2*inj_length !!!" << endl;
    exit(1);
  }

  // fix random choice setting seed=start+srcList_seed
  gRandom->SetSeed(int(start)+srcList_seed);

  source src;

  double timeStep = inj_rate>0 ? 1./inj_rate : 1./MDC_INJ_RATE; 

  int iStart=int(0.5+start/timeStep);
  int iStop=int(stop/timeStep);

  if(iStart==0) iStart+=1;
  vector<source> src_list;

  if(sky_distribution==MDC_LOGFILE) {

    if(inj_tree==NULL) {
      cout << "CWB::mdc::GetSourceList - Error : injection object is NULL" << endl;
      exit(1);
    }

    // gps are valid if are inside the range [start+inj_length, stop-inj_length]
    // inj_length is the scratch
    char cut[128];sprintf(cut,"time[0]>=%f && time[0]<%f",start+inj_length,stop-inj_length);
    inj_tree->Draw("Entry$",cut,"goff");
    int entries = inj_tree->GetSelectedRows();
    float   phi[4];              
    float   theta[4];            
    float   psi[2];              
    double  time[2*NIFO_MAX];   
    int     type[2];             

    inj_tree->SetBranchAddress("phi",phi);
    inj_tree->SetBranchAddress("theta",theta);
    inj_tree->SetBranchAddress("psi",psi);
    inj_tree->SetBranchAddress("time",time);
    inj_tree->SetBranchAddress("type",type);

    double* entry = inj_tree->GetV1();
    for(int i=0;i<entries;i++) {
      inj_tree->GetEntry(entry[i]);
      src.theta = theta[0];
      src.phi   = phi[0];
      src.psi   = psi[0];
      src.gps   = time[0]; 
      src.rho   = 10.; 
      src.iota  = 0.; 
      src.hrss  = 0.; 
      int TYPE = (type[0]-1)<mdcName.size() ? type[0]-1 : 0;
      int ID = GetWaveformID(mdcName[TYPE])!=-1 ? GetWaveformID(mdcName[TYPE]) : 0;
      int id = (int)gRandom->Uniform(0,1+wfList[ID].list.size());
      src.wf = GetWaveform(ID,id);
      src_list.push_back(src);
    }
  } else if((sky_distribution==MDC_CUSTOM)||(sky_distribution==MDC_XMLFILE)) {
    for(int i=0;i<(int)gpsList.size();i++) {
      // gps are valid if are inside the range [start+inj_length, stop-inj_length]
      // inj_length is the scratch
      src.gps   = gpsList[i]; 
      if(src.gps<=start+inj_length) continue;
      if(src.gps>=stop-inj_length)  continue;
      src.theta = thList[i];
      src.phi   = phList[i];
      src.psi   = psiList[i];
      src.rho   = rhoList[i]; 
      src.iota  = iotaList[i]; 
      src.hrss  = hrssList[i]; 
      src.ID    = IDList[i]; 
      src.id    = idList[i]; 
      if(src.ID<0) {              // select waveform by name (default)
        src.ID = GetWaveformID(nameList[i])!=-1 ? GetWaveformID(nameList[i]) : 0;
      } else {                    // select waveform by ID,id
        src.ID = (src.ID < wfList.size()) ? src.ID : 0;
      }
      if(src.id<0) {              // select waveform by name (default)
        src.id = (int)gRandom->Uniform(0,1+wfList[src.ID].list.size());
      } else {                    // select waveform by ID,id
        src.id = (src.id <= wfList[src.ID].list.size()) ? src.id : 0;
      }
      src.wf = GetWaveform(src.ID,src.id);
      src_list.push_back(src);
    }
  } else {
    for(int i=iStart;i<=iStop;i++) {
      if(gpsList.size()>0) src.gps=gpsList[0];  // fix gps
      else src.gps=inj_offset+i*timeStep+gRandom->Uniform(-inj_jitter,inj_jitter);
      // gps are valid if are inside the range [start+inj_length, stop-inj_length]
      // inj_length is the scratch
      if(src.gps<=start+inj_length) continue;
      if(src.gps>=stop-inj_length)  continue;
      GetSourceCoordinates(src.gps, src.theta, src.phi, src.psi, src.rho, src.iota, src.hrss, src.ID, src.id);
      src.wf = GetSourceWaveform(src.ID, src.id);

      if(src.wf.par.size() && src.wf.par[0].value==MDC_EBBH) {  // eBBH waveform
        int id; double m1,m2,rp0,e0,dist=0.;
        std::stringstream linestream(src.wf.par[0].name.Data());
        // get user defined distance
        if((linestream >> id >> m1 >> m2 >> rp0 >> e0 >> dist)) src.rho=dist;
      }

      src_list.push_back(src);
      if(gpsList.size()>0) break;  // fix gps, only one injection
    }
  }

  return src_list;
}

//______________________________________________________________________________
TString
CWB::mdc::GetBurstLog(source src, double FrameGPS, double SimHpHp, double SimHcHc, double SimHpHc) {
//
// Get Log string  
//
// Input: src      - source structure
//        FrameGPS - Frame gps time
//        SimHpHp  - Energy of hp component
//        SimHcHc  - Energy of hc component
//        SimHpHc  - Cross component of hp hc
//

  if(net==NULL) {
    cout << "CWB::mdc::GetBurstLog - Error : Dummy method : network is not initialized " << endl; 
    exit(1);
  }

  double deg2rad = TMath::Pi()/180.;

  // log burstMDC parameters
  char   logString[10000]="";

  char   GravEn_SimID[1024]; 
  double hrss = src.hrss>0 ? src.hrss : inj_hrss;
  TString hpPath = src.wf.hpPath.Sizeof()>1 ? src.wf.hpPath : src.wf.name; 
  sprintf(GravEn_SimID,"%s",hpPath.Data()); 
  if(src.wf.hxPath.Sizeof()>1) sprintf(GravEn_SimID,"%s;%s",GravEn_SimID,src.wf.hxPath.Data()); 
  // scale amplitude with the inverse of distance (standard candle @ 10Kpc)
  //double SimHrss      = src.rho>0 ? 10.*sqrt(SimHpHp+SimHcHc)/src.rho : sqrt(SimHpHp+SimHcHc);
  double SimHrss      = sqrt(SimHpHp+SimHcHc);
  double SimEgwR2     = 0.0;
  double GravEn_Ampl  = src.rho>0 ? 10.*hrss/src.rho : hrss;
  double Internal_x   = cos(src.iota*deg2rad);
  double Internal_phi = 0.0;
  double External_x   = cos(src.theta*deg2rad);     // DOM
  double External_phi = src.phi*deg2rad;
  double External_psi = src.psi*deg2rad;

  double EarthCtrGPS = src.gps;
  char   SimName[64]; strcpy(SimName,src.wf.name.Data());

  sprintf(logString,"%s",GravEn_SimID);
  sprintf(logString,"%s %e",logString,SimHrss);
  sprintf(logString,"%s %e",logString,SimEgwR2);
  sprintf(logString,"%s %e",logString,GravEn_Ampl);
  sprintf(logString,"%s %e",logString,Internal_x);
  sprintf(logString,"%s %e",logString,Internal_phi);
  sprintf(logString,"%s %e",logString,External_x);
  sprintf(logString,"%s %e",logString,External_phi);
  sprintf(logString,"%s %e",logString,External_psi);

  sprintf(logString,"%s %10.6f",logString,FrameGPS);
  sprintf(logString,"%s %10.6f",logString,EarthCtrGPS);
  sprintf(logString,"%s %s",logString,SimName);
  sprintf(logString,"%s %e",logString,SimHpHp);
  sprintf(logString,"%s %e",logString,SimHcHc);
  sprintf(logString,"%s %e",logString,SimHpHc);

  int nIFO=net->ifoListSize();
  TString refIFO=net->ifoName[0];
  for(int i=0;i<nIFO;i++) {
    TString ifo=net->ifoName[i];
    double IFOctrGPS = EarthCtrGPS;
    if(sky_distribution==MDC_XMLFILE) {
      // Time Delay is computed respect to geocenter
      // this is required to be LAL compliant
      IFOctrGPS += GetDelay(ifo,"",src.phi,src.theta);
    } else {
      // Time Delay is computed respect to the first detector in the network
      IFOctrGPS += GetDelay(ifo,refIFO,src.phi,src.theta);
    }
    double IFOfPlus  = GetAntennaPattern(ifo, src.phi, src.theta, src.psi, "hp");
    double IFOfCross = GetAntennaPattern(ifo, src.phi, src.theta, src.psi, "hx");
    if(src.wf.name=="eBBH") {  	// if MDC is eBBH add auxiliary infos
      int id; double m1,m2,rp0,e0,distance=0.;
      std::stringstream linestream(src.wf.par[0].name.Data());
      linestream >> id >> m1 >> m2 >> rp0 >> e0 >> distance;
      // compute effective distance
      double cosiota  = cos(src.iota*deg2rad);
      double eff_dist = distance / sqrt(pow(IFOfPlus*(1.+pow(cosiota,2)),2)/4.+pow(IFOfCross*cosiota,2));
      sprintf(logString,"%s %s %10.6f %e %e %g",logString,ifo.Data(),IFOctrGPS,IFOfPlus,IFOfCross,eff_dist/1000.);
    } else {
      sprintf(logString,"%s %s %10.6f %e %e",logString,ifo.Data(),IFOctrGPS,IFOfPlus,IFOfCross);
    }
  }

  if(src.wf.name=="eBBH") {  	// if MDC is eBBH add auxiliary infos
    int id; double m1,m2,rp0,e0,distance=0.,redshift=0.;
    std::stringstream linestream(src.wf.par[0].name.Data());
    linestream >> id >> m1 >> m2 >> rp0 >> e0 >> distance >> redshift;
    sprintf(logString, "%s ebbh ",logString);
    sprintf(logString, "%s mass1 %g ",logString, m1);
    sprintf(logString, "%s mass2 %g ",logString, m2);
    sprintf(logString, "%s rp0 %g ",logString, rp0);
    sprintf(logString, "%s e0 %g ",logString, e0);
    sprintf(logString, "%s distance %g ",logString, distance/1000.);	// converted in Mpc (compatible with LAL)
    sprintf(logString, "%s redshift %g ",logString, redshift);
    // chirp mass (solar units)
    double M = (m1+m2);
    double mu = (m1*m2)/(m1+m2);
    double eta = mu/M;
    double mchirp = pow(mu,3./5)*pow(M,2./5);
    sprintf(logString, "%s mchirp %g ",logString, mchirp);
  } else {
    // distance is converted in Mpc (compatible with LAL)
    if(src.rho>0) sprintf(logString, "%s distance %g ",logString, src.rho/1000.);
  }

  return logString;
}

//______________________________________________________________________________
watplot*
CWB::mdc::Draw(int ID, int id, TString polarization, MDC_DRAW type, TString options, Color_t color) {
//
// Draw waveform in time/frequency domain 
//
//
// Input: ID           - major id of waveform list
//        id           - minor id of waveform list (Ex : the list WNB with different random waveforms)
//        polarization - hp/hx
//        type         - MDC_TIME, MDC_FFT, MDC_TF
//        options      - graphic options
//        color        - option for MDC_TIME, MDC_FFT
//

  watplot* plot=NULL;
  polarization.ToUpper();
  waveform wf = GetWaveform(ID,id);
  if(wf.status) {
    if(polarization.Contains("HP")) plot=Draw(wf.hp,type,options,color);
    if(polarization.Contains("HX")) plot=Draw(wf.hx,type,options,color);
  }
  return plot;
}

//______________________________________________________________________________
watplot*
CWB::mdc::Draw(TString name, int id, TString polarization, MDC_DRAW type, TString options, Color_t color) {
//
// Draw waveform in time/frequency domain 
//
//
// Input: name         - waveform name
//        polarization - hp/hx
//        type         - MDC_TIME, MDC_FFT, MDC_TF
//        options      - graphic options
//        color        - option for MDC_TIME, MDC_FFT
//

  watplot* plot=NULL;
  polarization.ToUpper();
  waveform wf = GetWaveform(name,id);
  if(wf.status) {
    if(polarization.Contains("HP")) plot=Draw(wf.hp,type,options,color);
    if(polarization.Contains("HX")) plot=Draw(wf.hx,type,options,color);
  }
  return plot;
}

//______________________________________________________________________________
watplot*
CWB::mdc::Draw(TString ifo, double gpsStart, double gpsEnd, int id, 
               MDC_DRAW type, TString options, Color_t color) {
//
// Draw waveform in time/frequency domain contained in the interval [gpsStart,gpsEnd]
//
//
// Input: ifo      - detector name
//        gpsStart - start interval time 
//        gpsEnd   - stop interval time 
//        id       - sequential numer of the waveform contained in the interval
//        type     - MDC_TIME, MDC_FFT, MDC_TF
//        options  - graphic options
//        color    - option for MDC_TIME, MDC_FFT
//

  wavearray<double> y; 
  y.rate(MDC_SAMPLE_RATE);
  y.start(gpsStart);
  y.resize((gpsEnd-gpsStart)*y.rate());
  Get(y,ifo);

  if((int)mdcTime.size()==0) {
    cout << "CWB::mdc::Draw : Error - No events in the selected period" << endl;
    exit(1);
  }
  if(id>=(int)mdcTime.size()) id=(int)mdcTime.size()-1;
  double tOffset = mdcTime[id]-y.start();
  double tWindow = inj_length;
  int iStart = (tOffset-tWindow/2)*y.rate(); if(iStart<0) iStart=0;
  int iEnd   = (tOffset+tWindow/2)*y.rate(); if(iEnd>y.size()) iEnd=y.size();

  // make an array centered around the selected event 
  wavearray<double> x(tWindow*y.rate()); 
  x.rate(y.rate());
  for(int i=iStart;i<iEnd;i++) x[i-iStart]=y[i];

  watplot* plot=NULL;

  if(type==MDC_TIME) plot=DrawTime(x,options,color);
  if(type==MDC_FFT)  plot=DrawFFT(x,options,color);
  if(type==MDC_TF)   DrawTF(x);

  return plot;
}  

//______________________________________________________________________________
watplot*
CWB::mdc::Draw(wavearray<double>& x, MDC_DRAW type, TString options, Color_t color) {
//
// Draw waveform in time/frequency domain
//
//
// Input: x       - wavearray which contains waveform data
//        type    - MDC_TIME, MDC_FFT, MDC_TF
//        options - graphic options
//        color   - option for MDC_TIME, MDC_FFT
// 

  watplot* plot=NULL;

  if(type==MDC_TIME) plot=DrawTime(x,options,color);
  if(type==MDC_FFT)  plot=DrawFFT(x,options,color);
  if(type==MDC_TF)   DrawTF(x,options);

  return plot;
}  

//______________________________________________________________________________
watplot*
CWB::mdc::DrawTime(wavearray<double>& x, TString options, Color_t color) {
//
// Draw waveform in time domain
//
//
// Input: x       - wavearray which contains waveform data
//        options - draw options (same as TGraph)
//                  if contains ZOOM then the interval around the signal is showed
//

  if(x.rate()<=0) {
    cout << "CWB::mdc::DrawTime : Error - rate must be >0" << endl;
    exit(1);
  }

  options.ToUpper();
  options.ReplaceAll(" ","");
  if(stft!=NULL) delete stft;
  if(options.Contains("SAME")&&(pts!=NULL)) {
  } else {
    if(pts!=NULL) delete pts;
    char name[32];sprintf(name,"TIME-gID:%d",int(gRandom->Rndm(13)*1.e9));
    pts = new watplot(const_cast<char*>(name),200,20,800,500);
  }
  double tStart,tStop;
  if(options.Contains("ZOOM")) {
    options.ReplaceAll("ZOOM","");
    GetTimeRange(x, tStart, tStop, epzoom);
    tStart+=x.start();
    tStop+=x.start();
  } else {
    tStart=x.start();
    tStop=x.start()+x.size()/x.rate();
  }
  if(!options.Contains("SAME")) {
    if(!options.Contains("A")) options=options+" A";
    if(!options.Contains("L")) options=options+" L";
    if(!options.Contains("P")) options=options+" P";
  }
  pts->plot(x, const_cast<char*>(options.Data()), color, tStart, tStop);

  return pts;
}

//______________________________________________________________________________
watplot*
CWB::mdc::DrawFFT(wavearray<double>& x, TString options, Color_t color) {
//
// Draw waveform spectrum
//
//
// Input: x       - wavearray which contains waveform data
//        options - draw options (same as TGraph)
//                  if contains ZOOM then the interval around the signal is showed
//                  if contains NOLOGX/NOLOGY X/Y axis are linear   
//

  options.ToUpper();
  options.ReplaceAll(" ","");
  if(stft!=NULL) delete stft;
  if(options.Contains("SAME")&&(pts!=NULL)) {
  } else {
    if(pts!=NULL) delete pts;
    char name[32];sprintf(name,"FREQ-gID:%d",int(gRandom->Rndm(13)*1.e9));
    pts = new watplot(const_cast<char*>(name),200,20,800,500);
  }

  double fLow  = 32.;
  double fHigh = x.rate()/2.;
  double tStart,tStop;
  if(options.Contains("ZOOM")) {
    options.ReplaceAll("ZOOM","");
    GetTimeRange(x, tStart, tStop, epzoom);
    tStart+=x.start();
    tStop+=x.start();
  } else {
    tStart=x.start();
    tStop=x.start()+x.size()/x.rate();
  }
  bool logx = true;
  if(options.Contains("NOLOGX")) {logx=false;options.ReplaceAll("NOLOGX","");}
  bool logy = true;
  if(options.Contains("NOLOGY")) {logy=false;options.ReplaceAll("NOLOGY","");}
  pts->plot(x, const_cast<char*>(options.Data()), color, tStart, tStop, true, fLow, fHigh);
  pts->canvas->SetLogx(logx);
  pts->canvas->SetLogy(logy);

  return pts;
}

//______________________________________________________________________________
void
CWB::mdc::DrawTF(wavearray<double>& x, TString options) {
//
// Draw waveform spectrogram
//
//
// Input: x       - wavearray which contains waveform data
//        options - draw options (same as TH2D)
//

  int nfact=4;
  int nfft=nfact*512;
  int noverlap=nfft-10;
  double fparm=nfact*6;

  double tStart,tStop;
  if(options.Contains("ZOOM")) {
    options.ReplaceAll("ZOOM","");
    GetTimeRange(x, tStart, tStop, epzoom);
    tStart+=x.start();
    tStop+=x.start();
  } else {
    tStart=x.start();
    tStop=x.start()+x.size()/x.rate();
  }

  if(stft!=NULL) delete stft;
  if(pts!=NULL) delete pts;
  stft = new CWB::STFT(x,nfft,noverlap,"amplitude","gauss",fparm);

  double fLow  = 32.;
  double fHigh = x.rate()/2.;
  stft->Draw(tStart,tStop,fLow,fHigh,0,0,1);
  stft->GetCanvas()->SetLogy(true);

  return;
}

//______________________________________________________________________________
waveform
CWB::mdc::GetWaveform(int ID, int id) {
//
// Get waveform 
//
//
// Input: ID           - major id of waveform list
//        id           - minor id of waveform list (Ex : the list WNB with different random waveforms)
//
// Return waveform structure
//

  waveform wf;

  // check if waveform is declared in the wfList
  if(ID<0||ID>=(int)wfList.size()) {
    cout << "CWB::mdc::GetWaveform - Error : ID " << ID << " not in the list" << endl;
    wf.status=false;
    return wf;
  }

  id = (id<0||id>(int)wfList[ID].list.size()) ? 0 : id;
  if((int)wfList[ID].list.size()==0) {
    wf=wfList[ID];
  } else {
    wf = id==0 ? wf=wfList[ID] : wfList[ID].list[id-1];
  }

  GetWaveform(wf);   // if hp,hx are empty -> fill waveforms

  wf.status=true;
  return wf;
}

//______________________________________________________________________________
void
CWB::mdc::GetWaveform(waveform& wf) {
//
// if hp,hx are empty -> fill waveforms 
// only MDC_EBBH is implemented 
//

  if((wf.hp.size()==0)&&(wf.hx.size()==0)) {
#ifdef _USE_EBBH
    if(wf.par[0].value==MDC_EBBH) { // Get eBBH waveform 
      int id;
      double m1,m2,rp0,e0,dist=0.,redshift=0.;
      std::stringstream linestream(wf.par[0].name.Data());
      if(!(linestream >> id >> m1 >> m2 >> rp0 >> e0 >> dist >> redshift)) {
        linestream.str(wf.par[0].name.Data()); 
        linestream.clear(); 		// clear stringstream error status	
        if(!(linestream >> id >> m1 >> m2 >> rp0 >> e0 >> dist)) {
          linestream.str(wf.par[0].name.Data()); 
          linestream.clear(); 		// clear stringstream error status
          if(!(linestream >> id >> m1 >> m2 >> rp0 >> e0)) {
            cout << "CWB::mdc::GetWaveform - Wrong eBBH parameter format : " << endl;
            cout << "input line : " << endl;
            cout << wf.par[0].name.Data() << endl;
            cout << "must be : " << endl;
            cout << "event# " << "id " << " m1 " << " m2 " << " rp0 " << " e0 " << endl;
            cout << "or : " << endl;
            cout << "event# " << "id " << " m1 " << " m2 " << " rp0 " << " e0 " << " dist " << endl;
            cout << "or : " << endl;
            cout << "event# " << "id " << " m1 " << " m2 " << " rp0 " << " e0 " << " dist " << " redshift " << endl;
            exit(1);
          }
        }
      }

      if(wf.par.size()==1) {	// create eBBH on the fly

        getEBBH(m1,m2,rp0,e0,wf.hp,wf.hx);

        // the distance of source is G*M/c^2  meters
        double G  = watconstants::GravitationalConstant();
        double M  = watconstants::SolarMass();
        double c  = watconstants::SpeedOfLightInVacuo();
        double pc = watconstants::Parsec();
        double distance_source_Kpc = (m1+m2)*G*M/(c*c)/pc/1.e3;

        // rescale hp,hx to 10 Kpc
        wf.hp*=distance_source_Kpc/10.;
        wf.hx*=distance_source_Kpc/10.;

      } else 
      if(wf.par.size()==2) {	// read eBBH from file

        TString fName = wf.par[1].name;		// get root file name
        int entry = int(wf.par[1].value);	// get tree entry number

        TFile* efile = new TFile(fName);
        if(efile==NULL) {
          cout << "CWB::mdc::GetWaveform - Error opening root file : " << fName.Data() << endl;
          exit(1);
        }

        wavearray<double>* hp = new wavearray<double>;
        wavearray<double>* hx = new wavearray<double>;

        TTree* etree = (TTree *) efile->Get("ebbh");
        if(etree==NULL) {
          cout << "CWB::mdc::GetWaveform - file : " << fName.Data()  
               << " not contains tree ebbh" << endl;
          exit(1);
        }
  
        etree->SetBranchAddress("hp",&hp);
        etree->SetBranchAddress("hx",&hx);

        etree->GetEntry(entry);

        if(hp->size()!=hx->size()) {
          cout << "CWB::mdc::GetWaveform - Error : hp,hx size not equal !!!" << endl;
          exit(1);
        }
        if(hp->rate()!=hx->rate()) {
          cout << "CWB::mdc::GetWaveform - Error : hp,hx rate not equal  !!!" << endl;
          exit(1);
        }

        wf.hp = *hp;
        wf.hx = *hx;

        delete hp;
        delete hx;
        delete efile;

      } else {
        cout << "CWB::mdc::GetWaveform - Error : wrong number of parameters not !!!" << endl;
        exit(1);
      }
    }
#endif
    if((wf.hp.size()==0)&&(wf.hx.size()==0)) {
      cout << "CWB::mdc::GetWaveform - Error : hp,hp are empty !!!" << endl;
      exit(1);
    }
  }
}

//______________________________________________________________________________
waveform
CWB::mdc::GetWaveform(TString name, int id) {
//
// Get waveform 
//
//
// Input: name   - name of waveform
//        id     - minor id of waveform
//
// Return waveform structure
//

  waveform wf;

  // check if waveform is declared in the wfList
  int ID=-1;
  for(int i=0;i<(int)wfList.size();i++) if(wfList[i].name.CompareTo(name)==0) {ID=i;break;}
  if(ID<0) {
    wf.status=false;
    return wf;
  } else {
    id = (id<0||id>(int)wfList[ID].list.size()-1) ? 0 : id;
    if((int)wfList[id].list.size()==0) {
      wf=wfList[ID];
    } else {
      wf=wfList[ID].list[id];
    }

    GetWaveform(wf);   // if hp,hx are empty -> fill waveforms

    wf.status=true;
    return wf;
  }
}

//______________________________________________________________________________
int
CWB::mdc::GetWaveformID(TString name) {
//
// Get waveform 
//
//
// Input: name   - name of waveform
//
// Return waveform ID
//

  // check waveform name declared in the wfList
  int ID=-1;
  for(int i=0;i<(int)wfList.size();i++) if(wfList[i].name.CompareTo(name)==0) {ID=i;break;}
  return ID;
}

//______________________________________________________________________________
void
CWB::mdc::Print(int level) {
//
// Print list of waveforms
//
//
// Input: level   - if level=0 only the waveforms with major id are listed
//

  cout << endl;
  for(int i=0;i<(int)wfList.size();i++) {
    printf("ID : %02i (x% 3i) \t%s\n",i,(int)wfList[i].list.size()+1,wfList[i].name.Data());
    printf("% 3i - ",1);
    for(int k=0;k<(int)wfList[i].par.size();k++) 
      if(wfList[i].par[k].svalue!="") {
        printf("%s = %s ",wfList[i].par[k].name.Data(),wfList[i].par[k].svalue.Data());
      } else {
        printf("%s = %g ",wfList[i].par[k].name.Data(),wfList[i].par[k].value);
      }
    printf("\n");
    if(level>0) {
      for(int j=0;j<(int)wfList[i].list.size();j++) {
        printf("% 3i - ",j+2);
        for(int k=0;k<(int)wfList[i].par.size();k++) 
          if(wfList[i].list[j].par[k].svalue!="") {
            printf("%s = %s ",wfList[i].list[j].par[k].name.Data(),wfList[i].list[j].par[k].svalue.Data());
          } else {
            printf("%s = %g ",wfList[i].list[j].par[k].name.Data(),wfList[i].list[j].par[k].value);
          }
        printf("\n");
      }
    }
  }
  cout << endl;
  return;
}

//______________________________________________________________________________
double
CWB::mdc::GetCentralTime(waveform wf) {
//
// Get mdc central time
//
//
// Input: wf      - waveform structure which contains the waveform data
//
// Retun central time
//

  wavearray<double> z=wf.hp; z+=wf.hx;
  return GetCentralTime(z);
}

//______________________________________________________________________________
double
CWB::mdc::GetCentralTime(wavearray<double> x) {
//
// Get mdc central time
//
//
// Input: x      - wavearray which contains the waveform data
//
// Retun central time
//

  double a; 
  double E=0.,T=0.;
  int size=(int)x.size();
  double rate=x.rate(); 
  for(int j=0;j<size;j++) {
    a = x[j];
    T += a*a*j/rate;                   // central time
    E += a*a;                          // energy
  }
  T = E>0 ? T/E : 0.5*size/rate;

  return T;
}

//______________________________________________________________________________
double
CWB::mdc::GetCentralFrequency(waveform wf) {
//
// Get mdc central frequency
//
//
// Input: wf      - waveform structure which contains the waveform data
//
// Retun central frequency
//

  wavearray<double> z=wf.hp; z+=wf.hx;
  return GetCentralFrequency(z);
}

//______________________________________________________________________________
double
CWB::mdc::GetCentralFrequency(wavearray<double> x) {
//
// Get mdc central frequency
//
//
// Input: x      - wavearray which contains the waveform data
//
// Retun central frequency
//

  double a; 
  double E=0.,F=0.;
  int size=(int)x.size();
  double rate=x.rate(); 
  x.FFTW(1);
  double dF=(rate/(double)size)/2.;
  for(int j=0;j<size/2;j+=2) {
    a = x[j]*x[j]+x[j+1]*x[j+1];
    F += a*j*dF;                       // central frequency
    E += a;                            // energy
  }
  F = E>0 ? F/E : 0.5*rate;

  return F;
}

//______________________________________________________________________________
double
CWB::mdc::GetTimeRange(wavearray<double> x, double& tMin, double& tMax, double efraction) {
//
// Get mdc time interval 
//
//
// Input: x      - wavearray which contains the waveform data
//
// Output: tMin  - start time interval
//         tMax  - stop  time interval
//

  if(efraction<0) efraction=0;
  if(efraction>1) efraction=1;

  int N = x.size();

  double E = 0;                                                 // signal energy
  double avr = 0;                                               // average
  for(int i=0;i<N;i++) {avr+=i*x[i]*x[i]; E+=x[i]*x[i];}
  int M=int(avr/E);                                             // central index

  // search range which contains percentage P of the total energy E
  int jB=0;
  int jE=N-1;
  double a,b;
  double sum = ((M>=0)&&(M<N)) ? x[M]*x[M] : 0.;
  for(int j=1; j<N; j++) {
    a = ((M-j>=0)&&(M-j<N)) ? x[M-j] : 0.;
    b = ((M+j>=0)&&(M+j<N)) ? x[M+j] : 0.;
    if(a) jB=M-j;
    if(b) jE=M+j;
    sum += a*a+b*b;
    if(sum/E > efraction) break;
  }

  tMin = jB/x.rate();
  tMax = jE/x.rate();

  return tMax-tMin;
}

//______________________________________________________________________________
void
CWB::mdc::TimeShift(wavearray<double>& x, double tShift) {
//
// apply time shift
//
//
// Input: x      - wavearray which contains the waveform data
//        time   - time shift (sec)
//

  if(tShift==0) return;

  // search begin,end of non zero data
  int ibeg=0; int iend=0;
  for(int i=0;i<(int)x.size();i++) {
    if(x[i]!=0 && ibeg==0) ibeg=i;
    if(x[i]!=0) iend=i;
  }
  int ilen=iend-ibeg+1;
  // create temporary array for FFTW & add scratch buffer + tShift
  int ishift = fabs(tShift)*x.rate();
  int isize = 2*ilen+2*ishift;
  isize = isize + (isize%4 ? 4 - isize%4 : 0); // force to be multiple of 4
  wavearray<double> w(isize);
  w.rate(x.rate()); w=0;
  // copy x data !=0 in the middle of w array & set x=0
  for(int i=0;i<ilen;i++) {w[i+ishift+ilen/2]=x[ibeg+i];x[ibeg+i]=0;}

  double pi = TMath::Pi();
  // apply time shift to waveform vector
  w.FFTW(1);
  TComplex C;
  double df = w.rate()/w.size();
  //cout << "tShift : " << tShift << endl;
  for (int ii=0;ii<(int)w.size()/2;ii++) {
    TComplex X(w[2*ii],w[2*ii+1]);
    X=X*C.Exp(TComplex(0.,-2*pi*ii*df*tShift));  // Time Shift
    w[2*ii]=X.Re();
    w[2*ii+1]=X.Im();
  }
  w.FFTW(-1);

  // copy shifted data to input x array
  for(int i=0;i<(int)w.size();i++) {
    int j=ibeg-(ishift+ilen/2)+i;
    if((j>=0)&&(j<(int)x.size())) x[j]=w[i];
  }

  return;
}

//______________________________________________________________________________
void
CWB::mdc::PhaseShift(wavearray<double>& x, double pShift) {
//
// apply phase shift
//
//
// Input: x      - wavearray which contains the waveform data
//        pShift - phase shift (degrees)
//

  if(pShift==0) return;

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

  // apply phase shift to waveform vector
  w.FFTW(1);
  TComplex C;
  //cout << "pShift : " << pShift << endl;
  pShift*=TMath::Pi()/180.;   
  for (int ii=0;ii<(int)w.size()/2;ii++) {
    TComplex X(w[2*ii],w[2*ii+1]);
    X=X*C.Exp(TComplex(0.,-pShift));  // Phase Shift
    w[2*ii]=X.Re();
    w[2*ii+1]=X.Im();
  }
  w.FFTW(-1);

  // copy shifted data to input x array
  for(int i=0;i<(int)w.size();i++) {
    int j=ibeg-isize/4+i;
    if((j>=0)&&(j<(int)x.size())) x[j]=w[i];
  }

  return;
}

//______________________________________________________________________________
wavearray<double>
CWB::mdc::GetSGQ(double frequency, double Q) {
//
// Get SinGaussian 
//
//
// Input: frequency - SG frequency (Hz)
//        Q         - SG Q factor  
//
// Return wavearray data containing the waveform
//

  wavearray<double> x(int(inj_length*MDC_SAMPLE_RATE));
  x.rate(MDC_SAMPLE_RATE);
  x=0;
  double amplitude=1.;
  double duration = Q/(TMath::TwoPi()*frequency);
  AddSGBurst(x, amplitude, frequency, duration,0);

  return x;
}

//______________________________________________________________________________
wavearray<double>
CWB::mdc::GetCGQ(double frequency, double Q) {
//
// Get CosGaussian 
//
//
// Input: frequency - CG frequency (Hz)
//        Q         - CG Q factor  
//
// Return wavearray data containing the waveform
//

  wavearray<double> x(int(inj_length*MDC_SAMPLE_RATE));
  x.rate(MDC_SAMPLE_RATE);
  x=0;
  double amplitude=1.;
  double duration = Q/(TMath::TwoPi()*frequency);
  AddCGBurst(x, amplitude, frequency, duration,0);

  return x;
}

#define GA_FORMULA_WNB "[0]*TMath::Exp(-TMath::Power((x-[2])/[1],2)/2)"

//______________________________________________________________________________
wavearray<double>
CWB::mdc::GetWNB(double frequency, double bandwidth, double duration, int seed, bool mode) {
//
// Get White Band Gaussian Noise 
//
//
// Input: frequency - start white band frequency (Hz)
//        bandwidth - (Hz)
//        duration  - (sec)
//        mode      - if mode=1 -> apply Patrik Sutton Method
//                                 simmetric respect to the central frequency
//                    if mode=0 -> asymmetric respect to the central frequency   
//
// Return wavearray data containing the waveform
//

  gRandom->SetSeed(seed);

  // Generate white gaussian noise 1 sec
  wavearray<double> x(int(inj_length*MDC_SAMPLE_RATE));
  x.rate(MDC_SAMPLE_RATE);
  for (int i=0;i<(int)x.size();i++) x[i]=gRandom->Gaus(0,1);

  double dt = 1./x.rate();
  double df = x.rate()/x.size();

  // Apply a band limited cut in frequency
  x.FFTW(1);
  if(mode==1) {
    // set zero bins outside the range [0,bandwidth/2] 
    // Heterodyne up by frequency+bandwidth/2 (move zero frequency to center of desidered band)
    int bFrequency = frequency/df;
    int bBandwidth = (bandwidth/2.)/df;
    wavearray<double> y=x; y=0;
    int bin = 2*(bFrequency+bBandwidth);
    y[bin]=x[0];
    y[bin+1]=x[1];
    for (int i=1;i<bBandwidth;i++) {
      y[bin+2*i]=x[2*i];
      y[bin+2*i+1]=x[2*i+1];
  
      y[bin-2*i]=x[2*i];
      y[bin-2*i+1]=-x[2*i+1];
    }
    x=y;;
    x.FFTW(-1);
  } else {
    // set zero bins outside the range [frequency,frequency+bandwidth] 
    int bLow  = frequency/df;
    int bHigh = (frequency+bandwidth)/df;
    for (int i=0;i<bLow;i++) {x[2*i]=0;x[2*i+1]=0;}
    for (int i=bHigh;i<(int)x.size()/2;i++) {x[2*i]=0;x[2*i+1]=0;}
    x.FFTW(-1);
  }

  // Apply a gaussian shape in time
  double function_range = 1.;  // duration 1 sec
  TF1* ga_function = new TF1("Gaussian",GA_FORMULA_WNB,-function_range/2,function_range/2);
  ga_function->SetParameter(0,1);
  ga_function->SetParameter(1,duration);
  ga_function->SetParameter(2,function_range/2);
  for (int i=0;i<(int)x.size();i++) x[i]*=ga_function->Eval(dt*(i+1));

  // normalization 
  double hrss=0;
  for (int i=0;i<(int)x.size();i++) hrss+=x[i]*x[i];
  hrss=sqrt(hrss*dt);
  for (int i=0;i<(int)x.size();i++) x[i]*=(1./sqrt(2.))*1./hrss;

  delete ga_function;

  return x;
}

#define GA_FORMULA "[0]*TMath::Exp(-TMath::Power((x-[2])/[1],2))"

//______________________________________________________________________________
wavearray<double>  
CWB::mdc::GetGA(double duration) {
//
// Get Gaussian signal 
//
//
// Input: duration     - (sec)
//
// Return wavearray data containing the waveform

  // Apply a gaussian shape in time
  double function_range = 1.;  // duration 1 sec
  TF1* ga_function = new TF1("Gaussian",GA_FORMULA,-function_range/2,function_range/2);
  ga_function->SetParameter(0,1);
  ga_function->SetParameter(1,duration);
  ga_function->SetParameter(2,function_range/2);

  double dt = 1./MDC_SAMPLE_RATE;

  // plus component
  wavearray<double> x(int(inj_length*MDC_SAMPLE_RATE));
  x.rate(MDC_SAMPLE_RATE);
  for (int i=0;i<(int)x.size();i++) x[i]=ga_function->Eval(dt*(i+1));

  // normalization 
  double hrss=0;
  for (int i=0;i<(int)x.size();i++) hrss+=x[i]*x[i];
  hrss=sqrt(hrss*dt);
  for (int i=0;i<(int)x.size();i++) x[i]*=1./hrss;

  delete ga_function;

  return x;
}

#define RD_FORMULA "(((x-1./4./[2]+TMath::Abs(x-1./4./[2]))/2./(x-1./4./[2]))*[0]*TMath::Cos(TMath::TwoPi()*[2]*x)+[1]*TMath::Sin(TMath::TwoPi()*[2]*x))*TMath::Exp(-x/[3])" // Heaviside in cos like Andrea Vicere'

//______________________________________________________________________________
wavearray<double>  
CWB::mdc::GetRD(double frequency, double tau, double iota, bool polarization) {
//
// Get RingDown signal 
//
//
// Input: frequency    - (Hz)
//        tau          - (sec)
//        iota         - (degrees)
//        polarization - (degrees)
//
// Return wavearray data containing the waveform
//

  iota*=TMath::Pi()/180.;

  double c2 = TMath::Cos(iota);
  double c1 = (1+c2*c2)/2.;

  if (c1/c2<1e-10) c1=0;

  char rd_formula[256] = RD_FORMULA;
  double function_range = 1.;  // duration 1 sec
  TF1* rd_function = new TF1("RingDown",rd_formula,0.,function_range);
  rd_function->SetParameter(2,frequency);
  rd_function->SetParameter(3,tau);

  double time_offset=0.05;   // it is necessary to allows negative time shift 
  double dt = 1./MDC_SAMPLE_RATE;

  // plus component
  wavearray<double> x(int(inj_length*MDC_SAMPLE_RATE));
  x.rate(MDC_SAMPLE_RATE);
  rd_function->SetParameter(0,c1);
  rd_function->SetParameter(1,0);
  for (int i=0;i<(int)x.size();i++) {
    double time = dt*(i+1)-time_offset;
    if (time>=0) x[i]=rd_function->Eval(time); else x[i]=0;
  }

  // cross component
  wavearray<double> y(int(inj_length*MDC_SAMPLE_RATE));
  y.rate(MDC_SAMPLE_RATE);
  rd_function->SetParameter(0,0);
  rd_function->SetParameter(1,c2);
  for (int i=0;i<(int)y.size();i++) {
    double time = dt*(i+1)-time_offset;
    if (time>=0) y[i]=rd_function->Eval(time); else y[i]=0;
  }

  double hrss=0;
  for (int i=0;i<(int)x.size();i++) hrss+=x[i]*x[i]+y[i]*y[i];
  hrss=sqrt(hrss*dt);
  for (int i=0;i<(int)x.size();i++) x[i]*=1./hrss;
  for (int i=0;i<(int)y.size();i++) y[i]*=1./hrss;

  delete rd_function;

  return polarization==0 ? x : y;
}

//______________________________________________________________________________
void
CWB::mdc::AddGauss(wavearray<double> &td, double v, double u) {
//
// Add Gaussian Noise  
//
//
// Input: td       - wavearray with input data signal
//        v        - sigma of gaussian noise
//        u        - median of gaussian noise
//
// Output: td      - Return input signal + gaussian noise
//

  int n=td.size();
  for (int i=0; i < n; i++){
     td.data[i] += v*gRandom->Gaus(0.,1.)+u;
  }
}

//______________________________________________________________________________
void
CWB::mdc::AddExp(wavearray<double> &td, double v, int M) {
//
// Add exponential noise  
//
//
// Input: td       - wavearray with input data signal
//        v        - 1/tau of exponential
//        M        - number of random values added for each sample
//
// Output: td      - Return input signal + exponential noise
//

  int i,j;
  double x,y;
  int m=abs(M);
  int n=td.size();
  for (i=0; i<n; i++){
     y = 0.;
     for (j=0; j<m; j++){
        x = gRandom->Exp(v);
        if(M>0) y += x;
        else if(x>y) y=x;
     }
     td.data[i] += y;
  }
}

//______________________________________________________________________________
void
CWB::mdc::AddSGBurst(wavearray<double> &td, double a, double f, double s, double d) {
//
// Add SinGaussian burst  
//
//
// Input: td       - wavearray with input data signal
//        a        - amplitude 
//        f        - frequency (Hz)
//        s        - gaussian RMS
//        d        - duration (sec)
//
// Output: td      - Return input signal + SinGaussian burst
//

  int n=td.size();
  double r = td.rate();
  int m;
  double g, t;
  double sum = 0.;
  double delay = d;

  m = int(6*s*r);
  if(m > n/2-1) m = n/2-2;

  for (int i=0; i < m; i++){
     t = i/r;
     g = 2*TMath::Exp(-t*t/2/s/s)*TMath::Sin(2*PI*f*t);
     sum += g*g;
  }
  a *= TMath::Sqrt(r/sum);

  td.data[n/2+int(delay*r)] += 0;
  for (int i=1; i < m; i++){
     t = i/r;
     g = a*TMath::Exp(-t*t/2/s/s)*TMath::Sin(2*PI*f*t);
     td.data[n/2+i+int(delay*r)] += g;
     td.data[n/2-i+int(delay*r)] -= g;
  }
}

//______________________________________________________________________________
void
CWB::mdc::AddCGBurst(wavearray<double> &td, double a, double f, double s, double d) {
//
// Add CosGaussian burst  
//
//
// Input: td       - wavearray with input data signal
//        a        - amplitude 
//        f        - frequency (Hz)
//        s        - gaussian RMS
//        d        - duration (sec)
//
// Output: td      - Return input signal + CosGaussian burst
//

  int n=td.size();
  double r = td.rate();
  int m;
  double g, t;
  double sum = 0.;
  double delay = d;

  m = int(6*s*r);
  if(m > n/2-1) m = n/2-2;

  for (int i=0; i < m; i++){
     t = i/r;
     g = 2*TMath::Exp(-t*t/2/s/s)*TMath::Cos(2*PI*f*t);
     sum += g*g;
  }
  a *= TMath::Sqrt(r/sum);

  td.data[n/2+int(delay*r)] += a;
  for (int i=1; i < m; i++){
     t = i/r;
     g = a*TMath::Exp(-t*t/2/s/s)*TMath::Cos(2*PI*f*t);
     td.data[n/2+i+int(delay*r)] += g;
     td.data[n/2-i+int(delay*r)] += g;
  }
}

//______________________________________________________________________________
void
CWB::mdc::AddWGNoise(wavearray<double> &td, double a, double s) {
//
// Add windowed Gaussian noise 
//
//
// Input: td       - wavearray with input data signal
//        a        - amplitude 
//        s        - gaussian RMS
//
// Output: td      - Return input signal + windowed Gaussian noise
//

  int n=td.size();
  double r = td.rate();
  int m;
  double g, t;
  double sum = 0.;

  m = int(3*s*r);
  if(m > n/2-1) m = n/2-2;
  wavearray<double> gn(2*m);


  for (int i=0; i < m; i++){
     t = i/r;
     g   = TMath::Exp(-t*t/2/s/s);
     gn.data[m+i]   = g*gRandom->Gaus(0.,1.);
     gn.data[m-i-1] = g*gRandom->Gaus(0.,1.);
     sum += gn.data[m+i]*gn.data[m+i] + gn.data[m-i-1]*gn.data[m-i-1];
  }
  a *= TMath::Sqrt(r/sum);

  for (int i=0; i < m; i++){
     t = i/r;
     td.data[n/2+i] += a*gn.data[m+i];
     td.data[n/2-i] += a*gn.data[m-i-1];
  }
}

// --------------------------------------------------------
// Miyamoto-Nagai Galactic Disk Model
// www.astro.utu.fi/~cflynn/galdyn/lecture4.html
// a=[0]
// b=[1]  scale length
// R=x, z=y
// --------------------------------------------------------

#define MNGD_NUMERATOR   "([0]*((x-[3])*(x-[3])+y*y)+([0]+3*sqrt(z*z+[1]*[1]))*pow([0]+sqrt(z*z+[1]*[1]),2))"
#define MNGD_DENOMINATOR "(pow(((x-[3])*(x-[3])+y*y)+pow([0]+sqrt(z*z+[1]*[1]),2),5./2.)*pow(z*z+[1]*[1],3./2.))"

#define MNGD_B    0.3        // Kpc
#define MNGD_Md1  6.6e10     // Msol
#define MNGD_A1   5.81       // Kpc
#define MNGD_Md2  -2.9e10    // Msol
#define MNGD_A2   17.43      // Kpc
#define MNGD_Md3  3.3e9      // Msol
#define MNGD_A3   34.86      // Kpc

#define MNGD_XMAX 40
#define MNGD_YMAX 4

#define MNGD_SOLAR_SISTEM_DISTANCE_FROM_GC 7.62       // 7.62 [+/-0.32] Kpc it.wikipedia.org/wiki/Via_Lattea


//______________________________________________________________________________
void 
CWB::mdc::SetSkyDistribution(MDC_DISTRIBUTION sky_distribution, TString fName, int seed, bool add) {
//
// Set Sky Distribution
//
// see SetSkyDistribution(MDC_DISTRIBUTION sky_distribution, TString fName, vector<mdcpar> par, int seed)
//
//

  vector<mdcpar> par;
  SetSkyDistribution(sky_distribution, fName, par, seed, add);
  return;
}

//______________________________________________________________________________
void 
CWB::mdc::SetSkyDistribution(MDC_DISTRIBUTION sky_distribution, vector<mdcpar> par, int seed, bool add) {
//
// Set Sky Distribution
//
// see SetSkyDistribution(MDC_DISTRIBUTION sky_distribution, TString fName, vector<mdcpar> par, int seed)
//

  SetSkyDistribution(sky_distribution, "", par, seed, add);
  return;
}

//______________________________________________________________________________
void 
CWB::mdc::SetSkyDistribution(MDC_DISTRIBUTION sky_distribution, TString fName, vector<mdcpar> par, int seed, bool add) {
//
// Set Sky Distribution
//
//
// Input: sky_distribution - sky distribution type  
//        fName            - input file name 
//        par              - input sky distribution parameters
//        seed             - seed for random selection
//        add              - add (def=false) if true the distribution is added to the current one
//
//                   sky_distribution & par & fName can be one of the following choices :  
//
//                   MDC_RANDOM        = random sky distribution (fName not used)
//                                       fName can be used to define a custom rho distribution  
//                                             example : fName="pow(x,2)"  
//                                             do not use  par[n].name="rho_dist";                           
//                                       vector<mdcpar> par(1);
//                                       par[0].name="entries"; par[0].value=XXX; // number of random events 
//                                       par.resize(3);
//                                       par[1].name="rho_min"; par[1].value=MIN; // min distance in kPc
//                                       par[2].name="rho_max"; par[2].value=MAX; // max distance in kPc
//                                       // if(rho_min>0 && rho_max>=rho_min)
//                                       // the amplitude is rescaled to 10/rho (10 Kpc is the stantard candle distance)
//                                       // if rho_min,rho_max are not defined the amplitude is not rescaled
//                                       // 
//                                       par[n].name="rho_dist";     // define the rho distribution 
//                                       par[n].value=x;             // dist is rho^x
//                                       // if "rho_dist" is not defined then the distance of events 
//                                       // is randomly selected from a (rho*rho) distribution   
//                                       // so the sky distribution is uniform in volume
//
//                                       par[m].name="iota"; par[m].value=iota;   // ellipticity [0:180] deg
//                                                                                // if<0 || >180 -> iota=random 
//                                                                                // default : iota=0
//                                       par[p].name="psi";  par[p].value=psi;    // polarization (deg)
//                                                                                // if not defined -> random 
//
//                   MDC_EARTH_FIX     = earth fixed greographical coordinates (latitude, longitude)
//                   MDC_CELESTIAL_FIX = celestial fixed coordinates (DEC, RA)
//                                       vector<mdcpar> par(2);
//                                       par[0].name="theta"; par[0].value=XXX; // degrees
//                                       par[1].name="phi";   par[1].value=YYY; // degrees
//                                       par.resize(3);
//                                       par[2].name="psi";   par[2].value=ZZZ; // degrees
//                                       par.resize(4);
//                                       par[3].name="rho";   par[3].value=UUU; // kPc
//                                       par.resize(5);
//                                       par[4].name="gps";   par[4].value=VVV; // sec
//                                       par.resize(6);
//                                       par[5].name="entries"; par[5].value=TTT; // number of random events 
//                                       par.resize(7);
//                                       par[6].name="iota"; par[6].value=iota;   // ellipticity [0:180] deg
//                                                                                // if<0 || >180 -> random 
//
//                   MDC_MNGD          = Miyamoto-Nagai Galactic Disk Model (fName,par not used)
//                                       www.astro.utu.fi/~cflynn/galdyn/lecture4.html
//                                       vector<mdcpar> par(1);
//                                       par[0].name="entries"; par[0].value=XXX; // number of random events 
//                                       par.resize(2);
//                                       par[1].name="iota"; par[1].value=iota;   // ellipticity [0:180] deg
//                                                                                // if<0 || >180 -> random 
//
//                   MDC_GWGC          = Gravitation Wave Galaxy Catalog (par not used)
//                                       fName = path name of galaxy catalog    
//                                       vector<mdcpar> par(1);
//                                       par[0].name="distance_thr"; par[0].value=XXX; // distance max (Kpc) 
//                                       par.resize(2);
//                                       par[1].name="iota"; par[1].value=iota;   // ellipticity [0:180] deg
//                                                                                // if<0 || >180 -> random 
//
//                   MDC_XMLFILE       = mdc xml file
//                                       fName = path name of mdc xml file (LAL xml format)
//                                       par[0].name="start"; par[0].value=XXX; // start GPS time to search in XML file
//                                       par[1].name="stop";  par[1].value=YYY; // stop  GPS time to search in XML file
//                                       par[2].name="hrss_factor";  par[2].value=ZZZ; // used to rescale hrss (default=1)
//                                       par[3].name="decimals";  par[3].value=WWW; // used to format the MDC name (default=1) 
//                                       par[4].name="engine";par[4].value=V;   // mdc engine : 0->AAL, 1->CWB (default=1)
//                                       par[5].name="auto";par[5].value=T;     // 0/1->disaable/enable (default=0)
//                                       NOTE : if auto=0 the MDC names must be provided by user using the following setup: 
//                                       par[6].name="type";  par[6].svalue="MDC_NAME1"; // MDC name1 (with CWB Format)
//                                       par[X].name="type";  par[X].svalue="MDC_NAMEX"; // MDC nameX (with CWB Format)
//                                       par[N].name="type";  par[N].svalue="MDC_NAMEN"; // MDC nameN (with CWB Format)
//
//                   MDC_LOGFILE       = mdc log file
//                                       fName = path name of mdc log file 
//
//                   MDC_CUSTOM        = user define distribution
//                                       fName = path name of user file 
//					 file format :  gps  name  theta  phi  psi  rho  iota hrss ID id
//					 theta = [0:180] deg : phi = [0:360] deg : psi = [0:180] deg : iota = [0:180] deg
//					 ID/id are the major/minor id waveform numbers 
//					       if ID=-1 then 'name' is used, if id=1 then id is selected randomly
//
// Note1: the angle iota is the inclination of the system which originates the burst with respect to the line of sight
//        iota = 0  : the line of sight is perpendicular to the orbit plane 
//        iota = 90 : the line of sight has 0 deg inclination angle respect to the orbit plane 
//        the h+ ellipticity factor is : (1+cos[iota]*cos[iota])/2
//        the hx ellipticity factor is : cos[iota]
//
// Note2: rho defines the distance of the source in Kpc
//        if rho=0 then the hrss is the one defined with SetInjHrss (default)
//        if rho>0 then the hrss is rescaled by a factor 10/rho;
//                 the hrss defined with SetInjHrss is the hrss at 10Kpc  
//

  gRandom->SetSeed(seed);

  this->sky_distribution=sky_distribution;
  this->sky_file=fName;
  this->sky_parms=par;

  double rad2deg = 180./TMath::Pi();

  if(!add) {
   nameList.clear();
     thList.clear();
     phList.clear();
    psiList.clear();
    rhoList.clear();
   iotaList.clear();
   hrssList.clear();
    gpsList.clear();
     IDList.clear();
     idList.clear();
  }

  if(sky_distribution==MDC_RANDOM) {

    bool error=false;
    int entries = GetPar("entries",par,error);
    double rho_min = GetPar("rho_min",par,error);
    if(error) rho_min=0.;
    double rho_max = GetPar("rho_max",par,error);
    if(error) rho_max=0.;

    if(entries<=0) {
      cout << "CWB::mdc::SetSkyDistribution - "
           << "Error : entries must be positive" << endl;
      exit(1);
    }
    if(rho_max>0 && rho_min<=0) {
      cout << "CWB::mdc::SetSkyDistribution - " << "Error : rho_min must be > 0" << endl;
      exit(1);
    }
    if(rho_min>rho_max) {
      cout << "CWB::mdc::SetSkyDistribution - " << "Error : rho_min must be <= rho_max" << endl;
      exit(1);
    }

    cout << "CWB::mdc::SetSkyDistribution - All Sky random distribution" << endl;

    TF1 *rd = NULL;
    // define the rho distribution 
    if(rho_min!=rho_max) {
      double rho_dist = GetPar("rho_dist",par,error);
      char rho_dist_func[256]; 
           if((!error)&&(fName=="")) sprintf(rho_dist_func,"pow(x,%f)",rho_dist);  // custom pow(rho,rho_dist)
      else if(( error)&&(fName=="")) sprintf(rho_dist_func,"pow(x,2)");	           // default distribution
      else if(( error)&&(fName!="")) strcpy(rho_dist_func,fName.Data());	   // custom formula
      else if((!error)&&(fName!="")) {						   // ambiguoius definition
             cout << "CWB::mdc::SetSkyDistribution - " 
                  << "Error : ambiguous rho distribution - is defined twice " 
                  << " 1) as rho_dist params, 2) as formula" << endl;
             exit(1);
           }  
      rd = new TF1("rd",rho_dist_func,rho_min,rho_max);			      // sky distribution rho^rho_dist
      TF1* rdcheck = (TF1*)gROOT->GetListOfFunctions()->FindObject("rd");     // check if formula is correct
      if(rdcheck==NULL) {
         cout << "CWB::mdc::SetSkyDistribution - " 
              << "Error : wrong formula format : " << rho_dist_func << endl; 
         exit(1);
      }	
    }

    for(int n=0;n<entries;n++) {
      thList.push_back(rad2deg*acos(gRandom->Uniform(-1,1))-90.);
      phList.push_back(gRandom->Uniform(-180,180));
      if(rd!=NULL) rhoList.push_back(rd->GetRandom());
      else         rhoList.push_back(rho_min);

      error=false;double value=GetPar("psi",par,error);
      double psi = error ? gRandom->Uniform(0,180) : value;
      psiList.push_back(psi);

      error=false;double iota = GetPar("iota",par,error);
      if(error)                   iotaList.push_back(0);
      else if(iota>=0&&iota<=180) iotaList.push_back(iota);   
      else                        iotaList.push_back(rad2deg*acos(gRandom->Uniform(-1,1)));

      hrssList.push_back(0);
      IDList.push_back(-1);
      idList.push_back(-1);
    }

    if(rd!=NULL) delete rd;

  } else
  if(sky_distribution==MDC_MNGD) {

    bool error=false;
    double entries = GetPar("entries",par,error);

    if(error) {
      cout << "CWB::mdc::SetSkyDistribution - "
           << "Error : num par must be 1 [entries]" << endl;
      exit(1);
    }

    cout << "CWB::mdc::SetSkyDistribution - the Miyamoto-Nagai Galactic Disk Model" << endl;

    char formula[256];
    sprintf(formula,"[1]*[1]*[2]*%s/%s",MNGD_NUMERATOR,MNGD_DENOMINATOR);

    TF3 *gd1 = new TF3("gd1",formula,-MNGD_XMAX,MNGD_XMAX,-MNGD_XMAX,MNGD_XMAX,-MNGD_YMAX,MNGD_YMAX);
    gd1->SetParameter(0,MNGD_A1);      // a
    gd1->SetParameter(1,MNGD_B);       // b
    gd1->SetParameter(2,MNGD_Md1);     // b
    gd1->SetParameter(3,MNGD_SOLAR_SISTEM_DISTANCE_FROM_GC);

    TF3 *gd2 = new TF3("gd2",formula,-MNGD_XMAX,MNGD_XMAX,-MNGD_XMAX,MNGD_XMAX,-MNGD_YMAX,MNGD_YMAX);
    gd2->SetParameter(0,MNGD_A2);      // a
    gd2->SetParameter(1,MNGD_B);       // b
    gd2->SetParameter(2,MNGD_Md2);     // b
    gd2->SetParameter(3,MNGD_SOLAR_SISTEM_DISTANCE_FROM_GC);

    TF3 *gd3 = new TF3("gd3",formula,-MNGD_XMAX,MNGD_XMAX,-MNGD_XMAX,MNGD_XMAX,-MNGD_YMAX,MNGD_YMAX);
    gd3->SetParameter(0,MNGD_A3);      // a
    gd3->SetParameter(1,MNGD_B);       // b
    gd3->SetParameter(2,MNGD_Md3);     // b
    gd3->SetParameter(3,MNGD_SOLAR_SISTEM_DISTANCE_FROM_GC);

    TF3* gd = new TF3("gd","gd1+gd2+gd3",-MNGD_XMAX,MNGD_XMAX,-MNGD_XMAX,MNGD_XMAX,-MNGD_YMAX,MNGD_YMAX);

    gd->SetNpx(100);
    gd->SetNpy(100);
    gd->SetNpz(100);

    // Generate randomly sources from the Gatactic Disk

    XYZVector xyz;
    double xgc=MNGD_SOLAR_SISTEM_DISTANCE_FROM_GC,ygc=0.,zgc=0.;
    xyz.SetXYZ(xgc,ygc,zgc);

    double ilongitude = xyz.Phi()*rad2deg;
    double ilatitude = -(xyz.Theta()-TMath::Pi()/2.)*rad2deg;
    double olongitude,olatitude;
    GalacticToEquatorial(ilongitude, ilatitude, olongitude, olatitude);
    double gc_phi = olongitude;
    double gc_theta = olatitude;
    double gc_rho = sqrt(xyz.mag2());

    cout << "gc_phi : " << gc_phi << " gc_theta : " << gc_theta << " " << gc_rho << endl;

    double x,y,z;
    for(int n=0;n<entries;n++) {

      gd->GetRandom3(x,y,z);
      xyz.SetXYZ(x,y,z);

      double ilongitude  = xyz.Phi()*rad2deg;
      double ilatitude   = -(xyz.Theta()-TMath::Pi()/2.)*rad2deg;
      double olongitude,olatitude;
      GalacticToEquatorial(ilongitude, ilatitude, olongitude, olatitude);
      phList.push_back(olongitude);
      thList.push_back(olatitude);
      psiList.push_back(gRandom->Uniform(0,180));
      rhoList.push_back(sqrt(xyz.mag2()));

      error=false;double iota = GetPar("iota",par,error);
      if(error)                   iotaList.push_back(0);
      else if(iota>=0&&iota<=180) iotaList.push_back(iota);   
      else                        iotaList.push_back(rad2deg*acos(gRandom->Uniform(-1,1)));

      hrssList.push_back(0);
      IDList.push_back(-1);
      idList.push_back(-1);
    }

    delete gd1;
    delete gd2;
    delete gd3;
    delete gd;

  } else
  if(sky_distribution==MDC_GWGC) {

    bool error=false;
    double distance_thr = GetPar("distance_thr",par,error);

    if(error) {
      cout << "CWB::mdc::SetSkyDistribution - "
           << "Error : num par must be 1 [distance_thr]" << endl;
      exit(1);
    }

    cout << "CWB::mdc::SetSkyDistribution - Gravitational Wave Galaxy Catalog" << endl;
    cout << "CWB::mdc::SetSkyDistribution - Distance Threshold " << distance_thr << " Kpc" << endl;

    ifstream in;
    in.open(fName,ios::in);
    if (!in.good()) {cout << "CWB::mdc::SetSkyDistribution - Error Opening File : " << fName << endl;exit(1);}

    // get number of entries
    int entries=0;
    char str[1024];
    while(true) {
      in.getline(str,1024);
      if (!in.good()) break;
      if(str[0] != '#') entries++;
    }
    cout << "entries " << entries << endl;
    in.clear(ios::goodbit);
    in.seekg(0, ios::beg);

    char iline[1024];
    in.getline(iline,1024);  // skip first line (header)
    while (1) {

      in.getline(iline,1024);
      if (!in.good()) break;
      TObjArray* tok = TString(iline).Tokenize(TString('|'));

      TObjString* tra   = (TObjString*)tok->At(2);
      TObjString* tdec  = (TObjString*)tok->At(3);
      TObjString* tdist = (TObjString*)tok->At(14);

      delete tok;

      double ra   = tra->GetString().Atof();
      double dec  = tdec->GetString().Atof();
      double dist = tdist->GetString().Atof();  // Mpc

      dist *= 1000;  // Mpc -> Kpc 
      if (dist<distance_thr) {
        thList.push_back(dec);
        phList.push_back(ra*360./24.);
        psiList.push_back(gRandom->Uniform(0,180));
        rhoList.push_back(dist);

        error=false;double iota = GetPar("iota",par,error);
        if(error)                   iotaList.push_back(0);
        else if(iota>=0&&iota<=180) iotaList.push_back(iota);   
        else                        iotaList.push_back(rad2deg*acos(gRandom->Uniform(-1,1)));

        hrssList.push_back(0);
        IDList.push_back(-1);
        idList.push_back(-1);
      }
    }

  } else
  if(sky_distribution==MDC_CUSTOM) {

    cout << "CWB::mdc::SetSkyDistribution - User Custom Distribution" << endl;

    ifstream in;
    in.open(fName,ios::in);
    if (!in.good()) {cout << "CWB::mdc::SetSkyDistribution - Error Opening File : " << fName << endl;exit(1);}

    // get number of entries
    int entries=0;
    char str[1024];
    while(true) {
      in.getline(str,1024);
      if (!in.good()) break;
      if(str[0] != '#') entries++;
    }
    cout << "entries " << entries << endl;
    in.clear(ios::goodbit);
    in.seekg(0, ios::beg);

    char name[128];
    double gps,psi,rho,iota,hrss;  
    double theta,phi;         // geographic coordinates
    int ID,id;
    int fpos=0;
    while (1) {
      fpos=in.tellg();
      in.getline(str,1024);
      if(str[0] == '#') continue;
      if (!in.good()) break;

      std::stringstream linestream(str);
      if(!(linestream >> gps >> name >> theta >> phi >> psi >> rho >> iota >> hrss >> ID >> id)) {
         cout << "CWB::mdc::SetSkyDistribution - Wrong Format for File : " << fName << endl;
         cout << "input line : " << endl;
         cout << str << endl;
         cout << "must be : " << endl;
         cout << "gps " << "name " <<  "theta " <<  "phi " <<  "psi " 
              << "rho " << "iota " << "hrss " << "ID " << "id " << "..." << endl;
         exit(1);
      } 
      //cout << gps << " " <<  name << " " <<  theta << " " <<  phi << " " <<  psi 
      //     << " " <<  rho << " " << iota << " " << hrss << " " << ID << " " << id << endl;

       gpsList.push_back(gps);
      nameList.push_back(name);
        thList.push_back(theta);
        phList.push_back(phi);
       psiList.push_back(psi);
       rhoList.push_back(rho);
      iotaList.push_back(iota);
      hrssList.push_back(hrss);
        IDList.push_back(ID);
        idList.push_back(id);
    }

    in.close();

  } else
  if(sky_distribution==MDC_XMLFILE) {

#ifdef _USE_LAL
    cout << "CWB::mdc::SetSkyDistribution - User Distribution from XML" << endl;

#if LAL_VERSION_MAJOR >   6 || (LAL_VERSION_MAJOR ==  6 && \
   (LAL_VERSION_MINOR >  14 || (LAL_VERSION_MINOR == 14 && \
    LAL_VERSION_MICRO >=  0                             )))     // LAL_VERSION >= 6.14.0
#else
    cout << "CWB::mdc::SetSkyDistribution - LAL XML File can not be used with LAL ver < 6.14.0" << endl; 
    exit(1);
#endif

    Long_t id,fsize,flags,mt;
    int estat = gSystem->GetPathInfo(fName.Data(),&id,&fsize,&flags,&mt);
    if (estat!=0) {
      cout << "CWB::mdc::SetSkyDistribution - XML File : " << fName.Data() << " Not Exist" << endl;
      exit(1);
    }

    bool error;

    // get decimals
    error=false;
    int decimals = GetPar("decimals",par,error);
    if(error) decimals=1;

    // get hrss_factor : is a user factor used to rescale hrss
    error=false;
    double hrss_factor = GetPar("hrss_factor",par,error);
    if(error) hrss_factor=1;

    // get mdc engine : aal=0, cwb=1 (default)
    error=false;
    int engine = GetPar("engine",par,error);
    if(error) engine=1;

    // get auto mode for MDC type : if auto=0 the MDC names must be provided by user
    error=false;
    int AUTO = GetPar("auto",par,error);
    if(error) AUTO=0;

    // get time range for injections
    error=false;
    double start = GetPar("start",par,error);
    if(error) start=0;
    error=false;
    double stop = GetPar("stop",par,error);
    if(error) stop=std::numeric_limits<double>::max();;

    xmlType.clear();
    if(!AUTO) { 	// get user defined xml types
      for(int n=0;n<par.size();n++) {
        error=false;
        vector<mdcpar> xpar(1); xpar[0]=par[n];
        TString type = GetParString("type",xpar,error);
        if(!error) { 	// add to xmlType if type is not present in the xmlType list 	
          bool present=false;
          for(int m=0;m<xmlType.size();m++) {
            if(type==xmlType[m]) {present=true;break;}
          }
          if(!present) xmlType.push_back(type.Data());
        }
      }
    }

    // read sim_burst table from XML file (see LIGOLwXMLBurstRead.c)
    LIGOTimeGPS gpsStartTime = {(int)start, 0};
    LIGOTimeGPS gpsStopTime  = {(int)stop, 0};
    SimBurst* sim_burst = XLALSimBurstTableFromLIGOLw(fName.Data(),&gpsStartTime,&gpsStopTime); 
    if(sim_burst==NULL) { 
      cout << "CWB::mdc::SetSkyDistribution - XML File : " << fName.Data() << endl;
      cout << " No Events present in the range : " << (int)start << ":" << (int)stop <<endl;
    }
    skymap sm;		// used to convert ra -> phi
    int sim_burst_index=0;
    for(; sim_burst; sim_burst = sim_burst->next) {

      sim_burst_index++;

      // the gps time of the injection corresponds to the sample at t=0 
      // in the resultant time series. See XLALGenerateSimBurst in GenerateBurst.c
      double gps    = sim_burst->time_geocent_gps.gpsSeconds;
             gps   += 1.e-9*sim_burst->time_geocent_gps.gpsNanoSeconds;

      if(gps<start || gps>stop) continue;	// skip if injection is not in the time range

      char name[128]= ""; 

      // read waveform parameters (see XLALWriteLIGOLwXMLSimBurstTable in LIGOLwXML.c)
      double phi               = sm.RA2phi(rad2deg*sim_burst->ra,gps);
      double theta             = rad2deg*(TMath::PiOver2()-sim_burst->dec);
      double psi               = rad2deg*sim_burst->psi;
      double frequency         = sim_burst->frequency;
      double Q                 = sim_burst->q;
      double duration          = sim_burst->duration;
      double bandwidth         = sim_burst->bandwidth;
      double rho               = -1.;	// rho<0 -> no hrss rescaling in GetBurst
      double phase             = rad2deg*sim_burst->pol_ellipse_angle;	
      double eccentricity      = sim_burst->pol_ellipse_e;
      double hrss              = sim_burst->hrss * hrss_factor;
      double amplitude         = sim_burst->amplitude;
      double egw_over_rsquared = sim_burst->egw_over_rsquared;
      double waveform_number   = sim_burst->waveform_number;

      // NOTE : iota is set to eccentricity but is not used because LAL implicity 
      //        apply with the parameter sim_burst->pol_ellipse_e
      // NOTE : CWB uses ellipticity, instead LAL uses eccentricity : 
      //        see XLALSimBurstSineGaussian in LALSimBurst.c
      // NOTE : the parameter phase is used by LAL to generate SG,CG waveforms 
      //        when random polarization is used the phase parameter is not effective
      //        because the polarization angle is equivalent to the phase angle  

      // eccentricity is saved in iota parameter, it is used in CWB::mdc::GetBurst
      double iota = eccentricity;
      //double cosi = e2cosi(eccentricity);
      //double iota = acos(cosi)*180./TMath::Pi();

      mdcid waveid;
      vector<mdcpar> wpar; 

      // Convert LAL SimBurst defined in LALSimBurst.c into cWB MDC waveforms
      if(TString(sim_burst->waveform)=="Gaussian") {

        // add gaussian waveform
        wpar.resize(2);
        wpar[0].name="duration"; wpar[0].value=duration;
        wpar[1].name="decimals"; wpar[1].value=decimals;
        waveid = engine ? AddWaveform(MDC_GA, wpar) : 
                          AddWaveform(MDC_GA_LAL, sim_burst, wpar);
        sprintf(name,waveid.name.Data());

      } else 
      if(TString(sim_burst->waveform)=="SineGaussian") {

        // add sinegaussian waveform
        wpar.resize(3);
        wpar[0].name="frequency"; wpar[0].value=frequency;
        wpar[1].name="Q";         wpar[1].value=Q;
        wpar[2].name="decimals";  wpar[2].value=decimals;
        waveid = engine ? AddWaveform(MDC_SGE, wpar) : 
                          AddWaveform(MDC_SGE_LAL, sim_burst, wpar);
        sprintf(name,waveid.name.Data());

      } else 
      if(TString(sim_burst->waveform)=="BTLWNB") {

        // add whitenoiseburst waveform
        wpar.resize(7);
        wpar[0].name="frequency"; wpar[0].value=frequency;
        wpar[1].name="bandwidth"; wpar[1].value=bandwidth;
        wpar[2].name="duration";  wpar[2].value=duration;
        wpar[3].name="pseed";     wpar[3].value=seed+2*sim_burst_index;
        wpar[4].name="xseed";     wpar[4].value=seed+2*sim_burst_index+1;
        wpar[5].name="mode";      wpar[5].value=0;  // asymmetric
        wpar[6].name="decimals";  wpar[6].value=decimals;
        waveid = engine ? AddWaveform(MDC_WNB, wpar) : 
                          AddWaveform(MDC_WNB_LAL, sim_burst, wpar);
        sprintf(name,waveid.name.Data());

      } else 
      if(TString(sim_burst->waveform)=="StringCusp") {

        // add string cusp waveform
        wpar.resize(3);
        wpar[0].name="frequency"; wpar[0].value=frequency;
        wpar[1].name="amplitude"; wpar[1].value=amplitude;
        wpar[2].name="decimals";  wpar[2].value=decimals;
        waveid = AddWaveform(MDC_SC_LAL, sim_burst, wpar);
        sprintf(name,waveid.name.Data());

      }

      if(AUTO) {
	// fill xmlType
        for(int i=0;i<(int)wfList.size();i++) {
          bool save=true;
          for(int j=0; j<(int)xmlType.size(); j++){
            if(wfList[i].name.CompareTo(xmlType[j])==0) {save = false; break;}
          }
          if(save) {
            xmlType.push_back(wfList[i].name.Data());
          }
        }
      } else { 
        // check if waveid.name is in the xmlType list
        bool present=false;
        for(int n=0;n<xmlType.size();n++) {
          if(waveid.name==xmlType[n]) {present=true;break;}
        }
        if(!present) {
          cout << "CWB::mdc::SetSkyDistribution - "
               << "Error : mdc type " << waveid.name << " is not in the xmlType list" << endl;
          cout<<endl<<"xmlType list : "<<endl<<endl;
          for(int n=0;n<xmlType.size();n++) {
            cout << "xmlType[" << n << "] : " << xmlType[n] << endl;
          }
          cout<<endl;
          exit(1);
        }
      }

      //cout.precision(3); 
      //cout << (int)gps << " " <<  name << " th: " <<  theta << " ph: " <<  phi << " psi: " <<  psi   
      //     << " f: " <<  frequency << " Q: " << Q << " d: " << duration << " b: " << bandwidth 
      //     << " p:" << phase << " e: " << eccentricity << endl;
  
      // fill event parameter list   
       gpsList.push_back(gps);
      nameList.push_back(name);
        thList.push_back(theta);
        phList.push_back(phi);
       psiList.push_back(psi);
       rhoList.push_back(rho);
      iotaList.push_back(iota);
      hrssList.push_back(hrss);
        IDList.push_back(waveid.ID);
        idList.push_back(waveid.id);
    }
#else
    cout << "CWB::mdc::SetSkyDistribution - "
         << "Error : User Distribution from XML is enabled only with LAL" << endl;
    exit(1);
#endif

  } else
  if((sky_distribution==MDC_EARTH_FIX)||
     (sky_distribution==MDC_CELESTIAL_FIX)) {

    bool error=false;
    bool gerror=false;
    double theta = GetPar("theta",par,error);
    if(gerror) gerror=error;
    double phi   = GetPar("phi",par,error);
    if(gerror) gerror=error;

    if(gerror) {
      cout << "CWB::mdc::SetSkyDistribution - "
           << "Error : num par must be at least 2 [theta,phi,(psi,rho,gps)]" << endl;
      exit(1);
    }

    double value,psi,rho,iota,hrss;
    value=psi=rho=iota=0.;
    error=false;value=GetPar("entries",par,error);
    int entries = error ? 10000 : int(value);
    for(int n=0;n<entries;n++) {
       error=false;value=GetPar("psi",par,error);
       psi = error ? gRandom->Uniform(0,180) : value;
       error=false;value=GetPar("rho",par,error);
       rho = error ? 0. : value;
       error=false;value=GetPar("gps",par,error);
       if(!error) gpsList.push_back(value);

       error=false;iota = GetPar("iota",par,error);
       if(error)                   iotaList.push_back(0);
       else if(iota>=0&&iota<=180) iotaList.push_back(iota);   
       else                        iotaList.push_back(rad2deg*acos(gRandom->Uniform(-1,1)));

       thList.push_back(theta);
       phList.push_back(phi);
      psiList.push_back(psi);
      rhoList.push_back(rho);
     hrssList.push_back(0);
       IDList.push_back(-1);
       idList.push_back(-1);
    }
  } else
  if(sky_distribution==MDC_LOGFILE) {
    if(net==NULL) {
      cout << "CWB::mdc::SetSkyDistribution - Error : Dummy method : network is not initialized " << endl;
      exit(1);
    }
    int nInj=net->readMDClog(const_cast<char*>(fName.Data()),0.);
    printf("CWB::mdc::SetSkyDistribution - injections loaded : %d\n",nInj);
    for(int k=0;k<nInj;k++) net->mdc__ID.push_back(k); 
    size_t N = net->mdc__IDSize();
    size_t M = net->mdcTypeSize();
    cout << N << " " << M << endl;
    // fill mdcName vector with name of MDC (used in GetSourceList)
    mdcName.clear();
    for(int k=0;k<M;k++) mdcName.push_back(net->mdcType[k]);
    int nIFO=net->ifoListSize();
    inj = new injection(nIFO);
    inj_tree = inj->setTree();
    inj->output(inj_tree,net,1,false);
    delete inj;
    inj = new injection(inj_tree,nIFO);
    cout << "inj entries : " << inj->GetEntries() << endl;
  }

  return;
}

//______________________________________________________________________________
void 
CWB::mdc::DrawSkyDistribution(TString name, TString projection, TString coordinate, double resolution, bool background) {
//
// Draw Sky Distribution
//
//
// Input: name        - unique name for canvas plot 
//        projection  - hammer, sinusoidal, parabolic        
//        coordinate  - geographic, celestial
//        resolution  - default 2, define the resolution of pixels       
//        background  - false/true : if true draw a background color (blue)
//

  if(psp!=NULL) delete psp;
  psp = new gskymap(0.4,0,180,0,360);
  psp->SetOptions(projection,coordinate,resolution/2);

  if(background) {
    psp->SetGridxColor(kWhite);
    psp->SetGridyColor(kWhite);
  } else {
    psp->SetGridxColor(kBlack);
    psp->SetGridyColor(kBlack);
  }
//  psp->SetGalacticDisk(true);
//  psp->SetGalacticDiskColor(kYellow);

  int entries=0;
  int size=360*2*resolution*180*2*resolution;
  wavearray<double> x,y,z;

  double zmax=0.;
  if(sky_distribution==MDC_LOGFILE) {

    if(inj_tree==NULL) {
      cout << "CWB::mdc::DrawSkyDistribution - Error : injection object is NULL" << endl;
      exit(1);
    }

    inj_tree->Draw("theta[0]:phi[0]:distance[0]:time[0]","","goff");
    entries = inj_tree->GetSelectedRows();
    size += entries;
    x.resize(size); y.resize(size); z.resize(size);
    double* theta    = inj_tree->GetV1();
    double* phi      = inj_tree->GetV2();
    double* distance = inj_tree->GetV3();
    double* gps      = inj_tree->GetV4();
    coordinate.ToUpper();
    if(coordinate.CompareTo("CELESTIAL")==0) {
      for(int n=0;n<entries;n++) x[n] = sm.phi2RA(phi[n],gps[n]);  // celestial 2 earth coordinates
    } else {
      for(int n=0;n<entries;n++) x[n] = phi[n];  
    }
    for(int n=0;n<entries;n++) {
      y[n] = theta[n];
      z[n] = distance[n];
      if(zmax<z[n]) zmax=z[n];
    }

  } else {
    entries = (int)rhoList.size();
    size += entries;
    x.resize(size); y.resize(size); z.resize(size);
    for(int n=0;n<entries;n++) {
      x[n]=phList[n];
      y[n]=thList[n];
      z[n]=rhoList[n];
      if(zmax<z[n]) zmax=z[n];
    }
  }

  if(entries==0) {
    cout << "CWB::mdc::DrawSkyDistribution - Warning : no entries in sky distribution " << endl;
    exit(1);
  }

  size=entries;
  // add blue background
  if(background) {
    for (int i=0;i<360*2*resolution;i++) {
      for (int j=0;j<180*2*resolution;j++) {
        double ph = i/(double)(2*resolution);
        double th = j/(double)(2*resolution);
        x[size]=ph;
        y[size]=th;
        z[size]=zmax;
        size++;
      }
    }
  }

  // sort source respect to distance
  if(entries>1) {
    Int_t *index = new Int_t[size];
    TMath::Sort(size,z.data,index,true);
    wavearray<double> T(size);
    for (int i=0;i<size;i++) T[i]=x[index[i]];
    for (int i=0;i<size;i++) x[i]=T[i];
    for (int i=0;i<size;i++) T[i]=y[index[i]];
    for (int i=0;i<size;i++) y[i]=T[i];
    for (int i=0;i<size;i++) T[i]=z[index[i]];
    for (int i=0;i<size;i++) z[i]=T[i];
    T.resize(0);
    delete [] index;
  }

  psp->FillData(size, x.data, y.data, z.data);

  psp->SetZaxisTitle("Kpc");
  psp->Draw(-2);

  x.resize(0);
  y.resize(0);
  z.resize(0);

  return;
}

//______________________________________________________________________________
TString
CWB::mdc::WriteFrameFile(TString frDir, TString frLabel, size_t gps, size_t length, bool log, vector<TString> chName) {
//
// Write mdc to frame file
//
//
// Input: frDir     - output directory
//        frLabel   - label used for output file name
//                    file name path = frDir/network-frLabel-(gps/100000)/network-frLabel-gps.gwf
//        gps       - time of frame (sec - integer)     
//        length    - time length of frame (sec - integer)
//        log       - if 'true' then write log file
//        chName    - channel name list, list size must be 0 or nIFO
//                    if list size=0 then chName is IFO_NAME:GW-H 
//                    if list size=nIFO then is chName[n]="" then chName is IFO_NAME:GW-H 
//
// Return log string
//

  if(net==NULL) {
    cout << "CWB::mdc::WriteFrameFile - Error : Dummy method : network is not initialized " << endl;
    exit(1);
  }

  int nIFO=net->ifoListSize();
  char ifoLabel[64]="";
  for(int i=0;i<nIFO;i++) sprintf(ifoLabel,"%s%s",ifoLabel,net->ifoName[i]);

  // check input chName
  if(chName.size()!=0 && chName.size()!=nIFO) {
    cout << "CWB::mdc::WriteFrameFile - Error : chName list size must be equals to nIFO " << nIFO << endl;
    exit(1);
  }

  // make sub directory
  char sdir[64];
  sprintf(sdir,"%s-%s-%04d",ifoLabel,frLabel.Data(),int(gps/100000));
  char cmd[128];sprintf(cmd,"mkdir -p %s/%s",frDir.Data(),sdir);
  cout << cmd << endl;
  gSystem->Exec(cmd);

  char frFile[512];
  sprintf(frFile,"%s/%s/%s-%s-%lu-%lu.gwf",frDir.Data(),sdir,ifoLabel,frLabel.Data(),gps,length);
  cout << frFile << endl;
  FrFile *ofp = FrFileONew(frFile,1);   // gzip compression

  /*----------------------- Create a new frame ---------*/

  FrameH* simFrame = FrameNew(const_cast<char*>(ifoLabel));
  simFrame->frame = 0;
  simFrame->run = -1;
  simFrame->dt = length;
  simFrame->GTimeS = gps;
  simFrame->GTimeN = 0;

  wavearray<double> x(MDC_SAMPLE_RATE*length);
  x.rate(MDC_SAMPLE_RATE);
  x.start(gps);
  char chNAME[64];
  for(int i=0;i<nIFO;i++) {

    TString ifo=net->ifoName[i];
    if(chName.size()) {
      if(chName[i]!="") sprintf(chNAME,"%s",chName[i].Data());
      else              sprintf(chNAME,"%s:GW-H",net->ifoName[i]);
    } else {
      sprintf(chNAME,"%s:GW-H",net->ifoName[i]);
    }

    Get(x,ifo); 
    
    cout << "Size (sec) " << x.size()/x.rate() << endl;
    FrProcData* proc = FrProcDataNew(simFrame,chNAME,x.rate(),x.size(),-64);
    if(proc == NULL) {cout << "CWB::mdc::WriteFrameFile - Cannot create FrProcData" << endl; exit(-1);}
    proc->timeOffset = 0;
    proc->tRange = simFrame->dt;
    proc->type = 1;   // Time Serie

    for (int i=0;i<(int)proc->data->nData;i++) proc->data->dataD[i] = x[i];
  }

  int err=FrameWrite(simFrame,ofp);
  if (err) {cout << "CWB::mdc::WriteFrameFile - Error writing frame" << endl;exit(1);}
  FrameFree(simFrame);

  if (ofp!=NULL) FrFileOEnd(ofp);


  // log string
  TString logString = "";
  for(int i=0;i<(int)mdcList.size();i++) logString = logString+mdcList[i]+"\n";

  // write log file
  if(log) {
    TString logFile=frFile;
    logFile.ReplaceAll(".gwf","-Log.txt");
    DumpLog(logFile);
/*
    ofstream out;
    out.open(logFile.Data(),ios::out);
    if (!out.good()) {cout << "CWB::mdc::WriteFrameFile - Error Opening File : " << logFile.Data() << endl;exit(1);}
    out.precision(14);
    for(int i=0;i<(int)mdcList.size();i++) out << mdcList[i] << endl;
    out.close();
*/
  }

  return logString;
}

//______________________________________________________________________________
void
CWB::mdc::Dump(TString fname, int ID, int id, TString polarization) {
//
// Write waveform to file
//
//
// Input: fname        - output file name
//        ID           - major id of waveform list
//        id           - minor id of waveform list (Ex : the list WNB with different random waveforms)
//        polarization - hp/hx 
//

  polarization.ToUpper();
  waveform wf = GetWaveform(ID,id);
  if(wf.status) {
    if(polarization.Contains("HP")) Dump(fname, wf.hp);
    if(polarization.Contains("HX")) Dump(fname, wf.hx);
  }
  return;
}

//______________________________________________________________________________
void
CWB::mdc::Dump(TString fname, TString name, int id, TString polarization) {
//
// Write waveform to file
//
//
// Input: fname        - output file name
//        name         - name of waveform 
//        id           - minor id of waveform list (Ex : the list WNB with different random waveforms)
//        polarization - hp/hx 
//

  polarization.ToUpper();
  waveform wf = GetWaveform(name,id);
  if(wf.status) {
    if(polarization.Contains("HP")) Dump(fname, wf.hp);
    if(polarization.Contains("HX")) Dump(fname, wf.hx);
  }
  return;
}

//______________________________________________________________________________
void
CWB::mdc::Dump(TString fname, wavearray<double>& x) {
//
// Write waveform to file
//
//
// Input: fname        - output file name
//        x            - wavearray which contains the waveform data
//

  ofstream out;
  out.open(fname.Data(),ios::out);
  if (!out.good()) {cout << "CWB::mdc::WriteFrameFile - Error Opening File : " << fname.Data() << endl;exit(1);}
  out.precision(14);
  for(int i=0;i<(int)x.size();i++) out << x[i] << endl;
  out.close();

  return;
}  

//______________________________________________________________________________
void
CWB::mdc::Dump(TString dname) {
//
// Write all waveform to file
//
//
// Input: dname        - output directory name
//                       all file are written under dname
//

  // create dir
  char cmd[256];
  sprintf(cmd,"mkdir -p %s",dname.Data());
  gSystem->Exec(cmd);

  char sid[8];
  for(int i=0;i<(int)wfList.size();i++) {
    //cout << "ID : " << i << "\t" << wfList[i].name.Data() << endl;
    Dump(dname+"/"+wfList[i].name+"~01.txt", i, 0, "hp");
    Dump(dname+"/"+wfList[i].name+"~02.txt", i, 0, "hx");
    printf("%s\n",TString(dname+"/"+wfList[i].name+"~01.txt").Data());
    printf("%s\n",TString(dname+"/"+wfList[i].name+"~02.txt").Data());
    for(int j=0;j<(int)wfList[i].list.size();j++) {
      //cout << " id : " << j+1 << "\t" << wfList[i].list[j].name.Data() << endl;
      sprintf(sid,"~%02d",2*(j+1)+1);
      Dump(dname+"/"+wfList[i].name+sid+".txt", i, j+1, "hp");
      printf("%s\n",TString(dname+"/"+wfList[i].name+sid+".txt").Data());
      sprintf(sid,"~%02d",2*(j+1)+1);
      Dump(dname+"/"+wfList[i].name+sid+".txt", i, j+1, "hx");
      printf("%s\n",TString(dname+"/"+wfList[i].name+sid+".txt").Data());
    }
  }
  return;
}

//______________________________________________________________________________
double
CWB::mdc::GetAntennaPattern(TString ifo, double phi, double theta, double psi, TString polarization) {
//
// Get Antenna Pattern 
//
//
// Input: ifo          - name of detector
//        phi          - longitude    (degrees)
//        theta        - latitude     (degrees)
//        psi          - polarization (degrees) 
//        polarization - hp/hx         - 
//

  int nIFO = net ? net->ifoListSize() : 0;

  int ifoId=-1;
  detector* pD=NULL;
  for(int n=0; n<nIFO; n++) if(ifo.CompareTo(net->ifoName[n])==0) ifoId=n;
  if(ifoId>=0) pD = net->getifo(ifoId);
  else         pD = new detector(const_cast<char*>(ifo.Data()));
 
  wavecomplex F = pD->antenna(theta,phi,psi);
  if(ifoId<0) delete pD;
  polarization.ToUpper();
  if(polarization.Contains("HP")) return F.real(); else return F.imag();
}

//______________________________________________________________________________
double
CWB::mdc::GetDelay(TString ifo1, TString ifo2, double phi, double theta) {
//
// Get Delay Traveling Time between two detectors 
//
//
// Input: ifo1    - name of the first detector
//        ifo2    - name of the second detector
//                  if ifo2="" -> return ifo1 delay respect to geocenter
//        phi     - longitude (degrees)     
//        theta   - latitude  (degrees)
//

  if(ifo1.CompareTo(ifo2)==0) return 0.0;

  int nIFO = net ? net->ifoListSize() : 0;

  int ifoId1=-1;
  detector* pD1=NULL;
  for(int n=0; n<nIFO; n++) if(ifo1.CompareTo(net->ifoName[n])==0) ifoId1=n;
  if(ifoId1>=0) pD1 = net->getifo(ifoId1);
  else          pD1 = new detector(const_cast<char*>(ifo1.Data()));

  if(ifo2=="") {
    if(ifoId1<0) delete pD1;
    double delay = pD1->getTau(theta,phi);
    return delay;	// return Delay Traveling Time between ifo1 and geocenter
  }

  int ifoId2=-1;
  detector* pD2=NULL;
  for(int n=0; n<nIFO; n++) if(ifo2.CompareTo(net->ifoName[n])==0) ifoId2=n;
  if(ifoId2>=0) pD2 = net->getifo(ifoId2);
  else          pD2 = new detector(const_cast<char*>(ifo2.Data()));

  double delay = pD1->getTau(theta,phi)-pD2->getTau(theta,phi);

  if(ifoId1<0) delete pD1;
  if(ifoId2<0) delete pD2;

  return delay;		// return Delay Traveling Time between ifo1 and ifo2
}

//______________________________________________________________________________
void
CWB::mdc::DumpLogHeader(TString fName, TString label, int size) {
//
// Dump log header to file  
//
//
// Input: fName    - name of output file
//        label    - label of file
//         size    - number of mdc - if size<=0 use inj tree
//

  if(size<=0 && inj==NULL) {
    cout << "CWB::mdc::DumpLogHeader - Error : dump needs injection object " << endl;
    exit(1);
  }

  // WRITE field headings to log file

  ofstream out;
  out.open(fName.Data(),ios::out);
  if (!out.good()) {cout << "CWB::mdc::DumpLogHeader - Error Opening File : " << fName.Data() << endl;exit(1);}

  int nIFO=net->ifoListSize();

  out << "# Log File for Burst " << label.Data() << endl;
  out << "# Produced by CWB::mdc " << wat::Time("now").GetDateString().Data() << endl;
  out << "# The following detectors were simulated" << endl;
  for(int n=0;n<nIFO;n++) out << "# - " << net->ifoName[n] << endl;
  if(size>0) {
    out << "# There were " << size << " injections" << endl;
  } else {
    out << "# There were " << net->mdcList.size() << " injections" << endl;
//    for(int i=0;i<(int)wfList.size();i++) {
//      out << "# Burst MDC Sim type " << wfList[i].name.Data() << endl;
//    }
    for(int i=0;i<(int)net->mdcType.size();i++) {
      char cut[16];sprintf(cut,"type==%d",i+1);
      inj_tree->Draw("type",cut,"goff");
      int ninj = inj_tree->GetSelectedRows();
      out << "# Burst MDC Sim type " << net->mdcType[i].data() << "\t occurred " << ninj << " times" << endl;
    }
  }

  out << "#  GravEn_SimID                        SimHrss      SimEgwR2    GravEn_Ampl  ";
  out << "Internal_x   Internal_phi   ";
  out << "External_x   External_phi  External_psi  ";
  out << "FrameGPS     EarthCtrGPS     SimName               SimHpHp      SimHcHc      SimHpHc   ";
  for(int n=0;n<nIFO;n++) {
     out << net->ifoName[n] << "      "      << 
            net->ifoName[n] << "ctrGPS     " <<
            net->ifoName[n] << "fPlus      " <<
            net->ifoName[n] << "fCross     ";
  }
  out << endl;

  out.close();

  return;
}

#ifdef _USE_LAL
//______________________________________________________________________________
void
CWB::mdc::CreateBurstXML(TString fName, vector<mdcpar> xml_parms) {
//
// Create XML file  : Create XML LAL file + Add Header 
//                    see LAL LIGOLwXML.c
//
// Input: fName      - name of the xml output file
//        xml_parms  - contains custom auxiliary parameters to be added to the parameter table
//

  if(inspName!="") {
    cout << "CWB::mdc::CreateBurstXML : Error : This function is allowed only for Burst !!!" << endl;
    exit(1);
  }

  if(xml_filename!="") {
    cout << "CWB::mdc::CreateBurstXML : Error : file XML already opened !!!" << endl;
    exit(1);
  }

  bool error=false;
  // extract start_time - start GPS time of interval in which to synthesize injections 
  double start_time = GetPar("start_time",xml_parms,error); if(error) start_time = 0;       
  // extract end_time   - end GPS time of interval in which to synthesize injections 
  double end_time = GetPar("end_time",xml_parms,error); if(error) end_time = 0;       

  xml_filename = fName;
  xml_process_table_head = NULL;
  xml_process_params_table_head = NULL;
  xml_search_summary_table_head = NULL;
  xml_time_slide_table_head = NULL;
  xml_sim_burst_table_head = NULL;
  xml_sim_burst = &xml_sim_burst_table_head;

  TString ifos="";
  if(net!=NULL) for(int n=0;n<(int)net->ifoListSize();n++) ifos+=net->getifo(n)->Name+TString(",");
  if(ifos!="") ifos.Resize(ifos.Sizeof()-2);

  // Fill Process table 
  xml_process_table_head = xml_process = XLALCreateProcessTableRow();
  sprintf(xml_process->program, "cWB");			// program name entry
  XLALGPSTimeNow(&xml_process->start_time);		// start time
  sprintf(xml_process->version, watversion('r'));	// git version
  sprintf(xml_process->cvs_repository, "https://git.ligo.org/cWB/library/");
  char cmd_line[2048]="";
  for(int i=0;i<gApplication->Argc();i++) sprintf(cmd_line,"%s %s",cmd_line,gApplication->Argv(i));
  sprintf(xml_process->comment, cmd_line); 	// save command line used to start the application 
  sprintf(xml_process->domain, "cWB");		// online flag and domain
  // unix process id, username, host, and process_id 
  UserGroup_t* uinfo = gSystem->GetUserInfo();
  xml_process->unix_procid = gSystem->GetPid();
  sprintf(xml_process->node, gSystem->HostName());
  sprintf(xml_process->username, uinfo->fUser);
  sprintf(xml_process->ifos,ifos.Data());
  xml_process->process_id = 0;

  // Fill Process parameter table
  ProcessParamsTable **paramaddpoint = &xml_process_params_table_head;
  if(net!=NULL) for(int n=0;n<(int)net->ifoListSize();n++) {
     ADD_PROCESS_PARAM(paramaddpoint, xml_process, "string", "instrument", net->getifo(n)->Name); 
  }
  ADD_PROCESS_PARAM(paramaddpoint, xml_process, "double", "inj-hrss",  
                    TString::Format("%f",GetInjHrss()).Data()); 
  ADD_PROCESS_PARAM(paramaddpoint, xml_process, "double", "time-step",  
                    TString::Format("%f",1./GetInjRate()).Data()); 
  ADD_PROCESS_PARAM(paramaddpoint, xml_process, "double", "jitter",  
                    TString::Format("%f",GetInjJitter()).Data()); 
  ADD_PROCESS_PARAM(paramaddpoint, xml_process, "enum", "sky-distribution", 
                    DistributionToString(sky_distribution)); 
  ADD_PROCESS_PARAM(paramaddpoint, xml_process, "string", "sky-file", sky_file.Data()); 
  for(int i=0;i<sky_parms.size();i++) {
    ADD_PROCESS_PARAM(paramaddpoint, xml_process, "double", 
                      sky_parms[i].name.Data(), TString::Format("%f",sky_parms[i].value).Data()); 
  }
  for(int i=0;i<xml_parms.size();i++) {
    if(xml_parms[i].svalue!="") {
      ADD_PROCESS_PARAM(paramaddpoint, xml_process, "string", 
                        xml_parms[i].name.Data(), xml_parms[i].svalue.Data()); 
    } else {
      ADD_PROCESS_PARAM(paramaddpoint, xml_process, "double", 
                        xml_parms[i].name.Data(), TString::Format("%f",xml_parms[i].value).Data()); 
    }
  }

  // Fill Search summary table
  xml_search_summary_table_head = xml_search_summary = XLALCreateSearchSummaryTableRow(xml_process);
  sprintf(xml_search_summary->comment, "");
  sprintf(xml_search_summary->ifos,ifos.Data());
  xml_search_summary->nnodes = 1;
  wat::Time in_start_time(start_time);
  xml_search_summary->in_start_time.gpsSeconds = in_start_time.GetSec();
  xml_search_summary->in_start_time.gpsNanoSeconds = in_start_time.GetNSec();
  xml_search_summary->out_start_time.gpsSeconds = in_start_time.GetSec();
  xml_search_summary->out_start_time.gpsNanoSeconds = in_start_time.GetNSec();
  wat::Time in_end_time(end_time);
  xml_search_summary->in_end_time.gpsSeconds = in_end_time.GetSec();
  xml_search_summary->in_end_time.gpsNanoSeconds = in_end_time.GetNSec();
  xml_search_summary->out_end_time.gpsSeconds = in_end_time.GetSec();
  xml_search_summary->out_end_time.gpsNanoSeconds = in_end_time.GetNSec();

}
#endif

#ifdef _USE_LAL
//______________________________________________________________________________
void
CWB::mdc::FillBurstXML(bool verbose) {
//
// Fill XML LAL to file with event parameters
//
// NOTE : this function can be used only afer CreateBurstXML 
//

  if(xml_filename=="") {
    cout << "CWB::mdc::FillBurstXML : Error : file XML not declared !!! use CreateBurstXML" << endl;
    exit(1);
  }

  double NaN = std::numeric_limits<double>::quiet_NaN();
  const double deg2rad = TMath::Pi()/180.;

  if(verbose) {
    cout << "gps" << "\t\t" << "wf-name" << "\t\t\t" << "theta" << "\t" << "phi" << "\t" << "psi" 
         << "\t" << "rho" << "\t" << "iota" << "\t" << "hrss" << "\t" << "ID" << "\t" << "id" << endl;
  }

  bool error=false;
  for(int i=0;i<(int)srcList.size();i++) {
    waveform wf = GetWaveform(srcList[i].ID, srcList[i].id); 
    double hrss = srcList[i].hrss>0 ? srcList[i].hrss : inj_hrss;
    if(verbose) {
      cout << std::setprecision(13) << srcList[i].gps << "\t" << std::setprecision(6) 
           << wf.name.Data() << "\t" << srcList[i].theta << "\t" 
           << srcList[i].phi << "\t" << srcList[i].psi << "\t" << srcList[i].rho << "\t" 
           << srcList[i].iota << "\t" << hrss << "\t" << srcList[i].ID << "\t" << srcList[i].id << endl;
    }

    SimBurst *sim_burst = XLALCreateSimBurst();
    if(!sim_burst) {*xml_sim_burst=NULL; return;}

    // fill sim_burst structure
    strcpy(sim_burst->waveform, wf.name.Data());

    wat::Time time(srcList[i].gps);
    sim_burst->time_geocent_gps.gpsSeconds = time.GetSec();
    sim_burst->time_geocent_gps.gpsNanoSeconds = time.GetNSec();

    sim_burst->ra                = deg2rad*sm.phi2RA(srcList[i].phi,srcList[i].gps); 
    sim_burst->dec               = TMath::PiOver2()-deg2rad*srcList[i].theta;
    sim_burst->psi               = deg2rad*srcList[i].psi;
    error=false;
    sim_burst->frequency         = GetPar("frequency",wf.par,error);if(error) sim_burst->frequency=NaN;
    error=false;
    sim_burst->q                 = GetPar("Q",wf.par,error);        if(error) sim_burst->q=NaN;
    error=false;
    sim_burst->duration          = GetPar("duration",wf.par,error); if(error) sim_burst->duration=NaN;
    error=false;
    sim_burst->bandwidth         = GetPar("bandwidth",wf.par,error);if(error) sim_burst->bandwidth=NaN;
    sim_burst->pol_ellipse_angle = 0;
    sim_burst->pol_ellipse_e     = cos(deg2rad*srcList[i].iota);
    sim_burst->hrss              = hrss;
    sim_burst->amplitude         = 0;
    sim_burst->egw_over_rsquared = 0;
    sim_burst->waveform_number   = srcList[i].ID;

    *xml_sim_burst = sim_burst;
    xml_sim_burst = &(*xml_sim_burst)->next;
  }

}
#endif

#ifdef _USE_LAL
//______________________________________________________________________________
void
CWB::mdc::CloseBurstXML() {
//
// Close XML LAL file 
//
// NOTE : this function can be used only afer CreateBurstXML 
//

  if(xml_filename=="") {
    cout << "CWB::mdc::CloseBurstXML : Error : file XML not declared !!! use CreateBurstXML" << endl;
    exit(1);
  }

  XLALGPSTimeNow(&xml_process->end_time);

  // write xml file (see LAL binj.c)
  LIGOLwXMLStream *xml;

  xml = XLALOpenLIGOLwXMLFile(xml_filename.Data());

  /* process table */
  if(XLALWriteLIGOLwXMLProcessTable(xml, xml_process_table_head)) {
    cout << "CWB::mdc::CloseBurstXML : Error in XLALWriteLIGOLwXMLProcessTable" << endl;
    exit(1);
  }
  /* process params table */
  if(XLALWriteLIGOLwXMLProcessParamsTable(xml, xml_process_params_table_head)) {
    cout << "CWB::mdc::CloseBurstXML : Error in XLALWriteLIGOLwXMLProcessParamsTable" << endl;
    exit(1);
  }
  /* search summary table */
  if(XLALWriteLIGOLwXMLSearchSummaryTable(xml, xml_search_summary_table_head)) {
    cout << "CWB::mdc::CloseBurstXML : Error in XLALWriteLIGOLwXMLSearchSummaryTable" << endl;
    exit(1);
  }
  /* time slide table */
  if(XLALWriteLIGOLwXMLTimeSlideTable(xml, xml_time_slide_table_head)) {
    cout << "CWB::mdc::CloseBurstXML : Error in XLALWriteLIGOLwXMLTimeSlideTable" << endl;
    exit(1);
  }
  /* sim burst table */
  if(XLALWriteLIGOLwXMLSimBurstTable(xml, xml_sim_burst_table_head)) {
    cout << "CWB::mdc::CloseBurstXML : Error in XLALWriteLIGOLwXMLSimBurstTable" << endl;
    exit(1);
  }

  XLALCloseLIGOLwXMLFile(xml);

  xml_filename = "";
  XLALDestroyProcessTable(xml_process_table_head);
  XLALDestroyProcessParamsTable(xml_process_params_table_head);
  XLALDestroyTimeSlideTable(xml_time_slide_table_head);
  XLALDestroySearchSummaryTable(xml_search_summary_table_head);
  XLALDestroySimBurstTable(xml_sim_burst_table_head);
}
#endif

#ifdef _USE_LAL
//______________________________________________________________________________
TString 
CWB::mdc::GetBurstNameLAL(TString options) {
//
// convert LAL Burst options to the name used in mdc class 
//
// options : --name Gaussian     --duration value
// options : --name SineGaussian --frequency value --q value
// options : --name BTLWNB       --frequency value --bandwidth value --duration value
// options : --name StringCusp   --frequency value 
//
//           --decimals value    number of decimals used to format the waveform name
//                               if not defined ; used default format
//

  TString wfname = "";

  CWB::Toolbox TB;
  TString name = TB.getParameter(options,"--name");
  bool error=false;
  if(name=="") error=true; 
  if((name=="Gaussian")&&(name!="SineGaussian")&&(name=="BTLWNB")) error=true; 
  if(error) {
    cout << "CWB::mdc::GetBurstNameLAL : Error - Not valid --name option" << endl;
    cout << "available --name options are : Gaussian/SineGaussian/BTLWNB" << endl;
    return "";
  }

  // get number of decimals to be used to format the waveform name
  int decimals = -1;
  TString sdecimals = TB.getParameter(options,"--decimals");
  if(sdecimals.IsDigit()) decimals = sdecimals.Atoi();

  if(name=="Gaussian") {

    int d1 = decimals==-1 ? 1 : decimals;

    float duration;
    TString sduration = TB.getParameter(options,"--duration");
    if(sduration.IsFloat()) {
       duration = sduration.Atof();
    } else {
       cout << "CWB::mdc::GetBurstNameLAL : Error : --duration not defined or not a number!!! " << endl;
       return "";
    }

    wfname = TString::Format("LAL_GA%.*f",d1,1000.*duration);
    wfname.ReplaceAll(".","d");
    return wfname;

  } else if(name=="SineGaussian") { 

    int d1 = decimals==-1 ? 0 : decimals;
    int d2 = decimals==-1 ? 1 : decimals;

    float frequency;
    TString sfrequency = TB.getParameter(options,"--frequency");
    if(sfrequency.IsFloat()) {
       frequency = sfrequency.Atof();
    } else {
       cout << "CWB::mdc::GetBurstNameLAL : Error : --frequency not defined or not a number!!! " << endl;
       return "";
    }

    float q;
    TString sq = TB.getParameter(options,"--q");
    if(sq.IsFloat()) {
       q = sq.Atof();
    } else {
       cout << "CWB::mdc::GetBurstNameLAL : Error : --q not defined or not a number!!! " << endl;
       return "";
    }

    wfname = TString::Format("LAL_SGE%.*fQ%.*f",d1,frequency,d2,q);
    wfname.ReplaceAll(".","d");
    if(decimals==-1) wfname.ReplaceAll("d0","");
    return wfname;

  } else if(name=="BTLWNB") { 

    int d1 = decimals==-1 ? 0 : decimals;
    int d2 = decimals==-1 ? 0 : decimals;
    int d3 = decimals==-1 ? 3 : decimals;

    float frequency;
    TString sfrequency = TB.getParameter(options,"--frequency");
    if(sfrequency.IsFloat()) {
       frequency = sfrequency.Atof();
    } else {
       cout << "CWB::mdc::GetBurstNameLAL : Error : --frequency not defined or not a number!!! " << endl;
       return "";
    }

    float bandwidth;
    TString sbandwidth = TB.getParameter(options,"--bandwidth");
    if(sbandwidth.IsFloat()) {
       bandwidth = sbandwidth.Atof();
    } else {
       cout << "CWB::mdc::GetBurstNameLAL : Error : --bandwidth not defined or not a number!!! " << endl;
       return "";
    }

    float duration;
    TString sduration = TB.getParameter(options,"--duration");
    if(sduration.IsFloat()) {
       duration = sduration.Atof();
    } else {
       cout << "CWB::mdc::GetBurstNameLAL : Error : --duration not defined or not a number!!! " << endl;
       return "";
    }

    wfname = TString::Format("LAL_WNB%.*f_%.*f_%.*f",d1,frequency,d2,bandwidth,d3,duration);
    wfname.ReplaceAll(".","d");
    return wfname;

  } else if(name=="StringCusp") {

    int d1 = decimals==-1 ? 1 : decimals;

    float frequency;
    TString sfrequency = TB.getParameter(options,"--frequency");
    if(sfrequency.IsFloat()) {
       frequency = sfrequency.Atof();
    } else {
       cout << "CWB::mdc::GetBurstNameLAL : Error : --frequency not defined or not a number!!! " << endl;
       return "";
    }

    wfname = TString::Format("LAL_SC%.*f",d1,1000.*frequency);
    wfname.ReplaceAll(".","d");
    return wfname;

  }

  return wfname;
}
#endif

//______________________________________________________________________________
double 
CWB::mdc::e2cosi(double e) {
//
// convert the LAL eccentricity used to generate SineGaussian waveform into the
// equivalent iota value  
//
// e : eccentricity [0:1]
//
// return the cos(iota) value
//
// CWB uses the following formula :
// h_+ = 0.5*(1+cos(iota)^2) * A(t) cos(Phi(t))
// h_x = cos(iota) * A(t) sin(Phi(t))
// while LAL uses :
// h_+ = A * A(t) cos(Phi(t))
// h_x = B * A(t) sin(Phi(t))
// where :
// A = 1/sqrt(2-E*E) ; B = A*sqrt(1-E*E) 
// and E is the eccentricity
// it follows that A*A+B*B=1
//

  if(e<0 || e>1) {
    cout << "CWB::mdc::e2cosi - Error : eccentricity must be defined in the range [0:1] " << endl;
    return std::numeric_limits<double>::min();
  }

  double E  = e<0.9999999999 ? e : 0.9999999999;	// avoid x1,x2=nan 
  double A  = 1./sqrt(2-E*E);
  double B  = A*sqrt(1-E*E);
  double x1 = (A+sqrt(A*A-B*B))/(B*B);
  double x2 = (A-sqrt(A*A-B*B))/(B*B);
  double x  = fabs(x1*B)>1 ? x2 : x1;

  double cosi = -x*B;

  return cosi;
}

//______________________________________________________________________________
double 
CWB::mdc::cosi2e(double cosi) {
//
// convert cos(iota) into the equivalent LAL eccentricity used 
// to generate SineGaussian waveform 
//
// cosi : cos(iota) [-1:1]
//
// return the eccentricity value
//

  if(fabs(cosi)>1) {
    cout << "CWB::mdc::cosi2e - Error : cosi must be defined in the range [-1:1] " << endl;
    return std::numeric_limits<double>::min();
  }

  double A = (1+cosi*cosi)/2;
  double B = -cosi;

  double C = sqrt(A*A+B*B);
  A/=C; B/=C;

  double E = sqrt(1-(B*B)/(A*A));

  return E; 
}

//______________________________________________________________________________
void
CWB::mdc::DumpLog(TString fName, TString label, bool append) {
//
// Dump log to file  
//
//
// Input: fName    - name of output file
//        label    - label of file
//        append   - false/true -> create new file/append to the file
//

  if(fName.Contains(".root")) {

    if(inj==NULL) {
      cout << "CWB::mdc::DumpLog - Error : dump to root needs injection object " << endl;
      exit(1);
    }

    TFile froot(fName,"RECREATE");
    if(label.Sizeof()>1) this->Write(label); else this->Write("WaveMDC");
    inj_tree->Write();
    cout << "CWB::mdc::DumpLog - inj entries written : " << inj_tree->GetEntries() << endl;
    froot.Write();
    froot.Close();

  } else
  if(fName.Contains(".txt")) {

    // log string
    TString logString = "";
    for(int i=0;i<(int)mdcList.size();i++) logString = logString+mdcList[i]+"\n";

    // write log file
    ofstream out;
    if(append) out.open(fName.Data(),ios::app);
    else       out.open(fName.Data(),ios::out);
    if (!out.good()) {cout << "CWB::mdc::DumpLog - Error Opening File : " << fName.Data() << endl;exit(1);}
    out.precision(14);
    for(int i=0;i<(int)mdcList.size();i++) out << mdcList[i] << endl;
    out.close();

  } else 
  if(fName.Contains(".lst")) {

    // write source list file
    ofstream out;
    if(append) out.open(fName.Data(),ios::app);
    else       out.open(fName.Data(),ios::out);
    if (!out.good()) {cout << "CWB::mdc::DumpLog - Error Opening File : " << fName.Data() << endl;exit(1);}
    out.precision(14);
    out << "#gps\t\t" << "name\t\t" <<  "theta\t" <<  "phi\t" <<  "psi\t" 
        << "rho\t" << "iota\t" << "hrss\t" << "ID\t" << "id\t" << endl;
    for(int i=0;i<(int)srcList.size();i++) {
      waveform wf = GetWaveform(srcList[i].ID, srcList[i].id); 
      double hrss = srcList[i].hrss>0 ? srcList[i].hrss : inj_hrss;
      out << srcList[i].gps << "\t" << wf.name.Data() << "\t" << srcList[i].theta << "\t" 
          << srcList[i].phi << "\t" << srcList[i].psi << "\t" << srcList[i].rho << "\t" 
          << srcList[i].iota << "\t" << hrss << "\t" << srcList[i].ID << "\t" << srcList[i].id << endl;
    }
    out.close();

  } else {
    cout << "CWB::mdc::DumpLog - Error : file extention must be .root" << endl;
    exit(1);
  }  

  return;
}

#ifdef _USE_LAL
//______________________________________________________________________________
TString                                                                         
CWB::mdc::GetInspiral(wavearray<double>& x, TString ifo) {                      
//
// Get inspiral mdc data of the detector = ifo
//                                            
//                                            
// Input:  x       - the input start/stop gps time are obtained from the wavearray values
//                   start = x.start();  stop = x.start()+x.size()/x.rate()              
//         ifo     - name of the detector defined in the network                         
//                                                                                       
// Output: x       - x.data contains the mdc data                                        
//                                                                                       
// Return: log     - ascii string with mdc parameters                                    
//                                                                                       

  if(net==NULL) {
    cout << "CWB::mdc::GetInspiral - Error : Dummy method : network is not initialized " << endl;
    exit(1);                                                                                     
  }                                                                                              
  if(x.rate()!=MDC_SAMPLE_RATE) {                                                                
    cout << "CWB::mdc::GetInspiral - Error : x.rate() != " << MDC_SAMPLE_RATE << endl;           
    exit(1);                                                                                     
  }                                                                                              

  x=0.;

  // start/end times 
  int gpsStartSec = x.start();               
  int gpsEndSec = gpsStartSec + x.size()/x.rate(); 

  // injections 
  INT4 numInjections = 0;
  SimInspiralTable *injections = NULL;
  TString injectionFile = GenInspiralXML(0., 0., false);                           

  // set default debug level
  lal_errhandler = LAL_ERR_EXIT;

  XLALSetSilentErrorHandler();

  // read the injections 
  numInjections = SimInspiralTableFromLIGOLw(&injections, injectionFile.Data(), gpsStartSec, gpsEndSec);

  if ( numInjections == 0 ) {
    fprintf( stderr, "CWB::mdc::GetInspiral - Warning : No injections in specified time\n");
    return "";                                                     
  } else {                                                         
    fprintf( stdout, "CWB::mdc::GetInspiral : Read %d injection(s) from the file '%s'\n",  
        numInjections, injectionFile.Data() );                     
  }                                                                

  // reset lists
  mdcList.clear();
  mdcType.clear();
  mdcTime.clear();

  fprintf(stdout, "CWB::mdc::GetInspiral : Generating injection for: %s",ifo.Data());

  FindChirpInjectSignals(x, injections, ifo);                                         

  fprintf(stdout, "\n");

  // loop over injections 
  SimInspiralTable *thisInj = NULL;
  for (thisInj = injections; thisInj; thisInj = thisInj->next)
  {                                                           
    // add waveform log string to mdc log string              
    TString log = GetInspiralLog(inspName,x.start(),thisInj); 
    //cout << "LOG : " << log.Data() << endl;                 

    // add infos to lists
    mdcList.push_back(log.Data());
    double gps = thisInj->geocent_end_time.gpsSeconds+1.e-9*thisInj->geocent_end_time.gpsNanoSeconds;
    mdcTime.push_back(gps);                                                              
  }                                                                                      
  mdcType.push_back(inspName.Data());                                                    

  if(inspXML=="") gSystem->Exec(TString("rm ")+injectionFile);  // remove injectionFile file

  while ( injections ) {
    SimInspiralTable *thisInj = NULL;
    thisInj = injections;            
    injections = injections->next;   
    LALFree( thisInj );              
  }                                  

  LALCheckMemoryLeaks();
  return "";
}           
#endif      

#ifdef _USE_LAL
//______________________________________________________________________________
TString 
CWB::mdc::GetInspiralLog(TString inspName, double FrameGPS, SimInspiralTable *thisInj) {
//
// Get Log string for inspiral mdc 
//
// Input: inspName - name of inspiral
//        FrameGPS - Frame gps time
//        thisInj  - List of waveform parameters
//

  if(net==NULL) {
    cout << "CWB::mdc::GetInspiralLog - Error : Dummy method : network is not initialized " << endl;
    exit(1);
  }

  // network ifos
  INT4 nIFO=net->ifoListSize();

  // declare variables 
  float f_isco;
  REAL8 gmst;
  REAL8 longitude;
  char  output[2048]=""; 

  // GravEn_SimID 
  if (strncmp(thisInj->waveform, "NumRel", LIGOMETA_WAVEFORM_MAX) == 0)
    sprintf(output, "%s ", thisInj->numrel_data);
  else
    sprintf(output, "file ");
  // GravEn_Ampl 
  sprintf(output, "%s 1 ",output);
  // StartSamp1 
  sprintf(output, "%s 0 ",output);
  // StartSamp2 
  sprintf(output, "%s 0 ",output);
  // Internal_x 
  sprintf(output, "%s %g ",output, cos(thisInj->inclination));
  // Internal_phi 
  sprintf(output, "%s %g ",output, thisInj->coa_phase);
  // External_x 
  sprintf(output, "%s %g ",output, cos(thisInj->latitude - LAL_PI_2));
  // External_phi 
  gmst = fmod(XLALGreenwichMeanSiderealTime(&thisInj->geocent_end_time), LAL_TWOPI);
  longitude = fmod(thisInj->longitude - gmst, LAL_TWOPI);
  if (longitude < 0)
    longitude += LAL_TWOPI;
  sprintf(output, "%s %g ",output, longitude);
  // External_psi 
  sprintf(output, "%s %g ",output, thisInj->polarization);
  // FrameGPS 
  sprintf(output, "%s %d ",output, int(FrameGPS));
  // SimStartGPS 
  sprintf(output, "%s %d.%09d ",output, thisInj->geocent_end_time.gpsSeconds, thisInj->geocent_end_time.gpsNanoSeconds);
  // SimName 
  sprintf(output, "%s %s ",output, inspName.Data());
  // SimHpHp 
  sprintf(output, "%s 0 ",output);
  // SimHcHc 
  sprintf(output, "%s 0 ",output);
  // SimHpHc 
  sprintf(output, "%s 0 ",output);

  // IFO GPS F_+ F_x eff_dist 
  for(int n=0;n<nIFO;n++) {
    TString ifo;
    double latitude,longitude,polarization;
    double theta,phi,psi;
    double rad2deg = 180./TMath::Pi();

    ifo = net->ifoName[n];
    latitude=thisInj->latitude;
    longitude=thisInj->longitude;
    polarization=thisInj->polarization;
    theta = LAL_PI_2-latitude;  // LAL -> CWB
    theta*=rad2deg;
    longitude = fmod(longitude - gmst, LAL_TWOPI);
    if (longitude < 0) longitude += LAL_TWOPI;
    phi = longitude > 0 ? longitude : 2*TMath::Pi()+longitude;
    phi*= rad2deg;
    psi = polarization*rad2deg;

    double fp = GetAntennaPattern(ifo, phi, theta, psi, "hp");
    double fc = GetAntennaPattern(ifo, phi, theta, psi, "hx");

    // compute the gw delay between ifo and earth center
    double d2r = TMath::Pi()/180.;

    detector* D = net->getifo(n);
    detectorParams dP = D->getDetectorParams();	 // get ifo parameters
    double D_theta = dP.latitude;
    double D_phi   = dP.longitude;
    double D_R     = sqrt(D->Rv[0]*D->Rv[0]+D->Rv[1]*D->Rv[1]+D->Rv[2]*D->Rv[2]);  // earth radius
    GeographicToCwb(D_phi,D_theta,D_phi,D_theta);
    TVector3 dV(1,0,0);     // set vector to ifo direction
    dV.SetTheta(D_theta*d2r);
    dV.SetPhi(D_phi*d2r);

    TVector3 sV(1,0,0);     // set vector to source direction
    sV.SetTheta(theta*d2r);
    sV.SetPhi(phi*d2r);

    double cos_dVdS = dV.Dot(sV);   // get cos between ifo and source direction
    double D_geocent_delay=-cos_dVdS*D_R/speedlight;   // gw time delay ifo-geo_center
    int tShift_sec  = int(fabs(D_geocent_delay));
    int tShift_nsec = int(1e9*(fabs(D_geocent_delay)-tShift_sec));

    // add gw time delay ifo geo_center to geocent_end_time
    int end_time_gpsSeconds     = thisInj->geocent_end_time.gpsSeconds;
    int end_time_gpsNanoSeconds = thisInj->geocent_end_time.gpsNanoSeconds; 
    if(D_geocent_delay>=0) {
      end_time_gpsSeconds += tShift_sec;
      end_time_gpsNanoSeconds += tShift_nsec;
      if ( end_time_gpsNanoSeconds >= 1000000000 ) {
        end_time_gpsSeconds += INT_4S( end_time_gpsNanoSeconds / 1000000000 );
        end_time_gpsNanoSeconds = end_time_gpsNanoSeconds % 1000000000;
      }
    } else {
      end_time_gpsSeconds -= tShift_sec;
      if ( end_time_gpsNanoSeconds >= tShift_nsec ) {
        end_time_gpsNanoSeconds -= tShift_nsec;
      } else {
        --end_time_gpsSeconds;
        end_time_gpsNanoSeconds += 1000000000 - tShift_nsec;
      }
    }

    // compute effective distance
    double cosiota  = cos(thisInj->inclination);
    double eff_dist = thisInj->distance / sqrt(pow(fp*(1.+pow(cosiota,2)),2)/4.+pow(fc*cosiota,2));

    sprintf(output, "%s %s %d.%09d %e %e %g ",output, ifo.Data(), end_time_gpsSeconds,
            end_time_gpsNanoSeconds, fp, fc, eff_dist);
  }

  // calculate frequency of innermost stable circulat orbit 
  // taken from: gr-qc/0612100 
  f_isco = 205 * (20 / (thisInj->mass1 + thisInj->mass2));

  // numerical relativity specific parameters 
  sprintf(output, "%s insp ",output);
  sprintf(output, "%s distance %g ",output, thisInj->distance);
  sprintf(output, "%s mass1 %g ",output, thisInj->mass1);
  sprintf(output, "%s mass2 %g ",output, thisInj->mass2);
  sprintf(output, "%s mchirp %g ",output, thisInj->mchirp);
  sprintf(output, "%s spin1 %g %g %g ",output, thisInj->spin1x, thisInj->spin1y, thisInj->spin1z);
  sprintf(output, "%s spin2 %g %g %g ",output, thisInj->spin2x, thisInj->spin2y, thisInj->spin2z);
  sprintf(output, "%s freq %g %g",output, thisInj->f_lower, f_isco);
  sprintf(output, "%s redshift %g",output, thisInj->alpha3);
  sprintf(output, "%s alpha1 %g",output, thisInj->alpha1);	// used by cWB PE to store PE logl
  sprintf(output, "%s alpha2 %g\n",output, thisInj->alpha2);	// used by cWB PE to store PE prior

  return output;
}
#endif

#ifdef _USE_LAL
//______________________________________________________________________________
void 
CWB::mdc::SetInspiral(TString inspName, TString inspOptions) {
//
// Set inspiral mdc parameters
//
// Input: inspName    - name of inspiral
//        inspOptions - list of options (see lalapps_inspinj options)
//                      to dump all the available options do: 'lalapps_inspinj --help'
//                      for any details refer to the LAL documentation.
//                      There are some special options added only for the mdc class
//                      --xml    "file.xml" : read mdc injection's parameters from xml file
//                      --output "file.xml" : write mdc injection's parameters to xml file
//                      --dir "tmp dir"     : directory used to store the temporary xml file, 
//                                            default=/tmp
//
// Note : The LAL function used to generate the waveforms is defined in LALInspiralWave.c
//        for LAL version <  6.13.2 it is : XLALSimInspiralChooseWaveformFromSimInspiral
//        for LAL version >= 6.13.2 it is : XLALInspiralTDWaveformFromSimInspiral
//
//        If --waveform is used together with the --xml option the waveform name 
//        in the xml file is overwritten with the waveform user name  

//  const int nOpt = 12;
//  const int nOpt = 10;
  const int nOpt = 7;
  TString option[nOpt] = {"--help","--verbose","--user-tag","--write-compress",
                          "--ipn-gps-time","--t-distr",
                          "--write-sim-ring"};
//                          "--gps-start-time","--gps-end-time","--ipn-gps-time","--t-distr",
//                          "--time-step","--time-interval","--write-sim-ring"};

  if(inspName.Sizeof()<=1) {
    cout << "CWB::mdc::SetInspiral - Error : inspName not declared" << endl;
    exit(1);
  }

  for(int i=0;i<nOpt;i++) {
    if(inspOptions.Contains(option[i])) {
      cout << "CWB::mdc::SetInspiral - Error : option " << option[i].Data() 
           << " not allowed in CWB" << endl;
      exit(1);
    }   
  }

  // approximant is substitute with waveform 
  // approximant is maintained only for back compatibility
  inspOptions.ReplaceAll("--approximant","--waveform");

  TObjArray* token = inspOptions.Tokenize(TString(" "));

  // check options waveform and xml
  int bxml=0;
  int bwaveform=0;
  for(int i=0;i<token->GetEntries();i++) {
    TObjString* otoken = (TObjString*)token->At(i);
    TString stoken = otoken->GetString();
    if(stoken=="--xml") bxml++;
    if(stoken=="--waveform") bwaveform++;
  }
  if(!bxml) {  // xlm is build by mdc cwb
    if(bwaveform!=1) {
      cout << "CWB::mdc::SetInspiral - Error : ";
      if(bwaveform>1) cout << "imultiple declaration of waveform option" << endl;
      if(!bwaveform) cout << "missing waveform option" << endl;
      exit(1);
    }   
  }   

  TString inspDump="";
  TString inspOutput="";

  // exctract xml, dir, time-step & check waveform
  cout << "CWB::mdc::SetInspiral - Read options ..." << endl;
  for(int i=0;i<token->GetEntries();i++) {
    TObjString* otoken = (TObjString*)token->At(i);
    TString stoken = otoken->GetString();
    if(stoken=="--xml") {
      otoken = (TObjString*)token->At(i+1);
      stoken = otoken->GetString();
      cout << "--xml = " << stoken.Data() << endl;
      CWB::Toolbox TB;
      TB.checkFile(stoken.Data());
      this->inspXML=stoken;
    }
    if(stoken=="--dir") {
      // xml temporary dir 
      otoken = (TObjString*)token->At(i+1);
      stoken = otoken->GetString();
      cout << "--dir = " << stoken.Data() << endl;
      CWB::Toolbox TB;
      TB.checkFile(stoken.Data());
      this->inspDIR=stoken;
      // remove dir option from inspOptions
      inspOptions.ReplaceAll("--dir","");
      inspOptions.ReplaceAll(stoken,"");
    }
    if(stoken=="--output") {
      otoken = (TObjString*)token->At(i+1);
      stoken = otoken->GetString();
      cout << "--output = " << stoken.Data() << endl;
      inspOutput=stoken;
      // remove output option from inspOptions
      inspOptions.ReplaceAll("--output","");
      inspOptions.ReplaceAll(stoken,"");
    }
    if(stoken=="--dump") {
      otoken = (TObjString*)token->At(i+1);
      stoken = otoken->GetString();
      cout << "--dump = " << stoken.Data() << endl;
      inspDump=stoken;
      // remove dump option from inspOptions
      inspOptions.ReplaceAll("--dump","");
      inspOptions.ReplaceAll(stoken,"");
    }
    if(stoken=="--time-step") {
      otoken = (TObjString*)token->At(i+1);
      stoken = otoken->GetString();
      cout << "--time-step = " << stoken.Data() << endl;
      inj_rate = 1./stoken.Atof();
    }
    if(stoken=="--time-interval") {
      otoken = (TObjString*)token->At(i+1);
      stoken = otoken->GetString();
      cout << "--time-interval = " << stoken.Data() << endl;
      inj_jitter = stoken.Atof();
    }
    if(stoken=="--waveform") {
      otoken = (TObjString*)token->At(i+1);
      stoken = otoken->GetString();
      cout << "--waveform = " << stoken.Data() << endl;
      this->waveName=stoken;
      if(XLALGetApproximantFromString(stoken.Data()) == XLAL_FAILURE) {
        cout << "CWB::mdc::SetInspiral - waveform : \'" << stoken.Data() 
             << "\' not allowed !!!" << endl;
        exit(1);
      }
    }
  }

  delete token;

  if(this->inspDIR=="") {
    cout << endl;
    cout << "CWB::mdc::SetInspiral - Error : temporary dir not defined !!! Use --dir option" << endl;
    cout << endl;
    cout << "  if this option is called from a config plugin add the following line to the inspOptions : " << endl;
    cout << "    inspOptions+= \"--dir \"+TString(cfg->tmp_dir)+\" " << endl;
    cout << "  before to call : 'MDC.SetInspiral(inspOptions);'" << endl << endl; 
    cout << "  and these lines to the plugin : " << endl;
    cout << "    // export config object to the config plugin" << endl;
    cout << "    sprintf(cmd,\"CWB::config* cfg = (CWB::config*)%p;\",cfg);" << endl;
    cout << "    gROOT->ProcessLine(cmd);" << endl;
    cout << "  before to execute the config plugin : 'cfg->configPlugin.Exec(NULL,&error);'" << endl;
    cout << endl;
    exit(1);
  }

  this->inspName=inspName;
  this->inspOptions=inspOptions;

  // write xml file & exit
  if(inspDump!="") {
    TString injectionFile = GenInspiralXML(0, 0, false);
    cout << "Save inspiral injection file to " << inspDump.Data() << endl;
    if(inspXML=="") gSystem->Exec("mv "+injectionFile+" "+inspDump); 
    else            gSystem->Exec("cp "+injectionFile+" "+inspDump); 
    exit(0);
  }

  // write xml file 
  if(inspOutput!="") {
    // Check if file exist
    Long_t id,size=0,flags,mt;
    int estat = gSystem->GetPathInfo(inspOutput.Data(),&id,&size,&flags,&mt);
    if (estat!=0) { // skip if file already exist
      TString injectionFile = GenInspiralXML(0, 0, false);
      cout << "Save inspiral injection file to " << inspOutput.Data() << endl;
      gSystem->Exec("cp "+injectionFile+" "+inspOutput); 
    }
  }

  // dummy exec to check if parameters are correct
  if(inspXML=="") GenInspiralXML(900000000, 900000000, true);

  return; 
}
#endif

#ifdef _USE_LAL
//______________________________________________________________________________
TString 
CWB::mdc::GetInspiralOption(TString inspOption) {
//
// Get inspiral mdc parameters
//
// Input: inspOption    - name of inspiral xml option
//

  if(inspOptions=="") return "";	// inspOptions not defined

  if(inspXML!="") {  	// xml is provided by user

    // read XML
    ifstream in;
    in.open(inspXML.Data(),ios::in);
    if (!in.good()) {cout << "CWB::mdc::GetInspiralOption - Error Opening xml File : " 
                          << inspXML.Data() << endl;gSystem->Exit(1);}

    char str[2048];
    while(true) {
      in.getline(str,2048);
      if (!in.good()) break;
      if(TString(str).Contains("inspinj")&&TString(str).Contains(inspOption))  {
        TObjArray* token = TString(str).Tokenize(TString(","));
        if(token->GetEntries()!=5) {cout << "CWB::mdc::GetInspiralOption - bad line format : " 
                                         << str << endl;gSystem->Exit(1);}
        TObjString* otoken = (TObjString*)token->At(4);
        TString stoken = otoken->GetString();
        stoken.ReplaceAll("\"","");
        delete token;
        return stoken;	// return token
      }
    }

    in.close();
  } else {		// inspOptions provided by the user

    TObjArray* token = inspOptions.Tokenize(TString(" "));
    for(int i=0;i<token->GetEntries();i++) {
      TObjString* otoken = (TObjString*)token->At(i);
      TString stoken = otoken->GetString();
      if(stoken.Contains(inspOption)) {
        otoken = (TObjString*)token->At(i+1);
        stoken = otoken->GetString();
        stoken.ReplaceAll("\"","");
        return stoken;
      }
    }
    delete token;
  } 

  return "";
}
#endif

#ifdef _USE_LAL
//______________________________________________________________________________
TString 
CWB::mdc::GenInspiralXML(int gps_start_time, int gps_end_time, bool rmFile) {
//
// Get XML inspiral mdc parameters
//
// Input: gps_start_time  - start time
//        gps_end_time    - stop time
//        rmFile          - remove temporary file
//
// Return XML string
//

  if(inspXML!="") return inspXML;  	// xml is provided by user

  TString lalinspinj_exec;
  if(gSystem->Getenv("LALINSPINJ_EXEC")==NULL) {
    cout << "CWB::mdc::GenInspiralXML - Error : environment LALINSPINJ_EXEC is not defined!!!" << endl;exit(1);
  } else {
    lalinspinj_exec=TString(gSystem->Getenv("LALINSPINJ_EXEC"));
  }

  // test if command is well done (dummy exec)
  TString ofName = GetTemporaryFileName("inspiral","xml",inspDIR,true);

  TString cmdOptions = inspOptions;
  cmdOptions += " --output "+ofName;
  if(gps_start_time!=0) {
    cmdOptions += " --gps-start-time ";
    cmdOptions += gps_start_time;
  }
  if(gps_end_time!=0) {
    cmdOptions += " --gps-end-time ";
    cmdOptions += gps_end_time;
  }

  char cmd[1024];
  sprintf(cmd,"%s %s",lalinspinj_exec.Data(),cmdOptions.Data());
  cout << cmd << endl;
  int error = gSystem->Exec(cmd);
  if(rmFile) gSystem->Exec(TString("rm ")+ofName);  // remove temporary file

  if(error) {
    cout << endl;
    char help[1024];
    sprintf(help,"%s --help",lalinspinj_exec.Data());
    gSystem->Exec(help);
    cout << endl;
    cout << "CWB::mdc::GenInspiralXML - Error in exececuting :" << endl << endl;
    cout << cmd << endl << endl;
    exit(1);   
  }

  return ofName; 
}
#endif

//______________________________________________________________________________
TString 
CWB::mdc::GetTemporaryFileName(TString tag, TString ext, TString dir, bool mkdir) {
//
// Get a temporary file name 
//
// Input: tag   - tag string 
//        ext   - file extension tag
//        mkdir - if true create temporary dir
//


  gRandom->SetSeed(0);
  int rnID = int(gRandom->Rndm(13)*1.e9);

  if(dir=="/tmp") {		// use /tmp dir
    // create temporary dir
    UserGroup_t* uinfo = gSystem->GetUserInfo();
    TString uname = uinfo->fUser;
    dir="/tmp/"+uname;
    int error=0;
    if(mkdir) error=gSystem->Exec(TString("mkdir -p ")+dir);
    if(error) {
      cout << endl;
      cout << "CWB::mdc::GetTemporaryFileName - Error in exececuting :" << endl << endl;
      cout << TString(TString("mkdir -p ")+dir) << endl << endl;
      exit(1);   
    }
  }

  char fName[256];
  sprintf(fName,"%s/%s_%d.%s",dir.Data(),tag.Data(),rnID,ext.Data());

  return TString(fName); 
}

#ifdef _USE_LAL
//______________________________________________________________________________
void                                                                            
CWB::mdc::FindChirpInjectSignals(wavearray<double>& w, SimInspiralTable *injections, TString ifo) {
//                                                                                         
// Get inspiral for user define detectors                                                  
//                                                                                         
//                                                                                         
// Input:  injections  - structure which contains the parameters of waveforms                  
//                                                                                         
// Output: w           - wavearray which contains the waveform data                                
//                                                                                         

  double rad2deg = 180./TMath::Pi();
  double deg2rad = TMath::Pi()/180.;
  double fp,fx;
  double latitude,longitude,polarization;
  double theta,phi,psi;
  double gmst;

  wavearray<double> z;
  w=0; z=w;

  REAL8TimeSeries *hp=NULL;
  REAL8TimeSeries *hx=NULL;

  //TString ifoRef=net->ifoName[0];

  // get detector location
  int nIFO = net ? net->ifoListSize() : 0;
  int ifoId=-1;
  detector* pD=NULL;
  for(int n=0; n<nIFO; n++) if(ifo.CompareTo(net->ifoName[n])==0) ifoId=n;
  if(ifoId>=0) pD = net->getifo(ifoId);
  else         pD = new detector(const_cast<char*>(ifo.Data()));
  detectorParams dP = pD->getDetectorParams(); 
  REAL8 location[3] = {pD->Rv[0],pD->Rv[1],pD->Rv[2]};
  if(ifoId<0) delete pD;

  // loop over injections 
  SimInspiralTable *thisInj = NULL;
  for (thisInj = injections; thisInj; thisInj = thisInj->next)
  {                                                       
    REAL8 deltaT = 1./w.rate();
    if(waveName!="") strcpy(thisInj->waveform,waveName.Data());

    // get approximant ( see LALSimInspiral.c )
    Approximant approximant = (Approximant)XLALSimInspiralGetApproximantFromString(thisInj->waveform);
/*
    cout << endl << endl;
    cout << "---------------------------------------------------------------------------------" << endl;
    cout << "CWB::mdc::FindChirpInjectSignals" << endl;
    cout << "---------------------------------------------------------------------------------" << endl;
    cout << "thisSim->numrel_data : " << thisInj->numrel_data << endl;
    cout << "thisSim->waveform    : " << thisInj->waveform << endl;
    cout << "thisSim->amp_order   : " << thisInj->amp_order << endl;
    cout << "APPROXIMANT          : " << approximant << " (NR_hdf5=" << NR_hdf5 << ")" << endl;
    cout << endl;
*/
    if(approximant==NR_hdf5) {

#if LAL_VERSION_MAJOR >  6 || (LAL_VERSION_MAJOR ==  6 && \
   (LAL_VERSION_MINOR > 18 || (LAL_VERSION_MINOR == 18 && \
   (LAL_VERSION_MICRO >  0 || (LAL_VERSION_MICRO ==  0 && LAL_VERSION_DEVEL >= 1)))))     // LAL_VERSION >= 6.18.0.1

      // see LALInferenceReadData.c
      LALDict *LALpars=XLALCreateDict();
      XLALSimInspiralWaveformParamsInsertNumRelData(LALpars, thisInj->numrel_data);

      // add check if numrel_mode_min/max is empty then skip creating ModeArray
      LALValue* ModeArray = XLALSimInspiralCreateModeArray();

      int ell_min = thisInj->numrel_mode_min; // this should really only ever be numrel_mode_min=2
      int ell_max = thisInj->numrel_mode_max;

      if ((ell_min == 0) && (ell_max == 0)){ell_min=2;ell_max=LAL_SIM_L_MAX_MODE_ARRAY;}
      //loop over
      for (int ell=ell_min; ell<=ell_max; ell++){
          XLALSimInspiralModeArrayActivateAllModesAtL(ModeArray, ell);
      }

      XLALSimInspiralWaveformParamsInsertModeArray(LALpars, ModeArray);

      // see LALSimInspiral.c
      REAL8 f_min = XLALSimInspiralfLow2fStart(thisInj->f_lower, thisInj->amp_order, approximant);
/*
      cout << endl;
      printf("CWB::mdc::FindChirpInjectSignals - Injecting with approximant    = %d\n", approximant);
      printf("CWB::mdc::FindChirpInjectSignals - Injecting with numrel_data    = %s\n", thisInj->numrel_data);
      printf("CWB::mdc::FindChirpInjectSignals - Injecting with waveform       = %s\n", thisInj->waveform);
      printf("CWB::mdc::FindChirpInjectSignals - Injecting with f_min          = %f\n", f_min);
      printf("CWB::mdc::FindChirpInjectSignals - Injecting with f_lower        = %f\n", thisInj->f_lower);
      printf("CWB::mdc::FindChirpInjectSignals - Injecting with amp_order      = %d\n", thisInj->amp_order);
      printf("CWB::mdc::FindChirpInjectSignals - Injecting with mass1          = %f\n", thisInj->mass1);
      printf("CWB::mdc::FindChirpInjectSignals - Injecting with mass2          = %f\n", thisInj->mass2);
      printf("CWB::mdc::FindChirpInjectSignals - Injecting with spin1x         = %f\n", thisInj->spin1x);
      printf("CWB::mdc::FindChirpInjectSignals - Injecting with spin2x         = %f\n", thisInj->spin2x);
      printf("CWB::mdc::FindChirpInjectSignals - Injecting with spin1y         = %f\n", thisInj->spin1y);
      printf("CWB::mdc::FindChirpInjectSignals - Injecting with spin2y         = %f\n", thisInj->spin2y);
      printf("CWB::mdc::FindChirpInjectSignals - Injecting with spin1z         = %f\n", thisInj->spin1z);
      printf("CWB::mdc::FindChirpInjectSignals - Injecting with spin2z         = %f\n", thisInj->spin2z);
      printf("CWB::mdc::FindChirpInjectSignals - Injecting with distance       = %f\n", thisInj->distance);
      printf("CWB::mdc::FindChirpInjectSignals - Injecting with inclination    = %f\n", thisInj->inclination);
      printf("CWB::mdc::FindChirpInjectSignals - Injecting with coa_phase      = %f\n", thisInj->coa_phase);
      printf("CWB::mdc::FindChirpInjectSignals - Injecting with mode min      = %d\n", thisInj->numrel_mode_min);
      printf("CWB::mdc::FindChirpInjectSignals - Injecting with mode max      = %d\n", thisInj->numrel_mode_max);
      cout << endl;
*/
      if( thisInj->amp_order!=0 ) {
        cout << endl << "CWB::mdc::FindChirpInjectSignals - "
             << "Error : amp_order must be 0 for NR waveform  : "
             << "check inspiral parameters" << endl;
        exit(1);
      }

      // see LALSimInspiral.c
      int ret = XLALSimInspiralTD(&hp, &hx, thisInj->mass1*LAL_MSUN_SI, thisInj->mass2*LAL_MSUN_SI,
                                                thisInj->spin1x, thisInj->spin1y, thisInj->spin1z,
                                                thisInj->spin2x, thisInj->spin2y, thisInj->spin2z,
                                                thisInj->distance*LAL_PC_SI*1.0e6, thisInj->inclination,
                                                thisInj->coa_phase, 0., 0., 0.,
                                                deltaT, f_min, 0.,
                                                LALpars, approximant);

      XLALDestroyDict(LALpars);

      if( ret==XLAL_FAILURE ) {
        cout << endl << "CWB::mdc::GetInspiral - "
             << "Error in XLALInspiralTDWaveformFromSimInspiral(1) : "
             << "check inspiral parameters" << endl;
        exit(1);
      }
#else
      cout << endl << "CWB::mdc::GetInspiral - "
           << "Error - NR_hdf5 not supported for LAL_VERSION < 6.18.0.1" << endl;
      exit(1);
#endif

    } else {

#if LAL_VERSION_MAJOR >   6 || (LAL_VERSION_MAJOR ==  6 && \
   (LAL_VERSION_MINOR >  16 || (LAL_VERSION_MINOR == 16 && \
    LAL_VERSION_MICRO >=  1                             )))     // LAL_VERSION >= 6.16.1

      // check if source is produced by CWB::mdc::Posterior2XML
      TString source = TString(thisInj->source);
      if(source(0,5)=="#CWB:") {
        // extract version number
        TString version=source(source.Index(":",4)+1,source.Index(":",5)-source.Index(":",4)-1);
        // check version
        if(version.Atoi()<1) {
          cout << endl << "CWB::mdc::FindChirpInjectSignals - "
               << "Error: outdated version of XML produced by CWB::mdc::Posterior2XML, must be >= 1" << endl << endl;
          exit(1);
        }
      }

#if LAL_VERSION_MAJOR == 6 && LAL_VERSION_MINOR == 16 && LAL_VERSION_MICRO == 1	// used only for tagged version lalinference_o2 used for O2 posteriors 

      // NOTE: lambda1, lambda2 are provided with xml generated with CWB::mdc::Posterior2XML
      //       since there are not dedicated locations in the SimInspiralTable we use alpha,beta
      double lambda1 = source(0,5)=="#CWB:" ? thisInj->alpha : 0;
      double lambda2 = source(0,5)=="#CWB:" ? thisInj->beta  : 0;

      // see LALSimInspiral.c
      REAL8 f_min = XLALSimInspiralfLow2fStart(thisInj->f_lower, thisInj->amp_order, approximant);

      // get pn order from waveform
      int pn_order = XLALSimInspiralGetPNOrderFromString(thisInj->waveform);

      // see LALSimInspiral.c
      int ret = XLALSimInspiralTD(&hp, &hx, thisInj->coa_phase, deltaT,
                                  thisInj->mass1*LAL_MSUN_SI, thisInj->mass2*LAL_MSUN_SI,
                                  thisInj->spin1x, thisInj->spin1y, thisInj->spin1z,
                                  thisInj->spin2x, thisInj->spin2y, thisInj->spin2z,
                                  f_min, thisInj->f_final,
                                  thisInj->distance*LAL_PC_SI*1.0e6, 0., thisInj->inclination,
                                  lambda1, lambda2, NULL, NULL,
                                  thisInj->amp_order, pn_order, approximant);

#else
      // this function is reimplemented in the cWB code, see mdc::XLALInspiralTDWaveformFromSimInspiral
      int ret = this->XLALInspiralTDWaveformFromSimInspiral(&hp, &hx, thisInj, deltaT); 
#endif

      if( ret==XLAL_FAILURE ) {
        cout << endl << "CWB::mdc::GetInspiral - "
             << "Error in XLALInspiralTDWaveformFromSimInspiral(2) : "
             << "check inspiral parameters" << endl;
        exit(1);
      }
#else
      cout << endl << "CWB::mdc::GetInspiral - "
           << "Error - LAL_VERSION < 6.16.1 not supported " << endl;
      exit(1);
#endif

    }

    /* add the time of the injection at the geocentre to the
     * start times of the h+ and hx time series.  after this,
     * their epochs mark the start of those time series at the
     * geocentre. (see http://www.lsc-group.phys.uwm.edu/cgit/gstlal/tree/gstlal/gst/lal/gstlal_simulation.c) */

    XLALGPSAddGPS(&hp->epoch, &thisInj->geocent_end_time);
    XLALGPSAddGPS(&hx->epoch, &thisInj->geocent_end_time);

    wat::Time gpsbuf(w.start());
    wat::Time gpsinj(hp->epoch.gpsSeconds,hp->epoch.gpsNanoSeconds);
    wat::Time gpsoff=gpsinj;gpsoff-=gpsbuf;

    int    bininj = int(gpsoff.GetDouble()/deltaT);
    double binoff = fmod(gpsoff.GetDouble(),deltaT);

    bininj+=hp->data->length;
    if(bininj>w.size()) {
      cout << "CWB::mdc::FindChirpInjectSignals - Warning " << endl;
      cout << "          signal length " << hp->data->length*deltaT 
           << " (sec) is too high, signal is cutted to fit the buffer length" << endl; 
      cout << "          Decrease the injection rate and/or reduce the waveform length" << endl; 
      bininj=w.size()-1;
    }

    gmst = fmod(XLALGreenwichMeanSiderealTime(&thisInj->geocent_end_time), LAL_TWOPI);

    latitude=thisInj->latitude;
    longitude=thisInj->longitude;
    polarization=thisInj->polarization;
    theta = LAL_PI_2-latitude;  // LAL -> CWB
    theta*=rad2deg;
    longitude = fmod(longitude - gmst, LAL_TWOPI);
    if (longitude < 0) longitude += LAL_TWOPI;
    phi = longitude > 0 ? longitude : 2*TMath::Pi()+longitude;
    phi*= rad2deg;
    psi = polarization*rad2deg;

    fp = GetAntennaPattern(ifo, phi, theta, psi, "hp");
    fx = GetAntennaPattern(ifo, phi, theta, psi, "hx");

    // check if the injection length is greater of the buffer
    int length = hp->data->length<w.size() ? hp->data->length : w.size();
    z=0;
    for(int i=bininj;i>=0;i--) {
      int j=length-(bininj-i)-1;
      if(j<0) break; 
      z[i] = fp*hp->data->data[j]+fx*hx->data->data[j];  // detector wave projection
    }
    if(inspCLB!="") {
      int simulation_id = thisInj->bandpass;     // WARNING!!! used simTable->bandpass instead of simTable->simulation_id
      wat::Time wtime(thisInj->geocent_end_time.gpsSeconds, thisInj->geocent_end_time.gpsNanoSeconds);
      CalibrateInspiral(z, ifo, wtime.GetDouble(), simulation_id); // apply calibration if user has defined the calibration file inspCLB
    }
    // apply time shift respect to the ifoRef
    //double tShift = GetDelay(ifo,ifoRef,phi,theta);
    double tShift = XLALTimeDelayFromEarthCenter(location, thisInj->longitude, thisInj->latitude, &hp->epoch);
    // add to tShift the bin time offset
    tShift+=binoff; 
    TimeShift(z, tShift);
    for(int i=0;i<(int)w.size();i++) w[i]+=z[i];

    XLALDestroyREAL8TimeSeries( hp );
    XLALDestroyREAL8TimeSeries( hx );
    hp = NULL; hx = NULL;
  }
  return;
}
#endif

#ifdef _USE_LAL
//______________________________________________________________________________
wavearray<double>                                                                            
CWB::mdc::GetInspiral(TString pol, int gps_start_time, int gps_end_time) {
//                                                                                         
// Get inspiral waveform                                                  
//                                                                                         
//                                                                                         
// Input: pol             - wave polarization : hp, hx                
//        gps_start_time  - start time
//        gps_end_time    - stop time
//                                                                                         
// Output: w              - wavearray which contains the waveform data                                
//                                                                                         

  wavearray<double> w;
  w.rate(MDC_SAMPLE_RATE);

  if(pol!="hp" && pol!="hx") {
    fprintf( stderr, "CWB::mdc::GetInspiral - Warning : wrong polarization (enter hp or hx)\n");
    return w;                                                     
  }

  // start/end times 
  int gpsStartSec = gps_start_time;               
  int gpsEndSec = gps_end_time; 

  // injections 
  INT4 numInjections = 0;
  SimInspiralTable *injections = NULL;
  TString injectionFile = GenInspiralXML(0., 0., false);                           

  // set default debug level
  lal_errhandler = LAL_ERR_EXIT;

  XLALSetSilentErrorHandler();

  // read the injections 
  numInjections = SimInspiralTableFromLIGOLw(&injections, injectionFile.Data(), gpsStartSec, gpsEndSec);

  if ( numInjections == 0 ) {
    fprintf( stderr, "CWB::mdc::GetInspiral - Warning : No injections in specified time\n");
    return w;                                                     
  } else {                                                         
    fprintf( stdout, "CWB::mdc::GetInspiral : Read %d injection(s) from the file '%s'\n",  
        numInjections, injectionFile.Data() );                     
  }                                                                

  REAL8TimeSeries *hp=NULL;
  REAL8TimeSeries *hx=NULL;

  // loop over injections 
  SimInspiralTable *thisInj = NULL;
  for (thisInj = injections; thisInj; thisInj = thisInj->next)
  {                                                       
    REAL8 deltaT = 1./w.rate();
    if(waveName!="") strcpy(thisInj->waveform,waveName.Data());

    // get approximant ( see LALSimInspiral.c )
    Approximant approximant = (Approximant)XLALSimInspiralGetApproximantFromString(thisInj->waveform);
/*
    cout << endl << endl;
    cout << "---------------------------------------------------------------------------------" << endl;
    cout << "CWB::mdc::GetInspiral" << endl;
    cout << "---------------------------------------------------------------------------------" << endl;
    cout << "thisSim->numrel_data : " << thisInj->numrel_data << endl;
    cout << "thisSim->waveform    : " << thisInj->waveform << endl;
    cout << "thisSim->amp_order   : " << thisInj->amp_order << endl;
    cout << "APPROXIMANT          : " << approximant << " (NR_hdf5=" << NR_hdf5 << ")" << endl;
    cout << endl;
*/
    if(approximant==NR_hdf5) {

#if LAL_VERSION_MAJOR >  6 || (LAL_VERSION_MAJOR ==  6 && \
   (LAL_VERSION_MINOR > 18 || (LAL_VERSION_MINOR == 18 && \
   (LAL_VERSION_MICRO >  0 || (LAL_VERSION_MICRO ==  0 && LAL_VERSION_DEVEL >= 1)))))     // LAL_VERSION >= 6.18.0.1

      // see LALInferenceReadData.c
      LALDict *LALpars=XLALCreateDict();
      XLALSimInspiralWaveformParamsInsertNumRelData(LALpars, thisInj->numrel_data);

      // add check if numrel_mode_min/max is empty then skip creating ModeArray
      LALValue* ModeArray = XLALSimInspiralCreateModeArray();

      int ell_min = thisInj->numrel_mode_min; // this should really only ever be numrel_mode_min=2
      int ell_max = thisInj->numrel_mode_max;

      if ((ell_min == 0) && (ell_max == 0)){ell_min=2;ell_max=LAL_SIM_L_MAX_MODE_ARRAY;}
      //loop over
      for (int ell=ell_min; ell<=ell_max; ell++){
          XLALSimInspiralModeArrayActivateAllModesAtL(ModeArray, ell);

      }

      XLALSimInspiralWaveformParamsInsertModeArray(LALpars, ModeArray);

      // see LALSimInspiral.c
      REAL8 f_min = XLALSimInspiralfLow2fStart(thisInj->f_lower, thisInj->amp_order, approximant);
/*
      cout << endl;
      printf("CWB::mdc::GetInspiral - Injecting with approximant    = %d\n", approximant);
      printf("CWB::mdc::GetInspiral - Injecting with numrel_data    = %s\n", thisInj->numrel_data);
      printf("CWB::mdc::GetInspiral - Injecting with waveform       = %s\n", thisInj->waveform);
      printf("CWB::mdc::GetInspiral - Injecting with f_min          = %f\n", f_min);
      printf("CWB::mdc::GetInspiral - Injecting with f_lower        = %f\n", thisInj->f_lower);
      printf("CWB::mdc::GetInspiral - Injecting with amp_order      = %d\n", thisInj->amp_order);
      printf("CWB::mdc::GetInspiral - Injecting with mass1          = %f\n", thisInj->mass1);
      printf("CWB::mdc::GetInspiral - Injecting with mass2          = %f\n", thisInj->mass2);
      printf("CWB::mdc::GetInspiral - Injecting with spin1x         = %f\n", thisInj->spin1x);
      printf("CWB::mdc::GetInspiral - Injecting with spin2x         = %f\n", thisInj->spin2x);
      printf("CWB::mdc::GetInspiral - Injecting with spin1y         = %f\n", thisInj->spin1y);
      printf("CWB::mdc::GetInspiral - Injecting with spin2y         = %f\n", thisInj->spin2y);
      printf("CWB::mdc::GetInspiral - Injecting with spin1z         = %f\n", thisInj->spin1z);
      printf("CWB::mdc::GetInspiral - Injecting with spin2z         = %f\n", thisInj->spin2z);
      printf("CWB::mdc::GetInspiral - Injecting with distance       = %f\n", thisInj->distance);
      printf("CWB::mdc::GetInspiral - Injecting with inclination    = %f\n", thisInj->inclination);
      printf("CWB::mdc::GetInspiral - Injecting with coa_phase      = %f\n", thisInj->coa_phase);
      printf("CWB::mdc::GetInspiral - Injecting with mode min      = %d\n", thisInj->numrel_mode_min);
      printf("CWB::mdc::GetInspiral - Injecting with mode max      = %d\n", thisInj->numrel_mode_max);
      cout << endl;
*/
      if( thisInj->amp_order!=0 ) {
        cout << endl << "CWB::mdc::GetInspiral - "
             << "Error : amp_order must be 0 for NR waveform  : "
             << "check inspiral parameters" << endl;
        exit(1);
      }

      // see LALSimInspiral.c
      int ret = XLALSimInspiralTD(&hp, &hx, thisInj->mass1*LAL_MSUN_SI, thisInj->mass2*LAL_MSUN_SI,
                                                thisInj->spin1x, thisInj->spin1y, thisInj->spin1z,
                                                thisInj->spin2x, thisInj->spin2y, thisInj->spin2z,
                                                thisInj->distance*LAL_PC_SI*1.0e6, thisInj->inclination,
                                                thisInj->coa_phase, 0., 0., 0.,
                                                deltaT, f_min, 0.,
                                                LALpars, approximant);

      XLALDestroyDict(LALpars);

      if( ret==XLAL_FAILURE ) {
        cout << endl << "CWB::mdc::GetInspiral - "
             << "Error in XLALInspiralTDWaveformFromSimInspiral(3) : "
             << "check inspiral parameters" << endl;
        exit(1);
      }
#else
      cout << endl << "CWB::mdc::GetInspiral - "
           << "Error - NR_hdf5 not supported for LAL_VERSION < 6.18.0.1" << endl;
      exit(1);
#endif

    } else {

#if LAL_VERSION_MAJOR >   6 || (LAL_VERSION_MAJOR ==  6 && \
   (LAL_VERSION_MINOR >  16 || (LAL_VERSION_MINOR == 16 && \
    LAL_VERSION_MICRO >=  1                             )))     // LAL_VERSION >= 6.16.1

      // check if source is produced by CWB::mdc::Posterior2XML
      TString source = TString(thisInj->source);
      if(source(0,5)=="#CWB:") {
        // extract version number
        TString version=source(source.Index(":",4)+1,source.Index(":",5)-source.Index(":",4)-1);
        // check version
        if(version.Atoi()<1) {
          cout << endl << "CWB::mdc::FindChirpInjectSignals - "
               << "Error: outdated version of XML produced by CWB::mdc::Posterior2XML, must be >= 1" << endl << endl;
          exit(1);
        }
      }

#if LAL_VERSION_MAJOR == 6 && LAL_VERSION_MINOR == 16 && LAL_VERSION_MICRO == 1	// used only for tagged version lalinference_o2 used for O2 posteriors 

      // NOTE: lambda1, lambda2 are provided with xml generated with CWB::mdc::Posterior2XML
      //       since there are not dedicated locations in the SimInspiralTable we use alpha,beta
      double lambda1 = source(0,5)=="#CWB:" ? thisInj->alpha : 0;
      double lambda2 = source(0,5)=="#CWB:" ? thisInj->beta  : 0;

      // see LALSimInspiral.c
      REAL8 f_min = XLALSimInspiralfLow2fStart(thisInj->f_lower, thisInj->amp_order, approximant);

      // get pn order from waveform
      int pn_order = XLALSimInspiralGetPNOrderFromString(thisInj->waveform);

      // see LALSimInspiral.c
      int ret = XLALSimInspiralTD(&hp, &hx, thisInj->coa_phase, deltaT,
                                  thisInj->mass1*LAL_MSUN_SI, thisInj->mass2*LAL_MSUN_SI,
                                  thisInj->spin1x, thisInj->spin1y, thisInj->spin1z,
                                  thisInj->spin2x, thisInj->spin2y, thisInj->spin2z,
                                  f_min, thisInj->f_final,
                                  thisInj->distance*LAL_PC_SI*1.0e6, 0., thisInj->inclination,
                                  lambda1, lambda2, NULL, NULL,
                                  thisInj->amp_order, pn_order, approximant);

#else
      // this function is reimplemented in the cWB code, see mdc::XLALInspiralTDWaveformFromSimInspiral
      int ret = this->XLALInspiralTDWaveformFromSimInspiral(&hp, &hx, thisInj, deltaT); 
#endif

      if( ret==XLAL_FAILURE ) {
        cout << endl << "CWB::mdc::GetInspiral - "
             << "Error in XLALInspiralTDWaveformFromSimInspiral : "
             << "check inspiral parameters" << endl;
        exit(1);
      }
#else
      cout << endl << "CWB::mdc::GetInspiral - "
           << "Error - LAL_VERSION < 6.16.1 not supported " << endl;
      exit(1);
#endif

    }

    /* add the time of the injection at the geocentre to the
     * start times of the h+ and hx time series.  after this,
     * their epochs mark the start of those time series at the
     * geocentre. (see http://www.lsc-group.phys.uwm.edu/cgit/gstlal/tree/gstlal/gst/lal/gstlal_simulation.c) */

    XLALGPSAddGPS(&hp->epoch, &thisInj->geocent_end_time);
    XLALGPSAddGPS(&hx->epoch, &thisInj->geocent_end_time);

    // fill output wavearray
    if(pol=="hp") {
      w.resize(hp->data->length);
      wat::Time gps(hp->epoch.gpsSeconds,hp->epoch.gpsNanoSeconds);
      w.start(gps.GetDouble());
      for(int i=0;i<(int)w.size();i++) w[i] = hp->data->data[i];
    }
    if(pol=="hx") {
      w.resize(hx->data->length);
      wat::Time gps(hx->epoch.gpsSeconds,hx->epoch.gpsNanoSeconds);
      w.start(gps.GetDouble());
      for(int i=0;i<(int)w.size();i++) w[i] = hx->data->data[i];
    }

    XLALDestroyREAL8TimeSeries(hp);
    XLALDestroyREAL8TimeSeries(hx);
    hp = NULL; hx = NULL;

    break;	// break after the first waveform
  }

  if(inspXML=="") gSystem->Exec(TString("rm ")+injectionFile);  // remove injectionFile file

  while ( injections ) {
    SimInspiralTable *thisInj = NULL;
    thisInj = injections;
    injections = injections->next;
    LALFree( thisInj );
  }

  LALCheckMemoryLeaks();

  return w;
}
#endif

#ifdef _USE_LAL
// this function has been copied from LALInspiralWave.c (lalsuite-v6.54)
// we have chaged the line:
// REAL8 f_ref = 0.;
// with line:
// REAL8 f_ref = thisRow->f_final;
// to permits the use of the original f_ref stored in the posterior samples
// cWB stores f_ref into the thisRow->f_final
//
// new code has beed added to manage new approximants: SEOBNRv4P, SEOBNRv4PHM
// this is necessary because the LAL funcion XLALSimInspiralTD (LALSimInspiral.c)
// does a f_min correction which is not valid SEOBNRv4P, SEOBNRv4PHM
int 
CWB::mdc::XLALInspiralTDWaveformFromSimInspiral(
    REAL8TimeSeries **hplus,    /**< +-polarization waveform */
    REAL8TimeSeries **hcross,   /**< x-polarization waveform */
    SimInspiralTable *thisRow,  /**< row from the sim_inspiral table containing waveform parameters */
    REAL8 deltaT                /**< time step (s) */
    )
{
#if !(LAL_VERSION_MAJOR == 6 && LAL_VERSION_MINOR == 16 && LAL_VERSION_MICRO == 1) // not used for tagged version lalinference_o2 used for O2 posteriors
   int ret;
   LALPNOrder order;
   Approximant approximant;
   REAL8 phi0 = thisRow->coa_phase;
   REAL8 m1 = thisRow->mass1 * LAL_MSUN_SI;
   REAL8 m2 = thisRow->mass2 * LAL_MSUN_SI;
   REAL8 S1x = thisRow->spin1x;
   REAL8 S1y = thisRow->spin1y;
   REAL8 S1z = thisRow->spin1z;
   REAL8 S2x = thisRow->spin2x;
   REAL8 S2y = thisRow->spin2y;
   REAL8 S2z = thisRow->spin2z;
   REAL8 f_min = thisRow->f_lower;
   REAL8 f_ref = thisRow->f_final;	// new line 
//   REAL8 f_ref = 0.;			// old line
   REAL8 distance = thisRow->distance * LAL_PC_SI * 1e6;
   REAL8 inclination = thisRow->inclination;
   LALDict *params = XLALCreateDict();
   XLALSimInspiralWaveformParamsInsertPNAmplitudeOrder(params, thisRow->amp_order);

   /* get approximant */
   approximant = (Approximant)XLALSimInspiralGetApproximantFromString(thisRow->waveform);
   if ( (int) approximant == XLAL_FAILURE)
      XLAL_ERROR(XLAL_EFUNC);

   if (approximant == NR_hdf5) {
     char *filepath = thisRow->numrel_data;
     XLALSimInspiralWaveformParamsInsertNumRelData(params, filepath);
   }
   /* get phase PN order; this is an enum such that the value is twice the PN order */
   order = (LALPNOrder)XLALSimInspiralGetPNOrderFromString(thisRow->waveform);
   if ( (int) order == XLAL_FAILURE)
      XLAL_ERROR(XLAL_EFUNC);
   XLALSimInspiralWaveformParamsInsertPNPhaseOrder(params, order);

   /* note: the condition waveform already does tapering... ignore any request to do so get taper option */
   /* taper = XLALGetTaperFromString(thisRow->taper); */

   /* generate +,x waveforms */
#ifdef NEW_CODE_FOR_SEOBNRv4P_AND_SEOBNRv4PHM
   if (approximant == SEOBNRv4P || approximant == SEOBNRv4PHM) {
     ret = XLALSimInspiralChooseTDWaveform(hplus, hcross, m1, m2, S1x, S1y, S1z, S2x, S2y, S2z, distance, inclination, phi0, 0.0, 0.0, 0.0, deltaT, f_min, f_ref, params, approximant);
     XLALDestroyDict(params);
     if( ret == XLAL_FAILURE )
       XLAL_ERROR(XLAL_EFUNC);

     const double extra_time_fraction = 0.1; /* fraction of waveform duration to add as extra time for tapering */
     const double extra_cycles = 3.0; /* more extra time measured in cycles at the starting frequency */

     /* upper bound on the chirp time starting at f_min */
     double tchirp = XLALSimInspiralChirpTimeBound(f_min, m1, m2, S1z, S2z);

     /* extra time to include for all waveforms to take care of situations
      * where the frequency is close to merger (and is sweeping rapidly):
      * this is a few cycles at the low frequency */
     double textra = extra_cycles / f_min;

     /* condition the time domain waveform by tapering in the extra time
        * at the beginning and high-pass filtering above original f_min */
     XLALSimInspiralTDConditionStage1(*hplus, *hcross, extra_time_fraction * tchirp + textra, f_min);

     /* final tapering at the beginning and at the end to remove filter transients */

     /* waveform should terminate at a frequency >= Schwarzschild ISCO
        * so taper one cycle at this frequency at the end; should not make
        * any difference to IMR waveforms */
     double fisco = 1.0 / (pow(6.0, 1.5) * LAL_PI * (m1 + m2) * LAL_MTSUN_SI / LAL_MSUN_SI);
     XLALSimInspiralTDConditionStage2(*hplus, *hcross, f_min, fisco);
   } else {
     ret = XLALSimInspiralTD(hplus, hcross, m1, m2, S1x, S1y, S1z, S2x, S2y, S2z, distance, inclination, phi0, 0.0, 0.0, 0.0, deltaT, f_min, f_ref, params, approximant);
     XLALDestroyDict(params);
     if( ret == XLAL_FAILURE )
       XLAL_ERROR(XLAL_EFUNC);
   }
#else
/*
   // used to select the higher-order modes (or subleading modes) of the gravitational wave radiation
   // only for SEOBNRv4HM_ROM approximant
   int ell_min = thisRow->numrel_mode_min; // this should really only ever be numrel_mode_min=2
   int ell_max = thisRow->numrel_mode_max;
   if((ell_min>=1 && ell_max<=5)&&(ell_min<=ell_max)) {
     cout << endl << "CWB::mdc::XLALInspiralTDWaveformFromSimInspiral - select modes " << ell_min << " " << ell_max << endl;
     LALValue* ModeArray = XLALSimInspiralCreateModeArray();
     for(int ell=ell_min; ell<=ell_max; ell++) {
        if(ell==1) XLALSimInspiralModeArrayActivateMode(ModeArray, 2, -1);
        else       XLALSimInspiralModeArrayActivateMode(ModeArray, ell, -ell);
     }
     XLALSimInspiralWaveformParamsInsertModeArray(params, ModeArray);
   } else if(ell_min!=0 && ell_max!=0) {
     cout << endl << "CWB::mdc::XLALInspiralTDWaveformFromSimInspiral - error: select modes " << ell_min << " " << ell_max << endl;
     cout << "For SEOBNRv4_ROM the available modes are: 2,1,3,4,5 -> modes = [(2,-2),(2,-1),(3,-3),(4,-4),(5,-5)]" << endl << endl;
     exit(1);
   }
*/
   ret = XLALSimInspiralTD(hplus, hcross, m1, m2, S1x, S1y, S1z, S2x, S2y, S2z, distance, inclination, phi0, 0.0, 0.0, 0.0, deltaT, f_min, f_ref, params, approximant);
   XLALDestroyDict(params);
   if( ret == XLAL_FAILURE )
     XLAL_ERROR(XLAL_EFUNC);
#endif
#endif	// endif lalinference_o2

   return XLAL_SUCCESS;
}
#endif

#ifdef _USE_LAL
void 
CWB::mdc::Posterior2XML(TString sampleFile, TString xmlFile, TString options) {
//
// sampleFile         : input posterior sample file produced by PE
// xmlFile            : output xml file
// options
// --gps_start_time   : gps start time 
//                      if(gps_start_time<=0) the posterior time is used
//		        if gps_start_time>0 is integer (Ex: 123456789 without decimal point, NOT 123456789.) then 
//			   the GPS decimal part is taken from posterior time otherwise
//			   is the decimal part defined in gps_start_time (Ex: 123456789.123 -> 0.123) 	
// --gps_stop_time    : gps stop time 
// --time_step        : injection time step
// --seed             : if (<0) the parameters are read sequentially from posterior sample 
//                    : if (>0) the parameters are read randomly (uniformely with seed) from posterior sample 
// --waveform         : if declared it overwrite the approximant declared in the posteriors 
// --source           : if declared it is added to the SimInspiralTable->source (max length=20) 
// --ninjections      : number of xml entries (only for gps_stop_time<=0)
// --clb_file         : if defined then it is used to dump the calibration parameters (same format as for posteriors sample)
// --sample           : values [map,maxl]. If defined only the map/maxl sample is used to create the xml file
// --taper            : taper string to be included in the XML file. Valid values are:
//                      "TAPER_NONE","TAPER_START","TAPER_END","TAPER_STARTEND" -> see lalsimulation/src/LALSimInspiral.c 
// --f_lower          : if defined it overwrite the value defined in the sampleFile file
// --f_ref            : if defined it overwrite the value defined in the sampleFile file
// --pn_order         : if defined it overwrite the value defined in the sampleFile file
// --amp_order        : if defined it overwrite the value defined in the sampleFile file
//

  CWB::Toolbox TB;

  double gps_start_time=0, gps_stop_time=0, time_step = 0;
  int time_nsec = -1.;
  int seed = -1;
  int ninjections = 0;

  TString sgps_start_time = TB.getParameter(options,"--gps_start_time");
  if(sgps_start_time.IsDigit()) {
    gps_start_time = sgps_start_time.Atoi();
  } else if(sgps_start_time.IsFloat()) {
    gps_start_time = sgps_start_time.Atoi();
    sgps_start_time.Remove(0,sgps_start_time.Last('.'));
    time_nsec = int(1.e9*sgps_start_time.Atof());
  } else if(sgps_start_time!="") {
    cout << endl << "CWB::mdc::Posterior2XML - "
         << "Error : gps_start_time must be an interger number" << endl;
    return;
  }

  TString sgps_stop_time = TB.getParameter(options,"--gps_stop_time");
  if(sgps_stop_time.IsDigit()) {
    gps_stop_time = sgps_stop_time.Atoi();
  } else if(sgps_stop_time!="") {
    cout << endl << "CWB::mdc::Posterior2XML - "
         << "Error : gps_stop_time must be an interger number" << endl;
    return;
  }

  TString stime_step = TB.getParameter(options,"--time_step");
  if(stime_step.IsDigit()) {
    time_step = stime_step.Atoi();
  } else if(stime_step!="") {
    cout << endl << "CWB::mdc::Posterior2XML - "
         << "Error : time_step must be an interger number" << endl;
    return;
  }

  TString sseed = TB.getParameter(options,"--seed");
  if(sseed.IsDigit()) seed = sseed.Atoi();

  TString sninjections = TB.getParameter(options,"--ninjections");
  if(sninjections.IsDigit()) ninjections = sninjections.Atof();

  TString swaveform = TB.getParameter(options,"--waveform");

  TString ssource = TB.getParameter(options,"--source");
  if(ssource=="") ssource="NULL";
  if(ssource.Sizeof()>20) {
    cout << endl << "CWB::mdc::Posterior2XML - "
         << "Error : max source length is 20" << endl;
    return;
  }

  TString clb_file = TB.getParameter(options,"--clb_file");     // output calibration file
  TString extension = clb_file(clb_file.Last('.'),clb_file.Sizeof()-clb_file.Last('.')-1);
  if(extension!=".clb") {
    cout << "CWB::mdc::Posterior2XML - Error : Calibration file extension must be .clb" << endl;
    return;
  }  

  TString sample_opt = TB.getParameter(options,"--sample");
  if(sample_opt!="" && sample_opt!="map" && sample_opt!="maxl" && !sample_opt.IsDigit()) {
    cout << endl << "CWB::mdc::Posterior2XML - "
         << "Error : wrong --sample option type, must be map/maxl or an integer number >=0" << endl;
    return;
  }

  TString taper = TB.getParameter(options,"--taper");
  if(taper!="" && taper!="TAPER_NONE" && taper!="TAPER_START" && taper!="TAPER_END" && taper!="TAPER_STARTEND") {
    cout << endl << "CWB::mdc::Posterior2XML - "
         << "Error : wrong --taper option type, must be TAPER_NONE,TAPER_START,TAPER_END,TAPER_STARTEND" << endl;
    return;
  }

  TString sf_lower = TB.getParameter(options,"--f_lower");
  float f_lower=-1.;
  if(sf_lower.IsFloat()) {
    f_lower = sf_lower.Atof();
  } else if(sf_lower!="") {
    cout << endl << "CWB::mdc::Posterior2XML - "
         << "Error : f_lower must be a float number" << endl;
    return;
  }

  TString sf_ref = TB.getParameter(options,"--f_ref");
  float f_ref=-1.;
  if(sf_ref.IsFloat()) {
    f_ref = sf_ref.Atof();
  } else if(sf_ref!="") {
    cout << endl << "CWB::mdc::Posterior2XML - "
         << "Error : f_ref must be a float number" << endl;
    return;
  }

  TString spn_order = TB.getParameter(options,"--pn_order");
  int pn_order_opt=-1;
  if(spn_order.IsDigit()) {
    pn_order_opt = spn_order.Atoi();
  } else if(spn_order!="") {
    cout << endl << "CWB::mdc::Posterior2XML - "
         << "Error : pn_order must be an interger number" << endl;
    return;
  }
  if(pn_order_opt!=-1 && swaveform!="") {
    cout << endl << "CWB::mdc::Posterior2XML - "
         << "Error : only one of the following options can be defined: --pn_order,  --waveform" << endl;
    return;
  }

  TString samp_order = TB.getParameter(options,"--amp_order");
  int amp_order_opt=-1;
  if(samp_order.IsDigit()) {
    amp_order_opt = samp_order.Atoi();
  } else if(samp_order!="") {
    cout << endl << "CWB::mdc::Posterior2XML - "
         << "Error : amp_order must be an interger number" << endl;
    return;
  }

  if(gps_start_time>0) {
    if(gps_stop_time<gps_start_time) {
      cout << endl << "CWB::mdc::Posterior2XML - "
           << "Error : if gps_start_time>0 then gps_stop_time must be >= gps_start_time" << endl;
      return;
    }
    if(time_step<=0) {
      cout << endl << "CWB::mdc::Posterior2XML - "
           << "Error : if gps_start_time>0 then time_step must be > 0" << endl;
      return;
    }
    if(ninjections>0) {
      cout << endl << "CWB::mdc::Posterior2XML - "
           << "Error : if gps_start_time>0 then ninjections can not be declared" << endl;
      return;
    }
    // compute ninjections 
    ninjections = int((gps_stop_time-gps_start_time)/time_step);
  } else {
    if(ninjections<=0) {
      cout << endl << "CWB::mdc::Posterior2XML - "
           << "Error : if gps_start_time<=0 then --ninjections option must be > 0" << endl;
      return;
    }
  }

  ifstream in;
  in.open(sampleFile.Data(),ios::in);
  if(!in.good()) {
    cout << endl << "CWB::mdc::Posterior2XML - "
         << "Error Opening File : " << sampleFile.Data() << endl;
    return;
  }

  // get header line from posteriors files
  char hline[4*1024];
  in.getline(hline,4*1024);

  // converts new json header format into posterior_samples.dat format
  // the conversion is reported in:
  // https://git.ligo.org/lscsoft/pesummary/blob/master/pesummary/gw/file/standard_names.py
  // sline.ReplaceAll("json header item","posterior_samples.dat header item");
  TString sline = TString("\t")+hline+TString("\t");
  sline.ReplaceAll("#","\t");	// remove initial comment character presents in the O1/O2 versions
  sline.ReplaceAll(" ","\t");
  sline.ReplaceAll("\treference_frequency\t","\tf_ref\t");
  sline.ReplaceAll("\tlambda_1\t","\tlambda1\t");
  sline.ReplaceAll("\tlambda_2\t","\tlambda2\t");
  sline.ReplaceAll("\tluminosity_distance\t","\tdist\t");
  sline.ReplaceAll("\tmass_ratio\t","\tq\t");
  sline.ReplaceAll("\tchirp_mass\t","\tmc\t");
  sline.ReplaceAll("\trightascension\t","\tra\t");
  sline.ReplaceAll("\tdeclination\t","\tdec\t");
  sline.ReplaceAll("\tcos_theta_jn\t","\tcostheta_jn\t");
  sline.ReplaceAll("\ttilt_1\t","\ttilt1\t");
  sline.ReplaceAll("\ttilt_2\t","\ttilt2\t");
  sline.ReplaceAll("\tcos_tilt_1\t","\tcostilt1\t");
  sline.ReplaceAll("\tcos_tilt_2\t","\tcostilt2\t");
  sline.ReplaceAll("\tphi_12\t","\tphi12\t");
  sline.ReplaceAll("\ta_1\t","\ta1\t");
  sline.ReplaceAll("\ta_2\t","\ta2\t");
  sline.ReplaceAll("\tspin_1z\t","\ta1z\t");
  sline.ReplaceAll("\tspin_2z\t","\ta2z\t");
  sline.ReplaceAll("\ta_1\t","\ta1\t");
  sline.ReplaceAll("\ta_2\t","\ta2\t");
  sline.ReplaceAll("\tmarginalized_phase\t","\tphase_maxl\t");
  sline.ReplaceAll("\tgeocent_time\t","\ttime\t");
  sline.ReplaceAll("\tmarginalized_geocent_time\t","\ttime_maxl\t");
  sline.ReplaceAll("\tlog_likelihood\t","\tlogl\t");
  sline.ReplaceAll("\tlog_prior\t","\tlogprior\t");
  sline.ReplaceAll("\tH1_time\t","\th1_end_time\t");
  sline.ReplaceAll("\tL1_time\t","\tl1_end_time\t");
  sline.ReplaceAll("\tV1_time\t","\tv1_end_time\t");
  sline.ReplaceAll("\tG1_time\t","\tg1_end_time\t");
  sline=sline(1,sline.Sizeof()-3);	// remove initial and final \t
  if(sline!=hline) {
    cout << endl << "CWB::mdc::Posterior2XML - header converted " << endl;
    cout << endl << "before convertion:" << endl << hline << endl;
    cout << endl << "after  convertion:" << endl << sline << endl << endl;
  }
  strcpy(hline,sline.Data());

  std::map<TString, int> hsample;               // header sample
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
/*
  std::map<TString, int>::iterator iter;
  for (iter=hsample.begin(); iter!=hsample.end(); iter++) {
    cout << iter->first << "\t" << hsample[iter->first] << endl;
  }
*/

  // read samples from posteriors file
  std::vector<vector<double> > sample;
  int nsample=0;
  double value;
  int hsize=token->GetEntries();
  std::vector<double> vsample(hsample.size());
  while (1) {
    for(int i=0;i<hsize;i++) in >> vsample[i];
    vsample[hsize]=vsample[0]; 	// we move the first item to the end of the list (see comment above)
    if(!in.good()) break;
    sample.push_back(vsample); 
    nsample++;
  }
  in.close();
  cout << endl << "CWB::mdc::Posterior2XML - "
       << "Read\t" << nsample << "\tsamples from posterior samples file" << endl;

  // if sample_opt is map then get the map sample index
  int mapSampleID=-1;
  if(sample_opt=="map") {
    double post_max=-1e+20;
    for(int i=0;i<sample.size();i++) {
      double post = -1;
      if(hsample["logpost"]) {
        post = sample[i][hsample["logpost"]];
      } else if(hsample["logprior"] && hsample["logl"]) {
        post = sample[i][hsample["logprior"]] + sample[i][hsample["logl"]];
      } else if(hsample["post"]) {
        post = sample[i][hsample["post"]];
      } else {
        cout << endl << "CWB::mdc::Posterior2XML - "
             << "Error : none of values used by map option are present in the sample header" << endl;
        return;
      }
      if(post>post_max) {post_max=post;mapSampleID=i;}
    }
    cout << endl << "CWB::mdc::Posterior2XML - found map sample id = " << mapSampleID << endl << endl;
  }

  // if sample_opt is maxl then get the maxL sample index
  int maxlSampleID=-1;
  if(sample_opt=="maxl") {
    double maxl_max=-1e+20;
    for(int i=0;i<sample.size();i++) {
      double maxl = -1;
      if(hsample["logl"]) {
        maxl = sample[i][hsample["logl"]];
      } else {
        cout << endl << "CWB::mdc::Posterior2XML - "
             << "Error : logl is not present in the sample header" << endl;
        return;
      }
      if(maxl>maxl_max) {maxl_max=maxl;maxlSampleID=i;}
    }
    cout << endl << "CWB::mdc::Posterior2XML - found maxl sample id = " << maxlSampleID << endl << endl;
  }

  // if sample_opt is a number then get its sample index
  int SampleID=-1;
  if(sample_opt.IsDigit()) {
    SampleID=sample_opt.Atoi();
    if(SampleID<0 || SampleID>=sample.size()) {
      cout << endl << "CWB::mdc::Posterior2XML - "
           << "Error : sample index is not valid. Allowed valued for this posterior sample file are [0," << sample.size()-1 << "]" << endl;
      return;
    }
  }

  // write clb header to calibration file
  ofstream clb;
  if(clb_file!="") {
    if(!clb_file.Contains(".clb")) {
      cout << endl << "CWB::mdc::Posterior2XML - "
           << "Error : calibration file extention must be .clb" << endl;
      return;
    }
    clb.open(clb_file.Data(),ios::out);
    if(!clb.good()) {
      cout << endl << "CWB::mdc::Posterior2XML - "
           << "Error Opening Output Calibration File : " << clb_file.Data() << endl;
      return;
    }
    clb.precision(18);
    clb << "time\t";
    std::map<TString, int>::iterator iter;
    for (iter=hsample.begin(); iter!=hsample.end(); iter++) {
      if(TString(iter->first).Contains("spcal")) clb << iter->first << "\t";
    }
    clb << "simulation_id"; 
    clb << endl;
  }

  static LALStatus status;
  LIGOLwXMLStream  xmlfp;
  MetadataTable    injections;
  SimInspiralTable *simTable;                                   // see lalsuite include/lal/LIGOMetadataTables.h

  gps_start_time = int(gps_start_time/time_step)*time_step;     // force gps_start_time to be an multiple of time_step

  wat::Time gpsStartTime(gps_start_time);
  wat::Time gpsEndTime(gps_stop_time);
  wat::Time currentGpsTime(gps_start_time);
  wat::Time timeStep(time_step);

  // create the first injection 
  simTable = injections.simInspiralTable = (SimInspiralTable *) calloc( 1, sizeof(SimInspiralTable) );

//  if(seed>=0) gRandom->SetSeed(abs(seed));
  TRandom3 random(abs(seed));

  // loop over parameter generation until end time is reached 
  int ninj = 0;
  int simulation_id = 0;                        // used to syncronize the calibration entries with the xml entries
  int initialNanoSeconds=-1;
  while(1) {

    int id = (seed>=0) ? gRandom->Uniform(0,sample.size()-1) : ninj;
    if(sample_opt=="map")    id=mapSampleID;
    if(sample_opt=="maxl")   id=maxlSampleID;
    if(sample_opt.IsDigit()) id=SampleID;

    // get posterior sample parameters
    std::map<TString, double> psample = GetPsample(hsample, sample[id], f_lower, f_ref);
    if(psample.size()==0) {
      ninj++;
      if(seed>=0) continue;
      else {if(ninj<sample.size()) continue; else break;}
    }

    // store tag in table 
    char source_tag[30]; sprintf(source_tag,"#CWB:1:%s",ssource.Data());
    memcpy( simTable->source, source_tag, sizeof(CHAR) * LIGOMETA_SOURCE_MAX );

    // get time
    double time = psample["time"];

    // store time in table 
    wat::Time iwtime(time);
    if(gps_start_time>0) {
      simTable->geocent_end_time.gpsSeconds     = currentGpsTime.GetSec();
      if(time_nsec>=0) {
        simTable->geocent_end_time.gpsNanoSeconds = time_nsec;
      } else {
        simTable->geocent_end_time.gpsNanoSeconds = iwtime.GetNSec();
        // correct seconds when event times are gittering around the an integer (Ex: 1248242631.001, 1248242630.998, ...)
        if(initialNanoSeconds==-1) initialNanoSeconds=iwtime.GetNSec();
        int deltaNSec = int(iwtime.GetNSec())-initialNanoSeconds;
        if(deltaNSec> 500000000) simTable->geocent_end_time.gpsSeconds-=1;
        if(deltaNSec<-500000000) simTable->geocent_end_time.gpsSeconds+=1;
      }
    } else {
      simTable->geocent_end_time.gpsSeconds     = iwtime.GetSec();
      simTable->geocent_end_time.gpsNanoSeconds = iwtime.GetNSec();
    }
    wat::Time owtime(simTable->geocent_end_time.gpsSeconds, simTable->geocent_end_time.gpsNanoSeconds);

    // populate detector times (used by CBC pipelines)
    wat::Time h_time(psample["h1_end_time"]);
    simTable->h_end_time.gpsSeconds     = h_time.GetSec();
    simTable->h_end_time.gpsNanoSeconds = h_time.GetNSec();

    wat::Time l_time(psample["l1_end_time"]);
    simTable->l_end_time.gpsSeconds     = l_time.GetSec();
    simTable->l_end_time.gpsNanoSeconds = l_time.GetNSec();

    wat::Time v_time(psample["v1_end_time"]);
    simTable->v_end_time.gpsSeconds     = v_time.GetSec();
    simTable->v_end_time.gpsNanoSeconds = v_time.GetNSec();

    wat::Time g_time(psample["g1_end_time"]);
    simTable->g_end_time.gpsSeconds     = g_time.GetSec();
    simTable->g_end_time.gpsNanoSeconds = g_time.GetNSec();

    // populate location 
    simTable->longitude    = psample["ra"];
    simTable->latitude     = psample["dec"];
    simTable->polarization = psample["psi"];
    simTable->distance     = psample["dist"];

    if(gps_start_time>0) {
      skymap sm;
      double gps_source = time;
      double phi_source = sm.RA2phi(simTable->longitude*180./TMath::Pi(), gps_source);
      double ra_source  = sm.phi2RA(phi_source, owtime.GetDouble());
      if(ra_source>180) ra_source-=360.;
      simTable->longitude = ra_source*TMath::Pi()/180.;
    }

    // populate spins 
    simTable->spin1x = psample["s1x"];
    simTable->spin1y = psample["s1y"];
    simTable->spin1z = psample["s1z"];
    simTable->spin2x = psample["s2x"];
    simTable->spin2y = psample["s2y"];
    simTable->spin2z = psample["s2z"];

    // populate masses 
    simTable->mass1  = psample["m1"];
    simTable->mass2  = psample["m2"];
    simTable->mchirp = psample["mchirp"];

    // populate waveform and other parameters 
    int approx   = psample["approx"];
    int pn_order = psample["pn_order"];
    if(pn_order_opt>=0) pn_order = pn_order_opt; // it is overwritten by the value defined by user
    char waveform[256] = "";
    if(pn_order<0) {
      sprintf(waveform,"%s",XLALSimInspiralGetStringFromApproximant((Approximant)approx));
    } else {
      sprintf(waveform,"%s%s",XLALSimInspiralGetStringFromApproximant((Approximant)approx),XLALSimInspiralGetStringFromPNOrder((LALPNOrder)pn_order));
    }
    if(swaveform!="") sprintf(waveform,"%s",swaveform.Data());	// it is overwritten by the value defined by user
    memcpy( simTable->waveform, waveform, sizeof(CHAR) * LIGOMETA_WAVEFORM_MAX );
    //cout << "waveform : " << waveform << endl;

    simTable->f_lower = psample["flow"];
    if(f_lower>=0) simTable->f_lower = f_lower;	// it is overwritten by the value defined by user 
    simTable->f_final = psample["f_ref"];
    if(f_ref>=0) simTable->f_final = f_ref;	// it is overwritten by the value defined by user 
    simTable->amp_order = int(psample["amp_order"]);
    if(amp_order_opt>=0) simTable->amp_order = amp_order_opt; // it is overwritten by the value defined by user
    simTable->inclination = psample["iota"];
    simTable->coa_phase   = psample["phase"];

    // populate tidal parameters (since there are not dedicated locations in the simTable we use qmParameter1,qmParameter2)
    simTable->alpha = psample["lambda1"];
    simTable->beta  = psample["lambda2"];

    // populate logl and prior parameters. Such parameters are saved in simTable alpha1/2 unused locations. 
    // The values alpha1/2 are then stored in the output mdc log string and can be use as cross check
    // to verify the the distribution of reconstructed samples wrt the injected samples
    simTable->alpha1 = psample["logl"];
    simTable->alpha2 = psample["post"];

//    simTable->event_id = (EventIDColumn *)LALCalloc(1, sizeof(EventIDColumn) );
//    memcpy( simTable->event_id->textId, "CALIBRATION=1234.56789", sizeof(CHAR) * LIGOMETA_UNIQUE_MAX );
//    simTable->event_id->next = NULL;

    // simTable->taper is empty then it initialized to TAPER_START (used in FD)
    if(TString(simTable->taper)=="") memcpy( simTable->taper, "TAPER_START", sizeof(CHAR) * LIGOMETA_INSPIRALTAPER_MAX );
    // if taper is define by user then it overwrite the simTable->taper value
    if(taper!="") memcpy( simTable->taper, taper.Data(), sizeof(CHAR) * LIGOMETA_INSPIRALTAPER_MAX );

    ninj++;

    bool isBreak=false;
    if(gps_start_time>0) {
      currentGpsTime+=timeStep;
      if((currentGpsTime>gpsEndTime) || (seed<0 && ninj>=sample.size())) isBreak=true;
    } else {
      if(ninj>=ninjections) isBreak=true;
    }

    // dump calib parameters to calibration file 
    if(clb_file!="") {
      clb << std::scientific;
      clb << owtime.GetDouble() << "\t";
      std::map<TString, int>::iterator iter;
      for (iter=hsample.begin(); iter!=hsample.end(); iter++) {
        if(TString(iter->first).Contains("spcal")) clb << psample[iter->first] << "\t";
      }
      simTable->bandpass = simulation_id;       // WARNING!!! used simTable->bandpass instead of simTable->simulation_id 
                                                // because simulation_id is not defined in lalsuit branch lalinference_o2
                                                // simTable->bandpass is not used in CWB::mdc::FindChirpInjectSignals
      clb << simulation_id++; 
      clb << endl;
    }


    if(isBreak) break;

    // allocate and go to next SimInspiralTable 
    simTable = simTable->next = (SimInspiralTable *) calloc( 1, sizeof(SimInspiralTable) );
  }

  // close clb file
  if(clb_file!="") clb.close();

  // write xml file
  memset( &xmlfp, 0, sizeof(LIGOLwXMLStream) );
  LAL_CALL( LALOpenLIGOLwXMLFile( &status, &xmlfp, xmlFile.Data() ), &status );
  if(injections.simInspiralTable) {
    LAL_CALL( LALBeginLIGOLwXMLTable( &status, &xmlfp, sim_inspiral_table ), &status );
    LAL_CALL( LALWriteLIGOLwXMLTable( &status, &xmlfp, injections, sim_inspiral_table ), &status );
    LAL_CALL( LALEndLIGOLwXMLTable ( &status, &xmlfp ), &status );
  }
  LAL_CALL( LALCloseLIGOLwXMLFile ( &status, &xmlfp ), &status );
  LALCheckMemoryLeaks();
  cout << endl << "CWB::mdc::Posterior2XML - "
       << "Write\t" << ninj << "\tsamples to output xml file" << endl;

}
#endif

#ifdef _USE_LAL
std::map<TString, double> 
CWB::mdc::GetPsample(std::map<TString, int> hsample, vector<double> sample, float if_lower, float if_ref) {

  // ./lalapps/src/inspiral/inspinj.c

  std::map<TString, double> psample;
  std::map<TString, double> psample_null;

  bool check=false;

  int pn_order    = hsample["lal_pnorder"] ? int(sample[hsample["lal_pnorder"]]) : -1;
  int amp_order   = hsample["lal_amporder"] ? int(sample[hsample["lal_amporder"]]) : -1;
  int approx      = hsample["lal_approximant"] ? int(sample[hsample["lal_approximant"]]) : IMRPhenomPv2;  

  psample["pn_order"]    = pn_order;
  psample["amp_order"]   = amp_order;
  psample["approx"]      = approx;

  //cout << "TEST APPROXIMANT " << approx << " " << NumApproximants << " " << IMRPhenomPv2 << endl;

  double flow  = hsample["flow"] ? sample[hsample["flow"]] : 30.;
  double f_ref = hsample["f_ref"] ? sample[hsample["f_ref"]] : 30.;
#ifdef NEW_CODE_FOR_SEOBNRv4P_AND_SEOBNRv4PHM
  flow  = 2*flow /(2+psample["amp_order"]);
  f_ref = 2*f_ref/(2+psample["amp_order"]);
  psample["amp_order"] = 0;	// is set to 0 because flow and f_ref already renormalized ???
#endif
  if(if_lower>=0) flow = if_lower;	// it is overwritten by the value defined by user 
  if(if_ref>=0)   f_ref = if_ref;	// it is overwritten by the value defined by user

  psample["flow"]  = flow;
  psample["f_ref"] = f_ref;

  double lambda1  = hsample["lambda1"] ? sample[hsample["lambda1"]] : 0.;
  double lambda2  = hsample["lambda2"] ? sample[hsample["lambda2"]] : 0.;

  psample["lambda1"]  = lambda1;
  psample["lambda2"]  = lambda2;

  psample["logl"] = hsample["logl"] ? sample[hsample["logl"]] : 0.;
  if(hsample["logpost"]) {
    psample["post"] = sample[hsample["logpost"]];
  } else if(hsample["logprior"] && hsample["logl"]) {
    psample["post"] = sample[hsample["logprior"]] + sample[hsample["logl"]];
  } else if(hsample["post"]) {
    psample["post"] = sample[hsample["post"]];
  } else psample["post"] = 0.;

  double dist = 0;
  int    dist_check=0;
  if(hsample["distance"])    {dist = sample[hsample["distance"]];dist_check++;}
  if(hsample["dist"])        {dist = sample[hsample["dist"]];dist_check++;}
  if(hsample["logdistance"]) {dist = TMath::Exp(sample[hsample["logdistance"]]);dist_check++;}
  if(dist_check>1) {
    cout << endl << "CWB::mdc::Posterior2XML - "
         << "multiple definitions of distance in sample header " << endl;
    return psample_null;
  }
  if(dist==0) {
    cout << endl << "CWB::mdc::Posterior2XML - "
         << "distance not defined in sample header " << endl;
    return psample_null;
  }
  psample["dist"] = dist;

  double q        = sample[hsample["q"]];
  double factor   = sample[hsample["mc"]] * pow(1 + q, 1.0/5.0);

  double m1       = factor * pow(q, -3.0/5.0);
  double m2       = factor * pow(q, 2.0/5.0);

  double q1       = hsample["q1"] ? sample[hsample["q1"]] : 0.;
  double q2       = hsample["q2"] ? sample[hsample["q2"]] : 0.;

  double ra       = sample[hsample["ra"]];
  double dec      = sample[hsample["dec"]];

  double psi      = hsample["psi"] ? sample[hsample["psi"]] : sample[hsample["polarization"]];

  psample["ra"]  = ra;
  psample["dec"] = dec;
  psample["psi"] = psi;

  psample["q1"] = q1;
  psample["q2"] = q2;

  psample["m1"] = m1;
  psample["m2"] = m2;
  psample["mchirp"] = pow(m1*m2,3./5.)/pow(m1+m2,1./5.);

  // These are identical if non-precessing, in which case this will be overwritten
  double iota = -1;
  if(hsample["costheta_jn"]) {
    iota = TMath::ACos(sample[hsample["costheta_jn"]]);
  } else if(hsample["theta_jn"]) {
    iota = sample[hsample["theta_jn"]];
  } else {
    cout << endl << "CWB::mdc::Posterior2XML - "
         << "Error: costheta_jn or theta_jn not defined in posterior samples header file" << endl;
    return psample_null;
  }

  double phase;
  static bool phase_warning=true;
  if(hsample["phi_orb"]) {
    phase = sample[hsample["phi_orb"]];
  } else {
    if(hsample["phase"]) {
      phase = sample[hsample["phase"]];
    } else {
      phase = sample[hsample["phase_maxl"]];
      if(phase_warning) {
        cout << "WARNING: Samples already marginalized over phase." << endl;
        cout << "         Begrudgingly using max loglikelihood estimates." << endl;
        phase_warning=false;
      }
    }
  }
  psample["phase"] = phase;

  double s1x=0, s1y=0, s1z = 0.;
  double s2x=0, s2y=0, s2z = 0.;

  // Try to treat this as a precessing analysis
  double theta_jn, phi_jl, tilt1, tilt2, phi12, a1, a2;
  if((hsample["theta_jn"]||hsample["costheta_jn"]) && hsample["phi_jl"]) {

    theta_jn = hsample["theta_jn"] ? sample[hsample["theta_jn"]] : TMath::ACos(sample[hsample["costheta_jn"]]);
    phi_jl   = sample[hsample["phi_jl"]];

    if(hsample["tilt1"] && hsample["tilt2"]) {
      tilt1 = sample[hsample["tilt1"]];
      tilt2 = sample[hsample["tilt2"]];
    } else {
      tilt1 = TMath::ACos(sample[hsample["costilt1"]]);
      tilt2 = TMath::ACos(sample[hsample["costilt2"]]);
    }

    phi12 = sample[hsample["phi12"]];
    a1    = sample[hsample["a1"]];
    a2    = sample[hsample["a2"]];


    // routines for transforming initial conditions of precessing waveforms 
    // include/lal/LALSimInspiral.h
    double phiRef=phase;
#if LAL_VERSION_MAJOR == 6 && LAL_VERSION_MINOR == 16 && LAL_VERSION_MICRO == 1	// used only for tagged version lalinference_o2 used for O2 posteriors 
    int ret = XLALSimInspiralTransformPrecessingNewInitialConditions(&iota, &s1x, &s1y, &s1z, &s2x, &s2y, &s2z, theta_jn, 
                                                                     phi_jl, tilt1, tilt2, phi12, a1, a2, m1*LAL_MSUN_SI, m2*LAL_MSUN_SI, f_ref);      // version used by ben farr
#else
    int ret = XLALSimInspiralTransformPrecessingNewInitialConditions(&iota, &s1x, &s1y, &s1z, &s2x, &s2y, &s2z, theta_jn, 
                                                                     phi_jl, tilt1, tilt2, phi12, a1, a2, m1*LAL_MSUN_SI, m2*LAL_MSUN_SI, f_ref, phiRef);     // new version
#endif
    if( ret==XLAL_FAILURE ) {
      cout << endl << "CWB::mdc::Posterior2XML - "
           << "Error in XLALSimInspiralTransformPrecessingNewInitialConditions" << endl;
      return psample_null;
    }
    //cout << "SPINS : " << iota << " " << s1x << " " << s1y << " " << s1z << " " << s2x << " " << s2y << " " << s2z << endl;

  } else {

    if(hsample["a1z"] && hsample["a2z"]) {
      s1z = sample[hsample["a1z"]];
      s2z = sample[hsample["a2z"]];
    } else {
      if(hsample["a1"] && hsample["a2"]) {
        s1z = sample[hsample["a1"]];
        s2z = sample[hsample["a2"]];
      } else {
        cout << endl << "CWB::mdc::Posterior2XML - "
             << "Error: It must be non-spinning !!!" << endl;
      }
    }
  }

  psample["iota"] = iota;

  psample["s1x"] = s1x;
  psample["s1y"] = s1y;
  psample["s1z"] = s1z;
  psample["s2x"] = s2x;
  psample["s2y"] = s2y;
  psample["s2z"] = s2z;

  double time;
  static bool time_warning=true;
  if(hsample["time"]) {
    time = sample[hsample["time"]];
  } else {
    time = sample[hsample["time_maxl"]];
    if(time_warning) {
      cout << "WARNING: Samples already marginalized over time." << endl;
      cout << "         Begrudgingly using max loglikelihood estimates." << endl;
      time_warning=false;
    }
  }
  psample["time"] = time;

  if(hsample["h1_end_time"]) psample["h1_end_time"]=sample[hsample["h1_end_time"]]; else psample["h1_end_time"]=0.;
  if(hsample["l1_end_time"]) psample["l1_end_time"]=sample[hsample["l1_end_time"]]; else psample["l1_end_time"]=0.;
  if(hsample["v1_end_time"]) psample["v1_end_time"]=sample[hsample["v1_end_time"]]; else psample["v1_end_time"]=0.;
  if(hsample["g1_end_time"]) psample["g1_end_time"]=sample[hsample["g1_end_time"]]; else psample["g1_end_time"]=0.;

  // extract calibration parameters
  std::map<TString, int>::iterator iter;
  for (iter=hsample.begin(); iter!=hsample.end(); iter++) {
    if(TString(iter->first).Contains("spcal")) psample[iter->first] = sample[hsample[iter->first]];
  }

  return psample;
}
#endif

#ifdef _USE_LAL
//______________________________________________________________________________
void                                                                            
CWB::mdc::CalibrateInspiral(wavearray<double>& w, TString ifo, double time, int simulation_id) {
//                                                                                         
// Apply calibration to the detector waveforms                                                 
//                                                                                         
//                                                                                         
// Input:  ifo                 - ifo label             
//         time/simulation_id  - GPS time/simulation_id use to extract the calibration parameters 
//                                                                                         
// Output: w                   - wavearray which contains the calibrated waveform data                                
//                                                                                         

  if(inspCLB=="") return;

  TString extension = inspCLB(inspCLB.Last('.'),inspCLB.Sizeof()-inspCLB.Last('.')-1);
  if(extension!=".clb") {
    cout << "CWB::mdc::CalibrateInspiral - Error : Calibration file extension must be .clb" << endl;
    exit(1);
  }  

  ifstream in;
  in.open(inspCLB.Data(),ios::in);
  if(!in.good()) {
    cout << endl << "CWB::mdc::CalibrateInspiral - "
         << "Error Opening Calibration File : " << inspCLB.Data() << endl;
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
  while (1) {
    for(int i=0;i<hsize;i++) in >> sample[i];
    sample[hsize]=sample[0]; 	// we move the first item to the end of the list (see comment above)
    if(!in.good()) break;
    if((fabs(time-sample[hsample["time"]])<0.001)&&(simulation_id==sample[hsample["simulation_id"]])) {spcal=true;break;}
  }
  in.close();

  if(spcal) {
    printf("\nCWB::mdc::CalibrateInspiral - Calibration Found @ GPS time %10.6f\n",sample[hsample["time"]]);
  } else {
    return;
  }

  wavearray<double> logfreq;
  wavearray<double> amp;
  wavearray<double> phase;

  TString ifo_label=ifo;
  ifo_label.ToLower();
  if(sample[hsample["spcal_active"]]) {
    int spcal_npts = sample[hsample["spcal_npts"]];
    //cout << "spcal_npts : " << spcal_npts << endl;
    logfreq.resize(spcal_npts);
    amp.resize(spcal_npts);
    phase.resize(spcal_npts);
    for(int j=0;j<spcal_npts;j++) {
      logfreq[j]  = log(sample[hsample[TString::Format("%s_spcal_freq_%d",ifo_label.Data(),j)]]);
      amp[j]      = sample[hsample[TString::Format("%s_spcal_amp_%d",ifo_label.Data(),j)]];
      phase[j]    = sample[hsample[TString::Format("%s_spcal_phase_%d",ifo_label.Data(),j)]];
    }
  }

  //for(int j=0;j<logfreq.size();j++) {
  //  cout << ifo << "\t" << j << "\t" << exp(logfreq[j]) << "\t" << amp[j] << "\t" << phase[j] << endl;
  //}


  TSpline3 spamp("spamp",logfreq.data,amp.data,logfreq.size());
  TSpline3 spphs("spphs",logfreq.data,phase.data,logfreq.size());

  // apply calibration to waveform 
  w.FFTW(1);
  TComplex C;
  double df = w.rate()/w.size();
  for (int i=0;i<(int)w.size()/2;i++) {
    double logF = log(i*df);
    if(logF<logfreq[0] || logF>logfreq[logfreq.size()-1]) continue;     // exclude points outside the input freq range
    double dA = spamp.Eval(logF); 
    double dP = spphs.Eval(logF); 
    //cout << "1 -> " << exp(logF) << "\t" << dA << "\t" << dP << endl; 
    //TComplex cal_factor = (1.0+dA)*TComplex(2.0,dP)/TComplex(2.0,-dP);        // exp approximation used in gw_reconstruct package (Ben Farr))
    TComplex cal_factor = (1.0+dA)*C.Exp(TComplex(0.,dP));
    //cout << "2 -> " << exp(logF) << "\t" << cal_factor.Re() << "\t" << cal_factor.Im() << endl; 
    TComplex W(w[2*i],w[2*i+1]);
    W=W*cal_factor;  // calibrate
    w[2*i]=W.Re();
    w[2*i+1]=W.Im();
  }
  w.FFTW(-1);
}
#endif

//______________________________________________________________________________
wavearray<double>
CWB::mdc::GetXCorr(wavearray<double>& w1, wavearray<double>& w2) {
//
// Compute the cross correlation wavearray from w1,w2
//
//
// Input:  w1    - wavearray
//         w2    - wavearray
//
// return cross correlation wavearray

  if(w1.start()!=w2.start()) {
    cout << "CwB::mdc::GetXCorr error: start w1 != start w2" << endl;
    wavearray<double> w;
    return w;	// return empty wavearray
  }
  if(w1.size()!=w2.size()) {
    cout << "CwB::mdc::GetXCorr error: size w1 != size w2" << endl;
    wavearray<double> w;
    return w;	// return empty wavearray
  }

  int size = w1.size();

  wavearray<double> W1=w1; W1.FFTW(1);
  wavearray<double> W2=w2; W2.FFTW(1);

  wavearray<double> X(size);
  for(int i=0;i<size;i+=2) {
    complex<double> A(W1[i], W1[i+1]);
    complex<double> B(W2[i],-W2[i+1]); // conjugate
    complex<double> C = A*B;
    X[i]   = C.real();
    X[i+1] = C.imag();
  }
  X.FFTW(-1);

  return X;
}

//______________________________________________________________________________
wavearray<double>
CWB::mdc::GetAligned(wavearray<double>* w1, wavearray<double>* w2) {
//
// align w1 wrt w2
//
//
// Input:  w1    - wavearray
//         w2    - wavearray
//
// Return aligned w
//
// Input w1,w2:
// 
//                   |------- w1 ------|
// |---------- w2 ------------|
//
// Output w:
//
// |000000000000000000-- w1 --|

   wavearray<double> w = *w2;
   w=0;

   if(w1==NULL)      return w;
   if(w1->size()==0) return w;

   double R=w1->rate();

   double b_w1 = w1->start();
   double e_w1 = w1->start()+w1->size()/R;
   double b_w2 = w2->start();
   double e_w2 = w2->start()+w2->size()/R;

   int o_w1 = b_w1>b_w2 ? 0 : int((b_w2-b_w1)*R+0.5);
   int o_w2 = b_w1<b_w2 ? 0 : int((b_w1-b_w2)*R+0.5);

   double start = b_w1>b_w2 ? b_w1 : b_w2;
   double stop  = e_w1<e_w2 ? e_w1 : e_w2;
   int    size  = int((stop-start)*R+0.5);

   for(int i=0;i<size;i++) w[i+o_w2] = w1->data[i+o_w1];

   return w;
}

//______________________________________________________________________________
wavearray<double>
CWB::mdc::GetAdd(wavearray<double>* w1, wavearray<double>* w2) {
//
// compute w1+w2
//
//
// Input:  w1    - wavearray
//         w2    - wavearray
//
// Return  w1-w2
//
// Input w1,w2:
//
//                   |------- w1 ------|
// |---------- w2 ---------------|
//
// Output w:
//
//                   |__ w1+w2 __|

  wavearray<double> x = *w2;
  x*=-1;
  return GetDiff(w1, &x);
}

//______________________________________________________________________________
wavearray<double>
CWB::mdc::GetDiff(wavearray<double>* w1, wavearray<double>* w2) {
//
// compute w1-w2
//
//
// Input:  w1    - wavearray
//         w2    - wavearray
//
// Return  w1-w2
//
// Input w1,w2:
// 
//                   |------- w1 ------|
// |---------- w2 ---------------|
//
// Output w:
//
//                   |__ w1-w2 __|

   double R=w1->rate();

   double b_w1 = w1->start();
   double e_w1 = w1->start()+w1->size()/R;
   double b_w2 = w2->start();
   double e_w2 = w2->start()+w2->size()/R;

   int o_w1 = b_w1>b_w2 ? 0 : int((b_w2-b_w1)*R+0.5);
   int o_w2 = b_w1<b_w2 ? 0 : int((b_w1-b_w2)*R+0.5);

   double start = b_w1>b_w2 ? b_w1 : b_w2;
   double stop  = e_w1<e_w2 ? e_w1 : e_w2;
   int    size  = int((stop-start)*R+0.5);

   wavearray<double> w(size);
   w=0.;
   w.rate(R);
   w.start(b_w1+double(o_w1)/R);

   for(int i=0;i<size;i++) w[i] = w1->data[i+o_w1] - w2->data[i+o_w2];

   return w;
}

//______________________________________________________________________________
double
CWB::mdc::GetTimeBoundaries(wavearray<double> x, double P, double& bT, double& eT, double T, double Q) {
//
// compute the symmetric time interval including the energy fraction (1-P) 
//
//
// Input:  x     - wavearray
//         P     - energy fraction  (0:1)
//         T     - reference time 
//         Q     - left energy fraction  (0:1), if Q<0 then Q=P
//
// Output: bT    - interval begin time
//         eT    - interval end time
//
// Return eT-bT
//
// if T not belongs to [tstart,tstop] then
//
//    |------------ET-----------|
//    |xxx-------------------xxx|
//     bT <-               -> eT
//
//    where xxx is (1-P)/2 fraction of total energy ET
//
// if T belongs to [tstart,tstop] then
//
//    |---------EL------|--ER---|
//    |xxxx-------------|-----xx|
//      bT <-           T   -> eT
//
//    where xxx = EL*(1-P) : EL 
//    where xx  = ER*(1-Q) : ER
//    

  if(P<0) P=0;
  if(P>1) P=1;

  if(Q<0) Q=P;
  if(Q>1) Q=1;

  int    size  = x.size();
  double dt    = 1./x.rate();
  double start = x.start();
  double stop  = x.start()+size*dt;

  double E = 0;
  double ET = 0;		// total energy
  for(int i=0;i<size;i++) {ET+=x[i]*x[i];}

  double EL=0;
  double ER=0;
  if((T>=start)&&(T<=stop)) {
    for(int i=0;i<size;i++) {
      double t=i*dt+start;
      if(t<T) EL+=x[i]*x[i]; else ER+=x[i]*x[i];
    } 
  } else {
    EL=ER=ET/2;
  }

  double el = EL*(1-P);	// left  side energy
  double er = ER*(1-Q);	// right side energy

  // search jB sample index which contains EL energy in the left side of the array
  E=0;
  int jB=0;
  for(int i=0;i<size;i++) {
    E+=x[i]*x[i];
    if(E>el) break;
    jB=i;
  }

  // search jE sample index which contains ER energy in the right side of the array
  E=0;
  int jE=size-1;
  for(int i=size-1;i>0;i--) {
    E+=x[i]*x[i];
    if(E>er) break;
    jE=i;
  }

  bT = start+jB*dt;
  eT = start+jE*dt;

  return eT-bT;
}

//______________________________________________________________________________
double
CWB::mdc::PhaseSync(wavearray<double>& w1, wavearray<double>& w2, double& sync_phase) {
//
// wavearray w1 is synchronized in phase wrt w2
//
//
// Input:  w1    - wavearray to be synchronized wrt w2
//         w2    - reference wavearray
//
// Output: sync_phase - synchronization phase
//
// return xcor after synchronization

  int size = w1.size();

  // get 90 deg phase shift wavearray
  wavearray<double> W1=w1; CWB::mdc::PhaseShift(W1,90.);
  wavearray<double> W2=w2; CWB::mdc::PhaseShift(W2,90.);

  // compute sync_phase
  double num=0;
  double den=0;
  for(int i=0;i<(int)w1.size();i++) {
     num+=w1[i]*W2[i]-W1[i]*w2[i];
     den+=w1[i]*w2[i]+W1[i]*W2[i];
  }
  sync_phase = TMath::ATan2(num,den);
  sync_phase*=180./TMath::Pi();
  CWB::mdc::PhaseShift(w1,-sync_phase);

  // compute sync_xcor w1*w2
  double sync_xcor=0.0;
  for(int i=0;i<size;i++) sync_xcor+=w1[i]*w2[i];
  sync_xcor*=1./(w1.rms()*sqrt(size));
  sync_xcor*=1./(w2.rms()*sqrt(size));

  return sync_xcor;
}

//______________________________________________________________________________
void
CWB::mdc::PhaseSync(wavearray<double>& w, double sync_phase) {
//
// wavearray w1 is shifted in time and phase 
//
//
// Input:  w1         - wavearray 
//         sync_phase - phase shift
//
// Output: w1         - time/phase shifted wavearray

  CWB::mdc::PhaseShift(w,-sync_phase);
}

//______________________________________________________________________________
double
CWB::mdc::TimeSync(wavearray<double>& w1, wavearray<double>& w2, double& sync_time) {
//
// wavearray w1 is synchronized in time wrt w2
//
//
// Input:  w1    - wavearray to be synchronized wrt w2
//         w2    - reference wavearray
//
// Output: sync_time  - synchronization time
//
// return xcor after synchronization

  int size = w1.size();

  wavearray<double> X = GetXCorr(w1,w2);

  // compute max xcor -> sync_xcor
  double sync_xcor=-1.e-20;
  double isync_xcor=0;
  for(int i=0;i<size;i++) if(X[i]>sync_xcor) {isync_xcor=i;sync_xcor=X[i];}
  sync_xcor*=1./(w1.rms()*sqrt(size));
  sync_xcor*=1./(w2.rms()*sqrt(size));
  sync_xcor*=size;

  // compute sync_time
  int isync_time = isync_xcor<size/2 ? -isync_xcor : size-isync_xcor;
  sync_time = isync_time/w1.rate();
  CWB::mdc::TimeShift(w1,sync_time);

  return sync_xcor;
}

void
CWB::mdc::TimeSync(wavearray<double>& w, double sync_time) {
//
// wavearray w1 is shifted in time and phase 
//
//
// Input:  w1         - wavearray 
//         sync_phase - phase shift
//
// Output: w1         - time/phase shifted wavearray

  CWB::mdc::TimeShift(w,sync_time);
}

//______________________________________________________________________________
double
CWB::mdc::TimePhaseSync(wavearray<double>& w1, wavearray<double>& w2, double& sync_time, double& sync_phase) {
//
// wavearray w1 is synchronized in time and phase wrt w2
//
//
// Input:  w1    - wavearray to be synchronized wrt w2
//         w2    - reference wavearray
//
// Output: sync_time  - synchronization time
//         sync_phase - synchronization phase
//
// return xcor after synchronization

  int size = w1.size();

  // get 90 deg phase shift wavearray
  wavearray<double> W1=w1; CWB::mdc::PhaseShift(W1,90.);
  wavearray<double> W2=w2; CWB::mdc::PhaseShift(W2,90.);

  wavearray<double> w1W2 = GetXCorr(w1,W2);
  wavearray<double> W1w2 = GetXCorr(W1,w2);
  wavearray<double> w1w2 = GetXCorr(w1,w2);
  wavearray<double> W1W2 = GetXCorr(W1,W2);

  wavearray<double> phase(size);
  for(int i=0;i<size;i++) phase[i] = TMath::ATan2(w1W2[i]-W1w2[i],w1w2[i]+W1W2[i]);

  wavearray<double> z1(size);
  for(int i=0;i<size;i++) z1[i] = w1w2[i]*cos(phase[i])-W1w2[i]*sin(phase[i]);

  double sync_xcor=-1.e-20;
  double isync_xcor=0;
  for(int i=0;i<size;i++) if(z1[i]>sync_xcor) {isync_xcor=i;sync_xcor=z1[i];}

  int isync_time = isync_xcor<size/2 ? -isync_xcor : size-isync_xcor;
  sync_time = isync_time/w1.rate();
  CWB::mdc::TimeShift(w1,sync_time);

  int isync_phase = isync_xcor;
  sync_phase = phase[isync_phase]*180./TMath::Pi();
  CWB::mdc::PhaseShift(w1,-sync_phase);

  sync_xcor=0;
  for(int i=0;i<size;i++) sync_xcor+=w1[i]*w2[i];
  sync_xcor*=1./(w1.rms()*sqrt(size));
  sync_xcor*=1./(w2.rms()*sqrt(size));

  return sync_xcor;
}

//______________________________________________________________________________
void
CWB::mdc::TimePhaseSync(wavearray<double>& w, double sync_time, double sync_phase) {
//
// wavearray w1 is shifted in time and phase 
//
//
// Input:  w1         - wavearray 
//         sync_time  - time shift
//         sync_phase - phase shift
//
// Output: w1         - time/phase shifted wavearray

  CWB::mdc::TimeShift(w,sync_time);
  CWB::mdc::PhaseShift(w,-sync_phase);
}

//______________________________________________________________________________
int
CWB::mdc::Align(wavearray<double>& w1, wavearray<double>& w2) {
//
// align w2 wrt w1
//
//
// Input:  w1     - wavearray
//         w2     - wavearray
//
// Output: w1,w2 - aligned wavarrays
//
// Return no errors -> 0, errors -> -1 and output empty w1,w2
//
// Input:
// 
//                   |------- w1 ------|
// |---------- w2 ---------|
//
// Output:
//
// |000000000000000000------- w1 ------|
// |---------- w2 ----------00000000000|

   wavearray<double> w;

   if(w1.size()==0)          {cout << "CWB::mdc::Align - Error: w1 size is 0" << endl;w1=w;w2=w;return -1;}
   if(w2.size()==0)          {cout << "CWB::mdc::Align - Error: w2 size is 0" << endl;w1=w;w2=w;return -1;}
   if(w1.rate()!=w2.rate())  {cout << "CWB::mdc::Align - Error: w1 rate != w2 rate" << endl;w1=w;w2=w;return -1;}

   double R=w1.rate();

   double b_w1 = w1.start();
   double e_w1 = w1.start()+w1.size()/R;
   double b_w2 = w2.start();
   double e_w2 = w2.start()+w2.size()/R;

   double b_w = b_w1<b_w2 ? b_w1 : b_w2;;
   double e_w = e_w1>e_w2 ? e_w1 : e_w2;;

   int o_w1 = b_w1<b_w2 ? 0 : int((b_w1-b_w2)*R+0.5);
   int o_w2 = b_w2<b_w1 ? 0 : int((b_w2-b_w1)*R+0.5);

   w.rate(R);
   w.resize(int(R*((e_w-b_w)+0.5)));
   w.start(b_w);

   w=0.;
   for(int i=0;i<w1.size();i++) w[i+o_w1] = w1[i];
   w1=w;

   w=0.;
   for(int i=0;i<w2.size();i++) w[i+o_w2] = w2[i];
   w2=w;

   return 0;
}

//______________________________________________________________________________
double
CWB::mdc::GetMatchFactor(TString match, vector<wavearray<double> >& w1, vector<wavearray<double> >& w2, 
                         vector<double> tstart, vector<double> tstop) {
//
// compute match factor between vector of wavearrays w1 and w2
//		  - vector size is the number of detectors
//
// Input:  match  - ff/of/re: fitting-factor/overlap-factor/residual-energy 
//         w1     - wavearray vector 
//         w2     - wavearray vector - is the reference wavearray
//         tstart - vector of start times used in the match computation (if tstart.size()=0 then tstart is the start time of the array) 
//         tstop  - vector of stop  times used in the match computation (if tstop.size()=0  then tstop  is the stop  time of the array) 
//
// Output:        - return match value


    int vsize=w1.size();
    // vector wavearray size consistency
    if(w2.size()!=vsize) {
      cout << "CWB::mdc::GetMatchFactor - Error: input vectors wavearray are inconsistents, different size" << endl;
      return std::numeric_limits<double>::max();
    }
    // vector tstart size consistency
    if((tstart.size()>0) && (tstart.size()!=vsize)) {
      cout << "CWB::mdc::GetMatchFactor - Error: input vectors tstart are inconsistents, different size" << endl;
      return std::numeric_limits<double>::max();
    }
    // vector tstop size consistency
    if((tstop.size()>0) && (tstop.size()!=vsize)) {
      cout << "CWB::mdc::GetMatchFactor - Error: input vectors tstop are inconsistents, different size" << endl;
      return std::numeric_limits<double>::max();
    }
    // check if vsize>0
    if(vsize==0) {
      cout << "CWB::mdc::GetMatchFactor - Error: input vectors are empty" << endl;
      return std::numeric_limits<double>::max();
    }
    double rate=w1[0].rate();
    // wavearray rate consistency
    for(int n=0;n<vsize;n++) {
      if((w1[n].rate()!=rate)||(w2[n].rate()!=rate)) {
        cout << "CWB::mdc::GetMatchFactor - Error: input wavearray are inconsistents, different rate" << endl;
        return std::numeric_limits<double>::max();
      }
    }
    double size=w1[0].size();
    // wavearray size consistency
    for(int n=0;n<vsize;n++) {
      if((w1[n].size()!=size)||(w2[n].size()!=size)) {
        cout << "CWB::mdc::GetMatchFactor - Error: input wavearray are inconsistents, different size" << endl;
        return std::numeric_limits<double>::max();
      }
    }
    double start=w1[0].start();
    // wavearray start time consistency
    for(int n=0;n<vsize;n++) {
      if((w1[n].start()!=start)||(w2[n].start()!=start)) {
        cout << "CWB::mdc::GetMatchFactor - Error: input wavearray are inconsistents, different start time" << endl;
        return std::numeric_limits<double>::max();
      }
    }

    double dt = 1./rate;

    if(tstart.size()==0) {tstart.resize(vsize); for(int n=0;n<vsize;n++) tstart[n]=start;}
    if(tstop.size()==0)  {tstop.resize(vsize);  for(int n=0;n<vsize;n++) tstop[n] =start+size*dt;}

    double  ew1=0;                             // w1 energy
    double  ew2=0;                             // w2 energy
    double  ew12=0;                            // w1,w2 coherent energy
    double  re=0;                              // residual energy
    for(int n=0;n<vsize;n++) {
      for(int j=0;j<size;j++) {
         double t=j*dt+start;
         if((t<=tstart[n])||(t>tstop[n])) continue;
         ew1  += w1[n][j]*w1[n][j];
         ew2  += w2[n][j]*w2[n][j];
         ew12 += w1[n][j]*w2[n][j];
         re   += pow(w1[n][j]-w2[n][j],2);
      }
    }

    if((match=="ff")&&(ew1*ew2==0)) return std::numeric_limits<double>::max();
    if((match=="of")&&(ew2==0))     return std::numeric_limits<double>::max();

    if(match=="ff") return ew12/sqrt(ew1*ew2);  // fitting factor
    if(match=="of") return ew12/sqrt(ew2*ew2);  // overlap factor
    if(match=="re") return re;                  // residual energy 
}

//______________________________________________________________________________
wavearray<double>
CWB::mdc::GetEnvelope(wavearray<double>* x) {
//
// compute the envelope of wavearray x
//
// Input:
//         x      - input wavearray
//
// Output:
//         xq     - return the envelope of wavearray x

   wavearray<double> xq;	// quadrature

   if(x==NULL)      return xq;
   if(x->size()==0) return xq;

   // get quadrature
   xq = *x;
   CWB::mdc::PhaseShift(xq,90);

   // get envelope
   for(int i=0;i<x->size();i++) xq[i] = sqrt(pow(x->data[i],2)+pow(xq[i],2));

   return xq;
}

//______________________________________________________________________________
wavearray<double>
CWB::mdc::GetSpectrum(wavearray<double>* x, bool oneside) {
//
// compute the spectrum of wavearray x
//
// Input:
//         x      - input wavearray
//         oneside- trie/false(default) -> one-side/double-side
//
// Output:
//         xs     - return the spectrum of wavearray x

   wavearray<double> xs = *x;

   if(x==NULL)      return xs;
   if(x->size()==0) return xs;

   double dt   = xs.rate();
   int    size = xs.size()-xs.size()%2;  // size even 

   wavearray<double> xa(size);
   double df=(double)xs.rate()/(double)xs.size();
   xs.FFTW(1);
   for (int i=0;i<size;i+=2)  xa.data[i/2]=sqrt(pow(xs.data[i],2)+pow(xs.data[i+1],2));
   for (int i=0;i<size/2;i++) xs.data[i]=xa.data[i]*(1./df);     // double side spectra 1/Hz
   if(oneside) for(int i=0;i<size/2;i++) xs.data[i]*=sqrt(2.);   // one side spectra 1/Hz
   xs.resize(size/2.);
   xs.start(0);
   xs.rate(1./df);

   return xs;
}

//______________________________________________________________________________
wavearray<double>
CWB::mdc::GetBandpass(wavearray<double> x, double bF, double eF) {
//
// apply bandpass cut to wavearray x
//
// Input:
//         x      - input wavearray
//         bF,eF  - bandpass frequency range [bF,eF]
//
// Output:
//         y      - return bandpass x wavearray

  wavearray<double> y=x;

  // cut frequency range bF,eF
  double F=0.;
  double dF=(y.rate()/(double)y.size())/2.;
  y.FFTW(1);
  for(int j=0;j<y.size();j+=2) {
    F = j*dF;
    if(F<bF || F>eF) {y[j]=0;y[j+1]=0;}
  }
  y.FFTW(-1);

  return y;
}

//______________________________________________________________________________
double
CWB::mdc::GetFrequencyBoundaries(wavearray<double> x, double P, double& bF, double& eF) {
//
// compute the symmetric frequency interval including the energy fraction P
//
//
// Input:  x     - wavearray
//         P     - energy fraction  (0:1)
//
// Output: bF    - interval begin frequency
//         eF    - interval end frequency

  x.FFTW(1);

  if(P<0) P=0;
  if(P>1) P=1;

  int N = x.size();

  double E = 0;
  double ET = 0;                                                 // signal energy
  for(int i=0;i<N;i++) {ET+=x[i]*x[i];}

  double EF = ET*(1-P);

  // search jB sample index which contains EF/2 energy in the left side of the signal
  E=0;
  int jB=0;
  for(int i=0;i<N;i++) {
    E+=x[i]*x[i];
    if(E>EF/2) break;
    jB=i;
  }

  // search jE sample index which contains EF/2 energy in the right side of the signal
  E=0;
  int jE=0;
  for(int i=N-1;i>0;i--) {
    E+=x[i]*x[i];
    if(E>EF/2) break;
    jE=i;
  }

  double dF=(x.rate()/(double)x.size())/2.;

  bF = jB*dF;
  eF = jE*dF;

  return eF-bF;
}

#ifdef _USE_LAL
double
CWB::mdc::SimIMRSEOBNRv4ROMTimeOfFrequency(double freq, double m1, double m2, double chi1, double chi2) {
//
// wrapper to XLALSimIMRSEOBNRv4ROMTimeOfFrequency LAL function
// Compute the 'time' elapsed in the ROM waveform from a given starting frequency until the ringdown (model SEOBNRv4_ROM)
// lalsimulation/src/LALSimIMRSEOBNRv4ROM.c
// include/lal/LALSimIMR.h
// int XLALSimIMRSEOBNRv4ROMTimeOfFrequency(REAL8 *t, REAL8 frequency, REAL8 m1SI, REAL8 m2SI, REAL8 chi1, REAL8 chi2);
//
//
// Input:  freq      - frequency (Hz)
//         m1,m2     - mass  components in solar mass
//         chi1,chi2 - Dimensionless aligned spinz components
//
// Output: time     - time (sec) of 2,2 mode at frequency freq
//

  double Mo = watconstants::SolarMass();

  REAL8 time; 
  if(XLALSimIMRSEOBNRv4ROMTimeOfFrequency(&time, freq, m1*Mo, m2*Mo, chi1, chi2) == XLAL_FAILURE) {
    cout << "CWB::mdc::XLALSimIMRSEOBNRv4ROMTimeOfFrequency error" << endl;
    return -1;;
  }

  return time;
}

double
CWB::mdc::SimIMRSEOBNRv4ROMFrequencyOfTime(double time, double m1, double m2, double chi1, double chi2) {
//
// wrapper to XLALSimIMRSEOBNRv4ROMFrequencyOfTime LAL function
// Compute the starting frequency so that the given amount of 'time' elapses in the ROM waveform from the starting frequency until the ringdown.
// lalsimulation/src/LALSimIMRSEOBNRv4ROM.c
// include/lal/LALSimIMR.h
// int XLALSimIMRSEOBNRv4ROMFrequencyOfTime(REAL8 *frequency, REAL8 t, REAL8 m1SI, REAL8 m2SI, REAL8 chi1, REAL8 chi2);
//
//
// Input:  time      - time (sec)
//         m1,m2     - mass  components in solar mass
//         chi1,chi2 - Dimensionless aligned spinz components
//
// Output: freq      - freq (Hz) of 2,2 mode at time 'time'
//

  double Mo = watconstants::SolarMass();

  REAL8 freq; 

  if(XLALSimIMRSEOBNRv4ROMFrequencyOfTime(&freq, time, m1*Mo, m2*Mo, chi1, chi2) == XLAL_FAILURE) {
    cout << "CWB::mdc::XLALSimIMRSEOBNRv4ROMFrequencyOfTime error" << endl;
    return -1;
  }

  return freq;
}
#endif
