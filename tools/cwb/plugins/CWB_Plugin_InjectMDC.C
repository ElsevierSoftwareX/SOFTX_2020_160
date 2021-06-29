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
#include "TComplex.h"
#include "TMath.h"
#include "gnetwork.hh"

#define WAVEFORM_NAME "Waveforms/SG554Q8d9.txt"

void 
CWB_Plugin(TFile* jfile, CWB::config* cfg, network* net, WSeries<double>* x, TString ifo, int type)  {

// Plugin to inject 'on the fly' the MDC

/*
#  GravEn_SimID                        SimHrss      SimEgwR2    GravEn_Ampl  Internal_x   Internal_phi   External_x   Externa
l_phi  External_psi  FrameGPS     EarthCtrGPS     SimName               SimHpHp      SimHcHc      SimHpHc    
GEO      GEOctrG
PS       GEOfPlus       GEOfCross  H1       H1ctrGPS        H1fPlus        H1fCross  H2       H2ctrGPS        H2fPlus
H2fCross  L1       L1ctrGPS        L1fPlus        L1fCross  V1       V1ctrGPS        V1fPlus        V1fCross

Waveforms/SG554Q8d9.txt             1.213000e-21 1.011301e-47 1.213000e-21 +1.000000e+00 +0.000000e+00 +7.799320e-01 +7.720941e-01 +3.397006e+00 968653615   968654066.616913 SG554Q8d9            1.471369e-42 0.000000e+00 0.000000e+00 GEO   968654066.597115 +8.635698e-01 -3.501556e-01 H1   968654066.613762 +1.334841e-01 +7.593384e-02 H2   968654066.613762 +1.334841e-01 +7.593384e-02 L1   968654066.616641 -1.870505e-02 -1.359386e-02 V1   968654066.597488 -4.561425e-01 -7.932736e-01
*/

  if(type==CWB_PLUGIN_CONFIG) {
    cfg->mdcPlugin=true;  // disable read mdc from frames
  }

  if(type==CWB_PLUGIN_MDC) { // mdc
 
    double dt = 1./x->rate();

    cout << endl;
    cout << "-----> plugins/CWB_Plugin_InjectMDC.C" << endl;
    cout << "ifo " << ifo.Data() << endl;
    cout << "type " << type << endl;
    cout << endl;

    // log burstMDC parameters
    double burstMDC_theta = 7.799320e-01;
    double burstMDC_phi   = 7.720941e-01;
    double burstMDC_psi   = 3.397006;
  
    /* ------------------------------
      earthCenterTime = _EarthCtrGPS;
      phi = _External_phi;
      theta = _External_x;
      psi = _External_psi;
    --------------------------------- */
    double theta,phi,psi;
    double Pi = TMath::Pi();
    if(true) {
      theta = acos(burstMDC_theta);
      theta*= 180/Pi;
      phi = burstMDC_phi > 0 ? burstMDC_phi : 2*Pi+burstMDC_phi;
      phi*= 180/Pi;
//      phi = sm.phi2RA(phi,_EarthCtrGPS);
      psi = burstMDC_psi*180/Pi;
    }
    cout << "theta : " << theta << " phi : " << phi << " psi " << psi << endl;

    int nIFO=net->ifoListSize();
    TString ifos[NIFO_MAX];
    for(int n=0;n<nIFO;n++)ifos[n]= net->ifoName[n];
    gnetwork gNET(nIFO,ifos);
    for(int i=0;i<nIFO;i++) {
      for(int j=i+1;j<nIFO;j++) {
        cout << ifos[i].Data() << " " << ifos[j].Data() << " -> "
             << gNET.GetDelay(ifos[i].Data(),ifos[j].Data(),phi,theta) << " sec " << endl;
      }
    }
  

    char   logString[1024]="";

    char   GravEn_SimID[1024]="Waveforms/SG554Q8d9.txt";                        
    double SimHrss      = 1.213000e-21;      
    double SimEgwR2     = 1.011301e-47;
    double GravEn_Ampl  = 1.213000e-21;  
    double Internal_x   = 1.0;
    double Internal_phi = 0.0;
    double External_x   = burstMDC_theta;
    double External_phi = burstMDC_phi;
    double External_psi = burstMDC_psi;
  
    double FrameGPS = x->start();
    //double EarthCtrGPS = FrameGPS+dt*x->size()/2.;  // injected in the center of buffer 
    double EarthCtrGPS = 968654066.616913;  
    char   SimName[64] = "SG554Q8d9";               
    double SimHpHp = 1.471369e-42;     
    double SimHcHc = 0;     
    double SimHpHc = 0;
  
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
  
    double tShift=0;
    double fPlus=0;
    double fCross=0;
    for(int i=0;i<nIFO;i++) {
      double IFOctrGPS = EarthCtrGPS;
      if(i>0) IFOctrGPS += gNET.GetDelay(ifos[i].Data(),ifos[0].Data(),phi,theta);
      double IFOfPlus  = gNET.GetAntennaPattern(ifos[i], phi, theta, psi, true);
      double IFOfCross = gNET.GetAntennaPattern(ifos[i], phi, theta, psi, false);
      if(ifos[i].CompareTo(ifo)==0) {
        tShift=IFOctrGPS-EarthCtrGPS;
        fPlus=IFOfPlus;
        fCross=IFOfCross;
      }
      sprintf(logString,"%s %s %10.6f %e %e",logString,ifos[i].Data(),IFOctrGPS,IFOfPlus,IFOfCross);
    }
    cout << logString << endl;

    net->mdcList.clear();
    net->mdcList.push_back(logString);
    net->mdcTime.clear();
    net->mdcTime.push_back(EarthCtrGPS);
  
    cout << endl << endl;
    for (int nmdc=0; nmdc<(int)net->mdcListSize(); nmdc++) {
      TString mdcstring(net->getmdcList(nmdc));
//      if(nmdc==72) {
        cout << endl << endl;
        cout << "--------------------------------> " << nmdc << endl;
        cout << endl << endl;
        cout << mdcstring.Data() << endl;
//      }
    }
    //gSystem->Exit(0);

    // read waveform
    ifstream in;
    in.open(WAVEFORM_NAME, ios::in);
    if (!in.good()) {cout << "Error Opening File : " << WAVEFORM_NAME << endl;gSystem->Exit(1);}

    // fill mdc vector 
    (*x)=0;
    int offset = (EarthCtrGPS-FrameGPS)*x->rate();
    int size=0;
    double hrss=0;
    while (1) {
      in >> x->data[offset+size];
      hrss+=pow(x->data[offset+size],2);
      if (!in.good()) break;
      size++;
      if ((offset+size)>=(int)x->size()) break;
    }
    hrss=sqrt(hrss*dt);
    for(int i=0;i<(int)x->size();i++) x->data[i] = x->data[i]*fPlus*SimHrss/hrss;

    in.close();
 
    if(tShift==0) return;
 
    // apply time shift to mdc vector 
    x->FFTW(1);
    TComplex C;
    double df = x->rate()/x->size(); 
    //double tShift = 10./x->rate();
    cout << "tShift : " << tShift << endl;
    for (int ii=0;ii<(int)x->size()/2;ii++) {
      TComplex X(x->data[2*ii],x->data[2*ii+1]);
      X=X*C.Exp(TComplex(0.,-2*Pi*ii*df*tShift));  // Time Shift
      x->data[2*ii]=X.Re();
      x->data[2*ii+1]=X.Im();
    }
    x->FFTW(-1);
  }

  return;
}
