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
 * Package:      graphic wat Class Library
 * File name:    gnetwork.hh
 * Author:       Gabriele Vedovato (vedovato@lnl.infn.it)
 **********************************************************/


#ifndef GNETWORK_HH
#define GNETWORK_HH

#include "gskymap.hh"
#include "network.hh"
#include "TGraphPolar.h"

class gnetwork : public network 
{

public:
 
  // Constructors
  gnetwork() : network() {Init();}
  gnetwork(int nifo, TString* ifo, detectorParams* detParms=NULL);
  gnetwork(network*  net);
  virtual ~gnetwork();

  gskymap* GetGskymap() {return &gSM;}
  void     SetGskymap(gskymap& gSM) {this->gSM = gSM;}

  void setSkyMaps(double sms,double t1,double t2,double p1,double p2) {
                  network::setSkyMaps(sms,t1,t2,p1,p2);gSM=gskymap(sms,t1,t2,p1,p2);}
  void setSkyMaps(int healpix_order) {
                  network::setSkyMaps(int(healpix_order));gSM=gskymap(int(healpix_order));}

  void DrawAntennaPattern(int polarization=-1, int dpaletteId=0, bool btitle=true, int order=6);
  void DrawSites(Color_t mcolor=kBlack, Size_t msize=2.0, Style_t mstyle=20);
  void DrawSitesArms(double mlength=600000., Color_t lcolor=kBlack, Size_t lwidth=1.0, Style_t lstyle=1);
  void DrawSitesShortLabel(Color_t tcolor=kBlack, Size_t tsize=0.052, Font_t tfont=32);
  void DrawSitesLabel(Color_t tcolor=kBlack, Size_t tsize=0.045, Font_t tfont=32);
  void DrawDelay(TString ifo1, TString ifo2, 
                 double phi=1000., double theta=1000., double frequency=0.);
  double GetDelay(TString ifo1, TString ifo2, double phi, double theta);
  double GetDelay(TString ifo1, TString ifo2, TString mode="");
  double GetAntennaPattern(double phi, double theta, double psi=0., int polarization=1);
  double GetAntennaPattern(TString ifo, double phi, double theta, double psi=0., bool plus=true);

  int Delay2Coordinates(double t1, double t2, double t3, double& ph1, double& th1, double& ph2,double& th2);

  double GetSite(TString ifo, TString type="");

  void MakeNetworkDetectorIndex(double gamma, int loop, double snr=-1.0);
  void DrawNetworkDetectorIndex(double gamma, int loop, double snr=-1.0, 
                                int dpaletteId=DUMMY_PALETTE_ID, bool btitle=true);

  void DrawCircles(double phi, double theta, double gps, 
                   Color_t lcolor=kBlack, Width_t lwidth=1, Style_t lstyle=1, bool labels=false);
  void DrawCircles(double phi, double theta, 
                   Color_t lcolor=kBlack, Width_t lwidth=1, Style_t lstyle=1, bool labels=false);

  void Draw(char* smName, network* net=NULL) {Draw(TString(smName),net);}
  void Draw(TString smName, network* net=NULL);

  void SetScale(TString ifo, double scale);
  double GetScale(TString ifo);

  void ClearSites();
  void ClearSitesArms();
  void ClearSitesLabel();
  void ClearCircles();

  //: dump gnetwork into root file
  void DumpObject(char* file);
  //: load gnetwork object from root file
  void LoadObject(char* file);

  static inline double ndi(double*, double, int);

  TCanvas* DrawPolargram(int ptype, network* net=NULL);

private:

  void Init();
  void Draw(skymap sm, TString title="") {gSM=sm;gSM.SetTitle(TString("network ")+title);gSM.Draw();}
  void Draw(wavearray<int> w, TString title="") {
       skymap sm=nSensitivity;for(int l=0;l<(int)sm.size();l++) sm.set(l,double(w[l]));
       Draw(sm,title);return;}
  void Draw(wavearray<short> w, TString title="") {
       skymap sm=nSensitivity;for(int l=0;l<(int)sm.size();l++) sm.set(l,double(w[l]));
       Draw(sm,title);return;}
  void Draw(wavearray<double> w, TString title="") {
       skymap sm=nSensitivity;for(int l=0;l<(int)sm.size();l++) sm.set(l,double(w[l]));
       Draw(sm,title);return;}

  gskymap gSM;			

  double scale[NIFO_MAX];               //! default=1, used to scale the antenna pattern acccording the PSD contribution of each detector (at a fixed frequency) to the network

  unsigned int siteLabelType;

  std::vector<TMarker*>   siteM;	//!
  std::vector<TPolyLine*> siteA;	//!
  std::vector<TText*>     siteT;	//!
  std::vector<TLatex*>    siteL;	//!
  std::vector<TLatex*>    siteP;	//!
  std::vector<TPolyLine*> circleL;	//!
  std::vector<TMarker*>   daxisM;	//!

  detector* pD[NIFO_MAX]; 		//!

  TCanvas* canvasP;			//! canvas used for the polargrams
  TGraphPolar* grP[2];			//! polargraphs used for the polargrams

  ClassDef(gnetwork,3)
};  

inline double gnetwork::ndi(double* u, double um, int nIFO) {
  double d = 0.;
  for(int i=0;i<nIFO;i++) d+= (u[i]*u[i]*u[i]*u[i]);
  d = um > 0. ? 1./(d/(um*um)) : 0.;
  return d;
}

#endif
