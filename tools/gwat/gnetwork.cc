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



#include "gnetwork.hh"
#include "constants.hh"

using namespace std;
using namespace ROOT::Math;

ClassImp(gnetwork)


//______________________________________________________________________________
/* Begin_Html
<center><h2>gnetwork class</h2></center>
This class gnetwork is derived from the network class and uses the methods of the class gskymap
It is used to produce the skymap plots of the detector's network antenna pattern or the statistics produced in the cwb likelihood stage.

<p>
<h3><a name="example">Example</a></h3>
<p>
The macro <a href="./tutorials/gwat/DrawAntennaPattern.C.html">DrawAntennaPattern.C</a> is an example which shown how to use the gnetwork class.<br>
The pictures below gives the macro output plot.<br>
<p>

End_Html

Begin_Macro
DrawAntennaPattern.C
End_Macro */


//______________________________________________________________________________
gnetwork::gnetwork(int nifo, TString* ifo, detectorParams* detParms) : network() {
//
// constructor
//
// Input: nifo      - number of detectors
//        ifo       - array of built-in ifo names  
//                    Ex : TString ifo[3] = {"L1","H1","V1"} 
//        detParms  - array of usewr defined detector's parameters 
//                    it is used when ifo entry is empty -> ""
//

  for(int n=0; n<nifo; n++) {
    if((ifo!=NULL)&&(ifo[n]!="")) {
       pD[n] = new detector(const_cast<char*>(ifo[n].Data())); // built in detector
    } else {
       if(detParms!=NULL) {
         pD[n] = new detector(detParms[n]);                    // user define detector
       } else {
         cout << "gnetwork::gnetwork : Error - user define detector with no detParams" << endl;
         exit(1);
       }
    }
  }

  for(int n=0; n<nifo; n++) this->add(pD[n]);
  Init();
}

//______________________________________________________________________________
gnetwork::gnetwork(network* net) : network() {
//
// constructor
//
// Input: net       - pointer to network
//   

  int nifo = net->ifoListSize();
  for(int n=0; n<nifo; n++) pD[n] = new detector(*net->getifo(n));

  for(int n=0; n<nifo; n++) this->add(pD[n]);
  Init();
}

//______________________________________________________________________________
gnetwork::~gnetwork() {
//
// destructor
//

  ClearSites();
  ClearSitesArms();
  ClearSitesLabel();
  ClearCircles();
  for(int n=0;n<NIFO_MAX;n++) if(pD[n]!=NULL) delete pD[n];
  for(int i=0;i<2;i++) if(grP[i]!=NULL) delete grP[i];
  if(canvasP!=NULL) delete canvasP;
  for(size_t n=0; n<this->ifoList.size(); n++) if(this->ifoList[n]!=NULL) delete this->ifoList[n];
}

//______________________________________________________________________________
void
gnetwork::Init() { 
//
// initialize internal parameters
//

  siteLabelType=0;
  for(int n=0;n<NIFO_MAX;n++) pD[n]=NULL;
  for(int n=0;n<NIFO_MAX;n++) scale[n]=1;
  for(int i=0;i<2;i++) grP[i]=NULL;
  canvasP=NULL;
}

//______________________________________________________________________________
void 
gnetwork::Draw(TString smName, network* net) {
//
// Draw skymap network info
//
// Input: smName  - skymap name 
//
//                  nSensitivity   : network sensitivity
//                  nAlignment     : network alignment factor
//                  nCorrelation   : network correlation coefficient
//                  nLikelihood    : network likelihood
//                  nNullEnergy    : network null energy
//                  nPenalty       : signal * noise penalty factor
//                  nCorrEnergy    : reduced correlated energy
//                  nNetIndex      : network index
//                  nDisbalance    : energy disbalance
//                  nSkyStat       : sky optimization statistic
//                  nEllipticity   : waveform ellipticity
//                  nPolarisation  : polarisation angle
//                  nProbability   : probability skymap
//                     
//                  index          : theta, phi mask index array
//                  skyMask        : index array for setting sky mask
//                  skyMaskCC      : index array for setting sky mask Celestial Coordinates
//                  skyHole        : static sky mask describing "holes"
//                  veto           : veto array for pixel selection
//                  skyProb        : sky probability
//                  skyENRG        : energy skymap
//
//        net   - pointer to network object
//                if net!=NULL then net is used to retrive the data to be plotted
//

  network* NET = net!=NULL ? net : this;

  if(smName=="nSensitivity")  {Draw(NET->nSensitivity,smName);return;}  // network sensitivity
  if(smName=="nAlignment")    {Draw(NET->nAlignment,smName);return;}    // network alignment factor
  if(smName=="nCorrelation")  {Draw(NET->nCorrelation,smName);return;}  // network correlation coefficient
  if(smName=="nLikelihood")   {Draw(NET->nLikelihood,smName);return;}   // network likelihood
  if(smName=="nNullEnergy")   {Draw(NET->nNullEnergy,smName);return;}   // network null energy
  if(smName=="nPenalty")      {Draw(NET->nPenalty,smName);return;}      // signal * noise penalty factor
  if(smName=="nCorrEnergy")   {Draw(NET->nCorrEnergy,smName);return;}   // reduced correlated energy
  if(smName=="nNetIndex")     {Draw(NET->nNetIndex,smName);return;}     // network index
  if(smName=="nDisbalance")   {Draw(NET->nDisbalance,smName);return;}   // energy disbalance
  if(smName=="nSkyStat")      {Draw(NET->nSkyStat,smName);return;}      // sky optimization statistic
  if(smName=="nEllipticity")  {Draw(NET->nEllipticity,smName);return;}  // waveform ellipticity
  if(smName=="nPolarisation") {Draw(NET->nPolarisation,smName);return;} // polarisation angle
  if(smName=="nProbability")  {Draw(NET->nProbability,smName);return;}  // probability skymap

  if(smName=="index")         {Draw(NET->index,smName);return;}         // theta, phi mask index array
  if(smName=="skyMask")       {Draw(NET->skyMask,smName);return;}       // index array for setting sky mask
  if(smName=="skyMaskCC")     {Draw(NET->skyMaskCC,smName);return;}     // index array for setting sky mask Celestial Coordinates
  if(smName=="skyHole")       {Draw(NET->skyHole,smName);return;}       // static sky mask describing "holes"
  if(smName=="veto")          {Draw(NET->veto,smName);return;}          // veto array for pixel selection
  if(smName=="skyProb")       {Draw(NET->skyProb,smName);return;}       // sky probability
  if(smName=="skyENRG")       {Draw(NET->skyENRG,smName);return;}       // energy skymap

  cout << "gnetwork::Draw : Warning - skymap '" << smName.Data() << "' not valid !" << endl;
}

//______________________________________________________________________________
double
gnetwork::GetSite(TString ifo, TString type) {
//
// return detector infos
//
// Input: ifo     - detector name
//        type    - detector info type
//
//                  THETA 	: theta full 
//                  THETA_D	: theta degrees
//                  THETA_M	: theta minutes
//                  THETA_S	: theta seconds
//
//                  PHI         : phi full    
//                  PHI_D 	: phi degrees    
//                  PHI_M 	: phi minutes    
//                  PHI_S 	: phi seconds    
//
//                  ""          : print out theta,phi infos
//
// if type invalid return 1111
//

  int nIFO = this->ifoListSize();
  if(nIFO<1) {cout << "gnetwork::GetSite nIFO must be >= 1" << endl;return -1111;}

  detector *pD[NIFO];
  for(int n=0; n<nIFO; n++) pD[n] = this->getifo(n);

  int ifoId=0;
  int nifo=0; 
  for(int n=0; n<nIFO; n++) {
    if(ifo.CompareTo(this->getifo(n)->Name)==0) {nifo++;ifoId=n;}
  }                                                        
  if(nifo<1) {cout << "gnetwork::GetSite ifo not valid" << endl;return -1111;}

  detectorParams dP = pD[ifoId]->getDetectorParams();

  double phi=dP.longitude;
  double theta=dP.latitude;

  char LAT;
  double theta_t=theta;
  if(theta_t>0) LAT='N'; else {LAT='S';theta_t=-theta_t;}
  int theta_d = int(theta_t);                            
  int theta_m = int((theta_t-theta_d)*60);               
  float theta_s = (theta_t-theta_d-theta_m/60.)*3600.;   

  char LON;
  double phi_t=phi;
  if(phi_t>0) LON='E'; else {LON='W';phi_t=-phi_t;}
  int phi_d = int(phi_t);                          
  int phi_m = int((phi_t-phi_d)*60);               
  float phi_s = (phi_t-phi_d-phi_m/60.)*3600.;     

  type.ToUpper();
  if(type.CompareTo("THETA")==0) return theta;
  if(type.CompareTo("THETA_D")==0) return theta_d;
  if(type.CompareTo("THETA_M")==0) return theta_m;
  if(type.CompareTo("THETA_S")==0) return theta_s;
  if(type.CompareTo("PHI")==0) return phi;        
  if(type.CompareTo("PHI_D")==0) return phi_d;    
  if(type.CompareTo("PHI_M")==0) return phi_m;    
  if(type.CompareTo("PHI_S")==0) return phi_s;    

  if(type.CompareTo("")==0) pD[ifoId]->print();

  return  1111;
}              

//______________________________________________________________________________
void 
gnetwork::DrawSites(Color_t mcolor, Size_t msize, Style_t mstyle) {
//
// Draw a marker in the ifo site position (see ROOT TMarker input parameters)
//
// Input: mcolor  - color of marker
//        msize   - size of marker
//        mstyle  - style of marker
// 

  if(gSM.goff) return;
  if(gSM.canvas==NULL) return;

  int nIFO = this->ifoListSize();
  if(nIFO<1) return;            

  ClearSites();

  detector *pD[NIFO];
  for(int n=0; n<nIFO; n++) pD[n] = this->getifo(n);

  gSM.canvas->cd();

  double r2d = 180./TMath::Pi();
  for(int n=0;n<nIFO;n++) {     
    XYZVector Rv(pD[n]->Rv[0],pD[n]->Rv[1],pD[n]->Rv[2]);
    TVector3 vRv(Rv.X(),Rv.Y(),Rv.Z());                  

    double phi=vRv.Phi()*r2d;
    double theta=vRv.Theta()*r2d;

    theta=90-theta;

    if (gSM.projection.CompareTo("hammer")==0) {
      gSM.ProjectHammer(phi, theta, phi, theta);
    } else if (gSM.projection.CompareTo("sinusoidal")==0) {
      gSM.ProjectSinusoidal(phi, theta, phi, theta);       
    } else if (gSM.projection.CompareTo("parabolic")==0) { 
      gSM.ProjectParabolic(phi, theta, phi, theta);        
    }                                                  
    if (gSM.coordinate.CompareTo("CWB")==0) GeographicToCwb(phi,theta,phi,theta);
    TMarker* pM = new TMarker(phi,theta,mstyle);                             
    pM->SetMarkerSize(msize);                                                
    pM->SetMarkerColor(mcolor);                                              
    pM->Draw();                                                              
    siteM.push_back(pM);                                                     
    pM = new TMarker(phi,theta,mstyle);                                      
    pM->SetMarkerSize(msize/2);                                              
    pM->SetMarkerColor(kWhite);                                              
    pM->Draw();                                                              
    siteM.push_back(pM);                                                     
  }                                                                          

  return;
}        

//______________________________________________________________________________
void 
gnetwork::DrawSitesArms(double mlength, Color_t lcolor, Size_t lwidth, Style_t lstyle) {
//
// Draw detector's arms in the ifo site position (see ROOT TPolyLine input parameters)
//
// Input: mlength - arm's length in KM
//        mcolor  - color of line
//        msize   - size of line
//        mstyle  - style of line
// 

  
  if(gSM.goff) return;
  if(gSM.canvas==NULL) return;

  int nIFO = this->ifoListSize();
  if(nIFO<1) return;            

  ClearSitesArms();

  detector *pD[NIFO];
  for(int n=0; n<nIFO; n++) pD[n] = this->getifo(n);

  gSM.canvas->cd();

  double mlength_deg = 360.*mlength/(TMath::TwoPi()*watconstants::EarthEquatorialRadius());

  double x[nIFO][2][2];
  double y[nIFO][2][2];
  double phi,theta;    
  double r2d = 180./TMath::Pi();
  for(int n=0;n<nIFO;n++) {     
    XYZVector Rv(pD[n]->Rv[0],pD[n]->Rv[1],pD[n]->Rv[2]);
    XYZVector Ex(pD[n]->Ex[0],pD[n]->Ex[1],pD[n]->Ex[2]);
    XYZVector Ey(pD[n]->Ey[0],pD[n]->Ey[1],pD[n]->Ey[2]);
    TVector3 vRv(Rv.X(),Rv.Y(),Rv.Z());                  
    TVector3 vEx(Ex.X(),Ex.Y(),Ex.Z());                  
    TVector3 vEy(Ey.X(),Ey.Y(),Ey.Z());                  
    vEx=mlength*vEx+vRv;                                 
    vEy=mlength*vEy+vRv;                                 

    phi=vRv.Phi()*r2d;
    theta=vRv.Theta()*r2d;
    theta=90-theta;       
    if (gSM.projection.CompareTo("hammer")==0) {
      gSM.ProjectHammer(phi, theta, phi, theta);
    } else if (gSM.projection.CompareTo("sinusoidal")==0) {
      gSM.ProjectSinusoidal(phi, theta, phi, theta);       
    } else if (gSM.projection.CompareTo("parabolic")==0) { 
      gSM.ProjectParabolic(phi, theta, phi, theta);        
    }                                                  
    if (gSM.coordinate.CompareTo("CWB")==0) GeographicToCwb(phi,theta,phi,theta);
    x[n][0][0]=phi;y[n][0][0]=theta;                                         
    x[n][1][0]=phi;y[n][1][0]=theta;                                         

    phi=vEx.Phi()*r2d;
    theta=vEx.Theta()*r2d;
    theta=90-theta;       
    if (gSM.projection.CompareTo("hammer")==0) {
      gSM.ProjectHammer(phi, theta, phi, theta);
    } else if (gSM.projection.CompareTo("sinusoidal")==0) {
      gSM.ProjectSinusoidal(phi, theta, phi, theta);       
    } else if (gSM.projection.CompareTo("parabolic")==0) { 
      gSM.ProjectParabolic(phi, theta, phi, theta);        
    }                                                  
    if (gSM.coordinate.CompareTo("CWB")==0) GeographicToCwb(phi,theta,phi,theta);
    x[n][0][1]=phi;y[n][0][1]=theta;                                         
    // fix draw on the edge 
    if (gSM.coordinate.CompareTo("CWB")==0) {
       if((x[n][0][0]>360-mlength_deg)&&(x[n][0][0]+x[n][0][1]-360<mlength_deg)) x[n][0][1]=360;
    }
    if (gSM.coordinate.CompareTo("GEOGRAPHIC")==0) {
       if((x[n][0][0]>180-mlength_deg)&&(x[n][0][0]+x[n][0][1]<mlength_deg)) x[n][0][1]=180;
    }

    phi=vEy.Phi()*r2d;
    theta=vEy.Theta()*r2d;
    theta=90-theta;       
    if (gSM.projection.CompareTo("hammer")==0) {
      gSM.ProjectHammer(phi, theta, phi, theta);
    } else if (gSM.projection.CompareTo("sinusoidal")==0) {
      gSM.ProjectSinusoidal(phi, theta, phi, theta);       
    } else if (gSM.projection.CompareTo("parabolic")==0) { 
      gSM.ProjectParabolic(phi, theta, phi, theta);        
    }                                                  
    if (gSM.coordinate.CompareTo("CWB")==0) GeographicToCwb(phi,theta,phi,theta);
    x[n][1][1]=phi;y[n][1][1]=theta;                                         
    // fix draw on the edge 
    if (gSM.coordinate.CompareTo("CWB")==0) {
       if((x[n][1][0]>360-mlength_deg)&&(x[n][1][0]+x[n][1][1]-360<mlength_deg)) x[n][1][1]=360;
    }
    if (gSM.coordinate.CompareTo("GEOGRAPHIC")==0) {
       if((x[n][1][0]>180-mlength_deg)&&(x[n][1][0]+x[n][1][1]<mlength_deg)) x[n][1][1]=180;
    }
  }                                                                          

  // Draw X,Y arms
  for(int n=0;n<nIFO;n++) {
    for(int m=0;m<2;m++) { 
      TPolyLine *pL = new TPolyLine(2,x[n][m],y[n][m]);
      pL->SetLineColor(lcolor);                        
      pL->SetLineStyle(lstyle);                        
      pL->SetLineWidth(lwidth);                        
      pL->Draw();                                      
      siteA.push_back(pL);                             
    }                                                  
  }                                                    

  return;
}        

//______________________________________________________________________________
void 
gnetwork::DrawSitesShortLabel(Color_t tcolor, Size_t tsize, Font_t tfont) {
//
// Draw detector's short label in the ifo site position (see ROOT TLatex input parameters)
// Ex : L,H,V
//
// Input: tcolor  - color of text
//        tsize   - size of text
//        tfont   - font of text
//

  siteLabelType=1;                                                        
  DrawSitesLabel(tcolor, tsize, tfont);                                   
}                                                                         

//______________________________________________________________________________
void 
gnetwork::DrawSitesLabel(Color_t tcolor, Size_t tsize, Font_t tfont) {
//
// Draw detector's short label in the ifo site position (see ROOT TLatex input parameters)
// Ex : L1,H1,V1
//
// Input: tcolor  - color of text
//        tsize   - size of text
//        tfont   - font of text
//

  if(gSM.goff) return;
  if(gSM.canvas==NULL) return;

  int nIFO = this->ifoListSize();
  if(nIFO<1) return;            

  ClearSitesLabel();

  detector *pD[NIFO];
  for(int n=0; n<nIFO; n++) pD[n] = this->getifo(n);

  gSM.canvas->cd();

  double r2d = 180./TMath::Pi();
  for(int n=0;n<nIFO;n++) {     
    XYZVector Rv(pD[n]->Rv[0],pD[n]->Rv[1],pD[n]->Rv[2]);
    TVector3 vRv(Rv.X(),Rv.Y(),Rv.Z());                  

    double phi=vRv.Phi()*r2d;
    double theta=vRv.Theta()*r2d;

    theta=90-theta;

    if (gSM.projection.CompareTo("hammer")==0) {
      gSM.ProjectHammer(phi, theta, phi, theta);
    } else if (gSM.projection.CompareTo("sinusoidal")==0) {
      gSM.ProjectSinusoidal(phi, theta, phi, theta);       
    } else if (gSM.projection.CompareTo("parabolic")==0) { 
      gSM.ProjectParabolic(phi, theta, phi, theta);        
    }                                                  
    if (gSM.coordinate.CompareTo("CWB")==0) GeographicToCwb(phi,theta,phi,theta);
    TString ifoName=this->getifo(n)->Name;                                         
         if(ifoName.CompareTo("H2")==0) {phi+=2;theta-=10;}                  
    else if(ifoName.CompareTo("J1")==0) {phi+=6;theta-=1;}                   
    else if(ifoName.CompareTo("H1")==0) {phi+=6;theta-=1;}                   
    else if(ifoName.CompareTo("L1")==0) {phi+=3;theta+=5;}                   
    else if(ifoName.CompareTo("V1")==0) {phi+=4;theta+=2;}                   
    else {phi+=2;theta+=2;}                                                  
    if(siteLabelType==1) {                                                   
           if(ifoName.CompareTo("A1")==0) ifoName=TString("#tilde{A}");      
      else if(ifoName.CompareTo("H2")==0) ifoName=TString("#tilde{H}");      
      else if(ifoName.CompareTo("I2")==0) ifoName=TString("#tilde{I}");      
      else if(ifoName.CompareTo("R2")==0) ifoName=TString("#tilde{R}");      
      else ifoName.Resize(1);                                                
    }                                                                        
    TLatex *pL = new TLatex(phi+1.0,theta, ifoName);                         
    pL->SetTextFont(tfont);                                                  
    pL->SetTextSize(tsize);                                                  
    pL->SetTextColor(tcolor);                                                
    pL->Draw();                                                              
    siteL.push_back(pL);                                                     
  }                                                                          

  return;
}        

//______________________________________________________________________________
void
gnetwork::ClearCircles() {
//
// Clear all circles 
//

  for(int i=0; i<(int)circleL.size(); i++) {if(circleL[i]) delete circleL[i];}
  circleL.clear();
  for(int i=0; i<(int)daxisM.size(); i++) {if(daxisM[i]) delete daxisM[i];}
  daxisM.clear();
  return;
}

//______________________________________________________________________________
void 
gnetwork::ClearSites() {
//
// Clear all site markers 
//

  for(int i=0; i<(int)siteM.size(); i++) {if(siteM[i]) delete siteM[i];}
  siteM.clear();                                                        
  return;                                                               
}                                                                       

//______________________________________________________________________________
void 
gnetwork::ClearSitesArms() {
//
// Clear all site arms 
//

  for(int i=0; i<(int)siteA.size(); i++) {if(siteA[i]) delete siteA[i];}
  siteA.clear();                                                        
  return;                                                               
}                                                                       

//______________________________________________________________________________
void 
gnetwork::ClearSitesLabel() {
//
// Clear all site labels 
//

  for(int i=0; i<(int)siteT.size(); i++) {if(siteT[i]) delete siteT[i];}
  siteT.clear();                                                        
  for(int i=0; i<(int)siteL.size(); i++) {if(siteL[i]) delete siteL[i];}
  siteL.clear();                                                        
  for(int i=0; i<(int)siteP.size(); i++) {if(siteP[i]) delete siteP[i];}
  siteP.clear();                                                        
  return;                                                               
}                                                                       

//______________________________________________________________________________
double 
gnetwork::GetAntennaPattern(double phi, double theta, double psi, int polarization) {
//
// Get Network Antenna Pattern value in the DPF : Dominant Polarization Frame
//
// Input:  phi		- phi (degrees)
//         theta        - theta (degrees)
//         psi          - polarization angle (degrees)
//         polarization - DPF component
//
// 			  0 -> |Fx|                     	- DPF
// 			  1 -> |F+|                     	- DPF
// 			  2 -> |Fx|/|F+|                	- DPF
// 			  3 -> sqrt((|F+|^2+|Fx|^2)/nIFO)       - DPF
// 			  4 -> |Fx|^2                   	- DPF        
// 			  5 -> |F+|^2                   	- DPF        
// 			  6 -> Fx                       	- only with 1 detector
// 			  7 -> F+                       	- only with 1 detector
// 			  8 -> F1x/F2x                  	- only with 2 detectors
// 			  9 -> F1+/F2+                  	- only with 2 detectors
// 			  10 -> sqrt(|F1+|^2+|F1x|^2)/sqrt(|F2+|^2+|F2x|^2)  - only with 2 detectors
// 			  12 -> DPF Angle
//

  if((polarization<0||polarization>10)&&(polarization!=12)) 
    {cout << "polarization must be [0-10][12]" << endl;exit(1);}

  int nIFO = this->ifoListSize();
  if(nIFO==0)                   
    {cout << "detectors are no declared" << endl;exit(1);}
  if((nIFO!=1)&&(polarization==6||polarization==7))       
    {cout << "with polarization [6-7] works only with one detector" << endl;exit(1);}
  if((nIFO!=2)&&(polarization>=8&&polarization<=11))                                 
    {cout << "with polarization [8-11] works only two one detectors" << endl;exit(1);}

  detector *pD[NIFO];
  for(int n=0; n<nIFO; n++) pD[n] = this->getifo(n);

  double z=0.;
  std::vector<wavecomplex> F;
  F.resize(nIFO);

  for (int n=0;n<nIFO;n++) F[n] = pD[n]->antenna(theta,phi,psi);
  double gp=0,gx=0,gI=0;                                        
  for (int n=0;n<nIFO;n++) {                                    
    gp+=scale[n]*(F[n].real()*F[n].real());
    gx+=scale[n]*(F[n].imag()*F[n].imag());
    gI+=scale[n]*(F[n].real()*F[n].imag());
  }                                                             
  double gR = (gp-gx)/2.;                                       
  double gr = (gp+gx)/2.;                                       
  double gc = sqrt(gR*gR+gI*gI);                                

  if (polarization==0) z = sqrt(fabs(gr-gc)<1e-8?0.:gr-gc);
  if (polarization==1) z = sqrt(gr+gc);                    
  if (polarization==2&&sqrt(gr+gc)>0) z=sqrt((gr-gc)/(gr+gc));
  if (polarization==3) z = sqrt(2*gr/nIFO);                   
  if (polarization==4) z = (fabs(gr-gc)<1e-8?0.:gr-gc);       
  if (polarization==5) z = (gr+gc);                           
  if (polarization==6) z = F[0].imag();                       
  if (polarization==7) z = F[0].real();                       
  if (polarization==8) z = (F[1].imag()!=0) ? F[0].imag()/F[1].imag() : 0.; 
  if (polarization==9) z = (F[1].real()!=0) ? F[0].real()/F[1].real() : 0.; 
  if (polarization==10) {                                                   
    double F0=sqrt(pow(F[0].real(),2)+pow(F[0].imag(),2));                  
    double F1=sqrt(pow(F[1].real(),2)+pow(F[1].imag(),2));                  
    z = (F1!=0) ? F0/F1 : 0.;                                               
  }                                                                         
  if (polarization==12) {		// DPF Angle                                                   
    double cos2G = gc!=0 ? gR/gc : 0;	// (Fp2-Fc2)/Norm
    double sin2G = gc!=0 ? gI/gc : 0;	// 2*Fpc/Norm
    double r2d = 180./TMath::Pi();
    z = r2d*atan2(sin2G,cos2G)/2.; 
  }

  return z;
}          

//______________________________________________________________________________
double 
gnetwork::GetAntennaPattern(TString ifo, double phi, double theta, double psi, bool plus) {
//
// Get Detector Antenna Pattern value 
//
// Input:  ifo          - detector name (Ex: "L1","H1","V1")
//         phi	        - phi (degrees)
//         theta        - theta (degrees)
//         psi          - polarization angle (degrees)
//         plus         - true/false -> f+/fx
//

  int nIFO = this->ifoListSize();
  if(nIFO<1) {cout << "gnetwork::GetAntennaPattern nIFO must be >= 1" << endl;return -100;}

  detector *pD[NIFO];
  for(int n=0; n<nIFO; n++) pD[n] = this->getifo(n);

  int ifoId=0;
  int nifo=0; 
  for(int n=0; n<nIFO; n++) {
    if(ifo.CompareTo(this->getifo(n)->Name)==0) {nifo++;ifoId=n;}
  }                                                        
  if(nifo<1) {cout << "gnetwork::GetAntennaPattern ifo not valid" << endl;return -100;}

  wavecomplex F = pD[ifoId]->antenna(theta,phi,psi);
  if(plus) return F.real(); else return F.imag();   

}

//______________________________________________________________________________
void 
gnetwork::DrawAntennaPattern(int polarization, int dpaletteId, bool btitle, int order) {
//
// Draw Network Antenna Pattern sky map in the DPF : Dominant Polarization Frame
//
// Input:  polarization - DPF component
//
// 			  0 -> |Fx|                     	- DPF
// 			  1 -> |F+|                     	- DPF
// 			  2 -> |Fx|/|F+|                	- DPF
// 			  3 -> sqrt((|F+|^2+|Fx|^2)/nIFO)       - DPF
// 			  4 -> |Fx|^2                   	- DPF        
// 			  5 -> |F+|^2                   	- DPF        
// 			  6 -> Fx                       	- only with 1 detector
// 			  7 -> F+                       	- only with 1 detector
// 			  8 -> F1x/F2x                  	- only with 2 detectors
// 			  9 -> F1+/F2+                  	- only with 2 detectors
// 			  10 -> sqrt(|F1+|^2+|F1x|^2)/sqrt(|F2+|^2+|F2x|^2) - only with 2 detectors
// 			  11 -> The same as (10) but averaged over psi      - only with 2 detectors
// 			  12 -> DPF Angle
//
//         dpaletteId   - palette ID
//         btitle       - if true then a default title is provided
//         order        - healpix order used for Antenna Pattern
//

  if(gSM.goff) return;

  gSM=skymap(int(order));

  if(polarization<0||polarization>12) {
    cout << "---------------------------------------------------------------------------------------" << endl;
    cout << "polarization must be [0-12]" << endl;
    cout << "---------------------------------------------------------------------------------------" << endl;
    cout << "polarization=0  -> |Fx|                     	DPF" << endl;
    cout << "polarization=1  -> |F+|                     	DPF" << endl;
    cout << "polarization=2  -> |Fx|/|F+|                	DPF" << endl;
    cout << "polarization=3  -> sqrt((|F+|^2+|Fx|^2)/nIFO)      DPF" << endl;
    cout << "polarization=4  -> |Fx|^2                   	DPF" << endl;
    cout << "polarization=5  -> |F+|^2                   	DPF" << endl;
    cout << "polarization=6  -> Fx                       	only with 1 detector" << endl;
    cout << "polarization=7  -> F+                      	only with 1 detector" << endl;
    cout << "polarization=8  -> F1x/F2x                  	only with 2 detectors" << endl;
    cout << "polarization=9  -> F1+/F2+                  	only with 2 detectors" << endl;
    cout << "polarization=10 -> sqrt(|F1+|^2+|F1x|^2)/sqrt(|F2+|^2+|F2x|^2)    only with 2 detectors" << endl;
    cout << "polarization=11 -> The same as (10) but averaged over psi         only with 2 detectors" << endl;
    cout << "polarization=12 -> DPF Angle" << endl;
    cout << "---------------------------------------------------------------------------------------" << endl;
    return;
  }

  int nIFO = this->ifoListSize();
  if(nIFO==0)                   
    {cout << "detectors are no declared" << endl;return;}
  if((nIFO!=1)&&(polarization==6||polarization==7))      
    {cout << "with polarization [6-7] works only with one detector" << endl;return;}
  if((nIFO!=2)&&(polarization>=8&&polarization<=11))                                
    {cout << "with polarization [8-11] works only two one detectors" << endl;return;}

  detector *pD[NIFO];
  for(int n=0; n<nIFO; n++) pD[n] = this->getifo(n);

  TString pTitle="";
  switch(polarization) {
  case 0:               
    pTitle=TString("|F_{x}|");
    break;                    
  case 1:                     
    pTitle=TString("|F_{+}|");
    break;                    
  case 2:                     
    pTitle=TString("|F_{x}|/|F_{+}|");
    break;                            
  case 3:                             
    pTitle=TString("#sqrt{(|F_{+}|^{2}+|F_{x}|^{2})/nIFO}");
    break;                                                  
  case 4:                                                   
    pTitle=TString("|F_{x}|^{2}");                          
    break;                                                  
  case 5:                                                   
    pTitle=TString("|F_{+}|^{2}");                          
    break;                                                  
  case 6:                                                   
    pTitle=TString("F_{x}");                                
    break;                                                  
  case 7:                                                   
    pTitle=TString("F_{+}");                                
    break;                                                  
  case 8:                                                   
    pTitle=TString("F1_{x}/F2_{x}");                        
    break;                                                  
  case 9:                                                   
    pTitle=TString("F1_{+}/F2_{+}");                        
    break;                                                  
  case 10:                                                  
    pTitle=TString("#sqrt{|F1_{+}|^{2}+|F1_{x}|^{2})/sqrt(|F2_{+}|^{2}+|F2_{x}|^{2}}");
    break;                                                                             
  case 12:                                                  
    pTitle=TString("DPF Angle");
    break;                                                                             
  default:                                                                             
    break;                                                                             
  }                                                                                    
  TString netName="";                                                                  
  for(int n=0; n<nIFO; n++) {                                                          
    TString ifoName=this->getifo(n)->Name;                                                   
         if(ifoName.CompareTo("A1")==0) ifoName=TString("#tilde{A}");                  
    else if(ifoName.CompareTo("H2")==0) ifoName=TString("#tilde{H}");                  
    else if(ifoName.CompareTo("I2")==0) ifoName=TString("#tilde{I}");                  
    else if(ifoName.CompareTo("R2")==0) ifoName=TString("#tilde{R}");                  
    else ifoName.Resize(1);                                                            
    netName+=ifoName;                                                                  
  }                                                                                    
  if(btitle) gSM.h2->SetTitle("Network = "+netName+                                        
                          "                                    "+                      
                          "Antenna Pattern =  "+pTitle);                               

  std::vector<wavecomplex> F(nIFO);
  int size=int(180*2*gSM.resolution*360*2*gSM.resolution);
  double* x = new double[size];                   
  double* y = new double[size];                   
  double* z = new double[size];                   
  int cnt=0;                                      
  for (int i=0;i<360*2*gSM.resolution;i++) {          
    for (int j=0;j<180*2*gSM.resolution;j++) {        
      double phi = (double)(i+1)/(double)(2*gSM.resolution);
      double theta = (double)(j+1)/(double)(2*gSM.resolution);
      double psi=0;                                       
      x[cnt]=phi;                                         
      y[cnt]=theta;                                       
      if (polarization!=11) {                             
/*                                                        
        psi=0;                                            
        for (int n=0;n<nIFO;n++) F[n] = pD[n]->antenna(theta,phi,psi);
        double gp=0,gx=0,gI=0;                                        
        for (int n=0;n<nIFO;n++) {                                    
          gp+=F[n].real()*F[n].real();                                
          gx+=F[n].imag()*F[n].imag();                                
          gI+=F[n].real()*F[n].imag();                                
        }                                                             
        double gR = (gp-gx)/2.;                                       
        double gr = (gp+gx)/2.;                                       
        double gc = sqrt(gR*gR+gI*gI);                                

        if (polarization==0) z[cnt] = sqrt(fabs(gr-gc)<1e-8?0.:gr-gc);
        if (polarization==1) z[cnt] = sqrt(gr+gc);                    
        if (polarization==2&&sqrt(gr+gc)>0) z[cnt]=sqrt((gr-gc)/(gr+gc));
        if (polarization==3) z[cnt] = sqrt(2*gr/nIFO);                   
        //if (polarization==3) z[cnt] = sqrt(2*gr);                      
        if (polarization==4) z[cnt] = (fabs(gr-gc)<1e-8?0.:gr-gc);       
        if (polarization==5) z[cnt] = (gr+gc);                           
        if (polarization==6) z[cnt] = F[0].imag();                       
        if (polarization==7) z[cnt] = F[0].real();                       
        if (polarization==8) z[cnt] = (F[1].imag()!=0) ? F[0].imag()/F[1].imag() : 0.; 
        if (polarization==9) z[cnt] = (F[1].real()!=0) ? F[0].real()/F[1].real() : 0.; 
        //if (polarization==8) z[cnt] = F[0].imag()/F[1].imag()>0 ? 1. : -1.;          
        //if (polarization==9) z[cnt] = F[0].real()/F[1].real()>0 ? 1. : -1.;          
        if (polarization==10) {                                                        
          double F0=sqrt(pow(F[0].real(),2)+pow(F[0].imag(),2));                       
          double F1=sqrt(pow(F[1].real(),2)+pow(F[1].imag(),2));                       
          z[cnt] = (F1!=0) ? F0/F1 : 0.;                                               
        }                                                                              
*/                                                                                     
        z[cnt] = GetAntennaPattern(phi, theta, psi, polarization);                     
      } else {                                                                         
        int counts=0;                                                                  
        z[cnt] = 0.;                                                                   
        //for (int k=0;k<360;k+=8) {                                                   
        for (int k=0;k<1;k+=1) {                                                       
          psi=(double)k;                                                               
          for (int n=0;n<nIFO;n++) F[n] = pD[n]->antenna(theta,phi,psi);               
          //double F0=sqrt(pow(F[0].real(),2)+pow(F[0].imag(),2));                     
          //double F1=sqrt(pow(F[1].real(),2)+pow(F[1].imag(),2));                     
          double F0=F[0].real()+F[0].imag();                                           
          double F1=F[1].real()+F[1].imag();                                           
          z[cnt] += (F1!=0) ? fabs(F0/F1) : 0.;                                        
        }                                                                              
        z[cnt] /= (double)counts;                                                      
      }                   
      // fill skymap                                                             
      int ind = gSM.getSkyIndex(theta,phi);
      gSM.set(ind,z[cnt]);  
      cnt++;                                                                           
    }                                                                                  
  }                                                                                    

  if (gSM.coordinate.CompareTo("GEOGRAPHIC")==0) {
    for (int i=0;i<size;i++) CwbToGeographic(x[i],y[i],x[i],y[i]);
  }                                                               
  if (gSM.coordinate.CompareTo("CELESTIAL")==0) {                     
    for (int i=0;i<size;i++) CwbToCelestial(x[i],y[i],x[i],y[i]); 
  }                                                               
  gSM.FillData(size, x, y, z);                                        
  gSM.Draw(dpaletteId);                                               

  delete [] x;
  delete [] y;
  delete [] z;
}             

//______________________________________________________________________________
double 
gnetwork::GetDelay(TString ifo1, TString ifo2, TString mode) {
//
// Return the delay (sec) between two detector
//
// Input: ifo1	- detector name of the first detector  (Ex: "L1")
//        ifo2	- detector name of the second detector (Ex: "V1")
//        mode  - delay mde
//                ""  : (max-min)/2 - default
//                min : minum delay value
//                max : maximum delay value
//                MAX : maximum absolute delay value
//

  int nIFO = this->ifoListSize();
  if(nIFO<2) {cout << "nIFO must be >= 2" << endl;return -1.;}

  if(ifo1.CompareTo(ifo2)==0) return 0.0;

  detector *pD[NIFO];
  for(int n=0; n<nIFO; n++) pD[n] = this->getifo(n);

  int ifoId1=0;
  int ifoId2=0;
  int nifo=0;  
  for(int n=0; n<nIFO; n++) {
    if(ifo1.CompareTo(this->getifo(n)->Name)==0) {nifo++;ifoId1=n;}
    if(ifo2.CompareTo(this->getifo(n)->Name)==0) {nifo++;ifoId2=n;}
  }                                                          
  if(nifo<2) {cout << "gnetwork::GetDelay ifo not valid" << endl;return -1.;}

  skymap sm0, sTau;
  double maxTau = -1.;
  double minTau =  1.;
  double tmax, tmin;  

  tmax = tmin = 0.;
  for(int n=0; n<nifo; n++) {
    sTau = pD[n]->tau;       
    tmax = sTau.max();       
    tmin = sTau.min();       
    if(tmax > maxTau) maxTau = tmax;
    if(tmin < minTau) minTau = tmin;
  }                                 
  if(strstr(mode,"min")) return minTau;
  if(strstr(mode,"max")) return maxTau;
  if(strstr(mode,"MAX")) return fabs(maxTau)>fabs(minTau) ? fabs(maxTau) : fabs(minTau);
  return (maxTau-minTau)/2.;                                                            
}                                                                                       

//______________________________________________________________________________
double 
gnetwork::GetDelay(TString ifo1, TString ifo2, double phi, double theta) {
//
// Return the delay (sec) between two detector from the direction phi,theta
//
// Input: ifo1  - detector name of the first detector  (Ex: "L1")
//        ifo2  - detector name of the second detector (Ex: "V1")
//        phi   - phi direction (degrees)
//        theta - theta direction (degrees)
//

  int nIFO = this->ifoListSize();
  if(nIFO<2) {cout << "gnetwork::GetDelay nIFO must be >= 2" << endl;return -1.;}

  if(ifo1.CompareTo(ifo2)==0) return 0.0;

  detector *pD[NIFO];
  for(int n=0; n<nIFO; n++) pD[n] = this->getifo(n);

  int ifoId1=0;
  int ifoId2=0;
  int nifo=0;  
  for(int n=0; n<nIFO; n++) {
    if(ifo1.CompareTo(this->getifo(n)->Name)==0) {nifo++;ifoId1=n;}
    if(ifo2.CompareTo(this->getifo(n)->Name)==0) {nifo++;ifoId2=n;}
  }                                                          
  if(nifo<2) {cout << "gnetwork::GetDelay ifo not valid" << endl;return -1.;}

  return pD[ifoId1]->getTau(theta,phi)-pD[ifoId2]->getTau(theta,phi);
}                                                                    

//______________________________________________________________________________
void 
gnetwork::DrawDelay(TString ifo1, TString ifo2, double phi, double theta, double frequency) {

  if(gSM.goff) return;

  int nIFO = this->ifoListSize();
  if(nIFO<2) {cout << "gnetwork::DrawDelay nIFO must be >= 2" << endl;return;}

  if(ifo1.CompareTo(ifo2)==0) {cout << "gnetwork::DrawDelay ifo1 must != ifo2" << endl;return;}

  double ifo1d2_pos_delay=0.;
  // if ifo1d2_pos_delay!=0 draw delay respect to ifo1d2_pos_delay
  if(phi!=1000) ifo1d2_pos_delay=GetDelay(ifo1, ifo2, phi, theta);
                                                                  
  //double ifo1d2_max_delay=0.;                                   
  //if(frequency>0) ifo1d2_max_delay=GetDelay(ifo1, ifo2, "max"); 

  detector *pD[NIFO];
  for(int n=0; n<nIFO; n++) pD[n] = this->getifo(n);

  int ifoId1=0;
  int ifoId2=0;
  int nifo=0;  
  for(int n=0; n<nIFO; n++) {
    if(ifo1.CompareTo(this->getifo(n)->Name)==0) {nifo++;ifoId1=n;}
    if(ifo2.CompareTo(this->getifo(n)->Name)==0) {nifo++;ifoId2=n;}
  }                                                          
  if(nifo<2) {cout << "gnetwork::DrawDelay ifo not valid" << endl;return;}

  int size=int(180*2*gSM.resolution*360*2*gSM.resolution);
  double* x = new double[size];                   
  double* y = new double[size];                   
  double* z = new double[size];                   
  int cnt=0;                                      
  for (int i=0;i<360*2*gSM.resolution;i++) {          
    for (int j=0;j<180*2*gSM.resolution;j++) {        
      double phi = (double)(i+1)/(double)(2*gSM.resolution);
      double theta = (double)(j+1)/(double)(2*gSM.resolution);

      double delay=pD[ifoId1]->getTau(theta,phi)-pD[ifoId2]->getTau(theta,phi);

      delay-=ifo1d2_pos_delay;
      if(frequency>0) {       
        double fdelay=fabs(fmod(delay,0.5/frequency));
        int idelay = int((fabs(delay)-fdelay)/(0.5/frequency));
        delay = idelay%2==0 ? fdelay : (0.5/frequency)-fdelay; 
      }                                                        

      x[cnt]=phi;
      y[cnt]=theta;
      z[cnt]=delay;
      cnt++;       
    }              
  }                

  if (gSM.coordinate.CompareTo("GEOGRAPHIC")==0) {
    for (int i=0;i<size;i++) CwbToGeographic(x[i],y[i],x[i],y[i]);
  }                                                               
  if (gSM.coordinate.CompareTo("CELESTIAL")==0) {                     
    for (int i=0;i<size;i++) CwbToCelestial(x[i],y[i],x[i],y[i]); 
  }                                                               
  gSM.FillData(size, x, y, z);                                        
  gSM.Draw(0);                                                        

  delete [] x;
  delete [] y;
  delete [] z;
}             

//______________________________________________________________________________
int 
gnetwork::Delay2Coordinates(double t1, double t2, double t3, double& ph1, double& th1, double& ph2,double& th2) {
//
// Return the direction & mirror-direction uding the input arrival times of the 3 detectors
// NOTE: Works only with 3 detectors
//       With 3 detector there are 2 possible solutions
//
// Input:   t1,t2,t3   - arrival times of the 3 detectors (sec)
//
// Output:  ph1,th1    - phi,theta coordinates of the first direction
//          ph2,th2    - phi,theta coordinates of the second direction
//
// Return 0/1          - success/error
//	

  int nIFO = this->ifoListSize();
  if(nIFO!=3) {cout << "gnetwork::Delay2Source nIFO must be = 3" << endl;return 1;}

  TString ifoName[3];
  for(int n=0; n<nIFO; n++) ifoName[n]=TString(this->getifo(n)->Name);
  //for(int n=0; n<nIFO; n++) cout << n << " " << ifoName[n].Data() << " " << ifoId[n] << endl;
  for(int n=0; n<nIFO; n++) for(int m=n+1; m<nIFO; m++)                                        
    if(ifoName[n].CompareTo(ifoName[m])==0)                                                    
       {cout << "gnetwork::Delay2Source ifos must differents" << endl;return 2;}                
                                                                                               
  detector *pD[NIFO];                                                                          
  for(int n=0; n<nIFO; n++) pD[n] = this->getifo(n);                                            

  if(!gROOT->GetClass("XYZVector")) gSystem->Load("libMathCore");

  XYZVector vD[3];
  for(int n=0; n<nIFO; n++) vD[n].SetXYZ(pD[n]->Rv[0],pD[n]->Rv[1],pD[n]->Rv[2]);
  for(int n=0; n<nIFO; n++) vD[n]/=speedlight;                                   

  double a=sqrt((vD[1]-vD[0]).Mag2());
  double b=sqrt((vD[2]-vD[1]).Mag2());
  double c=sqrt((vD[2]-vD[0]).Mag2());

  double s=(a+b+c)/2;
  double b2=(2/a)*sqrt(s*(s-a)*(s-b)*(s-c));
  double b1=sqrt(c*c-b2*b2);                

  double delta=pow(b1*(t1-t2)-a*(t1-t3),2)+pow(b2*(t1-t2),2);

  double stheta = sqrt(delta)/(a*b2);
  double sphi = -b2*(t1-t2)/sqrt(delta);

  if (stheta>1) stheta=1;
  if (stheta<-1) stheta=-1;

  double theta = asin(stheta);
  double phi   = asin(sphi);  

  double ctheta = sqrt(pow(a*b2,2)-delta)/(a*b2);
  double cphi = (a*(t1-t3)-b1*(t1-t2))/sqrt(delta);

  if (TMath::IsNaN(ctheta)) ctheta=-1;
  if (ctheta>1) ctheta=-1;            
  if (ctheta<-1) ctheta=-1;           

  double thetacn = acos(-ctheta);

  XYZVector D1 = vD[1]-vD[0];
  XYZVector D3 = D1.Cross(vD[2]-vD[0]);
  XYZVector D2 = D1.Cross(D3);         

  XYZVector eD1 = D1.Unit();
  XYZVector eD2 = D2.Unit();
  XYZVector eD3 = D3.Unit();

  Rotation3D R(eD1,eD2,eD3);

  //XYZVector IS;
  //XYZVector S; 
  //S.SetR(1);   
  TVector3 IS;   
  TVector3 S(1,0,0); 

  if (cphi<0) {

    S.SetTheta(theta);
    S.SetPhi((TMath::PiOver2()+phi));
    IS = R*S;                        
    th1=IS.Theta()*180./TMath::Pi(); 
    ph1=IS.Phi()*180./TMath::Pi();   
    if (ph1<0) ph1+=360;             
    if (ph1>360) ph1-=360;           

    S.SetTheta(thetacn);
    S.SetPhi((TMath::PiOver2()+phi));
    IS = R*S;                        
    th2=IS.Theta()*180./TMath::Pi(); 
    ph2=IS.Phi()*180./TMath::Pi();   
    if (ph2<0) ph2+=360;             
    if (ph2>360) ph2-=360;           

  } else {

    S.SetTheta(theta);
    S.SetPhi(-(TMath::PiOver2()+phi));
    IS = R*S;                         
    th1=IS.Theta()*180./TMath::Pi();  
    ph1=IS.Phi()*180./TMath::Pi();    
    if (ph1<0) ph1+=360;              
    if (ph1>360) ph1-=360;            

    S.SetTheta(thetacn);
    S.SetPhi(-(TMath::PiOver2()+phi));
    IS = R*S;                         
    th2=IS.Theta()*180./TMath::Pi();  
    ph2=IS.Phi()*180./TMath::Pi();    
    if (ph2<0) ph2+=360;              
    if (ph2>360) ph2-=360;            
  }                                   

  return 0;
}          

//______________________________________________________________________________
void 
gnetwork::MakeNetworkDetectorIndex(double gamma, int loop, double snr) {
//
// Make Network Detector Index (1G analysis) and store it to the internal gSM sky map
//
// Input: gamma  - 1G gamma value
//        loop   - increase the statistic
//        snr    - add noise to signal -> snr          
//

  if(gamma>1) {cout << "gamma must be <=1" << endl;exit(1);}             
  if(loop<1) {cout << "gamma must be >=1" << endl;exit(1);}              

  int nIFO = this->ifoListSize();
  if(nIFO==0) {cout << "detectors are no declared" << endl;exit(1);}

  detector *pD[NIFO];
  for(int n=0; n<nIFO; n++) pD[n] = this->getifo(n);

  std::vector<wavecomplex> F(nIFO);
  double* X = new double[nIFO];
  double* u = new double[nIFO];
  double* v = new double[nIFO];

  int L = gSM.size();

  double h2max=0.;
  gRandom->SetSeed(0);
  for (int l=0;l<L;l++) {
    double ndi=0.;       
    for (int i=0;i<loop;i++) {

      double phi=gSM.getPhi(l);
      double theta=gSM.getTheta(l);
      double psi=gRandom->Uniform(0,180);
                                         
      double hp=gRandom->Uniform(-1,1);  
      double hx=gRandom->Uniform(-1,1);  

      double E=0.;
      for (int n=0;n<nIFO;n++) {
        F[n] = pD[n]->antenna(theta,phi,psi);
        X[n] = hp*F[n].real()+hx*F[n].imag();
        if(snr>=0.0) E+=X[n]*X[n];           
      }                                      
      // add noise to signal -> snr          
      E=sqrt(E);                             
      if(snr>=0.0) for (int n=0;n<nIFO;n++) X[n] = X[n]*snr/E + gRandom->Gaus(0.,1.);
                                                                                     
      double gp=0,gx=0,gI=0;                                                         
      for (int n=0;n<nIFO;n++) {                                                     
        gp+=F[n].real()*F[n].real();                                                 
        gx+=F[n].imag()*F[n].imag();                                                 
        gI+=F[n].real()*F[n].imag();                                                 
      }                                                                              
      gp+=1.e-12;                                                                    
      gx+=1.e-12;                                                                    
                                                                                     
      double Xp=0;                                                                   
      double Xx=0;                                                                   
      for (int k=0;k<nIFO;k++) {                                                     
        Xp+=X[k]*F[k].real();                                                        
        Xx+=X[k]*F[k].imag();                                                        
      }                                                                              
                                                                                     
      double uc = (Xp*gx - Xx*gI);            // u cos of rotation to PCF            
      double us = (Xx*gp - Xp*gI);            // u sin of rotation to PCF            
      double vc = (gp*uc + gI*us);            // (u*f+)/|u|^2 - 'cos' for v          
      double vs = (gx*us + gI*uc);            // (u*fx)/|u|^2 - 'sin' for v          
      double um=0;                                                                   
      double vm=0;                                                                   
      for (int k=0;k<nIFO;k++) {                                                     
        u[k]=0;                                                                      
        v[k]=0;                                                                      
        u[k]=F[k].real()*uc+F[k].imag()*us;                                          
        um+=u[k]*u[k];                                                               
        v[k]=-F[k].imag()*vc+F[k].real()*vs;                                         
        vm+=v[k]*v[k];                                                               
      }                                                                              
      vm+=1.e-24;                                                                    
      if(gamma>0) {                                                                  
        double GAMMA = 1.-gamma*gamma;                // network regulator           
        ndi+=this->netx(u,um,v,vm,GAMMA);                                             
      } else ndi+=this->ndi(u,um,nIFO);                                              
/*                                                                                   
      } else {                                                                       
        double xndi=this->ndi(u,um,nIFO);                                            
        if(xndi>=-gamma) ndi+=1;                                                     
      }                                                                              
*/                                                                                   
    }                                                                                
    double xndi=ndi/loop;                                                            
    if(h2max<xndi) h2max=xndi;                                                       
    gSM.set(l,xndi);                                                             
  }                                                                                  

  //if(gamma<0) h2->GetZaxis()->SetRangeUser(1,nIFO);
  if(gamma<0) gSM.h2->GetZaxis()->SetRangeUser(1,h2max); 

  delete [] X;
  delete [] u;
  delete [] v;
}             

//______________________________________________________________________________
void 
gnetwork::DrawNetworkDetectorIndex(double gamma, int loop, double snr, int dpaletteId, bool btitle) {
//
// Draw Network Detector Index sky map (1G analysis) 
//
// Input: gamma        - 1G gamma value
//        loop         - increase the statistic
//        snr          - add noise to signal -> snr          
//
//        dpaletteId   - palette ID
//        btitle       - if true then a default title is provided
//

  if(gSM.size()==0) {
    cout << "gnetwork::DrawNetworkDetectorIndex : Error - network skymap not initialized" << endl;
    exit(1);
  }

  MakeNetworkDetectorIndex(gamma, loop, snr);

  int nIFO = this->ifoListSize();

  TString netName="";
  for(int n=0; n<nIFO; n++) {
    TString ifoName=this->getifo(n)->Name;
         if(ifoName.CompareTo("A1")==0) ifoName=TString("#tilde{A}"); 
    else if(ifoName.CompareTo("H2")==0) ifoName=TString("#tilde{H}"); 
    else if(ifoName.CompareTo("I2")==0) ifoName=TString("#tilde{I}"); 
    else if(ifoName.CompareTo("R2")==0) ifoName=TString("#tilde{R}"); 
    else ifoName.Resize(1);                                           
    netName+=ifoName;                                                 
  }                                                                   
  if(btitle) gSM.h2->SetTitle("Network = "+netName+                       
                          "                                    "+     
                          "Network Detector Index");                  

  gSM.FillData();
  gSM.Draw(dpaletteId);                                                        
}

//______________________________________________________________________________
void 
gnetwork::DrawCircles(double phi, double theta, Color_t lcolor, Width_t lwidth, Style_t lstyle, bool labels) {
//
// Draw iso delay circles for each couple of detectors along the direction phi,theta
//
// Input: phi,theta - sky direction coordinates (degrees)
//        lcolor    - line color 
//        lwidth    - line width 
//        lstyle    - line style 
//

  DrawCircles(phi, theta, 0., lcolor, lwidth, lstyle, labels);                                          
}                                                                                               

//______________________________________________________________________________
void 
gnetwork::DrawCircles(double phi, double theta, double gps, Color_t lcolor, Width_t lwidth, Style_t lstyle, bool labels) {
//
// Draw iso delay circles for each couple of detectors along the direction phi,theta at time gps
//
// Input: phi,theta - sky direction coordinates (degrees)
//        gps       - gps time (sec)
//        lcolor    - line color 
//        lwidth    - line width 
//        lstyle    - line style 
//

  if(gSM.goff) return;
  if(gSM.canvas==NULL) return;

  int n,m,l;

  int nIFO = this->ifoListSize();
  if(nIFO<2) return;            

  detector *pD[NIFO];
  for(n=0; n<nIFO; n++) pD[n] = this->getifo(n);

  gSM.canvas->cd();

  double phi_off = gps>0 ? gSM.phi2RA(phi, gps) : phi;
  phi_off=fmod(phi_off-phi+360,360);                      

  double ph1,th1,ph2,th2;
  double r2d = 180./TMath::Pi();
  double d2r = TMath::Pi()/180.;
//Color_t colors[6] = {kBlue,kRed,kYellow,kGreen,kBlack,kWhite};
  Style_t marker[3] = {20,21,22}; // circle,square,triangle     
  int ncircle=0;                                                
  for(n=0;n<nIFO;n++) {                                         
    for(m=n+1;m<nIFO;m++) {                                     

      XYZVector D1(pD[n]->Rv[0],pD[n]->Rv[1],pD[n]->Rv[2]);
      XYZVector D2(pD[m]->Rv[0],pD[m]->Rv[1],pD[m]->Rv[2]);

      D1=D1/speedlight;
      D2=D2/speedlight;

      XYZVector D12 = D1-D2;
      if(D12.Mag2()==0) continue;      // the 2 detectors are located in the same site
      TVector3 vD12(D12.X(),D12.Y(),D12.Z());

      TVector3 vSD(1,0,0);
      // coordinates are transformed in the CWB frame
      if (gSM.coordinate.CompareTo("GEOGRAPHIC")==0) {   
        GeographicToCwb(phi,theta,ph1,th1);          
        vSD.SetTheta(th1*d2r);                       
        vSD.SetPhi(ph1*d2r);                         
      } else if (gSM.coordinate.CompareTo("CELESTIAL")==0) {
        CelestialToCwb(phi,theta,ph1,th1);              
        ph1=fmod(ph1-phi_off+360,360);                  
        vSD.SetTheta(th1*d2r);                          
        vSD.SetPhi(ph1*d2r);                            
      } else {                                          
        vSD.SetTheta(theta*d2r);                        
        vSD.SetPhi(phi*d2r);                            
      }                                                 

      TVector3 rvSD;
      double as;    
      int L=4*360;  
      wavearray<float> th(L);
      wavearray<float> ph(L);
      for (l=0;l<L;l++) {    
        as = l*d2r/4.;       
        TRotation vR;        
        vR.Rotate(as,vD12);  
        rvSD = vR*vSD;       
        th[l] = rvSD.Theta()*r2d;
        ph[l] = rvSD.Phi()*r2d+phi_off;
        ph[l]=fmod(ph[l]+360,360);
      }                                     

      wavearray<double> x(L);
      wavearray<double> y(L);
      int P=0;               
      for (l=0;l<L;l++) {    
        // the main circle is shown in sub-circles to avoid to draw fake lines from ph=0 to 360
        // if((fabs(ph[l]-ph[l==0?L-1:l-1])<180)&&(l<L-1)) {                          
        // TVector3 is in radians [phi=-Pi,+Pi] [theta=0,Pi]                                   
        ph1=ph[l];                                                                        
        ph2=ph[l==0?L-1:l-1];                                                             
        if (gSM.coordinate=="GEOGRAPHIC") {                                           
          CwbToGeographic(ph1,th1,ph1,th1);                                                    
          CwbToGeographic(ph2,th2,ph2,th2);                                                    
        }                                                                                      
        if((fabs(ph1-ph2)<180)&&(l<L-1)) {                                                     
          x[P]=ph[l]; y[P]=th[l];                                          

          double X=x[P];
          if (gSM.projection=="hammer") {
            if(X>180) {x[P]+=360; y[P]-=90;} 
            else      {y[P]-=90;}
            gSM.ProjectHammer(x[P], y[P], x[P], y[P]); 
            y[P]=y[P]+90;                     
          } else if (gSM.projection=="sinusoidal") {                    
            if(X>180) {x[P]-=360; y[P]+=90;} 
            else      {y[P]-=90;} 
            gSM.ProjectSinusoidal(x[P], y[P], x[P], y[P]); 
            if(X>180) {x[P]=-x[P]; y[P]-=90;}                         
            else      {y[P]+=90;}
          } else if (gSM.projection=="parabolic") {                         
            if(X>180) {x[P]-=360; y[P]-=90;} 
            else      {y[P]-=90;} 
            gSM.ProjectParabolic(x[P], y[P], x[P], y[P]); 
            y[P]+=90;                         
          } else {                                                                   
//            x[P]=x[P]-180; y[P]=y[P]-90;                       
          }                                                                          
          if (gSM.coordinate=="GEOGRAPHIC") CwbToGeographic(x[P],y[P],x[P],y[P]);                
          if (gSM.coordinate=="CELESTIAL")  CwbToCelestial(x[P],y[P],x[P],y[P]);                 
          P++;                                                                       
        } else {                                                                     
          TPolyLine *pL = new TPolyLine(P,x.data,y.data);                            
          //pL->SetLineColor(colors[ncircle%6]);                                     
          pL->SetLineColor(lcolor);                                                  
          pL->SetLineStyle(lstyle);                                                  
          pL->SetLineWidth(lwidth);                                                  
          pL->Draw();                                                                
          circleL.push_back(pL);                                                     
          P=0;                                                                       
        }                                                                            
      }                                                                              
      // draw detector axis                                                          
      double vph,vth;                                                                
      if(vSD.Dot(vD12)>0) {                                                          
        vph=D12.Phi()*r2d;                                                           
        vth=D12.Theta()*r2d;                                                         
      } else {                                                                       
        vph=(-D12).Phi()*r2d;                                                        
        vth=(-D12).Theta()*r2d;                                                      
      }                                                                              
      vph=fmod(vph+phi_off+360,360);                                                 

      double X=vph;
      if (gSM.projection=="hammer") {
        if(X>180) {vph+=360; vth-=90;}
        if(X<180) {vth-=90;}
        gSM.ProjectHammer(vph, vth, vph, vth); 
        vth=vth+90;                     
      } else if (gSM.projection=="sinusoidal") {
        if(X>180) {vph-=360; vth+=90;} 
        else      {vth-=90;} 
        gSM.ProjectSinusoidal(vph, vth, vph, vth);
        if(X>180) {vph=-vph; vth-=90;}                         
        else      {vth+=90;}
      } else if (gSM.projection=="parabolic") { 
        if(X>180) {vph-=360; vth-=90;} 
        else      {vth-=90;} 
        gSM.ProjectParabolic(vph, vth, vph, vth); 
        vth+=90;
      } else {                                           
//        vph=vph-180; vth=vth-90;                       
      }                                                  
      if (gSM.coordinate=="GEOGRAPHIC") CwbToGeographic(vph,vth,vph,vth);                
      if (gSM.coordinate=="CELESTIAL")  CwbToCelestial(vph,vth,vph,vth);                 
      //if (coordinate.CompareTo("CWB")==0) {vph+=180;vth+=90;}
      TMarker* pM = new TMarker(vph,vth,marker[ncircle%3]);    
      pM->SetMarkerSize(1.0);                                  
      pM->SetMarkerColor(kBlack);                              
      if(labels) pM->Draw();                                              
      daxisM.push_back(pM);                                    

      // draw ifos pair labels on the detector axis
      char ifoPair[16];sprintf(ifoPair,"%c%c",pD[n]->Name[0],pD[m]->Name[0]);
      TLatex *pL = new TLatex(vph+4.0,vth,ifoPair);
      pL->SetTextFont(32);
      pL->SetTextSize(0.045);
      pL->SetTextColor(kBlack);
      if(labels) pL->Draw();
      siteP.push_back(pL);

      ncircle++;
    }
  }  

  return;
}        

//______________________________________________________________________________
void gnetwork::DumpObject(char* file)
{
//
// dump gnetwork object into root file
//
// Input: file - root file name
//

  TFile* rfile = new TFile(file, "RECREATE");
  this->Write("gnetwork"); // write this object
  rfile->Close();
}

//______________________________________________________________________________
void gnetwork::LoadObject(char* file)
{
//
// load gnetwork object from root file
//
// Input: file - root file name
//

  TFile* rfile = new TFile(file);
  gnetwork* gnet = (gnetwork*)rfile->Get("gnetwork;1");
  if(gnet==NULL) {
    cout << "gnetwork::LoadObject : Error - input file don't contains object gnetwork" << endl;
    exit(1);
  }
  *this=*gnet;

  siteM.clear();
  siteA.clear();
  siteT.clear();
  siteL.clear();
  siteP.clear();
  circleL.clear();
  daxisM.clear();

  gSM.FillData();
  rfile->Close();
}

//______________________________________________________________________________
TCanvas* gnetwork::DrawPolargram(int ptype, network* net)
{
//
// Plot polargrams of the event pixels in the network plane
//
// Input: ptype - plot type [0/1]
//                0 : projection on network plane of 00/90 amplitudes
//                1 : standard response of 00/90 amplitudes
//        net   - pointer to network object
//                if net!=NULL then net is used to retrive the data pixels to be plotted
//
// Output: return canvas object
//

  if(ptype!=0 && ptype!=1) {
    cout << "gnetwork::DrawPolargram : Error - wrong input parameter, must be [0/1]" << endl;
    exit(1);
  }

  network* NET = net!=NULL ? net : this;

  // get size  
  int size = ptype==0 ? NET->p00_POL[0].size() : NET->r00_POL[0].size();    
  if(size==0) {
    cout << "gnetwork::DrawPolargram : Warning - event pixels size=0" << endl;
    return NULL;
  }

  // find max radius
  double rmax=0;     
  if(ptype==0) {
    for(int i=0;i<size;i++) if(NET->p00_POL[0][i]>rmax) rmax=NET->p00_POL[0][i];
    for(int i=0;i<size;i++) if(NET->p90_POL[0][i]>rmax) rmax=NET->p90_POL[0][i];
  }
  if(ptype==1) {
    for(int i=0;i<size;i++) if(NET->r00_POL[0][i]>rmax) rmax=NET->r00_POL[0][i];
    for(int i=0;i<size;i++) if(NET->r90_POL[0][i]>rmax) rmax=NET->r90_POL[0][i];
  }
  rmax=TMath::Nint(1.1*rmax+0.5);                                             

  for(int n=0;n<2;n++) if(grP[n]!=NULL) delete grP[n];
  if(canvasP!=NULL) delete canvasP;
  canvasP = new TCanvas("canvasP","polargram",800,800);
  gStyle->SetLineColor(kBlack);                            
  double* rpol;                                            
  double* apol;                                            
  for(int n=0;n<2;n++) {                                       
    // projection on network plane 00 amplitudes           
    if(ptype==0 && n==0) {rpol = NET->p00_POL[0].data; apol = NET->p00_POL[1].data;}
    // projection on network plane 90 amplitudes                                
    if(ptype==0 && n==1) {rpol = NET->p90_POL[0].data; apol = NET->p90_POL[1].data;}
    // standard response 00 amplitudes                                          
    if(ptype==1 && n==0) {rpol = NET->r00_POL[0].data; apol = NET->r00_POL[1].data;}
    // standard response 90 amplitudes                                          
    if(ptype==1 && n==1) {rpol = NET->r90_POL[0].data; apol = NET->r90_POL[1].data;}
    grP[n] = new TGraphPolar(size,apol,rpol);                                   

    grP[n]->SetMarkerStyle(20);
    grP[n]->SetMarkerSize(0.8);
    grP[n]->SetTitle(" ");     

    if(n==0) {
      grP[n]->SetMarkerColor(kBlue);
      grP[n]->SetLineColor(kBlue);  
      grP[n]->Draw("P");            
    }                               
                                    
    if(n==1) {                      
      grP[n]->SetMarkerColor(kRed); 
      grP[n]->SetLineColor(kRed);   
      grP[n]->Draw("PSAME");        
    }                               

    canvasP->Update();
    grP[n]->GetPolargram()->SetTextColor(8);
    grP[n]->GetPolargram()->SetRangePolar(-TMath::Pi(),TMath::Pi());
    grP[n]->SetMinRadial(0);                                        
    grP[n]->SetMaxRadial(rmax);                                     
    grP[n]->GetPolargram()->SetNdivPolar(612);                      
    grP[n]->GetPolargram()->SetNdivRadial(504);                     
    grP[n]->GetPolargram()->SetToDegree();                          
    grP[n]->GetPolargram()->SetPolarLabelSize(0.03);                
    grP[n]->GetPolargram()->SetRadialLabelSize(0.03);               
  }                                                                 

  canvasP->Update();

  return canvasP;
}                                  

//______________________________________________________________________________
void
gnetwork::SetScale(TString ifo, double scale) {
//
// set detector scale factor (PSD @ user given frequency)
// return detector scale factor 
//
// Input: ifo     - detector name
//
// if ifo invalid return -1
//

  int nIFO = this->ifoListSize();
  if(nIFO<1) {cout << "gnetwork::GetScale nIFO must be >= 1" << endl;exit(-1);}

  detector *pD[NIFO];
  for(int n=0; n<nIFO; n++) pD[n] = this->getifo(n);

  int ifoId=0;
  int nifo=0; 
  for(int n=0; n<nIFO; n++) {
    if(ifo.CompareTo(this->getifo(n)->Name)==0) {nifo++;ifoId=n;}
  }                                                        

  this->scale[ifoId]=scale;
}

//______________________________________________________________________________
double
gnetwork::GetScale(TString ifo) {
//
// return detector scale factor (PSD @ user given frequency)
//
// Input: ifo     - detector name
//
// if ifo invalid return -1
//

  int nIFO = this->ifoListSize();
  if(nIFO<1) {cout << "gnetwork::GetScale nIFO must be >= 1" << endl;return -1.;}

  detector *pD[NIFO];
  for(int n=0; n<nIFO; n++) pD[n] = this->getifo(n);

  int ifoId=0;
  int nifo=0; 
  for(int n=0; n<nIFO; n++) {
    if(ifo.CompareTo(this->getifo(n)->Name)==0) {nifo++;ifoId=n;}
  }                                                        

  return scale[ifoId];
}

