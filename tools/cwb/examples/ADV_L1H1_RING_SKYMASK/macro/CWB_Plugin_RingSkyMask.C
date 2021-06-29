#define XIFO 4

#pragma GCC system_header

#include "cwb.hh"
#include "cwb2G.hh"
#include "config.hh"
#include "network.hh"
#include "wavearray.hh"
#include "TString.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "TRandom.h"
#include "Toolbox.hh"
#include "skymap.hh"
#include "gskymap.hh"
#include "skycoord.hh"
#include "TVector3.h"
#include "TRotation.h"
#include "Math/Vector3Dfwd.h"
#include "Math/Rotation3D.h"

using namespace CWB;

// 
// Position of SN2007gr on the sky:

#define SN2007gr_NAME		"SN2007gr"
#define SN2007gr_GPS_START	870772910
#define SN2007gr_GPS_END	871215278
#define SN2007gr_RA		40.8666
#define SN2007gr_DEC		37.3458

#define MIN_RING_WIDTH		5		// define the min width of the ring (deg)

#define SAVE_SKYMASK_PLOT
#define APPLY_SKYMASK

static	bool first = true;

skymap GetRingSkymask(network* net, config* cfg, double gps_start, double gps_end);

void CWB_Plugin(TFile* jfile, CWB::config* cfg, network* net, WSeries<double>* x, TString ifo, int type) {
    
  cout << endl;
  cout << "-----> CWB_Plugin_RingSkyMask.C - type : " << type << endl;
  cout << endl;

  if(type==CWB_PLUGIN_ILIKELIHOOD) {

    if(first) {						// this code is executed only one time

      first=false;

      cwb2G* gCWB2G; IMPORT(cwb2G*,gCWB2G) 		// import cwb2G object pointer
 
      double seg_begin_gps = gCWB2G->GetSegBegin();	// segment begin gps
      double seg_end_gps   = gCWB2G->GetSegEnd();	// segment end   gps

      cout << "Segment Begin : " << seg_begin_gps << endl;
      cout << "Segment End   : " << seg_end_gps << endl;

      double gps = (seg_begin_gps+seg_end_gps)/2.;

      skymap sm = GetRingSkymask(net, cfg, seg_begin_gps, seg_end_gps);

#ifdef SAVE_SKYMASK_PLOT
      // save skymask in geograph coordinates
      TString fname = "plot/ring_skymask.png";
      gskymap gSM(sm);
      gSM.SetOptions("cartesian","geographic",1);
      TString title = Form("%s  gps: %9.0f ", SN2007gr_NAME, gps );
      gSM.Draw();
      double thSN = 90-SN2007gr_DEC;			// cwb coordinates
      double phSN = sm.RA2phi(SN2007gr_RA, gps);	// cwb coordinates
      CwbToGeographic(phSN, thSN, phSN, thSN);
      gSM.DrawMarker(phSN, thSN, 29, 2.0, kBlack);
      gSM.DrawText(phSN+10, thSN-10, SN2007gr_NAME, 0.04, kBlack);
      gSM.SetTitle(title);
      gSM.Print(fname);
      cout << "Write : " << fname << endl;
#endif 

#ifdef APPLY_SKYMASK
      // set ring skymask 
      char skycoord = 'e'; 			// e/c = earth/celestial coordinates 
      net->setSkyMask(sm,skycoord);
      if(skycoord=='e') net->setIndexMode(0);
#endif 
    }
  }
  return;
}

skymap
GetRingSkymask(network* net, config* cfg, double gps_start, double gps_end) {

  // -----------------------------------------------------------------
  // The position of HL detectors in geographical coordinates:    
  // H (46deg 72'18.53"N, 119deg24'27.56"W)
  // L (30deg 33'46.42"N,  90deg46'27.27"W)
  // -----------------------------------------------------------------

  double tgH,pgH;
  double tgL,pgL;
  for(int i=0;i<net->ifoListSize();i++) {
    detector* pD = net->getifo(i);	// get detector object
    //pD->print();			// print detector's parameters
    detectorParams dP = pD->getDetectorParams();
    cout << "detector name : " << dP.name << endl;
    cout << "latitude      : " << dP.latitude << endl;
    cout << "longitude     : " << dP.longitude << endl;
    if(TString(dP.name)=="L1") {tgL=dP.latitude;pgL=dP.longitude;}
    if(TString(dP.name)=="H1") {tgH=dP.latitude;pgH=dP.longitude;}
  }

  // cwb coordinates                    
  double thH =  90 - tgH;                  
  double phH = 360 + pgH;                  
  double thL =  90 - tgL;                  
  double phL = 360 + pgL;                  

  // Unit vectors
  TVector3 vX ( 1,0,0 );
  TVector3 vY ( 0,1,0 );
                        
  // Unit vectors along detectors and unit vectors in cWB coordinates
  TVector3 vH(sin(thH*PI/180) * cos(phH*PI/180), sin(thH*PI/180) * sin(phH*PI/180), cos(thH*PI/180));
  TVector3 vL(sin(thL*PI/180) * cos(phL*PI/180), sin(thL*PI/180) * sin(phL*PI/180), cos(thL*PI/180));
  TVector3 vHL = vH - vL;                                                                               
  // Position of HL detector on the sky:                                                                
  Double_t thHL = vHL.Theta() * 180 / PI;                                                               
  Double_t phHL = vHL.Phi()   * 180 / PI;                                                               
  cout << "The HL baseline postion on the sky is (thHL,phHL) = ( " << thHL << ", " << phHL << " )" << endl;

  // Define a skymask
  skymap* pSM=NULL;
  if(cfg->healpix) pSM = new skymap(int(cfg->healpix));
  else             pSM = new skymap(cfg->angle,cfg->Theta1,cfg->Theta2,cfg->Phi1,cfg->Phi2);
  skymap sm(*pSM); 				
  delete pSM;

  // define the width of the ring
  double phi_start  = sm.RA2phi(SN2007gr_RA, gps_start);
  double phi_end    = sm.RA2phi(SN2007gr_RA, gps_end);
  double phi_width  = phi_start>phi_end ? phi_start-phi_end : 360+phi_start-phi_end;
  double ring_width = phi_width>MIN_RING_WIDTH ? phi_width : MIN_RING_WIDTH;		
  //cout << "RING WIDTH : " << phi_start << " " << phi_end << " " << phi_width << " " << ring_width << endl;

  // Position of SN2007gr on the sky:
  double gps = (gps_start+gps_end)/2.;
  double thSN = 90-SN2007gr_DEC;		// cwb coordinates
  double phSN = sm.RA2phi(SN2007gr_RA, gps);	// cwb coordinates
  cout << "UTC of gps start: " << wat::Time(gps).GetDateString() << endl;
  TVector3 vSN(sin(thSN*PI/180)*cos(phSN*PI/180), 
               sin(thSN*PI/180)*sin(phSN*PI/180), cos(thSN*PI/180));
  Double_t angSN = vHL.Angle(vSN)*180/PI;  
  cout << "(thSN, phSN) = (" << thSN << ", " << phSN << ")" << " "<< "angSN = " << angSN <<endl;

  int L = sm.size();        			// skymap size
  for (int l=0;l<L;l++) {
    // get theta,phi in CWB coordinates
    double phi = sm.getPhi(l); 
    double theta = sm.getTheta(l); 
    TVector3 vpx(sin(theta*PI/180)*cos(phi*PI/180), 
                 sin(theta*PI/180)*sin(phi*PI/180), cos(theta*PI/180));
    double angpx = vHL.Angle(vpx)*180/PI;  

    // set 1 in the sky ring
    if((angpx < angSN+ring_width/2)&&(angpx > angSN-ring_width/2)) sm.set(l,1);
    else                                                           sm.set(l,0);
  } 

  return sm;
}

