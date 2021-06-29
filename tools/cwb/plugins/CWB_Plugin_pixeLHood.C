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
#include "mdc.hh"
#include <vector>
#include "watplot.hh"
#include <stdio.h>

#define OUTPUT_DIR 	"plots"
#define OUTPUT_TYPE	"png"

void Dump(WSeries<double> &w, double t1, double t2, const char* fname);

void 
CWB_Plugin(TFile* jfile, CWB::config* cfg, network* NET, WSeries<double>* x, TString ifo, int type)  {
//!EVENT_REPORT
// This plugin dump/plot the likelihood/null pixels of the detected/reconstructed event 

  cout << endl;
  cout << "-----> plugins/CWB_Plugin_pixeLHood.C" << endl;
  cout << "ifo " << ifo.Data() << endl;
  cout << "type " << type << endl;
  cout << endl;


  if(type==CWB_PLUGIN_OLIKELIHOOD) {  

    int i,j,k,n,m,l;
    int nIFO = NET->ifoListSize();
    int K = NET->nLag;            
    int M = NET->mdc__ID.size();  
    int ID;                       
    char search = NET->tYPe;                                   
  
    wavearray<double> id;
  
    bool batch = gROOT->IsBatch();
    gROOT->SetBatch(true);        
  
    watplot WTS(const_cast<char*>("wts"));
  
    char outdir[500] = OUTPUT_DIR;
    char ifostr[20]  = ""; 
    char strtime[1024]  = ""; 
    char fname[1024];   
    char gtype[32]   = OUTPUT_TYPE;
  
    int    minTimeDet=nIFO;
    double minTime=1e40;   
    double eventTime[NIFO];
    double lagTime[NIFO];  
    int ifoid[NIFO],sortid[NIFO];   
    double factor = 1;
 
    netevent* EVT = new netevent(nIFO);
 
    int rate = int(2*NET->getifo(0)->TFmap.resolution(0)+0.5); 
  
    for(k=0; k<K; k++) {  				// loop over the lags
  
      id = NET->getwc(k)->get(const_cast<char*>("ID"), 0, 'L', rate);
  
      for(j=0; j<(int)id.size(); j++) {  		// loop over cluster index
  
        ID = size_t(id.data[j]+0.5);			// cluster id
  
        EVT->output(NULL,NET,factor,ID,k); 		// get reconstructed parameters
        if(EVT->rho[1] < cfg->cedRHO) continue;    	// skip events under rho threshold
  
        int masterDet=0;
        int lagMin=2147483647;
        for(n=0; n<nIFO; n++) if(EVT->lag[n]<lagMin) {lagMin=int(EVT->lag[n]);masterDet=n;}
  
        NET->likelihood(search, NET->acor, ID, k);	// exec likelihood search
 
        // get event network time
        double gps_start = EVT->time[masterDet]-EVT->duration[masterDet];
        double gps_stop  = EVT->time[masterDet]+EVT->duration[masterDet];

        // create a unique label
        for(n=0; n<nIFO; n++) sprintf(strtime, "%s_%.3f",strtime,EVT->start[n]); 

        // dump event parameters
        sprintf(fname, "%s/eventDump%s.%s", outdir, strtime, "txt");
        EVT->dopen(fname,const_cast<char*>("w"));
        EVT->output(NULL,NET,factor,ID,k);
        EVT->dclose();
 
        // compute event time & lags time
        for(n=0; n<nIFO; n++) eventTime[n]=(EVT->start[n]+EVT->stop[n])/2.;
        for(n=0; n<nIFO; n++) minTime = (eventTime[n]<minTime) ? eventTime[n] : minTime;
        for(n=0; n<nIFO; n++) lagTime[n]=eventTime[n]-minTime;                          
        for(n=0; n<nIFO; n++) minTimeDet = (eventTime[n]<=minTime) ? n : minTimeDet;    

        // set likelihood maps to optimal resolution level
        double optRate  = (NET->getwc(k)->cRate[ID-1])[0];
        double optLayer = NET->getwc(k)->rate/optRate;
        double optLevel = int(log10(optLayer)/log10(2));
        for(int n=0;n<2;n++) {
          WSeries<double>* pTF;
          if(n==0) pTF = &NET->pixeLHood;
          if(n==1) pTF = &NET->pixeLNull;
          pTF->Inverse();				// transform to time domain
          pTF->Forward(optLevel);			// transform to optimal level
        }

        // plot likelihood map at optimal resolution level
        WTS.canvas->cd();                                 
        sprintf(fname, "%s/l_tfmap_scalogram%s.%s", outdir, strtime, gtype);
        WTS.plot(NET->pixeLHood, 0, gps_start-EVT->slag[masterDet],
                 gps_stop-EVT->slag[masterDet], const_cast<char*>("COLZ"));
        WTS.hist2D->GetYaxis()->SetRangeUser(NET->pixeLHood.getlow(),NET->pixeLHood.gethigh());
        WTS.hist2D->SetTitle("Scalogram");
        WTS.canvas->Print(fname); 
        WTS.clear();
  
        // plot null map at optimal resolution level
        sprintf(fname, "%s/n_tfmap_scalogram%s.%s", outdir, strtime, gtype);
        WTS.plot(NET->pixeLNull, 0, gps_start-EVT->slag[masterDet],
                 gps_stop-EVT->slag[masterDet], const_cast<char*>("COLZ"));
        WTS.hist2D->GetYaxis()->SetRangeUser(NET->pixeLNull.getlow(),NET->pixeLNull.gethigh());
        WTS.hist2D->SetTitle("Scalogram");
        WTS.canvas->Print(fname); 
        WTS.clear();

        // dump likelihood map at optimal resolution level
        sprintf(fname, "%s/l_tfmap_scalogram%s.%s", outdir, strtime, "txt");
        Dump(NET->pixeLHood, gps_start-EVT->slag[masterDet], gps_stop-EVT->slag[masterDet], fname);
        // dump null map at optimal resolution level
        sprintf(fname, "%s/n_tfmap_scalogram%s.%s", outdir, strtime, "txt");
        Dump(NET->pixeLNull, gps_start-EVT->slag[masterDet], gps_stop-EVT->slag[masterDet], fname);
  
      } 						// End loop on found events
    }   						// End loop on lags
  
    gROOT->SetBatch(batch);  // restore batch status

    delete EVT;  
  }

  return;
}

void Dump(WSeries<double> &w, double t1, double t2, const char* fname) {

  int ni = 1<<w.pWavelet->m_Level;
  int nb = int((t1-w.start())*w.rate()/ni);
  int nj = int((t2-t1)*w.rate())/ni;
  int ne = nb+nj;
  double rate=w.rate();
  double df=rate/2/ni;
  rate = rate/ni;
  double dt = 1./rate;

  wavearray<double> wl;

  FILE *fp;
  if ( (fp = fopen(fname, "w")) == NULL ) {
     cout << " Dump() error: cannot open file " << fname <<". \n";
     return;
  };

  for(int i=0;i<ni;i++) {
    w.getLayer(wl,i);
    for(int j=nb;j<=ne;j++) {
      float t = (j+0.5)*dt;
      float f = (i+0.5)*df;
      float x = wl.data[j];
      if(x>0.1) fprintf( fp,"%-3.6f   \t%-3.6f   \t%-3.3f   \t%-3.3f   \t%-3.3f\n", t, dt, f, df, x);
    }
  }

  fclose(fp);

  return;
}
