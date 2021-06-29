/*
# Copyright (C) 2019 Gabriele Vedovato, Sergey Klimenko
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


#include "ced.hh"
#include "watplot.hh"
#include "wseries.hh"

#include "THtml.h"
#include "TDocOutput.h"
#include "TH2.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TSystem.h"
#include "TMath.h"
#include "TMarker.h"
#include "TVector3.h"
#include "TRotation.h"
#include "TPolyLine.h"
#include "TMacro.h"
#include "Math/Rotation3D.h"
#include "Math/Vector3Dfwd.h"
#include "TLegend.h"

#include "Meyer.hh"  
#include "gskymap.hh"
#include "gnetwork.hh"
#include "STFT.hh"
#include <string.h>

#define REPLACE(STRING,DIR,EXT) TString(STRING).ReplaceAll("/","").ReplaceAll(DIR,"").ReplaceAll("."+TString(EXT),"")

using namespace ROOT::Math;

ClassImp(CWB::ced)	// used by THtml doc

void
CWB::ced::Write(double factor, size_t iID, int LAG, char* dirCED) {

  if(!NET || !NET->nLag) exit(1);

  int i,j,k,ind;
  int m,K,M,injID;
  detector* pd;
  netcluster* p;
  int nIFO = (int)NET->ifoListSize();   // number of detectors
  int bLag = LAG<0 ? 0 : LAG;
  int eLag = LAG<0 ? NET->nLag : LAG+1;
  wavecomplex Aa, gC;
  factor = factor<=0 ? 1 : fabs(factor);       // used to rescale injected events

  if(!nIFO) return;

  // check if it is one detector network  
  if(nIFO==2) if(strcmp(NET->ifoName[0],NET->ifoName[1])==0) nIFO=1;   

  watplot WTS(const_cast<char*>("wtswrc"));  
  watplot PTS(const_cast<char*>("ptswrc"),200,20,800,500);  
  gskymap gSM;

// injections

  injection INJ(nIFO);
  size_t mdcID, ID;
  double T, injTime;
  skymap* psm;
  std::vector<float>* vP;
  std::vector<int>*   vI;

  bool ellips = NET->tYPe=='i' || NET->tYPe=='I' || 
                NET->tYPe=='g' || NET->tYPe=='G' || 
                NET->tYPe=='s' || NET->tYPe=='S' ||
                NET->tYPe=='l' || NET->tYPe=='L' ||
                NET->tYPe=='c' || NET->tYPe=='C' ||
                NET->tYPe=='e' || NET->tYPe=='E' ||
                NET->tYPe=='p' || NET->tYPe=='P' ||
                NET->tYPe=='r' || NET->tYPe=='R';

  bool burst  = NET->tYPe=='b' || NET->tYPe=='B';  

  wavearray<double> skSNR(nIFO);
  wavearray<double> xkSNR(nIFO);

// arrays for cluster parameters
      
  wavearray<double> clusterID_net;

  EVT->run = NET->nRun;
  EVT->nevent = 0;

  for(int lag=bLag; lag<eLag; lag++){  		// loop on time lags

    p = NET->getwc(lag);               		// pointer to netcluster
    if(!p->size()) continue;

    if(NET->MRA) clusterID_net = p->get((char*)"ID",0,'S',0);
    else         clusterID_net = p->get((char*)"ID",0,'S',1);
    K = clusterID_net.size();

    if(!K) continue;
            
// read cluster parameters

    for(k=0; k<K; k++)                  	// loop on events
    {
      ID = size_t(clusterID_net.data[k]+0.1);
      if((iID)&&(ID!=iID)) continue; 
      if(ellips) { 
        if(NET->like()!='2') NET->SETNDM(ID,lag,true,NET->MRA ? 0 : 1);
      }
      else if(burst)
        NET->setndm(ID,lag,true,1);
      else
        {cout<<"netevent::output(): incorrect search option"<<endl; exit(1);}

      psm = &(NET->getifo(0)->tau);
      vI = &(NET->wc_List[lag].p_Ind[ID-1]);
      ind = (*vI)[0];               		// reconstructed sky index

// set injections if there are any

      M = NET->mdc__IDSize();
      if(!lag) {                                // only for zero lag
        injTime = 1.e12;
	injID   = -1;
	for(m=0; m<M; m++) {
	  mdcID = NET->getmdc__ID(m);
	  T = fabs(EVT->time[0] - NET->getmdcTime(mdcID));
     	  if(T<injTime && INJ.fill_in(NET,mdcID)) { 
	    injTime = T; 
	    injID = mdcID; 
	  } 
//	  printf("%d  %12.3f  %12.3f\n",mdcID,NET->getmdcTime(mdcID),T);
        }
	
	if(INJ.fill_in(NET,injID)) {            // set injections
 
//	printf("injection type: %d  %12.3f\n",INJ.type,INJ.time[0]);

          wavearray<double>** pwfINJ = new wavearray<double>*[nIFO];  
          wavearray<double>** pwfREC = new wavearray<double>*[nIFO];  
          pd = NET->getifo(0);  
          int idSize = pd->RWFID.size();  
          int wfIndex=-1;  
          for (int mm=0; mm<idSize; mm++) if (pd->RWFID[mm]==(int)ID) wfIndex=mm;  

	  for(j=0; j<int(nIFO); j++) {
	    pd = NET->getifo(j);
	    Aa = pd->antenna(EVT->theta[1],EVT->phi[1]);    // inj antenna pattern
	    EVT->hrss[j+nIFO] = INJ.hrss[j]*factor;
	    EVT->bp[j+nIFO]   = Aa.real();
	    EVT->bx[j+nIFO]   = Aa.imag();
	    EVT->time[j+nIFO] = INJ.time[j];
	    EVT->Deff[j]      = INJ.Deff[j]/factor;
            TString xtitle    = TString::Format("Time (sec) : GPS OFFSET = %.3f",EVT->gps[j]);
               
            pwfINJ[j]       = INJ.pwf[j];
            if (pwfINJ[j]==NULL) {
              cout << "Error : Injected waveform not saved !!! : detector " 
                   << NET->ifoName[j] << endl;
              continue;
            }
            if (wfIndex<0) {
              cout << "Error : Reconstructed waveform not saved !!! : ID -> " 
                   << ID << " : detector " << NET->ifoName[j] << endl;
              continue;
            }
            if (wfIndex>=0) pwfREC[j] = pd->RWFP[wfIndex];
            double R = pd->TFmap.rate();
            double rFactor = 1.;
            rFactor *= factor;
            wavearray<double>* wfINJ = pwfINJ[j];
            *wfINJ*=rFactor;
            wavearray<double>* wfREC = pwfREC[j];

            // rescale data when data are resampled (resample with wavelet change the amplitude)
            double rescale = 1./pow(sqrt(2.),TMath::Log2(inRate/wfINJ->rate())); 
            *wfINJ*=rescale; *wfREC*=rescale; 	// rescale data 

            double bINJ = wfINJ->start();
            double eINJ = wfINJ->start()+wfINJ->size()/R;
            double bREC = wfREC->start();
            double eREC = wfREC->start()+wfREC->size()/R;
            //cout.precision(14);
            //cout << "bINJ : " << bINJ << " eINJ : " << eINJ << endl;
            //cout << "bREC : " << bREC << " eREC : " << eREC << endl;

            int oINJ = bINJ>bREC ? 0 : int((bREC-bINJ)*R+0.5);
            int oREC = bINJ<bREC ? 0 : int((bINJ-bREC)*R+0.5);
            //cout << "oINJ : " << oINJ << " oREC : " << oREC << endl;

            double startXCOR = bINJ>bREC ? bINJ : bREC;
            double endXCOR   = eINJ<eREC ? eINJ : eREC;
            int sizeXCOR   = int((endXCOR-startXCOR)*R+0.5);
            //cout << "startXCOR : " << startXCOR << " endXCOR : " << endXCOR << " sizeXCOR :" << sizeXCOR << endl;

            if (sizeXCOR<=0) continue;   

            // the enINJ, enREC, xcorINJ_REC are computed in the INJ range

            double enINJ=0;
            for (int i=0;i<(int)wfINJ->size();i++) enINJ+=wfINJ->data[i]*wfINJ->data[i];
            //for (int i=0;i<sizeXCOR;i++) enINJ+=wfINJ->data[i+oINJ]*wfINJ->data[i+oINJ];

            double enREC=0;
            //for (int i=0;i<wfREC->size();i++) enREC+=wfREC->data[i]*wfREC->data[i];
            for (int i=0;i<sizeXCOR;i++) enREC+=wfREC->data[i+oREC]*wfREC->data[i+oREC];

            double xcorINJ_REC=0;
            for (int i=0;i<sizeXCOR;i++) xcorINJ_REC+=wfINJ->data[i+oINJ]*wfREC->data[i+oREC];

            WSeries<double> wfREC_SUB_INJ;  
            wfREC_SUB_INJ.resize(sizeXCOR);
            for (int i=0;i<sizeXCOR;i++) wfREC_SUB_INJ.data[i]=wfREC->data[i+oREC]-wfINJ->data[i+oINJ];
            wfREC_SUB_INJ.start(startXCOR);
            wfREC_SUB_INJ.rate(wfREC->rate());


            EVT->iSNR[j]    = enINJ;
            EVT->oSNR[j]    = enREC;
            EVT->ioSNR[j]   = xcorINJ_REC;

            //double erINJ_REC = enINJ+enREC-2*xcorINJ_REC;
            //cout << "enINJ : " << enINJ << " enREC : " << enREC << " xcorINJ_REC : " << xcorINJ_REC << endl;
            //cout << "erINJ_REC/enINJ : " << erINJ_REC/enINJ << endl;

            bool batch = gROOT->IsBatch();
            gROOT->SetBatch(true);

            // find TF maps to optimal resolution level
            double optRate  = (p->cRate[ID-1])[0];
            double optLayer = p->rate/optRate;
            double optLevel = int(log10(optLayer)/log10(2));

            // DUMP WRC !!!
            if (wfIndex>=0) {
              double gmin,gmax;
              char fname[1024];

              PTS.canvas->cd();

              sprintf(fname, "%s/%s_wf_white_inj.dat", dirCED, NET->ifoName[j]);
              if(sbasedirCED!=NULL) wfINJ->wavearray<double>::Dump(const_cast<char*>(fname),2);		// time, amp format
//              if(sbasedirCED!=NULL) wfINJ->Dump(fname);
              //cout << "Dump " << fname << endl;

              sprintf(fname, "%s/%s_wf_white_rec.dat", dirCED, NET->ifoName[j]);
              if(sbasedirCED!=NULL) wfREC->wavearray<double>::Dump(const_cast<char*>(fname),2);		// time, amp format
//              if(sbasedirCED!=NULL) wfREC->Dump(fname);
              //cout << "Dump " << fname << endl;

              sprintf(fname, "%s/%s_wf_white_inj_rec.png", dirCED, NET->ifoName[j]);
              //cout << "Print " << fname << endl;
              wfINJ->start(wfINJ->start()-EVT->gps[0]);
              wfREC->start(wfREC->start()-EVT->gps[0]);
              double bT,eT;
              GetBoundaries((wavearray<double>&)*wfINJ, 0.995, bT, eT);
              PTS.plot(wfINJ, const_cast<char*>("ALP"), 1, bT, eT);
              //         startXCOR-EVT->gps[0], endXCOR-EVT->gps[0]);
              PTS.graph[0]->SetLineWidth(1);
              PTS.graph[0]->GetXaxis()->SetTitle(xtitle);  
              if(bT<(startXCOR-EVT->gps[0])) bT=startXCOR-EVT->gps[0];
              if(eT>(endXCOR-EVT->gps[0]))   eT=endXCOR-EVT->gps[0];
              PTS.plot(wfREC, const_cast<char*>("SAME"), 2, bT, eT);
              //         startXCOR-EVT->gps[0], endXCOR-EVT->gps[0]);
              PTS.graph[1]->SetLineWidth(2);
              gmin = TMath::Min(PTS.graph[0]->GetHistogram()->GetMinimum(),PTS.graph[1]->GetHistogram()->GetMinimum());
              gmax = TMath::Max(PTS.graph[0]->GetHistogram()->GetMaximum(),PTS.graph[1]->GetHistogram()->GetMaximum());
              PTS.graph[0]->GetYaxis()->SetRangeUser(0.9*gmin,1.1*gmax);
              wfINJ->start(wfINJ->start()+EVT->gps[0]);
              wfREC->start(wfREC->start()+EVT->gps[0]);
              if(sbasedirCED!=NULL) PTS.canvas->Print(fname); else PTS.canvas->Write(REPLACE(fname,dirCED,gtype));
              //sprintf(fname, "%s/%s_wf_white_inj_rec.root", dirCED, NET->ifoName[j]);
              //PTS.canvas->Print(fname);
              PTS.clear();

              if((pd->TFmap.gethigh()-pd->TFmap.getlow())>200) PTS.canvas->SetLogx(true);  
              PTS.canvas->SetLogy(true);
              sprintf(fname, "%s/%s_wf_white_inj_rec_fft.png", dirCED, NET->ifoName[j]);
              //cout << "Print " << fname << endl;
              wfINJ->start(wfINJ->start()-EVT->gps[0]);
              wfREC->start(wfREC->start()-EVT->gps[0]);
              PTS.plot(wfINJ, const_cast<char*>("ALP"), 1, 
                       startXCOR-EVT->gps[0], endXCOR-EVT->gps[0],
                       true, pd->TFmap.getlow(), pd->TFmap.gethigh());
              PTS.graph[0]->SetLineWidth(1);
              PTS.plot(wfREC, const_cast<char*>("SAME"), 2,
                       startXCOR-EVT->gps[0], endXCOR-EVT->gps[0],
                       true, pd->TFmap.getlow(), pd->TFmap.gethigh());
              PTS.graph[1]->SetLineWidth(2);
              gmin = TMath::Min(PTS.graph[0]->GetHistogram()->GetMinimum(),PTS.graph[1]->GetHistogram()->GetMinimum());
              gmax = TMath::Max(PTS.graph[0]->GetHistogram()->GetMaximum(),PTS.graph[1]->GetHistogram()->GetMaximum());
              PTS.graph[0]->GetYaxis()->SetRangeUser(0.9*gmin,1.1*gmax);
              wfINJ->start(wfINJ->start()+EVT->gps[0]);
              wfREC->start(wfREC->start()+EVT->gps[0]);
              if(sbasedirCED!=NULL) PTS.canvas->Print(fname); else PTS.canvas->Write(REPLACE(fname,dirCED,gtype));
              PTS.clear();
              PTS.canvas->SetLogx(false);
              PTS.canvas->SetLogy(false);
              sprintf(fname, "%s/%s_wf_white_rec_sub_inj.png", dirCED, NET->ifoName[j]);
              //cout << "Print " << fname << endl;
              wfREC_SUB_INJ.start(wfREC_SUB_INJ.start()-EVT->gps[0]);
              wfREC->start(wfREC->start()-EVT->gps[0]);
              PTS.plot(wfREC, const_cast<char*>("ALP"), 2, 
                       startXCOR-EVT->gps[0], endXCOR-EVT->gps[0]);
              PTS.graph[0]->SetLineWidth(1);
              PTS.graph[0]->GetXaxis()->SetTitle(xtitle);  
              PTS.plot(wfREC_SUB_INJ, const_cast<char*>("SAME"), 1,
                       startXCOR-EVT->gps[0], endXCOR-EVT->gps[0]);
              PTS.graph[1]->SetLineWidth(2);
              gmin = TMath::Min(PTS.graph[0]->GetHistogram()->GetMinimum(),PTS.graph[1]->GetHistogram()->GetMinimum());
              gmax = TMath::Max(PTS.graph[0]->GetHistogram()->GetMaximum(),PTS.graph[1]->GetHistogram()->GetMaximum());
              PTS.graph[0]->GetYaxis()->SetRangeUser(0.9*gmin,1.1*gmax);
              wfREC_SUB_INJ.start(wfREC_SUB_INJ.start()+EVT->gps[0]);
              wfREC->start(wfREC->start()+EVT->gps[0]);
              if(sbasedirCED!=NULL) PTS.canvas->Print(fname); else PTS.canvas->Write(REPLACE(fname,dirCED,gtype));
              //sprintf(fname, "%s/%s_wf_white_rec_sub_inj.root", dirCED, NET->ifoName[j]);
              //PTS.canvas->Print(fname);
              PTS.clear();

              if((pd->TFmap.gethigh()-pd->TFmap.getlow())>200) PTS.canvas->SetLogx(true);  
              PTS.canvas->SetLogy(true);
              sprintf(fname, "%s/%s_wf_white_rec_sub_inj_fft.png", dirCED, NET->ifoName[j]);
              //cout << "Print " << fname << endl;
              wfREC_SUB_INJ.start(wfREC_SUB_INJ.start()-EVT->gps[0]);
              wfREC->start(wfREC->start()-EVT->gps[0]);
              PTS.plot(wfREC, const_cast<char*>("ALP"), 2, 
                       startXCOR-EVT->gps[0], endXCOR-EVT->gps[0],
                       true, pd->TFmap.getlow(), pd->TFmap.gethigh());
              PTS.graph[0]->SetLineWidth(1);
              PTS.plot(wfREC_SUB_INJ, const_cast<char*>("SAME"), 1,
                       startXCOR-EVT->gps[0], endXCOR-EVT->gps[0],
                       true, pd->TFmap.getlow(), pd->TFmap.gethigh());
              PTS.graph[1]->SetLineWidth(2);
              gmin = TMath::Min(PTS.graph[0]->GetHistogram()->GetMinimum(),PTS.graph[1]->GetHistogram()->GetMinimum());
              gmax = TMath::Max(PTS.graph[0]->GetHistogram()->GetMaximum(),PTS.graph[1]->GetHistogram()->GetMaximum());
              PTS.graph[0]->GetYaxis()->SetRangeUser(0.9*gmin,1.1*gmax);
              wfREC_SUB_INJ.start(wfREC_SUB_INJ.start()+EVT->gps[0]);
              wfREC->start(wfREC->start()+EVT->gps[0]);
              if(sbasedirCED!=NULL) PTS.canvas->Print(fname); else PTS.canvas->Write(REPLACE(fname,dirCED,gtype));
              PTS.clear();
              PTS.canvas->SetLogx(false);
              PTS.canvas->SetLogy(false);

              WTS.canvas->cd();
              //if((pd->TFmap.gethigh()-pd->TFmap.getlow())>200) WTS.canvas->SetLogy(true); 
              // length of temporary buffer for tf plots
              int xcor_length = sizeXCOR/wfINJ->rate();
	      int wts_size = xcor_length<8 ? 16 : 16*int(1+xcor_length/8.);
              wts_size*=wfINJ->rate();
/*
	      // find optimal level
              int oRate = int(wfINJ->rate()/EVT->rate[0]);
              int oLevel=0;
              while(oRate>1) {oRate/=2;oLevel++;}
              //cout << "EVT->rate[0] : " << EVT->rate[0] << " oLevel : " << oLevel << endl;
*/
              Meyer<double> S(1024,2);     // set wavelet for production

              wavearray<double> xINJ(wts_size);
              xINJ.start(wfINJ->start()-EVT->gps[0]+double(oINJ+wts_size/2)/wfINJ->rate());
              xINJ.rate(wfINJ->rate());
              xINJ=0.;
              for (int m=0;m<sizeXCOR;m++) 
	        if(m<(int)xINJ.size()/2) xINJ.data[m+wts_size/2]=wfINJ->data[m+oINJ];
              WSeries<double> wINJ(xINJ,S);
              xINJ.resize(0); 

	      double wts_start,wts_stop;

              double  rescale = 2.*pow(sqrt(2.),TMath::Log2(inRate/wINJ.rate())); 
              double  wts_offset = wfINJ->start()-double(oINJ+wts_size/2)/wfINJ->rate();
              TString xtitle = TString::Format("Time (sec) : GPS OFFSET = %.3f",wts_offset);

	      //scalogram maps
              sprintf(fname, "%s/%s_wf_white_inj_tf.png", dirCED, NET->ifoName[j]);
              //cout << "Print " << fname << endl;
//              wINJ.Forward(oLevel); 
              if(NET->wdm()) {if(NET->getwdm(optLayer+1)) wINJ.Forward(wINJ,*NET->getwdm(optLayer+1));} 
              else           wINJ.Forward(optLevel);
	      wts_start = wINJ.start()+(double)(wts_size/2)/wINJ.rate();
	      wts_stop  = sizeXCOR<wts_size/2 ? wts_start+sizeXCOR/wINJ.rate() : wts_start+(wts_size/2)/wINJ.rate();
              WTS.plot(&wINJ, 1, wts_start, wts_stop,const_cast<char*>("COLZ"));
              WTS.hist2D->GetYaxis()->SetRangeUser(pd->TFmap.getlow(), pd->TFmap.gethigh());
              WTS.hist2D->Scale(2*rescale);
//              WTS.hist2D->GetZaxis()->SetRangeUser(0.1, WTS.hist2D->GetMaximum());	// set to white color low Z values
              WTS.hist2D->GetXaxis()->SetTitle(xtitle);  
              WTS.hist2D->GetXaxis()->CenterTitle(false);
              if(sbasedirCED!=NULL) WTS.canvas->Print(fname); else WTS.canvas->Write(REPLACE(fname,dirCED,gtype));
              WTS.clear();

              wavearray<double> xREC(wts_size);
              xREC.start(wfREC->start()-EVT->gps[0]+double(oREC+wts_size/2)/wfREC->rate());
              xREC.rate(wfREC->rate());
              xREC=0.;
              for (int m=0;m<sizeXCOR;m++) 
		if(m<(int)xREC.size()/2) xREC.data[m+wts_size/2]=wfREC->data[m+oREC];
              WSeries<double> wREC(xREC,S);
              xREC.resize(0); 

	      //scalogram maps
              sprintf(fname, "%s/%s_wf_white_rec_tf.png", dirCED, NET->ifoName[j]);
              //cout << "Print " << fname << endl;
//              wREC.Forward(oLevel); 
              if(NET->wdm()) {if(NET->getwdm(optLayer+1)) wREC.Forward(wREC,*NET->getwdm(optLayer+1));} 
              else           wREC.Forward(optLevel);
	      wts_start = wREC.start()+(double)(wts_size/2)/wREC.rate();
	      wts_stop  = sizeXCOR<wts_size/2 ? wts_start+sizeXCOR/wREC.rate() : wts_start+(wts_size/2)/wREC.rate();
              WTS.plot(&wREC, 1, wts_start, wts_stop,const_cast<char*>("COLZ"));
              WTS.hist2D->GetYaxis()->SetRangeUser(pd->TFmap.getlow(), pd->TFmap.gethigh());
              WTS.hist2D->Scale(2*rescale);
//              WTS.hist2D->GetZaxis()->SetRangeUser(0.1, WTS.hist2D->GetMaximum());	// set to white color low Z values
              WTS.hist2D->GetXaxis()->SetTitle(xtitle);  
              WTS.hist2D->GetXaxis()->CenterTitle(false);
              if(sbasedirCED!=NULL) WTS.canvas->Print(fname); else WTS.canvas->Write(REPLACE(fname,dirCED,gtype));
              WTS.clear();

              wavearray<double> xDIF(wts_size);
              xDIF.start(wfREC->start()-EVT->gps[0]+double(oREC+wts_size/2)/wfREC->rate());
              xDIF.rate(wfREC->rate());
              xDIF=0.;
              for (int m=0;m<sizeXCOR;m++) 
                if(m<(int)xDIF.size()/2) xDIF.data[m+wts_size/2]=wfREC->data[m+oREC]-wfINJ->data[m+oINJ];
              WSeries<double> wDIF(xDIF,S);
              xDIF.resize(0); 

	      //scalogram maps
              sprintf(fname, "%s/%s_wf_white_dif_tf.png", dirCED, NET->ifoName[j]);
              //cout << "Print " << fname << endl;
//              wDIF.Forward(oLevel); 
              if(NET->wdm()) {if(NET->getwdm(optLayer+1)) wDIF.Forward(wDIF,*NET->getwdm(optLayer+1));} 
              else           wDIF.Forward(optLevel);
	      wts_start = wDIF.start()+(double)(wts_size/2)/wDIF.rate();
	      wts_stop  = sizeXCOR<wts_size/2 ? wts_start+sizeXCOR/wDIF.rate() : wts_start+(wts_size/2)/wDIF.rate();
              WTS.plot(&wDIF, 2, wts_start, wts_stop,const_cast<char*>("COLZ"));
              WTS.hist2D->GetYaxis()->SetRangeUser(pd->TFmap.getlow(), pd->TFmap.gethigh());
              if(sbasedirCED!=NULL) WTS.canvas->Print(fname); else WTS.canvas->Write(REPLACE(fname,dirCED,gtype));
              WTS.clear();

              double t1 = wDIF.start();
              double t2 = wDIF.start()+wDIF.size()/wDIF.rate();

              int ni = 1<<wDIF.pWavelet->m_Level;
              int nb = int((t1-wDIF.start())*wDIF.rate()/ni);
              int nj = int((t2-t1)*wDIF.rate())/ni;
              int ne = nb+nj;

              int nsize=0;
              double eeINJ=0;
              double eeREC=0;
              double eeDIF=0;
              wavearray<double> wlINJ;
              wavearray<double> wlREC;
              wavearray<double> wlDIF;
              for(int m=0;m<ni;m++) {
                wINJ.getLayer(wlINJ,m);
                wREC.getLayer(wlREC,m);
                wDIF.getLayer(wlDIF,m);
                for(int k=nb;k<ne;k++) {
                  if(fabs(wlREC.data[k])>0.1) {
                    eeINJ+=wlINJ.data[k]*wlINJ.data[k];
                    eeREC+=wlREC.data[k]*wlREC.data[k];
                    eeDIF+=wlDIF.data[k]*wlDIF.data[k];
                    nsize++;
                  }
                }
              }

              //cout << j << " tfNorm : " << eeINJ << " " << " " << eeREC << " " << eeDIF/eeINJ << " " << nsize << endl;

              wINJ.resize(0); 
              wREC.resize(0); 
              wDIF.resize(0); 

              WTS.canvas->SetLogy(false);
              gROOT->SetBatch(batch);
            }
            *wfINJ*=1./rFactor;
            *wfINJ*=1./rescale; *wfREC*=1./rescale; 	// restore original rescaled data 
	  }
        }
      }
/*
      if(EVT->fP!=NULL) {
        fprintf(EVT->fP,"\n# trigger %d in lag %d for \n",int(ID),int(n));
        EVT->Dump();
        vP = &(NET->wc_List[n].p_Map[ID-1]);
        vI = &(NET->wc_List[n].p_Ind[ID-1]);
        x = cos(psm->theta_1*PI/180.)-cos(psm->theta_2*PI/180.);
        x*= (psm->phi_2-psm->phi_1)*180/PI/psm->size();
        fprintf(EVT->fP,"sky_res:    %f\n",sqrt(fabs(x)));
        fprintf(EVT->fP,"map_lenght: %d\n",int(vP->size()));
        fprintf(EVT->fP,"#skyID  theta   DEC     step   phi     R.A    step  probability    cumulative\n");
        x = 0.;
        for(j=0; j<int(vP->size()); j++) {
          i = (*vI)[j];
          x+= (*vP)[j];
          fprintf(EVT->fP,"%6d  %5.1f  %5.1f  %6.2f  %5.1f  %5.1f  %6.2f  %e  %e\n",
                  int(i),psm->getTheta(i),psm->getDEC(i),psm->getThetaStep(i),
                  psm->getPhi(i),psm->getRA(i),psm->getPhiStep(i),(*vP)[j],x);
        }
      }
*/
      // dump to fits the probability skymap (only for healpix skymap)
      if(nIFO>1) {  
        bool batch = gROOT->IsBatch();
        gROOT->SetBatch(true);

        // Dump2fits probability skymap  (healpix)
        skymap skyprobcc = *psm;
        skyprobcc=0.;
        skymap skyprob = *psm;
        skyprob=1.e-12;

        vP = &(NET->wc_List[lag].p_Map[ID-1]);
        vI = &(NET->wc_List[lag].p_Ind[ID-1]);
        double th,ph,ra;
        int k;
	for(j=0; j<int(vP->size()); j++) {
	  i = (*vI)[j];
          th = skyprob.getTheta(i);
          ph = skyprob.getPhi(i);

          k=skyprob.getSkyIndex(th, ph);
          skyprob.set(k,(*vP)[j]);

          ra = skyprob.getRA(i);
          k=skyprob.getSkyIndex(th, ra);
          skyprobcc.set(k,(*vP)[j]);
        }

        char fname[1024];
#ifdef _USE_HEALPIX
        if(skyprobcc.getType()&&(sbasedirCED!=NULL)) { // healpix  in celestial coordinates
          sprintf(fname, "%s/skyprobcc.%s", dirCED, "fits");

          // build fits configur info
          TString analysis = "1G";
          if(NET->like()=='2') analysis="2G";
          if(NET->MRA) analysis+=":MRA";
          if(NET->pattern==0)                     analysis+=":Packet(0)";
          if((NET->pattern!=0 && NET->pattern<0)) analysis+=TString::Format(":Packet(%d)",NET->pattern);
          if((NET->pattern!=0 && NET->pattern>0)) analysis+=TString::Format(":Packet(+%d)",NET->pattern);

          char configur[64]="";
          char search = NET->tYPe;
          if (search=='r')                   sprintf(configur,"%s un-modeled",analysis.Data());
          if(analysis=="1G") {
            if((search=='i')||(search=='I')) sprintf(configur,"%s elliptical",analysis.Data());
            if((search=='s')||(search=='S')) sprintf(configur,"%s linear",analysis.Data());
            if((search=='g')||(search=='G')) sprintf(configur,"%s circular",analysis.Data());
          } else {
            if (search=='i')                 sprintf(configur,"%s iota-wave",analysis.Data());
            if (search=='p')                 sprintf(configur,"%s psi-wave",analysis.Data());
            if((search=='l')||(search=='s')) sprintf(configur,"%s linear",analysis.Data());
            if((search=='c')||(search=='g')) sprintf(configur,"%s circular",analysis.Data());
            if((search=='e')||(search=='b')) sprintf(configur,"%s elliptical",analysis.Data());
          }
          skyprobcc.Dump2fits(fname,EVT->time[0],configur,const_cast<char*>("PROB"),const_cast<char*>("pix-1"),'C'); 
        }
#endif

        // dump skyprob plot in cc coordinates
        gSM=skyprobcc;
        gSM.SetOptions("hammer","Celestial",2);

        int L     = gSM.size();               // number of pixels in the sphere
        double pi = TMath::Pi();
        double S  = 4*pi*pow(180/pi,2);       // solid angle of a sphere
        double dS = S/L;                      // solid angle of a pixel

        for(int l=0;l<L;l++) {                // loop over the sky grid
          double prob = gSM.get(l);           // get probability per pixel
          prob/=dS;                           // normalize probability to prob per deg^2
          if(prob==0) gSM.set(l,1e-40);	      // set !=0 to force dark blue background
        }

        gSM.SetZaxisTitle("prob. per deg^{2}");
        gSM.Draw(57);
        TH2D* hsm = gSM.GetHistogram();
        hsm->GetZaxis()->SetTitleOffset(0.85);
        hsm->GetZaxis()->SetTitleSize(0.03);
        sprintf(fname, "%s/skyprobcc.%s", dirCED, "png");
        if(sbasedirCED!=NULL) gSM.Print(fname); else gSM.Write(REPLACE(fname,dirCED,gtype));

        gROOT->SetBatch(batch);
      }
    }
  }
}

int
CWB::ced::Write(double factor, bool fullCED) {

  int i,j,k,n,m,l;
  int nIFO = NET->ifoListSize();
  int K = NET->nLag;
  int M = NET->mdc__ID.size();
  if(M==0&&simulation) {
    cout << "CWB::ced::Write : No injected events found -> force simulation=0" << endl;
    simulation=0;
  }
  int ID;
  char home_ced_path[1024]="~waveburst/WWW/LSC/waveburst/ced";
  char home_ced_www[1024]="";
  FILE* hP = NULL;  // html file
  char search = NET->tYPe;
  TString analysis = "1G";
  if(NET->like()=='2') analysis="2G";
  if(NET->optim) analysis+=":SRA";
  else           analysis+=":MRA";
  if(NET->pattern==0)                     analysis+=":Packet(0)";
  if((NET->pattern!=0 && NET->pattern<0)) analysis+=TString::Format(":Packet(%d)",NET->pattern);
  if((NET->pattern!=0 && NET->pattern>0)) analysis+=TString::Format(":Packet(+%d)",NET->pattern);
  int wmod = NET->optim ? 1 : 0;	// gwtMRAwave mode

  wavearray<double> id;
  wavearray<double> tm; 

  bool batch = gROOT->IsBatch();
  gROOT->SetBatch(true);

  watplot PCH(const_cast<char*>("chirp"));
  watplot WTS(const_cast<char*>("wts"));
  watplot PTS(const_cast<char*>("pts"),200,20,800,500);
  gskymap gSM; gSM.SetOptions("","Geographic",2);
  gnetwork gNET(NET); 			// used to plot circles

  int optRate,optLayer,optLevel;

  TMarker inj; inj.SetMarkerStyle(29);  // injected pos     (white star)
  TMarker INJ; INJ.SetMarkerStyle(30);  // injected pos     (empty black star)
  TMarker rec; rec.SetMarkerStyle(29);  // recostructed pos (black star) 
  TMarker REC; REC.SetMarkerStyle(30);  // recostructed pos (emply white star) 
  TMarker det; det.SetMarkerStyle(20);  // detected pos     (black dot ) 
  TMarker DET; DET.SetMarkerStyle(24);  // detected pos     (empty white dot ) 

  double r2d = 180./TMath::Pi();
  double d2r = TMath::Pi()/180.;

  char dirCED[1024]="";
  char command[1024];
  char ifostr[20]="";
  char fname[1024];

  int    minTimeDet=nIFO;
  double minTime=1e40;
  double eventTime[NIFO];
  double lagTime[NIFO];
  double startSegTime,stopSegTime;
  int ifoid[NIFO],sortid[NIFO];
  TString xtitle[NIFO];
  detector *pD[NIFO];

  int rate = 1;
  if(NET->like()=='2') rate = NET->MRA ? 0 : 1; 		// likelihood2G
  else rate = int(2*NET->getifo(0)->TFmap.resolution(0)+0.5);   // likelihoodI

  //Fill in all skymaps
  double old_cc = NET->netCC;
  double old_rho = NET->netRHO;
  NET->netCC = -1;
  NET->netRHO = 0;

  // check if it is one detector network  
  if(nIFO==2) if(strcmp(NET->ifoName[0],NET->ifoName[1])==0) nIFO=1;   

  for(n=0; n<nIFO; n++) {
    sprintf(ifostr,"%s%s",ifostr,NET->ifoName[n]);
    pD[n] = NET->getifo(n);
    ifoid[n]=n;
  }
  TMath::Sort(nIFO,ifoid,sortid,false);

  // create sbasedirCED if selected
  if(sbasedirCED!=NULL) {
    sprintf(command,"mkdir -p %s",sbasedirCED);
    if(gSystem->Exec(command)) return 0;
  } else {
    rbasedirCED->cd();
  }

  for(k=0; k<K; k++) {  // loop over the lags

    id = NET->getwc(k)->get(const_cast<char*>("ID"), 0, 'L', rate);
    tm = NET->getwc(k)->get(const_cast<char*>("TIME"), 0, 'L', 0);  


    for(j=0; j<(int)id.size(); j++) {  // loop over cluster index

      ID = size_t(id.data[j]+0.5);

      if(NET->wdm() && NET->getwc(k)->sCuts[ID-1]!=-1) continue;    // skip rejected/processed clusters

      if(NET->like()=='2') EVT->output2G(NULL,NET,ID,k,factor);
      else		   EVT->output(NULL,NET,factor,ID,k); 

      if(EVT->rho[analysis=="1G"?1:0] < rho) continue;

      int masterDet=0;
      int lagMin=2147483647;
      for(n=0; n<nIFO; n++) if(EVT->lag[n]<lagMin) {lagMin=int(EVT->lag[n]);masterDet=n;}

      // create ced directory & index.html link
      strcpy(ifostr,"");
      for(n=0; n<nIFO; n++) sprintf(ifostr,"%s%s",ifostr,NET->ifoName[n]);
      if(sbasedirCED!=NULL) {
        //if(nIFO==1&&strlen(chName)>1) sprintf(dirCED, "%s/%s_%s", sbasedirCED, ifostr,chName);
        //else sprintf(dirCED, "%s/%s", sbasedirCED, ifostr);
        sprintf(dirCED, "%s/%s", sbasedirCED, ifostr);
        for(n=0; n<nIFO; n++) sprintf(dirCED, "%s_%.3f",dirCED,EVT->start[n]); 
        sprintf(command,"mkdir -p %s",dirCED);
        if(gSystem->Exec(command)) return 0;

        if (getenv("HOME_CED_WWW")!=NULL) {
          sprintf(home_ced_www,"%s",getenv("HOME_CED_WWW"));
        }
        if (getenv("HOME_CED_PATH")!=NULL) {
          sprintf(home_ced_path,"%s",getenv("HOME_CED_PATH"));
        }
        strcpy(ifostr,"");
        for(n=0; n<nIFO; n++) sprintf(ifostr,"%s%s",ifostr,NET->ifoName[sortid[n]]);
        if(TString(home_ced_www).BeginsWith("/")) {	// -> CED is not published by WWW but local
          if(!simulation) sprintf(command,"cp %s/index/cedindex_%s.html  %s/index.html",home_ced_path,ifostr,dirCED);
          else sprintf(command,"cp %s/index/cedindex_sim_%s.html  %s/index.html",home_ced_path,ifostr,dirCED);
        } else {
          if(!simulation) sprintf(command,"ln -s %s/index/cedindex_%s.html  %s/index.html",home_ced_path,ifostr,dirCED);
          else sprintf(command,"ln -s %s/index/cedindex_sim_%s.html  %s/index.html",home_ced_path,ifostr,dirCED);
        }
        if(gSystem->Exec(command)) return 0;
      } else {
        //if(nIFO==1&&strlen(chName)>1) sprintf(dirCED, "%s_%s", ifostr,chName);
        //else sprintf(dirCED, "%s", ifostr);
        sprintf(dirCED, "%s", ifostr);
        for(n=0; n<nIFO; n++) sprintf(dirCED, "%s_%.3f",dirCED,EVT->start[n]); 
        rbasedirCED->cd(); rbasedirCED->mkdir(dirCED)->cd(); 
      }

      if(NET->like()=='2') EVT->output2G(NULL,NET,ID,k,factor);
      else		   EVT->output(NULL,NET,factor,ID,k); 

      //Write(factor,ID,k); 
      if(EVT->rho[analysis=="1G"?1:0] < rho) continue;

      double bPP =0.01;
      double TH = 0.2;
      TH2F* hist = NULL;
      if(NET->like()=='2') ;	// likelihood2G
      else                 NET->likelihood(search, NET->acor, ID, k);

      // dump event file & get event parameters
      if(sbasedirCED!=NULL) {
        sprintf(fname, "%s/eventDump.txt", dirCED);
      } else {
        gRandom->SetSeed(0);
        int rnID = int(gRandom->Rndm(13)*1.e9);   
        UserGroup_t* uinfo = gSystem->GetUserInfo();
        TString uname = uinfo->fUser;
        gSystem->Exec(TString("mkdir -p /dev/shm/")+uname);
        sprintf(fname,"/dev/shm/%s/eventDump-%d.txt",uname.Data(),rnID);
      }
      EVT->dopen(fname,const_cast<char*>("w"));

      if(NET->like()=='2') EVT->output2G(NULL,NET,ID,k,factor);
      else		   EVT->output(NULL,NET,factor,ID,k); 

      Write(factor,ID,k,dirCED); 
      EVT->dclose();
      if(sbasedirCED==NULL) {
        TMacro macro(fname);
        macro.Write("eventDump.txt");
        gSystem->Exec(TString("rm ")+fname);
      }

      // compute event time & lags time
      for(n=0; n<nIFO; n++) eventTime[n]=(EVT->start[n]+EVT->stop[n])/2.;
      for(n=0; n<nIFO; n++) minTime = ((eventTime[n]-EVT->gps[n])<minTime) ? eventTime[n]-EVT->gps[n] : minTime;
      for(n=0; n<nIFO; n++) minTimeDet = ((eventTime[n]-EVT->gps[n])<=minTime) ? n : minTimeDet;
      for(n=0; n<nIFO; n++) lagTime[n]=eventTime[n]-EVT->gps[minTimeDet]-minTime;
      startSegTime= EVT->gps[minTimeDet];  
      // NOTE : we use EVT->duration[1] because it caintains the event stop-start 
      //        while  EVT->duration[0] contains the duration estimated with rms
      stopSegTime = EVT->gps[minTimeDet]+EVT->left[minTimeDet]+EVT->right[minTimeDet]+EVT->duration[1];  
      for(n=0; n<nIFO; n++) xtitle[n] = TString::Format("Time (sec) : GPS OFFSET = %.3f",EVT->gps[n]);
 
      // create jobSummary.html & eventSummary.html
      if(sbasedirCED!=NULL) sprintf(fname, "%s/jobSummary.html", dirCED);
      else                  sprintf(fname, "/dev/null"); 
      if((hP = fopen(fname, "w")) != NULL) {
        fprintf(hP,"<table border=1 cellspacing=0 width=100%% height=100%%>\n");
        fprintf(hP,"<tr align=\"center\">\n");
        if(nIFO==1) {
          fprintf(hP,"<th>DETECTOR</th>\n");
        } else {
          fprintf(hP,"<th>NETWORK</th>\n");
        }
        fprintf(hP,"<td>%s</td>\n",ifostr);
        fprintf(hP,"</tr>\n");
        fprintf(hP,"<tr align=\"center\">\n");
        if(nIFO==1) {
          fprintf(hP,"<th>CHANNEL</th>\n");
          fprintf(hP,"<td>%s</td>\n",chName);
          fprintf(hP,"</tr>\n");
          fprintf(hP,"<tr align=\"center\">\n");
          fprintf(hP,"<th>SEARCH</th>\n");
          fprintf(hP,"<td>%s</td>\n",analysis.Data());
        } else {
          fprintf(hP,"<th>SEARCH</th>\n");
          if(analysis=="1G") {
            if(search=='E' || search=='E') fprintf(hP,"<td>%s un-modeled (%c)</td>\n",analysis.Data(),search);
            if(search=='b' || search=='B') fprintf(hP,"<td>%s un-modeled (%c)</td>\n",analysis.Data(),search);
            if(search=='r' || search=='R') fprintf(hP,"<td>%s un-modeled (%c)</td>\n",analysis.Data(),search);
            if(search=='i' || search=='I') fprintf(hP,"<td>%s elliptical (%c)</td>\n",analysis.Data(),search);
            if(search=='g' || search=='G') fprintf(hP,"<td>%s circular (%c)</td>\n",analysis.Data(),search);
            if(search=='s' || search=='S') fprintf(hP,"<td>%s linear (%c)</td>\n",analysis.Data(),search);
          } else {
            char _search = std::tolower(search);
            if(_search=='r')               fprintf(hP,"<td>%s un-modeled(%c)</td>\n",analysis.Data(),search);
            if(_search=='i')               fprintf(hP,"<td>%s iota-wave(%c)</td>\n",analysis.Data(),search);
            if(_search=='p')               fprintf(hP,"<td>%s psi-wave(%c)</td>\n",analysis.Data(),search);
            if(_search=='e'||_search=='b') fprintf(hP,"<td>%s elliptical(%c)</td>\n",analysis.Data(),search);
            if(_search=='c'||_search=='g') fprintf(hP,"<td>%s circular(%c)</td>\n",analysis.Data(),search);
            if(_search=='l'||_search=='s') fprintf(hP,"<td>%s linear(%c)</td>\n",analysis.Data(),search);
          }
        }
        fprintf(hP,"</tr>\n");
        fprintf(hP,"<tr align=\"center\">\n");
        fprintf(hP,"<th>START SEGMENT</th>");
        fprintf(hP,"<td>%9.3f</td>\n",startSegTime);
        fprintf(hP,"</tr>\n");
        fprintf(hP,"<tr align=\"center\">\n");
        fprintf(hP,"<th>STOP SEGMENT</th>");
        fprintf(hP,"<td>%9.3f</td>\n",stopSegTime);
        fprintf(hP,"</tr>\n");
        if(simulation) {
          fprintf(hP,"<tr align=\"center\">\n");
          fprintf(hP,"<th>MDC</th>\n");
          fprintf(hP,"<td>%s</td>\n",NET->getmdcType(EVT->type[1]-1).c_str());
          fprintf(hP,"</tr>\n");
        }
        if(!simulation&&nIFO>1) {
          fprintf(hP,"<tr align=\"center\">\n");
          fprintf(hP,"<th>LAG</th>\n");
          fprintf(hP,"<td>");
          for(n=0; n<nIFO-1; n++) fprintf(hP,"%3.3f / ",lagTime[n]);
          fprintf(hP,"%3.3f",lagTime[nIFO-1]);
          fprintf(hP,"</td>\n");
          fprintf(hP,"</tr>\n");
        }
        fprintf(hP,"</table>\n");
        fclose(hP);hP = NULL;
      } else {
        cout << "netevent::ced() error: cannot open file " << fname <<". \n";
      } 

      // convert user macro into html and copy to ced dir
      if(sbasedirCED!=NULL) {
        THtml html;
        html.SetEtcDir(gSystem->ExpandPathName("$HOME_WAT/html/etc/html"));
        html.SetProductName("CED");
        TString html_input_dir=dirCED;
        html.SetInputDir(html_input_dir.Data());

        char cmd[1024];
        sprintf(cmd,"cp %s/html/etc/html/ROOT.css %s/",gSystem->ExpandPathName("$HOME_WAT"),dirCED);
        gSystem->Exec(cmd);
        sprintf(cmd,"cp %s/html/etc/html/ROOT.js %s/",gSystem->ExpandPathName("$HOME_WAT"),dirCED);
        gSystem->Exec(cmd);

        TString cwb_parameters_file;
        if(gSystem->Getenv("CWB_PARAMETERS_FILE")==NULL) {
          cout << "Error : environment CWB_PARAMETERS_FILE is not defined!!!" << endl;exit(1);
        } else {
          cwb_parameters_file=TString(gSystem->Getenv("CWB_PARAMETERS_FILE"));
        }
        TString cwb_uparameters_file;
        if(gSystem->Getenv("CWB_UPARAMETERS_FILE")==NULL) {
          cout << "Error : environment CWB_UPARAMETERS_FILE is not defined!!!" << endl;exit(1);
        } else {
          cwb_uparameters_file=TString(gSystem->Getenv("CWB_UPARAMETERS_FILE"));
        }

        // redirect stderr to /dev/null to getrid of messages produced by html.Convert
        fpos_t poserr; fflush(stderr); fgetpos(stderr, &poserr);
        int fderr = dup(fileno(stderr)); freopen("/dev/null", "w", stderr);
        // redirect stdout to /dev/null to getrid of messages produced by html.Convert
        fpos_t posout; fflush(stdout); fgetpos(stdout, &posout);
        int fdout = dup(fileno(stdout)); freopen("/dev/null", "w", stdout);

        char title[1024];
        sprintf(title,"<h2 class=\"convert\" align=\"center\"> %s </h2>",cwb_parameters_file.Data());
        html.Convert(cwb_parameters_file.Data(),cwb_parameters_file.Data(),dirCED,"",0,title);
        if(analysis=="1G") sprintf(cmd,"mv %s/cwb1G_parameters.C.html %s/cwb_parameters.C.html",dirCED,dirCED);
        else               sprintf(cmd,"mv %s/cwb2G_parameters.C.html %s/cwb_parameters.C.html",dirCED,dirCED);
        gSystem->Exec(cmd);

        char work_dir[1024]; sprintf(work_dir,"%s",gSystem->WorkingDirectory()); // working dir
        TString cwb_uparameters_path = TString(work_dir)+"/"+cwb_uparameters_file;
        if(cwb_uparameters_file!="") {
          sprintf(title,"<h2 class=\"convert\" align=\"center\"> %s </h2>",cwb_uparameters_path.Data());
          html.Convert(cwb_uparameters_path.Data(),cwb_uparameters_file.Data(),dirCED,"",0,title);
          // rename file to the standard user_parameters.C name
          TString cwb_iparameters_path = TString(dirCED)+"/"+gSystem->BaseName(cwb_uparameters_path.Data())+".html";
          TString cwb_oparameters_path = TString(dirCED)+"/user_parameters.C.html";
          gSystem->Rename(cwb_iparameters_path.Data(),cwb_oparameters_path.Data());
        }

        // restore the stderr output
        fflush(stderr); dup2(fderr, fileno(stderr)); close(fderr);
        clearerr(stderr); fsetpos(stderr, &poserr);
        // restore the stdout output
        fflush(stdout); dup2(fdout, fileno(stdout)); close(fdout);
        clearerr(stdout); fsetpos(stdout, &posout);
      }

      // create jobSummary.html & eventSummary.html
      if(sbasedirCED!=NULL) sprintf(fname, "%s/eventSummary.html", dirCED);
      else                  sprintf(fname, "/dev/null"); 
      if((hP = fopen(fname, "w")) != NULL) {
        fprintf(hP,"<table border=1 cellspacing=0 width=100%% height=100%%>\n");
        fprintf(hP,"<tr align=\"center\">\n");
        fprintf(hP,"<th>GPS TIME</th>\n");
        fprintf(hP,"<th>SNR</th>\n");
        if(nIFO==1) {
          fprintf(hP,"<th>DURATION</th>\n");
          fprintf(hP,"<th>FREQUENCY</th>\n");
          fprintf(hP,"<th>BANDWIDTH</th>\n");
        } else {
          fprintf(hP,"<th>RHO[0/1]</th>\n");
          fprintf(hP,"<th>CC[0/1/2/3]</th>\n");
          fprintf(hP,"<th>ED</th>\n");
          fprintf(hP,"<th>PHI</th>\n");
          fprintf(hP,"<th>THETA</th>\n");
        }
        fprintf(hP,"</tr>\n");
        fprintf(hP,"<tr align=\"center\">\n");
        fprintf(hP,"<td>%9.3f</td>\n",EVT->time[minTimeDet]);
        if(nIFO==1) {
          fprintf(hP,"<td>%4.1f</td>\n",sqrt(EVT->likelihood/2.));
          fprintf(hP,"<td>%3.3f</td>\n",EVT->duration[0]);
          fprintf(hP,"<td>%3.1f</td>\n",EVT->frequency[0]);
          fprintf(hP,"<td>%3.1f</td>\n",EVT->bandwidth[0]);
        } else {
          fprintf(hP,"<td>%4.1f</td>\n",sqrt(EVT->likelihood));
          fprintf(hP,"<td>%3.1f / %3.1f</td>\n",EVT->rho[0],EVT->rho[1]);
          fprintf(hP,"<td>%1.2f / %1.2f / %1.2f / %1.2f</td>\n",
                  EVT->netcc[0],EVT->netcc[1],EVT->netcc[2],EVT->netcc[3]);
          fprintf(hP,"<td>%1.2f</td>\n",EVT->neted[0]/EVT->ecor);
          fprintf(hP,"<td>%3.1f</td>\n",EVT->phi[0]);
          fprintf(hP,"<td>%3.1f</td>\n",90.-EVT->theta[0]);
        }
        fprintf(hP,"</tr>\n");
        fprintf(hP,"</table>\n");
        fclose(hP);hP = NULL;
      } else {
        cout << "netevent::ced() error: cannot open file " << fname <<". \n";
      } 

      // plots nRMS
      int nIFO_RMS=0;
      std::vector<TGraph*> mgraph(nIFO); 
      Color_t psd_color[8] = {kBlue, kRed, kGreen, kBlack, 6, 3, 8, 43};
      for(n=0; n<nIFO; n++) {

        if(pD[n]->nRMS.size()==0) continue;	// it is zero when CED is produced starting from SUPERCLUSTER stage

        nIFO_RMS++; 

        WSeries<double> nRMS = pD[n]->nRMS;

        double fLow  = pD[n]->getTFmap()->getlow();
        double fHigh = pD[n]->getTFmap()->gethigh();
        double rate  = pD[n]->getTFmap()->rate();
        if(fLow<16) fLow=16; 
        if(fHigh>rate/2.) fHigh=rate/2.; 

        WDM<double>* wdm = (WDM<double>*) pD[n]->nRMS.pWavelet;
        int M = pD[n]->nRMS.getLevel();
        double* map00 = wdm->pWWS;
        double dF = rate/M/2./2.;			// psd frequency resolution

        wavearray<double> psd(2*M);			// PSD @ event time
        psd.start(0);
        psd.rate(1./dF);

        wavearray<double> PSD=psd;			// average PSD
        PSD=0;

        int nPSD = pD[n]->nRMS.size()/(M+1);     	// number of psd in the nRMS object
        double etime = eventTime[n]-EVT->gps[n];	// event time from the beginning of the buffer

        double segLen = stopSegTime-startSegTime;
        int I = (int)etime/(segLen/nPSD);		// psd event index 

        for(int i=0; i<nPSD; ++i){
          if(i==I) psd[0] = map00[0]; 			// first half-band
          PSD[0]+= pow(map00[0],2); 			// first half-band
          for(int j=1; j<M; j++){
             if(i==I) psd[2*j-1] = map00[j];
             if(i==I) psd[2*j]   = map00[j];
             PSD[2*j-1]+= pow(map00[j],2);
             PSD[2*j]  += pow(map00[j],2);
          }
          if(i==I) psd[2*M-1] = map00[M]; 		// last half-band
          PSD[2*M-1]+= pow(map00[M],2); 		// last half-band
          map00+=M+1; 
        }
        // NOTE : check if the following normalization is correct when cfg->fResample>0
        psd*=sqrt(2.)/sqrt(inRate);				// oneside psd normalization
        for(int j=0; j<PSD.size(); j++) PSD[j] = sqrt(PSD[j]/nPSD)*sqrt(2.)/sqrt(inRate);	// oneside PSD normalization

        PTS.canvas->cd();
        gStyle->SetTitleFont(12,"D");
        sprintf(fname, "%s/%s_psd.%s", dirCED, NET->ifoName[n], gtype);
	PTS.plot((wavearray<double>&)psd, const_cast<char*>("ALP"), 1, fLow+2*dF, fHigh-2*dF);
	//PTS.plot((wavearray<double>&)PSD, const_cast<char*>("ALP"), 1, fLow+2*dF, fHigh-2*dF);
        PTS.graph[0]->SetTitle(TString::Format("Power Spectral Density after lines' removal @ Time (sec) = %.0f",EVT->time[minTimeDet]));  
        PTS.graph[0]->GetXaxis()->SetTitle("Frequency (Hz)      ");  
        PTS.graph[0]->GetYaxis()->SetTitle("[strain / #sqrt{Hz}]          ");  
        PTS.graph[0]->SetMarkerStyle(1);  
        PTS.graph[0]->SetMinimum(1.e-24);
        PTS.graph[0]->SetMaximum(1.e-20);
        PTS.graph[0]->SetLineWidth(2);
        PTS.graph[0]->SetLineColor(psd_color[n%8]);
        if(TString(NET->ifoName[n])=="L1") PTS.graph[0]->SetLineColor(TColor::GetColor("#66ccff"));
        if(TString(NET->ifoName[n])=="L1") PTS.graph[0]->SetMarkerColor(TColor::GetColor("#66ccff"));
        if(TString(NET->ifoName[n])=="H1") PTS.graph[0]->SetLineColor(kRed);
        if(TString(NET->ifoName[n])=="H1") PTS.graph[0]->SetMarkerColor(kRed);
        if(TString(NET->ifoName[n])=="V1") PTS.graph[0]->SetLineColor(TColor::GetColor("#9900cc"));
        if(TString(NET->ifoName[n])=="V1") PTS.graph[0]->SetMarkerColor(TColor::GetColor("#9900cc"));
        PTS.canvas->SetLogx(true);  
        PTS.canvas->SetLogy(true);
        // draw the legend
        TLegend leg(0.753012,0.8-0.01,0.8885542,0.8793706,NULL,"brNDC");
        leg.SetTextAlign(22);
        leg.SetLineColor(kBlack);
        leg.AddEntry(PTS.graph[0],TString::Format("%s",NET->ifoName[n]),"lp");
        leg.Draw();
	//PTS.plot((wavearray<double>&)psd, const_cast<char*>("SAME"), 2, fLow+2*dF, fHigh-2*dF);
        if(sbasedirCED!=NULL) PTS.canvas->Print(fname); else PTS.canvas->Write(REPLACE(fname,dirCED,gtype));
        mgraph[n] = PTS.graph[0];
        PTS.graph.clear();
        //PTS.clear();

        sprintf(fname, "%s/%s_psd.dat", dirCED, NET->ifoName[n]);			
        if(sbasedirCED!=NULL) psd.wavearray<double>::Dump(const_cast<char*>(fname),2);	// save freq, psd format into ascii file
      }
      if(nIFO_RMS==nIFO) {		// nRMS is present for all detectors -> we can print the combined PSD
        PTS.graph=mgraph;;
        PTS.graph[0]->Draw("ALP"); for(n=1; n<nIFO; n++) PTS.graph[n]->Draw("SAME");
        // draw the legend
        TLegend leg(0.753012,0.8-nIFO*0.01,0.8885542,0.8793706,NULL,"brNDC");
        leg.SetTextAlign(22);
        leg.SetLineColor(kBlack);
        for(int n=0;n<nIFO;n++) leg.AddEntry(PTS.graph[n],TString::Format("%s",NET->ifoName[n]),"lp");
        leg.Draw();
        sprintf(fname, "%s/NET_psd.%s", dirCED, gtype);
        if(sbasedirCED!=NULL) PTS.canvas->Print(fname); else PTS.canvas->Write(REPLACE(fname,dirCED,gtype));
      }
      PTS.canvas->SetLogx(false);  
      PTS.canvas->SetLogy(false);
      PTS.clear();
      mgraph.clear();

      // get TF maps @ optimal resolution level
      optRate  = (NET->getwc(k)->cRate[ID-1])[0];
      optLayer = NET->getwc(k)->rate/optRate;
      optLevel = int(log10(optLayer)/log10(2));
      WSeries<double>* pTF[NIFO_MAX];
      for(int n=0; n<nIFO; n++) {
        pTF[n] = pD[n]->getTFmap();
        if(useSparse&&(analysis!="1G")) { // data not present in TF map (2G phase2), the sparse TF maps are used
          int optSparse = pD[n]->getSTFind(optRate);  // get optimal sparse map index
          pD[n]->vSS[optSparse].Expand(true); // rebuild wseries from sparse table
          pTF[n] = new WSeries<double>(*(WSeries<double>*)&(pD[n]->vSS[optSparse]));
          pTF[n]->setlow(pD[n]->getTFmap()->getlow());
          pTF[n]->sethigh(pD[n]->getTFmap()->gethigh());
        }
        // rescale data when data are resampled (resample with wavelet change the amplitude)
        double rescale = 1./pow(sqrt(2.),TMath::Log2(inRate/pTF[n]->rate())); 
        *pTF[n]*=rescale;	// rescale TF map
      }

      // plot spectrograms 
      for(int n=0; n<nIFO; n++) {
        if(pTF[n]->size()==0) continue;
        if(pTF[n]->getLevel()>0) pTF[n]->Inverse();

        int nfact=4;
        int nfft=nfact*512;
        int noverlap=nfft-10;
        double fparm=nfact*6;
        //int ystart = (EVT->start[n]-pTF[n]->start()-1)*pTF[n]->rate(); 
        //int ystop  = (EVT->stop[n]-pTF[n]->start()+1)*pTF[n]->rate(); 
        int ystart = int((EVT->start[n]-EVT->gps[n]-1)*pTF[n]->rate()); 
        int ystop  = int((EVT->stop[n]-EVT->gps[n]+1)*pTF[n]->rate()); 
        ystart-=nfft;
        ystop+=nfft;
        int ysize=ystop-ystart;
        wavearray<double> y;y.resize(ysize);y.rate(pTF[n]->rate());y.start(ystart/pTF[n]->rate());

        // stft use dt=y.rate() to normalize data but whitened data are already normalized by dt 
        // so before stft data must be divided by 1./sqrt(dt)
        for(int i=0;i<(int)y.size();i++) y.data[i]=pTF[n]->data[i+ystart]/sqrt(y.rate());

        CWB::STFT stft(y,nfft,noverlap,"energy","gauss",fparm);

        TCanvas* canvas;
        double tstart = nfft/pTF[n]->rate()+ystart/pTF[n]->rate();
        double tstop = (ysize-nfft)/pTF[n]->rate()+ystart/pTF[n]->rate();
        stft.Draw(tstart,tstop,pTF[n]->getlow(),pTF[n]->gethigh(),0,spectrogram_zmax,1);
        stft.GetHistogram()->GetXaxis()->SetTitle(xtitle[n]);  
        sprintf(fname, "%s/%s_spectrogram_1.%s", dirCED, NET->ifoName[n], gtype);
        canvas = stft.GetCanvas();
        if(sbasedirCED!=NULL) stft.Print(fname); else canvas->Write(REPLACE(fname,dirCED,gtype));
        canvas->SetLogy(true);
        stft.GetHistogram()->GetXaxis()->SetTitle(xtitle[n]);  
        sprintf(fname, "%s/%s_spectrogram_logy_1.%s", dirCED, NET->ifoName[n], gtype);
        if(sbasedirCED!=NULL) stft.Print(fname); else canvas->Write(REPLACE(fname,dirCED,gtype));

	tstart+=0.9;tstop-=0.9;
        stft.Draw(tstart,tstop,pTF[n]->getlow(),pTF[n]->gethigh(),0,spectrogram_zmax,1);
        stft.GetHistogram()->GetXaxis()->SetTitle(xtitle[n]);  
        sprintf(fname, "%s/%s_spectrogram_0.%s", dirCED, NET->ifoName[n], gtype);
        canvas = stft.GetCanvas();
        if(sbasedirCED!=NULL) stft.Print(fname); else canvas->Write(REPLACE(fname,dirCED,gtype));
        canvas->SetLogy(true);
        stft.GetHistogram()->GetXaxis()->SetTitle(xtitle[n]);  
        sprintf(fname, "%s/%s_spectrogram_logy_0.%s", dirCED, NET->ifoName[n], gtype);
        if(sbasedirCED!=NULL) stft.Print(fname); else canvas->Write(REPLACE(fname,dirCED,gtype));

        y.resize(0); 
      }

      // set TF maps to optimal resolution level
      for(int n=0; n<nIFO; n++) {
        float fLow = pTF[n]->getlow();
        float fHigh = pTF[n]->gethigh();
        if(useSparse&&(analysis!="1G")) { // data not present in TF map (2G phase2), use sparse TP map
          int optSparse = pD[n]->getSTFind(optRate);  // get optimal sparse map index
          delete pTF[n];
          pTF[n] = new WSeries<double>(*(WSeries<double>*)&(pD[n]->vSS[optSparse]));
        } else {
          if(pTF[n]->size()==0) continue;
          if(NET->wdm()) {if(NET->getwdm(optLayer+1)) pTF[n]->Forward(*pTF[n],*NET->getwdm(optLayer+1));} 
          else           pTF[n]->Forward(optLevel);
        }
        pTF[n]->setlow(fLow);
        pTF[n]->sethigh(fHigh);
      }

      WTS.canvas->cd();

      // plot scalogram maps at optimal resolution level
      for(int n=0; n<nIFO; n++) {
        if(pTF[n]->size()==0) continue;
        sprintf(fname, "%s/%s_scalogram_1.%s", dirCED, NET->ifoName[n], gtype);
        WTS.plot(pTF[n], 2, EVT->start[n]-EVT->slag[n]-1,
                 EVT->stop[n]-EVT->slag[n]+1,const_cast<char*>("COLZ"));
        WTS.hist2D->GetYaxis()->SetRangeUser(pTF[n]->getlow(),pTF[n]->gethigh());  
        WTS.hist2D->GetXaxis()->SetTitle(xtitle[n]);  
        if(sbasedirCED!=NULL) WTS.canvas->Print(fname); else WTS.canvas->Write(REPLACE(fname,dirCED,gtype));
        WTS.clear();
      }

      // plot zoomed scalogram maps at optimal resolution level
      for(int n=0; n<nIFO; n++) {
        if(pTF[n]->size()==0) continue;
        WTS.plot(pTF[n], 2, EVT->start[n]-EVT->slag[n]-0.1,
                 EVT->stop[n]-EVT->slag[n]+0.1, const_cast<char*>("COLZ"));
        WTS.hist2D->GetYaxis()->SetRangeUser(pTF[n]->getlow(),pTF[n]->gethigh());  
        WTS.hist2D->GetXaxis()->SetTitle(xtitle[n]);  
        sprintf(fname, "%s/%s_scalogram_0.%s", dirCED, NET->ifoName[n], gtype);
        if(sbasedirCED!=NULL) WTS.canvas->Print(fname); else WTS.canvas->Write(REPLACE(fname,dirCED,gtype));
        //sprintf(fname, "%s/%s_scalogram_0.%s", dirCED, NET->ifoName[n], "root");
        //WTS.canvas->Print(fname);
        WTS.clear();
      }

      // resize to 0 sseries
      for(int n=0; n<nIFO; n++) {
        if(useSparse&&(analysis!="1G")) { // data not present in TF map (2G phase2), the sparse TF maps are used
          int optSparse = pD[n]->getSTFind(optRate);  // get optimal sparse map index
          pD[n]->vSS[optSparse].Shrink(); // resize 0 wseries : leave only sparse table
          delete pTF[n];
        } else {
          // rescale data when data are resampled (resample with wavelet change the amplitude)
          double rescale = 1./pow(sqrt(2.),TMath::Log2(inRate/pTF[n]->rate())); 
          *pTF[n]*=1./rescale;	// restore original rescaled TF map
        }
      }

      // get reconstructed wave forms, signal

      int type = 1;

      PTS.canvas->cd();

      if(NET->like()=='2') {NET->getMRAwave(ID,k,'S',wmod,true);NET->getMRAwave(ID,k,'W',wmod,true);}
      else                 NET->getwave(ID, k, 'W');

      // rescale data when data are resampled (resample with wavelet change the amplitude)
      double rescale = 1./pow(sqrt(2.),TMath::Log2(inRate/pD[0]->waveForm.rate())); 
      for(int n=0; n<nIFO; n++) {
        pD[n]->waveForm*=rescale; pD[n]->waveBand*=rescale;	// rescale data
        double bT,eT;
        GetBoundaries((wavearray<double>&)pD[n]->waveForm, 0.999, bT, eT);
        sprintf(fname, "%s/%s_wf_signal.%s", dirCED, NET->ifoName[n], gtype);
        if(NET->wdm()) {
	  PTS.plot((wavearray<double>&)pD[n]->waveBand, const_cast<char*>("ALP"), kGray, bT, eT);
          PTS.plot((wavearray<double>&)pD[n]->waveForm, const_cast<char*>("SAME"), kRed, bT, eT);
        } else {
	  PTS.plot((wavearray<double>&)pD[n]->waveBand, const_cast<char*>("ALP"), kGray, bT, eT);
          PTS.plot((wavearray<double>&)pD[n]->waveForm, const_cast<char*>("SAME"), kRed, bT, eT);
        }
        PTS.graph[0]->GetXaxis()->SetTitle(xtitle[minTimeDet]);  
        if(sbasedirCED!=NULL) PTS.canvas->Print(fname); else PTS.canvas->Write(REPLACE(fname,dirCED,gtype));
        PTS.clear();

        sprintf(fname, "%s/%s_wf_signal_fft.%s", dirCED, NET->ifoName[n], gtype);
        double flow  = EVT->low[0];
        double fhigh = EVT->high[0];
        PTS.canvas->SetLogy(true);
        PTS.plot((wavearray<double>&)pD[n]->waveBand, const_cast<char*>("ALP"), kGray, 0., 0., true, flow, fhigh);
        PTS.plot((wavearray<double>&)pD[n]->waveForm, const_cast<char*>("SAME"), kRed, 0., 0., true, flow, fhigh);
        if(sbasedirCED!=NULL) PTS.canvas->Print(fname); else PTS.canvas->Write(REPLACE(fname,dirCED,gtype));
        PTS.canvas->SetLogy(false);
        PTS.clear();
      }

      //save whitenet waveforms (SNR^2 = the sum of the squares of the samples)

      for(int n=0; n<nIFO; n++) {
        pD[n]->waveForm*=1./rescale; pD[n]->waveBand*=1./rescale;	// restore original rescaled data
        //pD[n]->waveForm*=sqrt(pD[n]->waveForm.rate()); 			// normalize data to get snr=sqrt(sum(amp*amp*dt))
        sprintf(fname, "%s/%s_wf_signal.dat", dirCED, NET->ifoName[n]);
        double wstart = pD[n]->waveForm.wavearray<double>::start();						// save relative start time
        pD[n]->waveForm.wavearray<double>::start(EVT->gps[n]+pD[n]->waveForm.wavearray<double>::start());	// set absolute time
        if(sbasedirCED!=NULL) pD[n]->waveForm.wavearray<double>::Dump(const_cast<char*>(fname),2);		// time, strain format
        //if(sbasedirCED!=NULL) pD[n]->waveForm.Dump(const_cast<char*>(fname));
        //pD[n]->waveForm*=1./sqrt(pD[n]->waveForm.rate()); 		// restore original rescaled data
      }

      // get reconstructed waveforms, signal + noise

      for(int n=0; n<nIFO; n++) pD[n]->waveBand.sethigh(0);

      if(NET->like()=='2') {NET->getMRAwave(ID,k,'s',wmod,true);NET->getMRAwave(ID,k,'w',wmod,true);}
      else                 NET->getwave(ID, k, 'S');
      for(int n=0; n<nIFO; n++) {
        pD[n]->waveForm*=rescale; pD[n]->waveBand*=rescale;	// rescale data
        sprintf(fname, "%s/%s_wf_noise.%s", dirCED, NET->ifoName[n], gtype);
        if(NET->wdm()) {
          PTS.plot((wavearray<double>&)pD[n]->waveBand, const_cast<char*>("ALP"), kGray);
          PTS.plot((wavearray<double>&)pD[n]->waveForm, const_cast<char*>("SAME"), kRed);
        } else {
          PTS.plot((wavearray<double>&)pD[n]->waveBand, const_cast<char*>("ALP"), kGray);
          PTS.plot((wavearray<double>&)pD[n]->waveForm, const_cast<char*>("SAME"), kRed);
        }
        PTS.graph[0]->GetXaxis()->SetTitle(xtitle[minTimeDet]);  
        if(sbasedirCED!=NULL) PTS.canvas->Print(fname); else PTS.canvas->Write(REPLACE(fname,dirCED,gtype));
        PTS.clear();
      }

      // get reconstructed waveforms, strain
      double bT,eT;
      for(int n=0; n<nIFO; n++) pD[n]->waveBand.sethigh(0);
      if(NET->like()=='2') NET->getMRAwave(ID,k,'s',wmod,true);
      else                 NET->getwave(ID, k, 'w');
      for(int n=0; n<nIFO; n++) {
        pD[n]->waveForm*=rescale; 	// rescale data
        GetBoundaries((wavearray<double>&)pD[n]->waveForm, 0.999, bT, eT);
        sprintf(fname, "%s/%s_wf_strain.%s", dirCED, NET->ifoName[n], gtype);
        PTS.plot((wavearray<double>&)pD[n]->waveForm, const_cast<char*>("ALP"), kRed, bT, eT);
        PTS.graph[0]->GetXaxis()->SetTitle(xtitle[minTimeDet]);  
        if(sbasedirCED!=NULL) PTS.canvas->Print(fname); else PTS.canvas->Write(REPLACE(fname,dirCED,gtype));
        PTS.clear();

        GetBoundaries((wavearray<double>&)pD[n]->waveForm, 0.95, bT, eT);
        sprintf(fname, "%s/%s_wf_strain_zoom.%s", dirCED, NET->ifoName[n], gtype);
        PTS.plot((wavearray<double>&)pD[n]->waveForm, const_cast<char*>("ALP"), kRed, bT, eT);
        PTS.graph[0]->GetXaxis()->SetTitle(xtitle[minTimeDet]);  
        if(sbasedirCED!=NULL) PTS.canvas->Print(fname); else PTS.canvas->Write(REPLACE(fname,dirCED,gtype));
        PTS.clear();

        double flow  = EVT->low[0];
        double fhigh = EVT->high[0];
        sprintf(fname, "%s/%s_wf_strain_fft.%s", dirCED, NET->ifoName[n], gtype);
        PTS.plot((wavearray<double>&)pD[n]->waveForm, const_cast<char*>("ALP"), kRed, 0., 0., true, flow, fhigh);
        if(sbasedirCED!=NULL) PTS.canvas->Print(fname); else PTS.canvas->Write(REPLACE(fname,dirCED,gtype));
        PTS.clear();
      }

      //save waveforms

      for(int n=0; n<nIFO; n++) {
        sprintf(fname, "%s/%s_wf_strain.dat", dirCED, NET->ifoName[n]);
        double wstart = pD[n]->waveForm.wavearray<double>::start();						// save relative start time
        pD[n]->waveForm.wavearray<double>::start(EVT->gps[n]+pD[n]->waveForm.wavearray<double>::start());	// set absolute time
        if(sbasedirCED!=NULL) pD[n]->waveForm.wavearray<double>::Dump(const_cast<char*>(fname),2);		// time, strain format
        pD[n]->waveForm.wavearray<double>::start(wstart);							// restore relative start time
        //if(sbasedirCED!=NULL) pD[n]->waveForm.Dump(const_cast<char*>(fname));
      }

      if(nIFO==1 || !fullCED) {
        if(NET->wdm()) NET->getwc(k)->sCuts[ID-1]=1;     // mark clusters as processed (rejected)
        goto skip_skymap;   // skip skymaps  
      }

      // plot polargrams

      if(analysis.Contains("2G")) {
        for(m=0;m<2;m++) {
          gStyle->SetLineColor(kBlack);
          TCanvas* CPol = gNET.DrawPolargram(m,NET);
          if(CPol!=NULL) {
            sprintf(fname, "%s/polargram_%d.%s", dirCED, m+1, gtype);
            if(sbasedirCED!=NULL) CPol->Print(fname); else CPol->Write(REPLACE(fname,dirCED,gtype));
          }
          gStyle->SetLineColor(kWhite);
        }
      } 

      //likelihood skymaps

      inj.SetX(EVT->phi[1]<180?EVT->phi[1]:EVT->phi[1]-360); inj.SetY(90.-EVT->theta[1]); // injected pos (white star)
      INJ.SetX(EVT->phi[1]<180?EVT->phi[1]:EVT->phi[1]-360); INJ.SetY(90.-EVT->theta[1]); // injected pos (ewhite star)
      rec.SetX(EVT->phi[0]<180?EVT->phi[0]:EVT->phi[0]-360); rec.SetY(90.-EVT->theta[0]); // recostr. pos (black star)
      REC.SetX(EVT->phi[0]<180?EVT->phi[0]:EVT->phi[0]-360); REC.SetY(90.-EVT->theta[0]); // recostr. pos (ewhite star)
      det.SetX(EVT->phi[3]<180?EVT->phi[3]:EVT->phi[3]-360); det.SetY(90.-EVT->theta[3]); // detected pos (black dot )
      DET.SetX(EVT->phi[3]<180?EVT->phi[3]:EVT->phi[3]-360); DET.SetY(90.-EVT->theta[3]); // detected pos (ewhite dot )

      inj.SetMarkerSize(2.5); inj.SetMarkerColor(kWhite);
      INJ.SetMarkerSize(2.5); INJ.SetMarkerColor(kBlack);
      rec.SetMarkerSize(2.5); rec.SetMarkerColor(kBlack);
      REC.SetMarkerSize(2.5); REC.SetMarkerColor(kWhite);
      det.SetMarkerSize(1.5); det.SetMarkerColor(kBlack);
      DET.SetMarkerSize(1.5); DET.SetMarkerColor(kWhite);

      gSM=NET->nSensitivity; 
      sprintf(fname, "%s/sensitivity_plus.%s", dirCED, gtype);
      gSM.Draw(); rec.Draw();REC.Draw(); det.Draw();DET.Draw(); if(simulation) {inj.Draw();INJ.Draw();}
      if(sbasedirCED!=NULL) gSM.Print(fname); else gSM.Write(REPLACE(fname,dirCED,gtype));
      gSM=NET->nAlignment; 
      sprintf(fname, "%s/sensitivity_cross.%s", dirCED, gtype);
      gSM.Draw(); rec.Draw();REC.Draw(); det.Draw();DET.Draw(); if(simulation) {inj.Draw();INJ.Draw();}
      if(sbasedirCED!=NULL) gSM.Print(fname); else gSM.Write(REPLACE(fname,dirCED,gtype));
      gSM=NET->nSkyStat;; 
      sprintf(fname, "%s/skystat.%s", dirCED, gtype);
      gSM.Draw(); rec.Draw();REC.Draw(); det.Draw();DET.Draw(); if(simulation) {inj.Draw();INJ.Draw();}
      if(sbasedirCED!=NULL) gSM.Print(fname); else gSM.Write(REPLACE(fname,dirCED,gtype));
      gSM=NET->nLikelihood;; 
      sprintf(fname, "%s/likelihood.%s", dirCED, gtype);
      gSM.Draw(); rec.Draw();REC.Draw(); det.Draw();DET.Draw(); if(simulation) {inj.Draw();INJ.Draw();}
      if(sbasedirCED!=NULL) gSM.Print(fname); else gSM.Write(REPLACE(fname,dirCED,gtype));
      gSM=NET->nNullEnergy;; 
      sprintf(fname, "%s/null_energy.%s", dirCED, gtype);
      gSM.Draw(); rec.Draw();REC.Draw(); det.Draw();DET.Draw(); if(simulation) {inj.Draw();INJ.Draw();}
      if(sbasedirCED!=NULL) gSM.Print(fname); else gSM.Write(REPLACE(fname,dirCED,gtype));
      gSM=NET->nCorrEnergy;; 
      sprintf(fname, "%s/corr_energy.%s", dirCED, gtype);
      gSM.Draw(); rec.Draw();REC.Draw(); det.Draw();DET.Draw(); if(simulation) {inj.Draw();INJ.Draw();}
      if(sbasedirCED!=NULL) gSM.Print(fname); else gSM.Write(REPLACE(fname,dirCED,gtype));
      gSM=NET->nPenalty;; 
      sprintf(fname, "%s/penalty.%s", dirCED, gtype);
      gSM.Draw(); rec.Draw();REC.Draw(); det.Draw();DET.Draw(); if(simulation) {inj.Draw();INJ.Draw();}
      if(sbasedirCED!=NULL) gSM.Print(fname); else gSM.Write(REPLACE(fname,dirCED,gtype));
      gSM=NET->nDisbalance;; 
      sprintf(fname, "%s/disbalance.%s", dirCED, gtype);
      gSM.Draw(); rec.Draw();REC.Draw(); det.Draw();DET.Draw(); if(simulation) {inj.Draw();INJ.Draw();}
      if(sbasedirCED!=NULL) gSM.Print(fname); else gSM.Write(REPLACE(fname,dirCED,gtype));
      gSM=NET->nCorrelation;; 
      sprintf(fname, "%s/correlation.%s", dirCED, gtype);
      gSM.Draw(); rec.Draw();REC.Draw(); det.Draw();DET.Draw(); if(simulation) {inj.Draw();INJ.Draw();}
      if(sbasedirCED!=NULL) gSM.Print(fname); else gSM.Write(REPLACE(fname,dirCED,gtype));
      gSM=NET->nNetIndex;; 
      gSM.GetHistogram()->GetZaxis()->SetRangeUser(1,nIFO);
      sprintf(fname, "%s/netindex.%s", dirCED, gtype);
      gSM.Draw(); rec.Draw();REC.Draw(); det.Draw();DET.Draw(); if(simulation) {inj.Draw();INJ.Draw();}
      if(sbasedirCED!=NULL) gSM.Print(fname); else gSM.Write(REPLACE(fname,dirCED,gtype));
      gSM.GetHistogram()->GetZaxis()->UnZoom();
      gSM=NET->nEllipticity;; 
      sprintf(fname, "%s/ellipticity.%s", dirCED, gtype);
      gSM.Draw(); rec.Draw();REC.Draw(); det.Draw();DET.Draw(); if(simulation) {inj.Draw();INJ.Draw();}
      if(sbasedirCED!=NULL) gSM.Print(fname); else gSM.Write(REPLACE(fname,dirCED,gtype));
      gSM=NET->nPolarisation;; 
      sprintf(fname, "%s/polarisation.%s", dirCED, gtype);
      gSM.Draw(); rec.Draw();REC.Draw(); det.Draw();DET.Draw(); if(simulation) {inj.Draw();INJ.Draw();}
      if(sbasedirCED!=NULL) gSM.Print(fname); else gSM.Write(REPLACE(fname,dirCED,gtype));
      gSM=NET->nProbability;; 
      gSM.GetHistogram()->GetZaxis()->SetRangeUser(0,gSM.max());
      sprintf(fname, "%s/probability.%s", dirCED, gtype);
      gSM.Draw(); rec.Draw();REC.Draw(); det.Draw();DET.Draw(); if(simulation) {inj.Draw();INJ.Draw();}
      if(sbasedirCED!=NULL) gSM.Print(fname); else gSM.Write(REPLACE(fname,dirCED,gtype));

      // Dump probability skymap 
      //int L = NET->nProbability.size();
      //wavearray<float> xprob(L);
      //for(int ll=0;ll<L;ll++) xprob.data[ll]=NET->nProbability.get(ll);
      sprintf(fname, "%s/probability.%s", dirCED, "root");
      if(sbasedirCED!=NULL) {
        gskymap gSM(NET->nProbability);
        gSM.SetOptions("","Geographic");
        gSM.DumpObject(fname);
      }
      //xprob.resize(0);

#ifdef _USE_HEALPIX
      // Dump2fits probability skymap  (healpix) 
      if(sbasedirCED!=NULL) { 
        sprintf(fname, "%s/probability.%s", dirCED, "fits");
        if(NET->nProbability.getType()) 
           NET->nProbability.Dump2fits(fname,EVT->time[0],const_cast<char*>(""),
                                       const_cast<char*>("PROB"),const_cast<char*>("pix-1"),'C');
      }
#endif

      // add great circles (rec pos) to probability skymap
      sprintf(fname, "%s/probability_circles.%s", dirCED, gtype);
      gNET.SetGskymap(gSM); 
      gNET.GetGskymap()->Draw();
      double phi,theta;
      CwbToGeographic(EVT->phi[0],EVT->theta[0],phi,theta);
      gNET.DrawCircles(phi,theta,(Color_t)kBlack,1,1,true);
      rec.Draw(); det.Draw();
      if(sbasedirCED!=NULL) gNET.GetGskymap()->Print(fname); 
      else gNET.GetGskymap()->Write(REPLACE(fname,dirCED,gtype));

      skip_skymap:	

      // get network time 
      double gps_start = EVT->time[masterDet]-EVT->duration[1];
      double gps_stop  = EVT->time[masterDet]+EVT->duration[1];

      // draw chirp : f^(-8/3) vs time
      if(NET->like()=='2') {	
        PCH.canvas->cd();
        clusterdata* pcd = &(NET->getwc(k)->cData[ID-1]);
        sprintf(fname, "%s/mchirp.%s", dirCED, gtype);
        if(pcd->chirp.GetN()>0) {	// number of points > 0
          PCH.plot(pcd, simulation ? EVT->chirp[0] : 0);
          if(sbasedirCED!=NULL) PCH.canvas->Print(fname); 
          else PCH.canvas->Write(REPLACE(fname,dirCED,gtype));
        }
        PCH.clear();
      }

      // in likelihood2G pixeLHood,pixeLNull are not defined 
      // use monster event display (multi resolution analysis)
      if(NET->like()=='2') {	
        bool isPCs = !(NET->optim&&std::isupper(search));	// are Principal Components ?
        WTS.canvas->cd();
        netcluster* pwc = NET->getwc(k);
        sprintf(fname, "%s/l_tfmap_scalogram.%s", dirCED, gtype);
        WTS.plot(pwc, ID, nIFO, isPCs?'L':'l', 0, const_cast<char*>("COLZ"),256,NET->pattern>0);
        WTS.hist2D->GetXaxis()->SetTitle(xtitle[minTimeDet]);  
        if(sbasedirCED!=NULL) WTS.canvas->Print(fname); 
        else WTS.canvas->Write(REPLACE(fname,dirCED,gtype));
        WTS.clear();
        sprintf(fname, "%s/n_tfmap_scalogram.%s", dirCED, gtype);
        WTS.plot(pwc, ID, nIFO, isPCs?'N':'n', 0, const_cast<char*>("COLZ"),256,NET->pattern>0);
        WTS.hist2D->GetXaxis()->SetTitle(xtitle[minTimeDet]);  
        if(sbasedirCED!=NULL) WTS.canvas->Print(fname); 
        else WTS.canvas->Write(REPLACE(fname,dirCED,gtype));
        WTS.clear();
      } else {
        // plot likelihood map 
        WTS.canvas->cd();
        sprintf(fname, "%s/l_tfmap_scalogram.%s", dirCED, gtype);
        WTS.plot(NET->pixeLHood, 0, gps_start-EVT->slag[masterDet], 
                 gps_stop-EVT->slag[masterDet], const_cast<char*>("COLZ"));
        WTS.hist2D->GetYaxis()->SetRangeUser(EVT->low[0],EVT->high[0]);  
        WTS.hist2D->SetTitle("Scalogram");
        WTS.hist2D->GetXaxis()->SetTitle(xtitle[minTimeDet]);  
        if(sbasedirCED!=NULL) WTS.canvas->Print(fname); else WTS.canvas->Write(REPLACE(fname,dirCED,gtype));
        WTS.clear();
  
        // plot null map 
        sprintf(fname, "%s/n_tfmap_scalogram.%s", dirCED, gtype);
        WTS.plot(NET->pixeLNull, 0, gps_start-EVT->slag[masterDet], 
                 gps_stop-EVT->slag[masterDet], const_cast<char*>("COLZ"));
        WTS.hist2D->GetYaxis()->SetRangeUser(EVT->low[0],EVT->high[0]);  
        WTS.hist2D->SetTitle("Scalogram");
        WTS.hist2D->GetXaxis()->SetTitle(xtitle[minTimeDet]);  
        if(sbasedirCED!=NULL) WTS.canvas->Print(fname); else WTS.canvas->Write(REPLACE(fname,dirCED,gtype));
        WTS.clear();
      }

      // mark clusters as processed (rejected)
      if(NET->wdm()) NET->getwc(k)->sCuts[ID-1]=1;     

    } // End loop on found events
  }   // End loop on lags

  // restore NET thresholds
  NET->netCC = old_cc;
  NET->netRHO = old_rho;

  gROOT->SetBatch(batch);  // restore batch status

  return 1;
}

double
CWB::ced::GetBoundaries(wavearray<double> x, double P, double& bT, double& eT) {

  if(P<0) P=0;
  if(P>1) P=1;

  int N = x.size();

  double E = 0;							// signal energy
  double avr = 0;						// average
  for(int i=0;i<N;i++) {avr+=i*x[i]*x[i]; E+=x[i]*x[i];}
  int M=int(avr/E);						// central index

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
    if(sum/E > P) break;
  }

  bT = x.start()+jB/x.rate();
  eT = x.start()+jE/x.rate();

  return eT-bT;
}

void 
CWB::ced::SetOptions(TString options) {

  if(TString(options)!="") {

    TObjArray* token = TString(options).Tokenize(TString(' '));
    for(int j=0;j<token->GetEntries();j++) {

      TObjString* tok = (TObjString*)token->At(j);
      TString stok = tok->GetString();

      // set max z value used for spectrograms
      if(stok.Contains("spectrogram_zmax=")) {
        TString _spectrogram_zmax=stok;
        _spectrogram_zmax.Remove(0,_spectrogram_zmax.Last('=')+1);
        if(_spectrogram_zmax.IsFloat()) spectrogram_zmax=_spectrogram_zmax.Atof();
      }
    }
  }
}

