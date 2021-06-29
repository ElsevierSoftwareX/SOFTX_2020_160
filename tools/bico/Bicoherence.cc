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


/***************************************************************************
                          CWBBicoherence.cpp  -  description
                             -------------------
    begin                : Sat Dec 21 2011
    copyright            : (C) 2011 by Gabriele Vedovato
    email                : vedovato@lnl.infn.it
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#include "Bicoherence.hh"

#define BLC_WARNING_SIZE (64*pow(2.,24.))    

using namespace CWB;

CWB::Bicoherence::Bicoherence(TString ch1name, TString ch2name, double srate, int bsize, 
                              int x_num_slices, int y_num_slices,
                              int x_slice_index, int y_slice_index,
                              int order, char* window_type) {

  // check parameters
  if (x_num_slices<=0) {
    cout << "DAQ_BCM_X_NUM_SLICES : " << x_num_slices << endl;
    {cout << "CWB::Bicoherence::Bicoherence : bad input parameters" << endl;exit(1);}
  }
  if (y_num_slices<=0) {
    cout << "DAQ_BCM_Y_NUM_SLICES : " << y_num_slices << endl;
    {cout << "CWB::Bicoherence::Bicoherence : bad input parameters" << endl;exit(1);}
  }
  if (ch1name.Sizeof()==1) {
    {cout << "CWB::Bicoherence::Bicoherence : ch1name is not defined" << endl;exit(1);}
  }
  if (fmod(srate,x_num_slices)) {
    cout << "DAQ_BCM_X_NUM_SLICES : " << x_num_slices << endl;
    cout << "CWB::Bicoherence::Bicoherence : bad input parameters" << endl;
    cout << "RATE/DAQ_BCM_X_NUM_SLICES must be an integer !!!" << endl;exit(1);
  }
  if (fmod(srate,y_num_slices)) {
    cout << "DAQ_BCM_Y_NUM_SLICES : " << y_num_slices << endl;
    cout << "CWB::Bicoherence::Bicoherence : bad input parameters" << endl;
    cout << "RATE/DAQ_BCM_Y_NUM_SLICES must be an integer !!!" << endl;exit(1);
  }

  canvas = NULL;
  blhist = NULL;
  for (int i=0;i<3;i++) sgraph[i]=NULL;
  pgraph_id = -1;
  canvas_number = 1;

  this->chnames[0] = strdup(ch1name.Data()); 
  if (ch2name.Sizeof()!=1) this->chnames[1] = strdup(ch2name.Data()); 
  else this->chnames[1] = NULL;

  fxmin = -1;
  fxmax = -1;
  fymin = -1;
  fymax = -1;

  this->x_slice_index  = x_slice_index;    // x axis slice index
  this->y_slice_index  = y_slice_index;    // y axis slice index
  this->order = order;                     // order

  this->bstart = 0;
  this->srate = srate;
  this->bsize = bsize;
  blength = bsize/srate;	
  xbsize = bsize/x_num_slices;
  ybsize = bsize/y_num_slices;
  df = srate/bsize;

  // check parameters
  if ((x_slice_index+1)*xbsize/2 > bsize/2) {
    cout << "DAQ_BCM_X_SLICE_INDEX = " << x_slice_index << endl;
    {cout << "CWB::Bicoherence::Bicoherence : bad DAQ_BCM_X_SLICE_INDEX" << endl;exit(1);}
  }
  if ((y_slice_index+1)*ybsize/2 > bsize/2) {
    cout << "DAQ_BCM_Y_SLICE_INDEX = " << y_slice_index << endl;
    {cout << "CWB::Bicoherence::Bicoherence : bad DAQ_BCM_Y_SLICE_INDEX" << endl;exit(1);}
  }
  //if ((order==2) && (((x_slice_index+1)*xbsize/2)+((y_slice_index+1)*ybsize/2) > bsize/2)) {
  if ((order==2) && (((x_slice_index)*xbsize/2)+((y_slice_index)*ybsize/2) > bsize/2)) {
    cout << "DAQ_BCM_X_SLICE_INDEX = " << x_slice_index << endl;
    cout << "DAQ_BCM_Y_SLICE_INDEX = " << y_slice_index << endl;
    {cout << "CWB::Bicoherence::Bicoherence : bad DAQ_BCM_*_SLICE_INDEX" << endl;exit(1);}
  }
  if (order<1 || order>2)
    {cout << "CWB::Bicoherence::Bicoherence : DAQ_BCM_ORDER must be 1 or 2" << endl;exit(1);}

  unsigned int mSize = (xbsize/2)*(ybsize/2)*sizeof(complex<double>*);
  if (mSize > BLC_WARNING_SIZE) {
  char answer[256];
  strcpy(answer,"");			
  do {
    cout << endl;
    cout << "Requested memory size : " << mSize/pow(2.,20.) << " MB" 
         << " Warning : " << BLC_WARNING_SIZE/pow(2.,20.) << " MB" << endl; 
    cout << "Warning !!! Requested memory too big, possible troubles, do you want to continue ? (y/n) ";
    cin >> answer;
    cout << endl << endl;
  } while ((strcmp(answer,"y")!=0)&&(strcmp(answer,"n")!=0));    			
    if (strcmp(answer,"n")==0) exit(0);
  }

  B = (complex<double>**)malloc(xbsize/2*sizeof(complex<double>*));
  for (int i=0;i<xbsize/2;i++) B[i] = new complex<double>[ybsize/2];
  f = new double[bsize/2];
  SX = new double[xbsize/2];
  SY = new double[ybsize/2];
  xout  = new complex<double>[xbsize/2];
  yout  = new complex<double>[ybsize/2];

  for (int i=0;i<3;i++) S[i] = new double[bsize/2];

  for (int i=0;i<bsize/2;i++) f[i] = i*df;

  // ---------------------------------------
  // Compute Window
  // ---------------------------------------
	
  CWB::Window wnd(window_type,bsize);
  window = new double[bsize];
  for (int i=0;i<bsize;i++) window[i] = wnd.GetValue(i);

  Reset();

}

CWB::Bicoherence::~Bicoherence() {

  delete [] window;
  for (int i=0;i<xbsize/2;i++) delete [] B[i];
  free(B);
  delete [] SX;
  delete [] SY;
  delete [] f;
  for (int i=0;i<3;i++) delete [] S[i];
  delete [] xout;
  delete [] yout;
  for (int i=0;i<3;i++) delete sgraph[i];
}

bool
CWB::Bicoherence::MakeBicoherence(wavearray<double> x, wavearray<double> y) {

  if (x.rate()!=srate)
    {cout << "CWB::Bicoherence::MakeBicoherence - rate is not equals to the one declared in the setup" << endl;exit(1);}
  if ((int)x.size()!=bsize)
    {cout << "CWB::Bicoherence::MakeBicoherence - data size is not equals to the one declared in the setup" << endl;exit(1);}
  if(y.size()>0) { // 2 channels
    if (y.start()!=x.start())
      {cout << "CWB::Bicoherence::MakeBicoherence - the 2 channels have different start time" << endl;exit(1);}
    if (y.size()/y.rate()!=x.size()/x.rate())
      {cout << "CWB::Bicoherence::MakeBicoherence - the 2 channels have different length" << endl;exit(1);}
    if (((y_slice_index+1)*ybsize)>(int)y.size()) {
      double yfband = df*y.size()/2;
      double max_yfreq = df*((y_slice_index+1)*ybsize);
      cout << "CWB::Bicoherence::MakeBicoherence - the max yfreq " << max_yfreq << " is > yband " << yfband << endl;
      exit(1);
    }
  }

  if(!segListCheck(x.start(),x.start()+x.size()/x.rate())) return false;

  if(bstart==0) bstart=x.start();

  for (int j=0;j<bsize;j++) x[j]*=window[j];
  x.FFTW(1);
  if(y.size()>1) {
    for (int j=0;j<bsize;j++) y[j]*=window[j];
    y.FFTW(1);
  }

  for (int i=0;i<xbsize/2;i++) {
    int n = (x_slice_index*xbsize/2)+i;
    xout[i] = complex<double>(x[2*n],x[2*n+1]);
  }
  if (y.size()>0) {                                 // cross bcm
    for (int i=0;i<ybsize/2;i++) {
      int n = (y_slice_index*ybsize/2)+i;
      yout[i] = complex<double>(y[2*n],y[2*n+1]);
    }
  } else {
    for (int i=0;i<ybsize/2;i++) {
      int n = (y_slice_index*ybsize/2)+i;
      yout[i] = complex<double>(x[2*n],x[2*n+1]);
    }
  }

  switch (order) {
  case 1:
    for (int i=0;i<xbsize/2;i++) {
      for (int j=0;j<ybsize/2;j++) {
        B[i][j] += xout[i]*conj(yout[j]);
      }
    }
    break;
  case 2:
    int max_index = xbsize>ybsize ? xbsize/2 : ybsize/2;
    for (int i=0;i<xbsize/2;i++) {
      for (int j=0;j<ybsize/2;j++) {
        if(i+j>max_index) continue; 
        B[i][j] += xout[i]*yout[j]*conj(xout[i+j]);
      }
    }
    break;
  }

  for (int i=0;i<xbsize/2;i++) SX[i] += pow(abs(xout[i]),2);
  for (int i=0;i<ybsize/2;i++) SY[i] += pow(abs(yout[i]),2);

  for (int i=0;i<bsize/2;i++) {
    complex<double> X = complex<double>(x[2*i],x[2*i+1]);
    S[0][i] += pow(abs(X),2);
  }
  if (y.size()>0) {
    for (int i=0;i<(int)y.size()/2;i++) {
      complex<double> X = complex<double>(x[2*i],x[2*i+1]);
      complex<double> Y = complex<double>(y[2*i],y[2*i+1]);
      S[1][i] += pow(abs(Y),2);
      S[2][i] += abs(X*conj(Y));
    }
  } else {
    for (int i=0;i<bsize/2;i++) {
      complex<double> X = complex<double>(x[2*i],x[2*i+1]);
      S[1][i] += pow(abs(X),2);
      S[2][i] = S[1][i];
    }
  }

  bic_averages++;

  return true;
}

void
CWB::Bicoherence::Reset() {
  bic_averages = 0;
  for (int i=0;i<xbsize/2;i++) bzero(B[i],sizeof(complex<double>)*ybsize/2);	
  for (int i=0;i<3;i++) bzero(S[i],sizeof(double)*bsize/2);		
  bzero(SX,sizeof(double)*xbsize/2);		
  bzero(SY,sizeof(double)*ybsize/2);		
  bzero(xout,sizeof(complex<double>)*xbsize/2);		
  bzero(yout,sizeof(complex<double>)*ybsize/2);		
}

vector<bico>
CWB::Bicoherence::GetBicoherence(float threshold, int rebin) {

  if(threshold<0) threshold=0.;
  if(threshold>1) threshold=1.;

  rebin = xbsize%rebin==0 ? rebin : rebin-1; 

  int rxbsize = xbsize/2/rebin;
  int rybsize = ybsize/2/rebin;

  float** RB = (float**)malloc(rxbsize*sizeof(float*));
  for (int i=0;i<rxbsize;i++) RB[i] = new float[rybsize];
  for (int i=0;i<rxbsize;i++) bzero(RB[i],sizeof(float)*rybsize);	

  float** xRB = (float**)malloc(rxbsize*sizeof(float*));
  for (int i=0;i<rxbsize;i++) xRB[i] = new float[rybsize];
  for (int i=0;i<rxbsize;i++) bzero(xRB[i],sizeof(float)*rybsize);	

  float** yRB = (float**)malloc(rxbsize*sizeof(float*));
  for (int i=0;i<rxbsize;i++) yRB[i] = new float[rybsize];
  for (int i=0;i<rxbsize;i++) bzero(yRB[i],sizeof(float)*rybsize);	

  double xoff = srate*((x_slice_index)*xbsize/2)/bsize;	// slice starting freq in x axis
  double yoff = srate*((y_slice_index)*ybsize/2)/bsize;	// slice starting freq in y axis
 
  switch (order) {
  case 1:
    for (int i=0;i<xbsize/2;i++) {
      for (int j=0;j<ybsize/2;j++) {
        int I = TMath::FloorNint(i/rebin);
        int J = TMath::FloorNint(j/rebin);
        double binc = RB[I][J];
        double val = abs(B[i][j])/sqrt(SX[i]*SY[j]);  
        if(val>binc) {
          RB[I][J]=val;
          xRB[I][J]=i*df+xoff;
          yRB[I][J]=j*df+yoff;
        }
      }
    }
    break;
  case 2:
    int max_index = xbsize>ybsize ? xbsize/2 : ybsize/2;
    for (int i=1;i<xbsize/2;i++) {
      for (int j=1;j<ybsize/2;j++) {
        if((i+j)>max_index) continue; 
        int I = TMath::FloorNint(i/rebin);
        int J = TMath::FloorNint(j/rebin);
        double binc = RB[I][J];
        double val = abs(B[i][j])/sqrt(SX[i]*SY[j]*SX[i+j]);  
        if(val>binc) {
          RB[I][J]=val;
          xRB[I][J]=i*df+xoff; 
          yRB[I][J]=j*df+yoff;
        }
      }
    }
    break;
  }

  bico BICO; 
  vector<bico> olist;
  switch (order) {
  case 1:
    for (int i=0;i<rxbsize;i++) {
      for (int j=0;j<rybsize;j++) {
        double binc = RB[i][j];
        if(binc>threshold) {
          BICO.x=xRB[i][j];
          BICO.y=yRB[i][j];
          BICO.c=binc;
          olist.push_back(BICO);
        }
      }
    }
    break;
  case 2:
    for (int i=0;i<rxbsize;i++) {
      for (int j=0;j<rybsize;j++) {
        double binc = RB[i][j];
        if(binc>threshold) {
          BICO.x=xRB[i][j];
          BICO.y=yRB[i][j];
          BICO.c=binc;
          olist.push_back(BICO);
        }
      }
    }
    break;
  }

  for (int i=0;i<rxbsize;i++) delete [] RB[i];
  for (int i=0;i<rxbsize;i++) delete [] xRB[i];
  for (int i=0;i<rxbsize;i++) delete [] yRB[i];
  free(RB);
  free(xRB);
  free(yRB);

  return olist;
}

void
CWB::Bicoherence::DrawBicoherence(int rebin, int graph_id, TString ofname, bool batch) {

  if((graph_id<0)||(graph_id>2)) 
    {cout << "CWB::Bicoherence::DrawBicoherence - wrong graph id {0,1,2}" << endl;exit(1);}

  rebin = xbsize%rebin==0 ? rebin : rebin-1; 

  gROOT->SetBatch(batch);  

  // ---------------------------------------
  // Canvas
  // ---------------------------------------

  char  canvas_name[64];
  sprintf(canvas_name,"BicoherenceCanvas");
  if(canvas!=NULL) delete canvas;
  canvas = new TCanvas(const_cast<char*>(canvas_name),"Bicoherence",10,10,800,800);
  canvas->SetFillColor(10);
  canvas->Divide(1,2);
  canvas->ToggleEventStatus();
  canvas->SetCursor(kWatch);
  canvas->Update();
  canvas->SetGridx(true);
  canvas->SetGridy(true);

  // remove the red box around canvas
  gStyle->SetFrameBorderMode(0);
  gROOT->ForceStyle();

  gStyle->SetTitleH(0.12);
  gStyle->SetTitleW(0.95);
  gStyle->SetTitleY(0.98);
  gStyle->SetTitleFont(12,"D");
  gStyle->SetTitleColor(kBlue,"D");
  gStyle->SetTextFont(12);
  gStyle->SetTitleFillColor(kWhite);
//  gStyle->SetLineColor(kWhite);
  gStyle->SetNumberContours(256);
//  gStyle->SetMarkerStyle(7);
//  gStyle->SetMarkerSize(2);
  gStyle->SetTitleBorderSize(0);
  gStyle->SetCanvasColor(kWhite);
  gStyle->SetStatBorderSize(1);
  gStyle->SetPalette(1,0);

//  TStyle* myStyle = new TStyle("myStyle","my style");
//  myStyle->SetPalette(1,0);
//  myStyle->cd();

  double xmin = srate*((x_slice_index)*xbsize/2)/bsize;
  double xmax = srate*((x_slice_index+1)*xbsize/2)/bsize;

  double ymin = srate*((y_slice_index)*ybsize/2)/bsize;
  double ymax = srate*((y_slice_index+1)*ybsize/2)/bsize;
//  ymax = ymin+(ymax-ymin)/2;

  if(blhist!=NULL) delete blhist;
  if(graph_id==0 || graph_id==2) {
    blhist = new TH2F("BilinearCoupling","Bicoherence",xbsize/2/rebin,xmin,xmax,ybsize/2/rebin,ymin,ymax);
  } else {
    blhist = new TH2F("BilinearCoupling","Bicoherence",ybsize/2/rebin,ymin,ymax,xbsize/2/rebin,xmin,xmax);
  }
  blhist->SetStats(kFALSE);
  blhist->SetFillColor(46);

  char titles[3][256];
  if (chnames[1]!=NULL) {
    sprintf(titles[0],"%s   -   Frequency (Hz)",chnames[0]);
    sprintf(titles[1],"%s   -   Frequency (Hz)",chnames[1]);
    sprintf(titles[2],"%s %s   -   Frequency (Hz)",chnames[0],chnames[1]);
  } else {
    for (int i=0;i<2;i++) sprintf(titles[i],"%s   -   Frequency (Hz)",chnames[0]);
    sprintf(titles[2],"%s %s Frequency (Hz)",chnames[0],chnames[0]);
  }

  blhist->GetXaxis()->SetLabelSize(0.04);
  blhist->GetYaxis()->SetLabelSize(0.04);
  blhist->GetXaxis()->SetTitleSize(0.04);
  blhist->GetYaxis()->SetTitleSize(0.04);
  blhist->GetXaxis()->SetTitleFont(42);
  blhist->GetYaxis()->SetTitleFont(42);
  blhist->GetXaxis()->SetLabelFont(42);
  blhist->GetYaxis()->SetLabelFont(42);
  blhist->GetYaxis()->SetLabelOffset(0.01);
  blhist->GetYaxis()->SetTitleOffset(1.3);
  blhist->GetZaxis()->SetLabelSize(0.04);
  if(graph_id==0 || graph_id==2) {
    blhist->SetXTitle(titles[0]);
    blhist->SetYTitle(titles[1]);
  } else {
    blhist->SetXTitle(titles[1]);
    blhist->SetYTitle(titles[0]);
  }

  char	time_offset_str[128];
  char	time_start[128];
  char  time_range[128];
  sprintf(time_start," gps %d",int(bstart));
  sprintf(time_range," (duration = %4.3f)",bic_averages*blength);
  sprintf(time_offset_str,"avr %d ",bic_averages);
  strncpy(time_offset_str+strlen(time_offset_str),time_start,strlen(time_start)+1); 
  strcat(time_offset_str, time_range);

  char	blhist_title[128] = "";
  switch (order) {
  case 1:
    strcpy(blhist_title,"Coherence   ");
    break;
  case 2:
    strcpy(blhist_title,"Bicoherence   ");
    break;
  }
  strcat(blhist_title,time_offset_str);
  blhist->SetTitle(blhist_title);

  gStyle->SetTitleH(0.05);
  gStyle->SetTitleW(0.98);
  gStyle->SetTitleFont(72);
  gStyle->SetPalette(1,0);

  blhist->SetEntries((xbsize/2*ybsize/2)/(rebin*rebin));
  switch (order) {
  case 1:
    for (int i=0;i<xbsize/2;i++) {
      for (int j=0;j<ybsize/2;j++) {
        int I = TMath::FloorNint(i/rebin);
        int J = TMath::FloorNint(j/rebin);
        double val = abs(B[i][j])/sqrt(SX[i]*SY[j]);  
        if(graph_id==0 || graph_id==2) {
          double binc = blhist->GetBinContent(I+1,J+1);
          if(val>binc) blhist->SetBinContent(I+1,J+1,val);
        } else {
          double binc = blhist->GetBinContent(J+1,I+1);
          if(val>binc) blhist->SetBinContent(J+1,I+1,val);
        }
      }
    }
    break;
  case 2:
    int max_index = xbsize>ybsize ? xbsize/2 : ybsize/2;
    for (int i=1;i<xbsize/2;i++) {
      for (int j=1;j<ybsize/2;j++) {
        //blhist->SetBinContent(i/rebin+1,j/rebin+1,0.);
        if((i+j)>max_index) continue; 
        int I = TMath::FloorNint(i/rebin);
        int J = TMath::FloorNint(j/rebin);
        double val = abs(B[i][j])/sqrt(SX[i]*SY[j]*SX[i+j]);  
        if(graph_id==0 || graph_id==2) {
          double binc = blhist->GetBinContent(I+1,J+1);
          if(val>binc) blhist->SetBinContent(I+1,J+1,val);
        } else {
          double binc = blhist->GetBinContent(J+1,I+1);
          if(val>binc) blhist->SetBinContent(J+1,I+1,val);
        }
      }
    }
    break;
  }

  double* psd[3];
  for (int i=0;i<3;i++) psd[i] = new double[bsize/2];
  double norm = bic_averages*df;	
  for (int i=0;i<2;i++) for (int j=0;j<bsize/2;j++) psd[i][j]=sqrt(S[i][j]/norm);  // one side spectrum
  for (int j=0;j<bsize/2;j++) psd[2][j]=S[2][j]/sqrt(S[0][j]*S[1][j]);
  double max_psd[3]={0,0,0};
  for (int i=0;i<3;i++) {
    for (int j=0;j<bsize/2;j++) {
      double f=j*df;
      if(i==0 || i==2) {
        if(f>xmin&&f<xmax) if(psd[i][j]>max_psd[i]) max_psd[i]=psd[i][j];
      } else {
        if(f>ymin&&f<ymax) if(psd[i][j]>max_psd[i]) max_psd[i]=psd[i][j];
      }
    }
  }
  double min_psd[3];
  for (int i=0;i<3;i++) {
    min_psd[i]=max_psd[i];
    for (int j=0;j<bsize/2;j++) {
      double freq=j*df;
      if(i==0 || i==2) {
        if(freq>xmin&&freq<xmax) if(psd[i][j]<min_psd[i]) min_psd[i]=psd[i][j];
      } else {
        if(freq>ymin&&freq<ymax) if(psd[i][j]<min_psd[i]) min_psd[i]=psd[i][j];
      }
    }
  }

  canvas->cd(1);
  blhist->Draw("colz");

  for (int i=0;i<3;i++) {
    if(sgraph[i]!=NULL) delete sgraph[i];
    sgraph[i] = new TGraph(bsize/2, f, psd[i]);
    char sgraph_title[128] = "Power Spectrum   ";
    strcat(sgraph_title,time_offset_str);
    sgraph[i]->SetTitle(sgraph_title);

    char sgraph_name[128];
    if (i==0) sprintf(sgraph_name,"%s-psd",chnames[i]);
    if ((i==1) && (chnames[1]!=NULL)) sprintf(sgraph_name,"%s-psd",chnames[1]);
    if ((i==1) && (chnames[1]==NULL)) sprintf(sgraph_name,"%s-psd",chnames[0]);
    if ((i==2) && (chnames[1]!=NULL)) sprintf(sgraph_name,"%s-%s-psd",chnames[0],chnames[1]);
    if ((i==2) && (chnames[1]==NULL)) sprintf(sgraph_name,"%s-%s-psd",chnames[0],chnames[0]);
    sgraph[i]->SetName(sgraph_name);

    canvas->cd(2);
    sgraph[i]->Draw("APL");
    gPad->SetLogy(kTRUE);

    sgraph[i]->GetHistogram()->SetXTitle(titles[i]);
    sgraph[i]->GetHistogram()->SetYTitle("PSD");
    sgraph[i]->GetHistogram()->GetXaxis()->SetLabelSize(0.04);
    sgraph[i]->GetHistogram()->GetYaxis()->SetLabelSize(0.03);
    sgraph[i]->GetHistogram()->GetXaxis()->SetTitleSize(0.04);
    sgraph[i]->GetHistogram()->GetYaxis()->SetTitleSize(0.04);
    sgraph[i]->GetHistogram()->GetXaxis()->SetTitleFont(42);
    sgraph[i]->GetHistogram()->GetYaxis()->SetTitleFont(42);
    sgraph[i]->GetHistogram()->GetXaxis()->SetLabelFont(42);
    sgraph[i]->GetHistogram()->GetYaxis()->SetLabelFont(42);
    sgraph[i]->GetHistogram()->GetYaxis()->SetLabelOffset(0.01);
    sgraph[i]->GetHistogram()->GetYaxis()->SetTitleOffset(1.3);
    sgraph[i]->GetHistogram()->GetYaxis()->SetLabelSize(0.04);
    if ((fxmin!=-1)&&(fxmax!=-1)) {
      sgraph[i]->GetHistogram()->GetXaxis()->SetRange(fxmin,fxmax);
    } else {
      if(i==0 || i==2)
        sgraph[i]->GetHistogram()->GetXaxis()->SetRangeUser(xmin, xmax);
      else
        sgraph[i]->GetHistogram()->GetXaxis()->SetRangeUser(ymin, ymax);
    }
    sgraph[i]->GetHistogram()->SetMinimum(min_psd[i]);
    sgraph[i]->GetHistogram()->SetMaximum(max_psd[i]);
  }
  gPad->Clear();
  if (graph_id==2) gPad->SetLogy(kFALSE); else gPad->SetLogy(kTRUE);
  sgraph[graph_id]->Draw("APL");

  canvas->cd();
  canvas->Update();

  canvas->SetTitle(time_offset_str);

  if(ofname.Sizeof()>1) {
    TString plname=ofname; 
    if(plname.Contains(".png")) {
      plname.ReplaceAll(".png",".gif");
      cout << plname.Data() << endl;
      canvas->Print(plname.Data());
      char cmd[256];
      sprintf(cmd,"convert %s %s",plname.Data(),ofname.Data());
      cout << cmd << endl;
      gSystem->Exec(cmd);
      sprintf(cmd,"rm %s",plname.Data());
      cout << cmd << endl;
      gSystem->Exec(cmd);
    } else {
      canvas->Print(ofname.Data());
    }
  }

// START INTERACTIVE
/*
//  try {
    bool running = true;
    while (running) {
//      if (interactive) canvas->Run(true);
      double hxmin,hxmax;
      double hymin,hymax;
      EEventType event = canvas->GetInputEvent();
      if (event == kKeyPress) {
        Char_t chCharCodeX = (Char_t) canvas->GetInputX();
        Char_t chCharCodeY = (Char_t) canvas->GetInputY();
        unsigned  int keyPressed = 0x100*(int)chCharCodeY+(int)chCharCodeX;
        //printf("keyPressed %x\n",keyPressed);
        switch (keyPressed) {
        case  0x5151: // q
        case  0x7171: // Q
          exit(0);
          break;
        case  0x5858: // x   zoom x
        case  0x7878: // X
            sgraph[0]->GetHistogram()->GetXaxis()->SetRangeUser(xmin, xmax);
            canvas->cd(2);
            gPad->Clear();
            gPad->SetLogy(kTRUE);
            sgraph[0]->Draw("APL");
            canvas->cd();
            canvas->Update();
            graph_id = 0;
          break;
        case  0x5959: // y   zoom y
        case  0x7979: // Y
            sgraph[1]->GetHistogram()->GetXaxis()->SetRangeUser(ymin, ymin+(ymax-ymin));
            canvas->cd(2);
            gPad->Clear();
            gPad->SetLogy(kTRUE);
            sgraph[1]->Draw("APL");
            canvas->cd();
            canvas->Update();
            graph_id = 1;
          break;
        case  0x5a5a: // z   zoom xy
        case  0x7a7a: // Z
            sgraph[2]->GetHistogram()->GetXaxis()->SetRangeUser(ymin, ymin+(ymax-ymin));
            canvas->cd(2);
            gPad->Clear();
            gPad->SetLogy(kFALSE);
            sgraph[2]->Draw("APL");
            canvas->cd();
            canvas->Update();
            graph_id = 2;
          break;
        case  0x4646: // f   full
        case  0x6666: // F
            sgraph[graph_id]->GetHistogram()->GetXaxis()->SetRangeUser(0, df*bsize/2);
            canvas->cd(2);
            gPad->Clear();
            if (graph_id==2) gPad->SetLogy(kFALSE); else gPad->SetLogy(kTRUE);
            sgraph[graph_id]->Draw("APL");
            canvas->cd();
            canvas->Update();
          break;
        default:
          running = false;
        }
      } else {
        running = false;
      }
    }
//  } catch (AalCanvasException e) {
//    canvas = NULL;
//    throw e;
//  }

// END INTERACTIVE
*/
  canvas->SetCursor(kWatch);
  canvas->Update();

  fxmin = sgraph[graph_id]->GetHistogram()->GetXaxis()->GetFirst();
  fxmax = sgraph[graph_id]->GetHistogram()->GetXaxis()->GetLast();
  if (graph_id == pgraph_id) {
    fymin = sgraph[graph_id]->GetHistogram()->GetMinimum();
    fymax = sgraph[graph_id]->GetHistogram()->GetMaximum();
  } else {
    fymin=-1;fymax=-1;
    pgraph_id=graph_id;
  }

  for (int i=0;i<3;i++) delete [] psd[i];
}

int   
CWB::Bicoherence::readSegList(char* seglist, double shift, bool invert, bool c4) {

  CWB::Toolbox TB;

  dqfile DQF[1]={ {"L1", "", CWB_CAT0, shift, invert, c4} };
  strcpy(DQF[0].file, seglist);

  segList = TB.readSegList(1, DQF, CWB_CAT1);

  return segList.size();
}

bool
CWB::Bicoherence::segListCheck(double start, double stop) {
  
  int size = segList.size();
  if(size==0) return true;
  waveSegment SEG;
  bool check=false;
  for(int i=0;i<size;i++) {
     SEG = segList[i];
     if(start>=SEG.start && stop<=SEG.stop) {check=true;break;}
  }

  return check;
}
