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


#include "STFT.hh"
#include "TPaletteAxis.h"

CWB::STFT::STFT(wavearray<double> x, int nfft, int noverlap, TString ztype,
                TString fwindow, double fparam, TString name) : 
                isLogz(false), title("") {

  this->name = name;
  this->canvas = NULL;
  this->h2 = NULL;
  this->ztype = ztype;
  int iztype=0; 
  if(ztype.CompareTo("amplitude")==0) iztype=1;
  if(ztype.CompareTo("energy")==0) iztype=2;
  if(iztype==0) {cout << "CWB::STFT::STFT: Error wrong ztype (amplitude/energy)" << endl;exit(1);} 

  if (nfft<=0) {cout << "CWB::STFT::STFT: Error nfft must be positive" << endl;exit(1);} 
  if (noverlap<=0||noverlap>=nfft) {cout << "CWB::STFT::STFT: Error noverlap must be > 0 & < nfft" << endl;exit(1);} 
  if (x.rate()<=0) {cout << "CWB::STFT::STFT: Error sample rate must be > 0" << endl;exit(1);} 
  if (x.start()<0) {cout << "CWB::STFT::STFT: Error start time must be positive" << endl;exit(1);} 
  if (int(x.size())<nfft) {cout << "CWB::STFT::STFT: Error size must be > nfft " << endl;exit(1);} 

  //cout << "nfft        : " << nfft << endl;
  //cout << "sample rate : " << x.rate() << endl;
  //cout << "start time  : " << x.start() << endl;

  // Compute Window
  CWB::Window wnd(const_cast<char*>(fwindow.Data()),nfft,fparam);
  window = new double[nfft];
  for (int i=0;i<nfft;i++) window[i] = wnd.GetValue(i);

  double df=(double)x.rate()/(double)(nfft);
  int loops = x.size()/nfft;

  int nshift=nfft-noverlap;
  int N=int(nfft/nshift);  
  double Tmin=x.start();
  double Tmax=x.size()/x.rate()+x.start();
  int nT=int(x.size()/nshift);
  double Fmin=0.;
  double Fmax=x.rate()/2.;
  int nF=nfft/2;
  h2 = new TH2D(name.Data(),"Spectrogram", nT, Tmin, Tmax, nF, Fmin, Fmax);
  wavearray<double> y(nfft);
//double tEN=0;
//double fEN=0;
  double dt=1./x.rate();
  double norm=sqrt(x.rate());
  for (int n=0;n<N*loops;n++) {
    int shift=n*nshift;
    if(shift+nfft>int(x.size())) break;
    for (int i=0;i<nfft;i++) y.data[i]=norm*x.data[i+shift]*window[i];
//for (int i=0;i<nfft;i++) tEN+=pow(y.data[i],2);
    y.FFT(1);
    double time = (shift+nfft/2)*dt;
    for (int i=0;i<nfft;i+=2) {
      double psd=sqrt(pow(y.data[i],2)+pow(y.data[i+1],2))*sqrt(1/df);
//fEN+=pow(psd,2);
      double frequency=(double)i*df/2.;
      if(iztype==2) psd*=psd;  // energy
      h2->SetBinContent(int(nT*time/(Tmax-x.start())),int(nF*frequency/Fmax),psd);
    }
//cout << tEN/x.rate() << " " << fEN*df << " " << df << endl;exit(0);
  }
  y.resize(0);
}

CWB::STFT::~STFT() {

  if(h2!=NULL) delete h2;
  if(canvas!=NULL) delete canvas;
  if(window!=NULL) delete [] window;
}

void
CWB::STFT::Draw(double t1, double t2, double f1, double f2, double z1, double z2,
                int dpaletteId, Option_t* goption) {

//  TCanvas* tcanvas = (TCanvas*) gROOT->FindObject(name);
//  if (tcanvas!=NULL) {cout << "CWB::STFT::STFT: Error Canvas " << name.Data() << " already exist" << endl;exit(1);} 
  if(canvas!=NULL) delete canvas;
  //canvas = new TCanvas(name, "LVC experiment", 300,40, 1000, 600);
  canvas = new TCanvas(name, "LVC experiment", 300,40, 800, 600);
  canvas->Clear();
  canvas->ToggleEventStatus();
  canvas->SetGridx(false);
  canvas->SetGridy(false);
  canvas->SetLogz(isLogz);
  canvas->SetFillColor(kWhite);
  canvas->SetRightMargin(0.10);
  canvas->SetLeftMargin(0.10);
  canvas->SetBottomMargin(0.13);
  canvas->SetBorderMode(0);
  //canvas->SetWindowSize(1200,600);
  //canvas->SetGrayscale();

  // remove the red box around canvas 
  gStyle->SetFrameBorderMode(0);
  gROOT->ForceStyle();

  gStyle->SetTitleH(0.050);
  gStyle->SetTitleW(0.95);
  gStyle->SetTitleY(0.98);
  gStyle->SetTitleFont(12,"D");
  gStyle->SetTitleColor(kBlue,"D");
  gStyle->SetTextFont(12);
  gStyle->SetTitleFillColor(kWhite);
  gStyle->SetLineColor(kWhite);
  gStyle->SetNumberContours(256);
//  gStyle->SetMarkerStyle(7);
//  gStyle->SetMarkerSize(2);
  gStyle->SetCanvasColor(kWhite);
  gStyle->SetStatBorderSize(1);

  if (dpaletteId==DUMMY_PALETTE_ID) {
    if (paletteId!=0) {
      SetPlotStyle(paletteId);
    } else {
      gStyle->SetPalette(1,0);
    }
  } else {
    if (dpaletteId!=0) {
      SetPlotStyle(dpaletteId);
    } else {
      gStyle->SetPalette(1,0);
    }
  }

  canvas->cd();
  canvas->SetLogz(isLogz);

  h2->SetStats(kFALSE);
  h2->SetTitleFont(12);
  h2->SetTitle(title);
  h2->SetFillColor(kWhite);

  h2->GetXaxis()->SetNdivisions(506);
  h2->GetXaxis()->SetLabelFont(42);
  h2->GetXaxis()->SetLabelOffset(0.014);
  h2->GetXaxis()->SetTitleOffset(1.4);
  h2->GetYaxis()->SetTitleOffset(1.2);
  h2->GetYaxis()->SetNdivisions(506);
  h2->GetYaxis()->SetLabelFont(42);
  h2->GetYaxis()->SetLabelOffset(0.01);
  h2->GetZaxis()->SetLabelFont(42);
  h2->GetZaxis()->SetNoExponent(false);
  h2->GetZaxis()->SetNdivisions(506);

  h2->GetXaxis()->SetTitleFont(42);
  h2->GetXaxis()->SetTitle("Time (sec)");
  h2->GetXaxis()->CenterTitle(true);

  h2->GetYaxis()->SetTitleFont(42);
  h2->GetYaxis()->SetTitle("Frequency (Hz)");
  h2->GetYaxis()->CenterTitle(true);

  //char ztitle[256];
  //sprintf(ztitle,"Normalized tile %s",ztype.Data());

  h2->GetZaxis()->SetTitleOffset(0.6);
  h2->GetZaxis()->SetTitleFont(42);
  //h2->GetZaxis()->SetTitle(ztitle);
  h2->GetZaxis()->CenterTitle(true);

  h2->GetXaxis()->SetLabelSize(0.03);
  h2->GetYaxis()->SetLabelSize(0.03);
  h2->GetZaxis()->SetLabelSize(0.03);

  if(title.Sizeof()==1) {
    char stitle[256];
    sprintf(stitle,"Spectrogram (Normalized tile %s)",ztype.Data());
    h2->SetTitle(stitle);
  }
  double dt=h2->GetXaxis()->GetBinWidth(0);
  double df=h2->GetYaxis()->GetBinWidth(0);

  int nt1=0;
  int nt2=h2->GetNbinsX();
  int nf1=0;
  int nf2=h2->GetNbinsY();

  if(t2>t1) nt1=int((t1-h2->GetXaxis()->GetXmin())/dt);
  if(t2>t1) nt2=int((t2-h2->GetXaxis()->GetXmin())/dt);
  if(f2>f1) nf1=int(f1/df);
  if(f2>f1) nf2=int(f2/df);

  double h2max=0.0;
  for (int i=nt1;i<nt2;i++) {
    for (int j=nf1;j<nf2;j++) {
      double binc=h2->GetBinContent(i,j);
      if(binc>h2max) h2max=binc;
    }
  }
  //h2->GetZaxis()->SetRangeUser(0,int(h2max)+1);
  h2->GetZaxis()->SetRangeUser(0,h2max);

  if(t2>t1) h2->GetXaxis()->SetRangeUser(t1,t2);
  if(f2>f1) h2->GetYaxis()->SetRangeUser(f1,f2);
  if(z2>z1) h2->GetZaxis()->SetRangeUser(z1,z2);

  h2->Draw(goption);

  // change palette's width
  canvas->Update();
  TPaletteAxis *palette = (TPaletteAxis*)h2->GetListOfFunctions()->FindObject("palette");
  palette->SetX1NDC(0.91);
  palette->SetX2NDC(0.933);
  palette->SetTitleOffset(0.92);
  palette->GetAxis()->SetTickSize(0.01);
  canvas->Modified();

}

// -----------------------------------------------------------------
// http://ultrahigh.org/2007/08/20/making-pretty-root-color-palettes/
// -----------------------------------------------------------------
void
CWB::STFT::SetPlotStyle(int paletteId) {

  const Int_t NRGBs = 5;
  const Int_t NCont = 255;

  if (fabs(paletteId)==1) {
    Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
    Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
    Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
    Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
    if (paletteId<0) {
      TColor::CreateGradientColorTable(NRGBs, stops, blue, green, red, NCont);
    } else {
      TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
    }
  } else
  if (fabs(paletteId)==2) {
    Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
    Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
    //Double_t red[NRGBs]   = { 0.00, 0.00, 0.00, 1.00, 1.00 };
    Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
    //Double_t green[NRGBs] = { 0.00, 1.00, 1.00, 1.00, 0.00 };
    Double_t blue[NRGBs]  = { 1.00, 1.00, 0.00, 0.00, 0.00 };
    if (paletteId<0) {
      TColor::CreateGradientColorTable(NRGBs, stops, blue, green, red, NCont);
    } else {
      TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
    }
  } else
  if (fabs(paletteId)==3) {
    Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
    Double_t red[NRGBs]   = { 0.00, 0.09, 0.18, 0.09, 0.00 };
    Double_t green[NRGBs] = { 0.01, 0.02, 0.39, 0.68, 0.97 };
    Double_t blue[NRGBs]  = { 0.17, 0.39, 0.62, 0.79, 0.97 };
    if (paletteId<0) {
      TColor::CreateGradientColorTable(NRGBs, stops, blue, green, red, NCont);
    } else {
      TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
    }
  } else
  if (fabs(paletteId)==4) {
    Double_t stops[NRGBs] = { 0.00, 0.50, 0.75, 0.875, 1.00 };
    Double_t red[NRGBs]   = { 1.00, 1.00, 1.00, 1.00, 1.00 };
    Double_t green[NRGBs] = { 1.00, 0.75, 0.50, 0.25, 0.00 };
    Double_t blue[NRGBs]  = { 0.00, 0.00, 0.00, 0.00, 0.00 };
    if (paletteId<0) {
      TColor::CreateGradientColorTable(NRGBs, stops, blue, green, red, NCont);
    } else {
      TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
    }
  } else
  if (fabs(paletteId)==5) {  // Greyscale palette
    Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
    Double_t red[NRGBs]   = { 1.00, 0.84, 0.61, 0.34, 0.00 };
    Double_t green[NRGBs] = { 1.00, 0.84, 0.61, 0.34, 0.00 };
    Double_t blue[NRGBs]  = { 1.00, 0.84, 0.61, 0.34, 0.00 };
    if (paletteId<0) {
      TColor::CreateGradientColorTable(NRGBs, stops, blue, green, red, NCont);
    } else {
      TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
    }
  } else
  if (fabs(paletteId)==57) {  // kBird palette (is defined only in ROOT6)
    Double_t stops[9] = { 0., 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1.0};
    Double_t red[9]   = { 0.2082, 0.0592, 0.0780, 0.0232, 0.1802, 0.5301, 0.8186, 0.9956, 0.9764};
    Double_t green[9] = { 0.1664, 0.3599, 0.5041, 0.6419, 0.7178, 0.7492, 0.7328, 0.7862, 0.9832};
    Double_t blue[9]  = { 0.5293, 0.8684, 0.8385, 0.7914, 0.6425, 0.4662, 0.3499, 0.1968, 0.0539};
    if (paletteId<0) {
      TColor::CreateGradientColorTable(9, stops, blue, green, red, NCont);
    } else {
      TColor::CreateGradientColorTable(9, stops, red, green, blue, NCont);
    }
  }
  gStyle->SetNumberContours(NCont);

  return;
}

void 
CWB::STFT::Print(TString pname) {

  canvas->Print(pname);

/*
  if(TString(pname).Contains(".png")!=0) { // fix gray background for png plots
    TString gname=pname;
    gname.ReplaceAll(".png",".gif");
    canvas->Print(gname);
    char cmd[1024];
    sprintf(cmd,"convert %s %s",gname.Data(),pname.Data());
    //cout << cmd << endl;
    gSystem->Exec(cmd);
    sprintf(cmd,"rm %s",gname.Data());
    //cout << cmd << endl;
    gSystem->Exec(cmd);
  } else {
    canvas->Print(pname);
  }
*/
}
