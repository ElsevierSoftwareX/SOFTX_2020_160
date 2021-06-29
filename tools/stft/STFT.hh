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
 * Package:      STFT Class Library
 * File name:    STFT.hh
 * Author:       Gabriele Vedovato (vedovato@lnl.infn.it)
 **********************************************************/


#ifndef STFT_HH
#define STFT_HH

#include "TCanvas.h"
#include "TH2D.h"
#include "TPolyLine.h"
#include "TStyle.h"
#include "TColor.h"
#include "TObjString.h"
#include "TObjArray.h"
#include "TString.h"
#include "TMarker.h"
#include "TROOT.h"
#include "TMath.h"
#include "TPaveText.h"
#include "TSystem.h"
#include "TLatex.h"

#include <string>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>

#include "wavecomplex.hh"
#include "wavearray.hh"

#include "Window.hh"

#define DUMMY_PALETTE_ID 1000000000

using namespace std;

namespace CWB {

class STFT {

public:
  
  STFT(wavearray<double> x, int nfft, int noverlap, TString ztype="amplitude", 
       TString fwindow="hann", double fparam=0.0, TString name="stft");  
  ~STFT();  

  TCanvas* GetCanvas() {return canvas;}
  TH2D* GetHistogram() {return h2;}

  void SetLogz(bool isLogz=true) {this->isLogz=isLogz;}
  bool GetLogz() {return isLogz;}
  TString GetZtype() {return ztype;}
  void SetTitle(TString title) {h2->SetTitle(title);this->title=title;}
  TString GetTitle() {return title;}
  void SetPalette(int paletteId=1) {this->paletteId=paletteId;}
  bool GetPaletteId() {return paletteId;}

  void Draw(double t1=0.0, double t2=0.0, double f1=0.0, double f2=0.0, double z1=0.0, double z2=0.0,
            int dpaletteId = DUMMY_PALETTE_ID, Option_t* option = "colfz");
  void Print(TString pname);

private:

  void SetPlotStyle(int paletteId = 1);

  TCanvas* canvas;
  TH2D* h2;

  bool isLogz;
  TString title;
  TString name;
  int paletteId;
  TString ztype;

  double* window;
};  

} // end namespace

#endif
