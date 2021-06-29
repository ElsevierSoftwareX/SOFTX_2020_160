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
                          Bicoherence.hh  -  description
                             -------------------
    begin                : Wen Dec 21 2011
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

#ifndef CWBBICOHERENCE_H
#define CWBBICOHERENCE_H


/**
  *@author Gabriele Vedovato
  */

#include "Window.hh"
#include "TGraph.h"
#include <TH2F.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TPaletteAxis.h>
#include <TCutG.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TString.h>
#include "wavearray.hh"
#include <time.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <string>
#include <iostream>
#include <fstream>
#include <complex>
#include <vector>
#include "Toolbox.hh"


struct bico {
  float x;
  float y;
  float c;
};

namespace CWB {

class Bicoherence {
public: 
  Bicoherence(TString ch1name, TString ch2name, double srate, int bsize,
              int x_num_slices=1, int y_num_slices=1, 
              int x_slice_index=0, int y_slice_index=0,
              int order=4, char* window_type = const_cast<char*>("welch"));	
  ~Bicoherence();

  bool  MakeBicoherence(wavearray<double> x, wavearray<double> y);
  vector<bico> GetBicoherence(float threshold, int rebin=1);
//  void  GetBicoherence();
  void  DrawBicoherence(int rebin=1, int graph_id=0, TString ofname="", bool batch=false);
  int   GetAverages() {return bic_averages;}
  void  Reset();
  int   readSegList(char* seglist, double shift=0, bool invert=false, bool c4=false);
  bool  segListCheck(double start, double stop);

private:

  TCanvas* canvas;

  int nbl;
  double bstart;
  double blength;
  int bsize;
  int xbsize;
  int ybsize;
  complex<double>** B;
  char* chnames[2];
  double* SX;
  double* SY;
  double* S[3];
  double* f;
  complex<double>* xout;
  complex<double>* yout;
  double* window;
  int bic_averages;
  int order;
  int x_slice_index;
  int y_slice_index;
  double df;
  double fxmin,fxmax;
  double fymin,fymax;
  int pgraph_id;
  double srate;
  int canvas_number;
  TH2F* blhist;
  TGraph* sgraph[3];
  vector<waveSegment> segList;

};

} // end namespace

#endif
