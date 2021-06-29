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


// macro used to generate the Overlap Catalog
{    

  #define NMAX_RES 32

  // get low resolution level
  int cwb_xtalk_low_level=0;
  if(gSystem->Getenv("CWB_XTALK_LOW_LEVEL")==NULL) {
    cout << "cwb_xtalk.C - Error : environment CWB_XTALK_LOW_LEVEL is not defined!!!" << endl;exit(1);
  } else {
    if(TString(gSystem->Getenv("CWB_XTALK_LOW_LEVEL")).IsDigit()) {
      cwb_xtalk_low_level=TString(gSystem->Getenv("CWB_XTALK_LOW_LEVEL")).Atoi();
    } else {
      cout << "cwb_xtalk.C - Error : environment CWB_XTALK_LOW_LEVEL is not defined!!!" << endl;exit(1);
    }
  }

  // get high resolution level
  int cwb_xtalk_high_level=0;
  if(gSystem->Getenv("CWB_XTALK_HIGH_LEVEL")==NULL) {
    cout << "cwb_xtalk.C - Error : environment CWB_XTALK_HIGH_LEVEL is not defined!!!" << endl;exit(1);
  } else {
    if(TString(gSystem->Getenv("CWB_XTALK_HIGH_LEVEL")).IsDigit()) {
      cwb_xtalk_high_level=TString(gSystem->Getenv("CWB_XTALK_HIGH_LEVEL")).Atoi();
    } else {
      cout << "cwb_xtalk.C - Error : environment CWB_XTALK_HIGH_LEVEL is not defined!!!" << endl;exit(1);
    }
  }

  int nRes = cwb_xtalk_high_level-cwb_xtalk_low_level+1; 
  if(nRes<=0) {
    cout << "cwb_xtalk.C - Error : low res Level must be < high res level" << endl;exit(1);
  }
  if(nRes>NMAX_RES) {
    cout << "cwb_xtalk.C - Error : number of max resolutions is : " << NMAX_RES << endl;exit(1);
  }

  // get iNu : defines the sharpness of the 'edge' of the basis function in Fourier domain 
  int cwb_xtalk_iNu=4;
  if(gSystem->Getenv("CWB_XTALK_INU")!=NULL) {
    if(TString(gSystem->Getenv("CWB_XTALK_INU")).IsDigit()) {
      cwb_xtalk_iNu=TString(gSystem->Getenv("CWB_XTALK_INU")).Atoi();
    } else {
      cout << "cwb_xtalk.C - Error : environment CWB_XTALK_INU is not defined!!!" << endl;exit(1);
    }
  }

  // get precison : defines filter length by truncation error quantified by P = -log10(1 - norm_of_filter)
  int cwb_xtalk_precision=10;
  if(gSystem->Getenv("CWB_XTALK_PRECISION")!=NULL) {
    if(TString(gSystem->Getenv("CWB_XTALK_HIGH_LEVEL")).IsDigit()) {
      cwb_xtalk_precision=TString(gSystem->Getenv("CWB_XTALK_PRECISION")).Atoi();
    } else {
      cout << "cwb_xtalk.C - Error : environment CWB_XTALK_PRECISION is not defined!!!" << endl;exit(1);
    }
  }

  // define wdm
  cout << "cwb_xtalk.C - define wdm ..." << endl;
  WDM<double>* wdm[NMAX_RES];
  for(int level=cwb_xtalk_low_level; level<=cwb_xtalk_high_level; level++) {
    int layers = level>0 ? 1<<level : 0;
    wdm[level-cwb_xtalk_low_level] = new WDM<double>(layers, layers, cwb_xtalk_iNu, cwb_xtalk_precision);
  }

  // generate xtalk catalog
  cout << "cwb_xtalk.C - generate xtalk catalog : be patient, it takes a while ..." << endl;
  monster x(wdm, nRes);

  // write xtalk catalog
  cout << "cwb_xtalk.C - write xtalk catalog ..." << endl;
  char fName[1024];
  sprintf(fName,"OverlapCatalog-ilLev%d-hLev%d-iNu%d-P%d.xbin",
          cwb_xtalk_low_level,cwb_xtalk_high_level,cwb_xtalk_iNu,cwb_xtalk_precision);
  x.write(fName);
  cout << "OverlapCatalog Name : " << fName << endl;
  gSystem->Exit(1);
} 
