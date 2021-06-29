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
 * File name:    gskymap.hh
 * Author:       Gabriele Vedovato (vedovato@lnl.infn.it)
 **********************************************************/


#ifndef GSKYMAP_HH
#define GSKYMAP_HH

#include "TCanvas.h"
#include "TH2F.h"
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
#include "TVector3.h"
#include "TRotation.h"
#include "TPolyLine.h"
#include "Math/Rotation3D.h"
#include "Math/Vector3Dfwd.h"
#include "TSystem.h"
#include "TLatex.h"
#include "TRandom.h"
#include "TGaxis.h"

#include <iostream>  
#include <fstream>
#include <stdlib.h>
#include "skymap.hh"
#include "network.hh"
#include "detector.hh"
#include "wavecomplex.hh"
#include "watversion.hh"
#include "skycoord.hh"

#include "gwat.hh"

#define WorlMapCoastLine "WorlMapCoastLine.txt"
#define WM_ENTRIES 426345
#define DUMMY_PALETTE_ID 1000000000

#define GSKYMAP_INIT canvas(NULL), h2(NULL), changed(false), isGridx(true), colorGridx(kBlack),    \
                     isGridy(true), colorGridy(kBlack), isLogz(false), title(""),                  \
                     drawWorldMap(false), goff(0), wm_size(0), paletteId(1), gpsGalacticDisk(-1.), \
                     colorGalacticDisk(kBlack), zAxisTitle("")

using namespace ROOT::Math;

class gskymap : public skymap 
{

public:
 
  // Constructors
  gskymap() : skymap(),GSKYMAP_INIT {SetOptions();}
  gskymap(double sms,double t1=0.,double t2=180.,double p1=0.,double p2=360.) : 
         skymap(sms,t1,t2,p1,p2),GSKYMAP_INIT {SetOptions();}
  gskymap(char* ifile) : skymap(ifile),GSKYMAP_INIT {SetOptions();changed=true;}
  gskymap(int healpix_order) : skymap(healpix_order),GSKYMAP_INIT {SetOptions();}
  gskymap(TString ifile, TString name="gskymap") : 
         skymap(ifile,name),GSKYMAP_INIT {SetOptions();changed=true;}
  gskymap(const skymap& sm) : skymap(sm),GSKYMAP_INIT {SetOptions();changed=true;}
  virtual ~gskymap();  

  // operators
  gskymap& operator  = (const gskymap& sm) {
    skymap::operator=(sm);
    isGridx=sm.isGridx;
    colorGridx=sm.colorGridx;
    isGridy=sm.isGridy;
    colorGridy=sm.colorGridy;
    isLogz=sm.isLogz;
    title=sm.title;
    drawWorldMap=sm.drawWorldMap;
    resolution=sm.resolution;
    goff=sm.goff;
    wm_size=0;
    worldMapPath=sm.worldMapPath;
    paletteId=sm.paletteId;
    coordinate=sm.coordinate;
    projection=sm.projection;
    gpsGalacticDisk=sm.gpsGalacticDisk;
    colorGalacticDisk=sm.colorGalacticDisk;
    zAxisTitle=sm.zAxisTitle;
    name="";
    wtopx=sm.wtopx;
    wtopy=sm.wtopy;
    ww=sm.ww;
    wh=sm.wh;
    canvas = NULL;
    h2 = new TH2D(*sm.h2);
    return *this;
  }
  gskymap& operator  = (const skymap& sm) {skymap::operator=(sm);changed=true;return *this;}
  gskymap& operator += (const skymap& sm) {skymap::operator+=(sm);changed=true;return *this;}
  gskymap& operator -= (const skymap& sm) {skymap::operator-=(sm);changed=true;return *this;}
  gskymap& operator *= (const skymap& sm) {skymap::operator*=(sm);changed=true;return *this;}
  gskymap& operator /= (const skymap& sm) {skymap::operator/=(sm);changed=true;return *this;}
  gskymap& operator  = (const double a)   {skymap::operator =(a);changed=true;return *this;}
  gskymap& operator *= (const double a)   {skymap::operator*=(a);changed=true;return *this;}
  gskymap& operator += (const double a)   {skymap::operator+=(a);changed=true;return *this;}

  void set(size_t i, double a) {skymap::set(i,a);changed=true;}

  void SetOptions(TString projection = "hammer", TString coordinate = "Geographic",
                  double resolution = 1, bool goff = false);
  void SetOptions(const char* name, Int_t wtopx, Int_t wtopy, Int_t ww, Int_t wh) {
                  char sname[64];sprintf(sname,"%s-%d",name,int(gRandom->Rndm(13)*1.e9));
                  this->name=sname;this->wtopx=wtopx;this->wtopy=wtopy;this->ww=ww;this->wh=wh;
                  if(canvas!=NULL) Draw();} 

  TCanvas* GetCanvas() {return canvas;}
  TH2D* GetHistogram() {return h2;}

  void SetGridx(bool isGridx=true) {this->isGridx=isGridx;}
  bool GetGridx() {return isGridx;}
  void SetGridxColor(Color_t colorGridx = kBlack) {this->colorGridx=colorGridx;}
  bool GetGridxColor() {return colorGridx;}
  void SetGridy(bool isGridy=true) {this->isGridy=isGridy;}
  bool GetGridy() {return isGridy;}
  void SetGridyColor(Color_t colorGridy = kBlack) {this->colorGridy=colorGridy;}
  bool GetGridyColor() {return colorGridy;}
  void SetLogz(bool isLogz=true) {this->isLogz=isLogz;}
  bool GetLogz() {return isLogz;}
  void SetZaxisTitle(TString zAxisTitle) {h2->GetZaxis()->SetTitle(zAxisTitle);this->zAxisTitle=zAxisTitle;}
  TString GetZaxisTitle() {return zAxisTitle;}
  void SetTitle(TString title) {h2->SetTitle(title);this->title=title;}
  TString GetTitle() {return title;}
  void SetWorldMap(bool drawWorldMap=true) {this->drawWorldMap=drawWorldMap;}
  bool GetWorldMap() {return drawWorldMap;}
  void SetWorldMapPath(TString worldMapPath) {this->worldMapPath=worldMapPath;}
  TString GetWorldMapPath() {return worldMapPath;}
  double GetResolution() {return resolution;}
  void SetPalette(int paletteId=kBird) {this->paletteId=paletteId;}
  int GetPaletteId() {return paletteId;}
  TString GetProjection() {return projection;}
  TString GetCoordinate() {return coordinate;}
  void SetGalacticDisk(double gpsGalacticDisk=0.0) {this->gpsGalacticDisk=gpsGalacticDisk;}
  double GetGalacticDisk() {return gpsGalacticDisk;}
  void SetGalacticDiskColor(Color_t colorGalacticDisk = kBlack) {this->colorGalacticDisk=colorGalacticDisk;}
  bool GetGalacticDiskColor() {return colorGalacticDisk;}

  void FillData(int size, double* phi, double* theta, double* binc);
  void FillData(char* fname);

  void Plot();							// *MENU*
  void Draw(int dpaletteId = 1, Option_t* option = "colfz");
  void Print(TString pname);					// *MENU*

  void ProjectHammer(Double_t l,Double_t b,Double_t &Al,Double_t &Ab);
  void ProjectSinusoidal(Double_t l, Double_t b, Double_t &Al, Double_t &Ab);
  void ProjectParabolic(Double_t l, Double_t b, Double_t &Al, Double_t &Ab);

  void DrawMarker(double phi, double theta, int marker, Size_t msize = 1, Color_t tcolor = 1);
  void DrawText(double phi, double theta, TString text, double tsize = 0.04, Color_t tcolor = 1);

  void DrawMarker(double ra, double dec, double gps, int marker, Size_t msize = 1, Color_t tcolor = 1);
  void DrawText(double ra, double dec, double gps, TString text, double tsize = 0.04, Color_t tcolor = 1);

  virtual void Browse(TBrowser *b) {Plot();}

  void ClearAxisLabel();
  void ClearGalacticDisk();
  void ClearWorldMap();
  void ClearGridx();
  void ClearGridy();

  //: dump to ascii file
  void DumpSkyMap(char* fname);

  //: dump skymap into root file
  void DumpObject(const char* file, const char* name="gskymap");
  //: load gnetwork object from root file
  void LoadObject(const char* file, const char* name="gskymap");

  friend class gnetwork;

private:

  void SetPlotStyle(int paletteId = 1);
  int  ReadWorlMapCoastLine(double*& wm_lon, double*& wm_lat);
  void HeapSort( double* data, double lenght);
  void ReverseXAxis(TH2D *h2);
  void FillData();
  void CreateCanvas(); 

  TCanvas* canvas;		//!`
  TH2D* h2;			//!

  bool changed;

  bool isGridx;
  Color_t colorGridx;
  bool isGridy;
  Color_t colorGridy;
  bool isLogz;
  TString title;
  bool drawWorldMap;
  double resolution;
  bool goff;
  int wm_size;
  TString worldMapPath;
  int paletteId;
  TString coordinate;
  TString projection;
  double gpsGalacticDisk;
  Color_t colorGalacticDisk;
  TString zAxisTitle;

  // canvas parameters
  TString name;
  Int_t wtopx, wtopy, ww, wh;

  std::vector<TPolyLine*> gridxL;	//!
  std::vector<TPolyLine*> gridyL;	//!
  std::vector<TPolyLine*> gdL;		//!
  std::vector<TMarker*>   wmM;		//!
  std::vector<TText*>     axisT;	//!

  double *wm_lon;			//!
  double *wm_lat;			//!

  ClassDef(gskymap,1)
};  

#endif
