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
 * File name:    gwavearray.hh
 * Author:       Gabriele Vedovato (vedovato@lnl.infn.it)
 **********************************************************/

#ifndef GWAVEARRAY_HH
#define GWAVEARRAY_HH

#ifndef WAVEARRAY_HH
#include "wavearray.hh"
#endif

#include "gwat.hh"
#include "watplot.hh"
#include "time.hh"
#include "STFT.hh"
#include "TComplex.h"


template<class DataType_t> 
class gwavearray : public wavearray<DataType_t>
{
public:

  gwavearray() : wavearray<DataType_t>() {Init();wextern=false;}   // Default constructor
  gwavearray(const wavearray<DataType_t>&);                        // copy Constructor
  gwavearray(const wavearray<DataType_t>*);        

  virtual ~gwavearray();                                           // Destructor

  void   Draw(GWAT_DRAW type=GWAT_TIME,
              TString options="ALP", Color_t color=kBlack);

  void   PlotTime(double edge=0) {
           SetTimeRange(this->start()+edge,this->start()+this->size()/this->rate()-edge);
           DrawTime("CUSTOM ALP");} 		// *MENU*
  void   PlotFFT(double edge=0) {
           SetTimeRange(this->start()+edge,this->start()+this->size()/this->rate()-edge);
           DrawFFT("CUSTOM ALP");}  	     	// *MENU*
  void   PlotPSD(double edge=0) {
           SetTimeRange(this->start()+edge,this->start()+this->size()/this->rate()-edge);
           DrawFFT("CUSTOM ALP PSD");} 		// *MENU*
  void   PlotTF(double edge=0)  {DrawTF();;}  					   

  watplot*   DrawTime(TString options="ALP", Color_t color=kBlack);  
  watplot*   DrawFFT(TString options="ALP", Color_t color=kBlack);   
  CWB::STFT* DrawTF(TString options="");                             

  void   Draw(wavearray<DataType_t>* x, GWAT_DRAW type=GWAT_TIME,
              TString options="ALP", Color_t color=kBlack);

  watplot*   DrawTime(wavearray<DataType_t>* x, TString options="ALP", Color_t color=kBlack);  
  watplot*   DrawFFT(wavearray<DataType_t>* x, TString options="ALP", Color_t color=kBlack);   
  CWB::STFT* DrawTF(wavearray<DataType_t>* x, TString options="");                             

  void SetTimeRange(double tMin, double tMax) {
    if(tMin>=tMax) return;
    double tStart=this->start();
    double tStop=this->start()+this->size()/this->rate();
    this->tMin = tMin<tStart||tMin>tStop ? tStart : tMin;
    this->tMax = tMax<tStart||tMax>tStop ? tStop  : tMax;
  }
  double GetTimeRange(double& tMin, double& tMax, double efraction = 0.9999999);
  double GetCentralTime();
  void   TimeShift(double tShift=0.);
  void   PhaseShift(double pShift=0.);

  CWB::STFT* GetSTFT()  {return this->stft;}
  watplot* GetWATPLOT() {return this->pts;}

  // print wseries parameters
  void PrintParameters() {wavearray<DataType_t>::print();}         // *MENU*
  virtual void Browse(TBrowser *b) {PlotTime();}
  void SetWATPLOT(watplot* pts) {this->pts=pts;}

  void SetComment(TString comment) {this->comment=comment;}
  TString GetComment() {return comment;}
  void PrintComment() {cout << comment << endl;}		   // *MENU*

  void DumpToFile(char* fname);   			   	   // *MENU*

private:

  CWB::STFT*    stft; //!
  watplot*      pts;  //!

  void Init();

  bool wextern;	      	// external wavearray

  TString comment;	// store comment

  double tMin;
  double tMax;
  double tSave;         // start time : used to manage DrawTime with option "SAME"

  ClassDef(gwavearray,3)
};

#endif
