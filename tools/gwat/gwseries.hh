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
 * File name:    gwseries.hh
 * Author:       Gabriele Vedovato (vedovato@lnl.infn.it)
 **********************************************************/

#ifndef GWSERIES_HH
#define GWSERIES_HH

#ifndef WSERIES_HH
#include "wseries.hh"
#endif

#include "gwavearray.hh"

template<class DataType_t> 
class gWSeries : public WSeries<DataType_t>
{
public:

  gWSeries() : WSeries<DataType_t>() {this->Init();wsextern=false;}    // Default constructor
  gWSeries(const WSeries<DataType_t>&);                                // copy Constructor
  gWSeries(const WSeries<DataType_t>*);        

  void Forward(int n = -1) 
              {WSeries<DataType_t>::Forward(n);this->Init();}
  void Forward(wavearray<DataType_t> &w, int n = -1) 
              {WSeries<DataType_t>::Forward(w,n);this->Init();}
  void Forward(wavearray<DataType_t> &w, Wavelet &s, int n = -1) 
              {WSeries<DataType_t>::Forward(w,s,n);this->Init();}
                       
  virtual ~gWSeries();                                                 // Destructor

  void   Draw(GWAT_DRAW type=GWAT_TIME, TString options="ALP", Color_t color=kBlack)
              {if(type==GWAT_SG) {DrawSG(NULL,options);return;}
               if(this->getLevel()) return; else return gw->Draw(type,options,color);}

  watplot*   DrawTime(TString options="ALP", Color_t color=kBlack)  // *MENU*
              {if(this->getLevel()) return NULL; else return gw->DrawTime(options,color);}
  watplot*   DrawFFT(TString options="ALP", Color_t color=kBlack)   // *MENU*
              {if(this->getLevel()) return NULL; else return gw->DrawFFT(options,color);}
  CWB::STFT* DrawTF(TString options="")                             // *MENU*
              {if(this->getLevel()) return NULL; else return gw->DrawTF(options);}
  watplot*   DrawSG(TString options="")                             // *MENU*
              {return DrawSG(NULL,options);} 

  void   Draw(wavearray<DataType_t>* x, GWAT_DRAW type=GWAT_TIME,
              TString options="ALP", Color_t color=kBlack)
              {if(this->getLevel()) return; else return gw->Draw(x,type,options,color);}

  void   Draw(WSeries<DataType_t>* x, TString options="")
              {DrawSG(x,options);return;}

  watplot*   DrawTime(wavearray<DataType_t>* x, TString options="ALP", Color_t color=kBlack)
              {if(this->getLevel()) return NULL; else return gw->DrawTime(x,options,color);}
  watplot*   DrawFFT(wavearray<DataType_t>* x, TString options="ALP", Color_t color=kBlack)   
              {if(this->getLevel()) return NULL; else return gw->DrawFFT(x,options,color);}
  CWB::STFT* DrawTF(wavearray<DataType_t>* x, TString options="")                             
              {if(this->getLevel()) return NULL; else return gw->DrawTF(x,options);}
  watplot*   DrawSG(WSeries<DataType_t>* x, TString options=""); 

  double GetTimeRange(double& tMin, double& tMax) {return gw->GetTimeRange(tMin,tMax);}
  double GetCentralTime()                         {return gw->GetCentralTime();}
  void   TimeShift(double tShift=0.)              {return gw->TimeShift(tShift);}
  void   PhaseShift(double pShift=0.)             {return gw->PhaseShift(pShift);}

  CWB::STFT* GetSTFT()  {return gw->GetSTFT();}
  watplot* GetWATPLOT() {return gw->GetWATPLOT();}

  // print wseries parameters
  virtual void Browse(TBrowser *b) {Draw();}

private:

//  watplot*      wts;  //!

  gwavearray<DataType_t>* gw;

  void Init();

  int  rnID;
  bool wsextern;	      // external WSeries

  ClassDef(gWSeries,1)
};

#endif
