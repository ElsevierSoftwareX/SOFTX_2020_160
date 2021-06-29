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


#include "gwseries.hh"

ClassImp(gWSeries<double>)


//______________________________________________________________________________
/* Begin_Html
<center><h2>gWSeries class</h2></center>
This class gWSeries is derived from the WSeries class
It uses the methods implemented in gwavearray and the scalogram plot

<p>
<h3><a name="example">Example</a></h3>
<p>
The macro <a href="./tutorials/gwat/DrawGWSeries.C.html">DrawGWSeries.C</a> is an example which shown how to use the gWSeries class.<br>
The pictures below gives the macro output plots.<br>
<p>

End_Html

Begin_Macro
DrawGWSeries.C
End_Macro */


using namespace std;


//______________________________________________________________________________
// destructor
template<class DataType_t>
gWSeries<DataType_t>::~gWSeries() 
{
  if(gw!=NULL) delete gw;

  if(wsextern) {
    this->data = NULL;
    this->pWavelet = NULL;
  }

//  if(wts!=NULL)  delete wts;
}

//______________________________________________________________________________
template<class DataType_t>
gWSeries<DataType_t>::gWSeries(const WSeries<DataType_t> &w) :
WSeries<DataType_t>(w) {this->Init();wsextern=false;}

//______________________________________________________________________________
template<class DataType_t>
gWSeries<DataType_t>::gWSeries(const WSeries<DataType_t> *w) 
{
  wsextern = true;

  this->data  = w->data;
  this->Rate  = w->rate();
  this->Size  = w->size();
  this->Start = w->start();
  this->Stop  = w->start() + w->size()/w->rate();
  this->Slice = std::slice(0,w->size(),1);

  this->pWavelet = w->pWavelet;
  this->bpp = w->getbpp();
  this->wRate = w->wRate;
  this->f_low = w->getlow();
  this->f_high = w->gethigh();
  this->w_mode = w->w_mode;

  this->Init();
}

//______________________________________________________________________________
template<class DataType_t>
void gWSeries<DataType_t>::Init() {

  wavearray<DataType_t>* p = this;
  gw = new gwavearray<DataType_t>(p);

  gRandom->SetSeed(0);
  rnID = int(gRandom->Rndm(13)*1.e9);   // random name ID

//  wts=NULL;
}

//______________________________________________________________________________
template<class DataType_t>                                                      
watplot* gWSeries<DataType_t>::DrawSG(WSeries<DataType_t>* x, TString options) {
//                                                                                    
// Draw waveform scalogram                                                          
//                                                                                    
// Input:                                                                             
//        options - draw options (same as TH2D)                                       
//                                                                                    
// return pointer to watplot object                                                 
//                                                                                    

  double tStart,tStop;
  if(options.Contains("FULL")) {
    options.ReplaceAll("FULL","");
    tStart=this->start();         
    tStop=this->start()+this->size()/this->rate();
  } else {                                        
    GetTimeRange(tStart, tStop);                  
  }                                               

  watplot* wts = GetWATPLOT();
  if(wts!=NULL) delete wts;  

  char name[32];sprintf(name,"sg-%d",rnID);
  wts = new watplot(const_cast<char*>(name));
  WSeries<DataType_t>* ws = x!=NULL ? x : this;
  wts->plot(ws, 2, tStart, tStop,const_cast<char*>("COLZ"));

  double fLow  = 32.;
  double fHigh = this->rate()/2.;
  wts->hist2D->GetYaxis()->SetRangeUser(fLow,fHigh);

  gw->SetWATPLOT(wts);

  return wts;
}             

// instantiations

#define CLASS_INSTANTIATION(class_) template class gWSeries< class_ >;

//CLASS_INSTANTIATION(short)
//CLASS_INSTANTIATION(int)
//CLASS_INSTANTIATION(unsigned int)
//CLASS_INSTANTIATION(long)
//CLASS_INSTANTIATION(long long)
//CLASS_INSTANTIATION(float)
CLASS_INSTANTIATION(double)

