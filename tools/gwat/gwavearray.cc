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


#include "gwavearray.hh"

ClassImp(gwavearray<double>)

//______________________________________________________________________________
/* Begin_Html
<center><h2>gwavearray class</h2></center>
This class gwavearray is derived from the wavearray class
It is used to produce time,freq and time-freq plots of a signal stored in a wavearray structure.

<p>
<h3><a name="example">Example</a></h3>
<p>
The macro <a href="./tutorials/gwat/DrawGWaveArray.C.html">DrawGWaveArray.C</a> is an example which shown how to use the gWaveArray class.<br>
The pictures below gives the macro output plots.<br>
<p>

End_Html

Begin_Macro
DrawGWaveArray.C
End_Macro */

using namespace std;


//______________________________________________________________________________
// destructor
template<class DataType_t>
gwavearray<DataType_t>::~gwavearray() 
{
  if(wextern) this->data = NULL;

  if(stft!=NULL) delete stft;
  if(pts!=NULL)  delete pts;
}

//______________________________________________________________________________
template<class DataType_t>
gwavearray<DataType_t>::gwavearray(const wavearray<DataType_t> &w) :
wavearray<DataType_t>(w) {Init();wextern=false;}

//______________________________________________________________________________
template<class DataType_t>
gwavearray<DataType_t>::gwavearray(const wavearray<DataType_t> *w) 
{
  Init();
  wextern = true;

  this->data  = w->data;
  this->Rate  = w->rate();
  this->Size  = w->size();
  this->Start = w->start();
  this->Stop  = w->start() + w->size()/w->rate();
  this->Slice = std::slice(0,w->size(),1);
}

//______________________________________________________________________________
template<class DataType_t>
void gwavearray<DataType_t>::Init() {

  gRandom->SetSeed(0);
  tMin=0;
  tMax=0;

  stft=NULL;pts=NULL;
}

//______________________________________________________________________________
template<class DataType_t>
void gwavearray<DataType_t>::Draw(GWAT_DRAW type, TString options, Color_t color) {
//
// Draw waveform in time/frequency domain
//
// Input:        
//        type    - GWAT_TIME, GWAT_FFT, GWAT_TF
//        options - graphic options
//        color   - option for GWAT_TIME, GWAT_FFT
//

  return Draw(NULL,type,options,color);
}

//______________________________________________________________________________
template<class DataType_t>
void gwavearray<DataType_t>::Draw(wavearray<DataType_t>* x, GWAT_DRAW type, TString options, Color_t color) {
//
// Draw waveform in time/frequency domain
//
// Input:        
//      wavearray - data array
//      type      - GWAT_TIME, GWAT_FFT, GWAT_TF
//      options   - graphic options
//      color     - option for GWAT_TIME, GWAT_FFT
//

  if(type==GWAT_TIME) DrawTime(x,options,color);
  if(type==GWAT_FFT)  DrawFFT(x,options,color);
  if(type==GWAT_TF)   DrawTF(x,options);
  if(type==GWAT_SG)   cout << "gwavearray<DataType_t>::Draw - Error : draw type " 
                           << "GWAT_SG" << " not allowed" << endl;
}

//______________________________________________________________________________
template<class DataType_t>
watplot* gwavearray<DataType_t>::DrawTime(TString options, Color_t color) {
//
// Draw waveform in time domain
//
// Input:        
//        options - draw options (same as TGraph)
//                  if contains FULL the full interval is showed
//
// return pointer to watplot object
//

   return DrawTime(NULL,options,color);
}

//______________________________________________________________________________
template<class DataType_t>
watplot* gwavearray<DataType_t>::DrawTime(wavearray<DataType_t>* x, TString options, Color_t color) {
//
// Draw waveform in time domain
//
// Input:        
//      wavearray - data array
//        options - draw options (same as TGraph)
//                  if contains FULL the full interval is showed
//
// return pointer to watplot object
//

  if(this->rate()<=0) {
    cout << "gwavearray::DrawTime : Error - rate must be >0" << endl;
    exit(1);
  }

  options.ToUpper();
  if(stft!=NULL) delete stft;
  if(options.Contains("SAME")&&(pts!=NULL)) {
  } else {
    if(pts!=NULL) delete pts;
    char name[32];sprintf(name,"TIME-gID:%d",int(gRandom->Rndm(13)*1.e9));
    pts = new watplot(const_cast<char*>(name),200,20,800,500);
  }
  double tStart,tStop;
  if(options.Contains("FULL")) {
    options.ReplaceAll("FULL","");
    tStart=this->start();
    tStop=this->start()+this->size()/this->rate();
    tStop-=this->start();tStart-=this->start();
  } else if(options.Contains("CUSTOM")) {
    options.ReplaceAll("CUSTOM","");
    tStart=this->start();
    tStop=this->start()+this->size()/this->rate();
    tStart = tMin>tStart&&tMin<tStop ? tMin : tStart;
    tStop  = tMax>tStart&&tMax<tStop ? tMax : tStop;
    tStop-=this->start();tStart-=this->start();
  } else {
    GetTimeRange(tStart, tStop);
  }
  wavearray<double> tf;
  if(x!=NULL) waveAssign(tf,*x);
  else        waveAssign(tf,*this);
  char title[256];
  sprintf(title,"START : %10.3f (gps) - LENGHT : %g (sec) - RATE : %g (hz)",
          tf.start()+tStart,tStop-tStart,tf.rate());
  if(!options.Contains("SAME")) {
    tSave=tf.start();		// save start time of first plot
    tf.start(0);
  } else {
    tf.start(tf.start()-tSave);
  }
  if(!options.Contains("SAME")) options=options+TString(" ALP");
  options.ReplaceAll(" ","");
  pts->plot(tf, const_cast<char*>(options.Data()), color, tStart, tStop);
  pts->graph[pts->graph.size()-1]->SetTitle(title);

  return pts;
}

//______________________________________________________________________________
template<class DataType_t>
watplot* gwavearray<DataType_t>::DrawFFT(TString options, Color_t color) {
//
// Draw waveform spectrum
//
// Input:        
//      wavearray - data array
//        options - draw options (same as TGraph)
//                  if contains FULL the full interval is showed
//                  if contains NOLOGX/NOLOGY X/Y axis are linear
//
// return pointer to watplot object
//

   return DrawFFT(NULL,options,color);
}

//______________________________________________________________________________
template<class DataType_t>
watplot* gwavearray<DataType_t>::DrawFFT(wavearray<DataType_t>* x, TString options, Color_t color) {
//
// Draw waveform spectrum
//
// Input: 
//      wavearray - data array
//        options - draw options (same as TGraph)
//                  if contains FULL the full interval is showed
//                  if contains NOLOGX/NOLOGY X/Y axis are linear
//
// return pointer to watplot object
//

  options.ToUpper();
  if(stft!=NULL) delete stft;
  if(options.Contains("SAME")&&(pts!=NULL)) {
  } else {
    if(pts!=NULL) delete pts;
    char name[32];sprintf(name,"FFT-gID:%d",int(gRandom->Rndm(13)*1.e9));
    pts = new watplot(const_cast<char*>(name),200,20,800,500);
  }

  double fLow  = 4.;
  double fHigh = this->rate()/2.;
  double tStart,tStop;
  bool psd = options.Contains("PSD") ? true : false;
  options.ReplaceAll("PSD","");
  if(options.Contains("FULL")) {
    options.ReplaceAll("FULL","");
    tStart=this->start();
    tStop=this->start()+this->size()/this->rate();
    tStop-=this->start();tStart-=this->start();
  } else if(options.Contains("CUSTOM")) {
    options.ReplaceAll("CUSTOM","");
    tStart=this->start();
    tStop=this->start()+this->size()/this->rate();
    tStart = tMin>tStart&&tMin<tStop ? tMin : tStart;
    tStop  = tMax>tStart&&tMax<tStop ? tMax : tStop;
    tStop-=this->start();tStart-=this->start();
  } else {
    GetTimeRange(tStart, tStop);
  }
  bool logx = true;
  if(options.Contains("NOLOGX")) {logx=false;options.ReplaceAll("NOLOGX","");}
  bool logy = true;
  if(options.Contains("NOLOGY")) {logy=false;options.ReplaceAll("NOLOGY","");}
  wavearray<double> tf;
  if(x!=NULL) waveAssign(tf,*x);
  else        waveAssign(tf,*this);
  // use only the selected range
  int nb=tStart*tf.rate();
  int ne=tStop*tf.rate();
  for(int i=nb;i<ne;i++) tf[i-nb] = tf[i];
  tf.resize(ne-nb);
  char title[256];
  sprintf(title,"START : %10.3f (gps) - LENGHT : %g (sec) - RATE : %g (hz)",
          tf.start()+tStart,tStop-tStart,tf.rate());
  tf.start(0);
  if(!options.Contains("SAME")) options=options+TString(" ALP");
  options.ReplaceAll(" ","");
  pts->plot(tf, const_cast<char*>(options.Data()), color, 0, tStop-tStart, true, fLow, fHigh, psd, 8);
  pts->graph[pts->graph.size()-1]->SetTitle(title);
  if(logx) pts->canvas->SetLogx();
  if(logy) pts->canvas->SetLogy();

  return pts;
}

//______________________________________________________________________________
template<class DataType_t>
CWB::STFT* gwavearray<DataType_t>::DrawTF(TString options) {
//
// Draw waveform spectrogram
//
// Input: 
//      wavearray - data array
//        options - draw options (same as TH2D)
//
// return pointer to CWB::STFT object
//

   return DrawTF(NULL,options);
}

//______________________________________________________________________________
template<class DataType_t>
CWB::STFT* gwavearray<DataType_t>::DrawTF(wavearray<DataType_t>* x, TString options) {
//
// Draw waveform spectrogram
//
// Input: 
//        options - draw options (same as TH2D)
//
// return pointer to CWB::STFT object
//

  int nfact=4;
  int nfft=nfact*512;
  int noverlap=nfft-10;
  double fparm=nfact*6;

  double tStart,tStop;
  if(options.Contains("FULL")) {
    options.ReplaceAll("FULL","");
    tStart=this->start();
    tStop=this->start()+this->size()/this->rate();
    tStop-=this->start();tStart-=this->start();
  } else if(options.Contains("CUSTOM")) {
    options.ReplaceAll("CUSTOM","");
    tStart=this->start();
    tStop=this->start()+this->size()/this->rate();
    tStart = tMin>tStart&&tMin<tStop ? tMin : tStart;
    tStop  = tMax>tStart&&tMax<tStop ? tMax : tStop;
    tStop-=this->start();tStart-=this->start();
  } else {
    GetTimeRange(tStart, tStop);
  }

  if(stft!=NULL) delete stft;
  if(pts!=NULL) delete pts;

  wavearray<double> tf;
  if(x!=NULL) waveAssign(tf,*x);
  else        waveAssign(tf,*this);
  tf.start(0);
  char name[32];sprintf(name,"TF-gID:%d",int(gRandom->Rndm(13)*1.e9));
  stft = new CWB::STFT(tf,nfft,noverlap,"amplitude","gauss",fparm,name);

  double fLow  = 4.;
  double fHigh = this->rate()/2.;
  stft->Draw(tStart,tStop,fLow,fHigh,0,0,1);
  stft->GetCanvas()->SetLogy(true);

  return stft;
}

//______________________________________________________________________________
template<class DataType_t>
double gwavearray<DataType_t>::GetTimeRange(double& tMin, double& tMax, double efraction) {
//
// Get mdc time interval
//
// Input:  efraction - energy fraction used to evaluate the range which contains efraction
//                     default = 0.9999999
//
// Output: tMin  - start time interval
//         tMax  - stop  time interval
//

  double a,b;

  double T = GetCentralTime();

  double E=0.;
  int size=(int)this->size();
  double rate=this->rate();
  for(int j=0;j<size;j++) E += this->data[j]*this->data[j];

  int OS = int(T*rate);
  int M  = OS<size/2 ? OS : size-OS;

  int I=0;
  double sum = this->data[OS]*this->data[OS];
  for(int j=1; j<int(M); j++) {
    a = this->data[OS-j];
    b = this->data[OS+j];
    sum += a*a+b*b;
    if(sum/E > efraction) break;
    I=j;
  }
  I++;
  tMin = tMin>=0 ? double((OS-I)/rate) : 0.;
  tMax = tMax<size/rate ? double((OS+I)/rate) : size/rate;

  if(sum/E<efraction) {tMin=0.;tMax=size/rate;}

  return T;
}

//______________________________________________________________________________
template<class DataType_t>
double gwavearray<DataType_t>::GetCentralTime() {
//
// Get mdc central time
//
// Input: 
//
// Retun central time
//

  double a;
  double E=0.,T=0.;
  int size=(int)this->size();
  double rate=this->rate();
  for(int j=0;j<size;j++) {
    a = this->data[j];
    T += a*a*j/rate;                   // central time
    E += a*a;                          // energy
  }
  T = E>0 ? T/E : 0.5*size/rate;

  return T;
}

//______________________________________________________________________________
template<class DataType_t>
void gwavearray<DataType_t>::TimeShift(double tShift) {
//
// apply time shift
//
//
// Input: 
//        time   - time shift (sec)
//

  if(tShift==0) return;

  // search begin,end of non zero data
  int ibeg=0; int iend=0;
  for(int i=0;i<(int)this->size();i++) {
    if(this->data[i]!=0 && ibeg==0) ibeg=i;
    if(this->data[i]!=0) iend=i;
  }
  int ilen=iend-ibeg+1;
  // create temporary array for FFTW & add scratch buffer + tShift
  int isize = 2*ilen+2*fabs(tShift)*this->rate();
  isize = isize + (isize%4 ? 4 - isize%4 : 0); // force to be multiple of 4
  wavearray<double> w(isize);
  w.rate(this->rate()); w=0;
  // copy this->data data !=0 in the middle of w array & set this->data=0
  for(int i=0;i<ilen;i++) {w[i+isize/4]=this->data[ibeg+i];this->data[ibeg+i]=0;}

  double pi = TMath::Pi();
  // apply time shift to waveform vector
  w.FFTW(1);
  TComplex C;
  double df = w.rate()/w.size();
  //cout << "tShift : " << tShift << endl;
  for (int ii=0;ii<(int)w.size()/2;ii++) {
    TComplex X(w[2*ii],w[2*ii+1]);
    X=X*C.Exp(TComplex(0.,-2*pi*ii*df*tShift));  // Time Shift
    w[2*ii]=X.Re();
    w[2*ii+1]=X.Im();
  }
  w.FFTW(-1);

  // copy shifted data to input this->data array
  for(int i=0;i<(int)w.size();i++) {
    int j=ibeg-isize/4+i;
    if((j>=0)&&(j<(int)this->size())) this->data[j]=w[i];
  }

  return;
}

//______________________________________________________________________________
template<class DataType_t>
void gwavearray<DataType_t>::PhaseShift(double pShift) {
//
// apply phase shift
//
// Input: 
//        pShift - phase shift (degrees)
//

  if(pShift==0) return;

  // search begin,end of non zero data
  int ibeg=0; int iend=0;
  for(int i=0;i<(int)this->size();i++) {
    if(this->data[i]!=0 && ibeg==0) ibeg=i;
    if(this->data[i]!=0) iend=i;
  }
  int ilen=iend-ibeg+1;
  // create temporary array for FFTW & add scratch buffer + tShift
  int isize = 2*ilen;
  isize = isize + (isize%4 ? 4 - isize%4 : 0); // force to be multiple of 4
  wavearray<double> w(isize);
  w.rate(this->rate()); w=0;
  // copy this->data data !=0 in the middle of w array & set this->data=0
  for(int i=0;i<ilen;i++) {w[i+isize/4]=this->data[ibeg+i];this->data[ibeg+i]=0;}

  // apply phase shift to waveform vector
  w.FFTW(1);
  TComplex C;
  //cout << "pShift : " << pShift << endl;
  pShift*=TMath::Pi()/180.;
  for (int ii=0;ii<(int)w.size()/2;ii++) {
    TComplex X(w[2*ii],w[2*ii+1]);
    X=X*C.Exp(TComplex(0.,-pShift));  // Phase Shift
    w[2*ii]=X.Re();
    w[2*ii+1]=X.Im();
  }
  w.FFTW(-1);

  // copy shifted data to input this->data array
  for(int i=0;i<(int)w.size();i++) {
    int j=ibeg-isize/4+i;
    if((j>=0)&&(j<(int)this->size())) this->data[j]=w[i];
  }

  return;
}

template <class DataType_t>
void gwavearray<DataType_t>::DumpToFile(char* fname)
{
  int n=this->size();

  FILE *fp;

  if ( (fp = fopen(fname, "w")) == NULL ) {
     cout << " gwavearray::DumpToFile() error: cannot open file " << fname <<". \n";
     return;
  };

  double dt= this->rate() ? 1./this->rate() : 1.;
  for (int i = 0; i < n; i++) fprintf( fp,"%e\t%e\n", i*dt, this->data[i]);
  fclose(fp);
}

/*
//______________________________________________________________________________
template <class DataType_t>
void gwavearray<DataType_t>::Streamer(TBuffer &R__b)
{
   // Stream an object of class gwavearray<DataType_t>.

   UInt_t R__s, R__c;
   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c); if (R__v) { }
      wavearray<DataType_t>::Streamer(R__b);
      R__b >> wextern;
      R__b.CheckByteCount(R__s, R__c, gwavearray<DataType_t>::IsA());
   } else {
      R__c = R__b.WriteVersion(gwavearray<DataType_t>::IsA(), kTRUE);
      wavearray<DataType_t>::Streamer(R__b);
      R__b << wextern;
      R__b.SetByteCount(R__c, kTRUE);
   }
}
*/

// instantiations

#define CLASS_INSTANTIATION(class_) template class gwavearray< class_ >;

CLASS_INSTANTIATION(short)
CLASS_INSTANTIATION(int)
CLASS_INSTANTIATION(unsigned int)
CLASS_INSTANTIATION(long)
CLASS_INSTANTIATION(long long)
CLASS_INSTANTIATION(float)
CLASS_INSTANTIATION(double)

