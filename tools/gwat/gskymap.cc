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



#include "gskymap.hh"
#include "TPaletteAxis.h"


ClassImp(gskymap)

//______________________________________________________________________________
/* Begin_Html
<center><h2>gskymap class</h2></center>
This class gskymap is derived from the skymap class
 It is used to produce the skymap plots.

<p>

<h3><a name="usage">Usage</a></h3>

create gskymap with HEALPix order=7
<pre>
   <a class="funcname" href="#gskymap:gskymap">gskymap</a> gSM((int)7);
</pre>
set gskymap options
<pre>
  gSM.<a class="funcname" href="#gskymap:SetOptions">SetOptions</a>("cartesian","geographic",1);
</pre>
set title
<pre>
  gSM.<a class="funcname" href="#gskymap:SetTitle">SetTitle</a>("Projection : cartesian  -  Coordinates : geographic");
</pre>
set world map
<pre>
  gSM.<a class="funcname" href="#gskymap:SetWorldMap">SetWorldMap</a>();
</pre>
draw skymap
<pre>
  gSM.<a class="funcname" href="#gskymap:Draw">Draw</a>(0,"col");
</pre>

<p>
<h3><a name="example">Example</a></h3>
<p>
The macro <a href="./tutorials/gwat/DrawGskymap.C.html">DrawGskymap.C</a> is an example which shown how to use the gskymap class.<br>
The pictures below gives the macro output plots.<br>
<p>

End_Html

Begin_Macro
DrawGskymap.C("cartesian")
End_Macro 
Begin_Macro
DrawGskymap.C("hammer")
End_Macro 
Begin_Macro
DrawGskymap.C("parabolic")
End_Macro 
Begin_Macro
DrawGskymap.C("sinusoidal")
End_Macro */


using namespace std;

//______________________________________________________________________________
void
gskymap::SetOptions(TString projection, TString coordinate, double resolution, bool goff) { 
//
// Set graphical options
//
// Input: projection    - set projection type 
//                        - "hammer"            (default)
//                        - "sinusoidal"
//                        - "parabolic"
//                        - "cartesian" or ""
//
//        coordinate    - set the coordinate type
//                        - "geographic"        (default) 
//                        - "celestial" 
//                        - "cwb"
//
//        resolution    - factor used to increase the resolution of the skymap (default=1)
//                        phi = 2*360*resolution points, theta = 2*180*resolution	
//
//        goff          - true -> turn off the graphical output (default = false)
//  

  coordinate.ToUpper();

  if (resolution<=0) {cout << "gskymap::gskymap: Error sky resolution must be >0" << endl;exit(1);}

  if(gSystem->Getenv("CWB_GWAT")!=NULL) {
    worldMapPath=TString(gSystem->Getenv("CWB_GWAT"))+"/data";
  } else {
    worldMapPath=".";
  }

  // generate a unique base name for canvas and hist2d
  gRandom->SetSeed(0);
  char sname[64];sprintf(sname,"gskymap_%d",int(gRandom->Rndm(13)*1.e9));
  name = sname;
  this->projection = projection;
  this->coordinate = coordinate;
  this->resolution = resolution;
  this->goff = goff;

  wtopx=35; wtopy=46; ww=817; wh=472;
 
  double phMin = -180;
  double phMax =  180;
  double thMin = -90;
  double thMax =  90;
  if (coordinate.CompareTo("CWB")==0) {phMin=0;phMax=360;thMin=0;thMax=180;}
  if (coordinate.CompareTo("CELESTIAL")==0) {phMin=0;phMax=360;thMin=-90;thMax=90;}

  if(h2!=NULL) delete h2;
  char nameh[64];sprintf(nameh,"h_%s",name.Data());
  h2 = new TH2D(nameh,"gskymap", int(360*resolution), phMin, phMax, int(180*resolution), thMin, thMax);
  h2->SetStats(kFALSE);
  h2->SetTitleFont(12);
  h2->SetTitle(title);
  h2->SetFillColor(kWhite);

  h2->GetXaxis()->SetNdivisions(-604);
  h2->GetXaxis()->SetLabelFont(42);
  h2->GetXaxis()->SetLabelOffset(0.014);
  h2->GetXaxis()->SetTitleOffset(1.2);
  h2->GetYaxis()->SetTitleOffset(0.8);
  h2->GetYaxis()->SetNdivisions(-602);
  h2->GetYaxis()->SetLabelFont(42);
  h2->GetYaxis()->SetLabelOffset(0.01);
  h2->GetZaxis()->SetLabelFont(42);
  h2->GetZaxis()->SetNoExponent(false);

  h2->GetXaxis()->SetTitleFont(42);
  if ((coordinate.CompareTo("CWB")==0)||(coordinate.CompareTo("GEOGRAPHIC")==0)) {
    h2->GetXaxis()->SetTitle("#phi (deg)");
  } else if (coordinate.CompareTo("CELESTIAL")==0) {
    h2->GetXaxis()->SetTitle("RA (deg)");
  } else {
    {cout << "gskymap::gskymap: Error : coordinate system not declared" << endl;exit(1);} 
  }
  h2->GetXaxis()->CenterTitle(true);

  h2->GetYaxis()->SetTitleFont(42);
  if ((coordinate.CompareTo("CWB")==0)||(coordinate.CompareTo("GEOGRAPHIC")==0)) {
    h2->GetYaxis()->SetTitle("#theta (deg)");
  } else if (coordinate.CompareTo("CELESTIAL")==0) {
    h2->GetYaxis()->SetTitle("DEC (deg)");
  } else {
    {cout << "gskymap::gskymap: Error : coordinate system not declared" << endl;exit(1);} 
  }
  h2->GetYaxis()->CenterTitle(true);

  h2->GetZaxis()->SetTitleOffset(0.6);
  h2->GetZaxis()->SetTitleFont(42);
  h2->GetZaxis()->SetTitle(zAxisTitle);
  h2->GetZaxis()->CenterTitle(true);

  if(this->size()) FillData();
}

//______________________________________________________________________________
void
gskymap::CreateCanvas() { 
//
// create canvas object used for skymap plot
//

//  TCanvas* tcanvas = (TCanvas*) gROOT->FindObject(name);
//  if (tcanvas!=NULL) {cout << "gskymap::gskymap: Error Canvas Name already exist" << endl;exit(1);} 
  if(canvas!=NULL) delete canvas;
  if(!goff) {
    char namec[64];sprintf(namec,"c_%s",name.Data());
    canvas = new TCanvas(namec, "gskymap", wtopx, wtopy, ww, wh);
    canvas->Clear();
    canvas->ToggleEventStatus();
    canvas->SetGridx(false);
    canvas->SetGridy(false);
    canvas->SetLogz(isLogz);
    canvas->SetFillColor(kWhite);
    canvas->SetLeftMargin(0.08);
    canvas->SetBottomMargin(0.13);
    canvas->SetBorderMode(0);
    //canvas->SetWindowSize(1200,600);
    //canvas->SetGrayscale();
  } else {
    canvas=NULL;
  }

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
  gStyle->SetMarkerStyle(7);
  gStyle->SetMarkerSize(2);
  gStyle->SetCanvasColor(kWhite);
  gStyle->SetStatBorderSize(1);
}

//______________________________________________________________________________
gskymap::~gskymap() {
//
// destructor
//

  if(h2!=NULL) delete h2;
  if(canvas!=NULL) delete canvas;

  if (wm_size) delete [] wm_lon; 
  if (wm_size) delete [] wm_lat; 

  ClearGalacticDisk();
  ClearWorldMap();
  ClearGridx();
  ClearGridy();
  ClearAxisLabel();
}

//______________________________________________________________________________
void
gskymap::FillData(char* fname) {
//
// fill the internal skymap reading from ascii file
//
// Input: fname  - input file name
//
// format :  
//
// header :         sky_res  -> area of a pixel (degrees) 
//                  theta_1  -> theta min
//                  theta_2  -> theta max
//                  phi_1    -> phi min
//                  phi_2    -> phi max
//                  lenght   -> number of entries 
//
// list of data:    sky_index value
// 

  char iline[1024];
  TObjArray* tok;

  ifstream in;
  in.open(fname,ios::in);
  if (!in.good()) {cout << "gskymap::FillData: Error Opening File : " << fname << endl;exit(1);}

  in.getline(iline,1024);
  in.getline(iline,1024);
  if(TString(iline).Contains("skymap")==0) 
    {cout << "gskymap::FillData: " << fname << " is not a skymap file" << endl;exit(1);}
  in.getline(iline,1024);
  tok = TString(iline).Tokenize(TString(' '));
  TObjString* tsms  = (TObjString*)tok->At(1);
  if(tsms->GetString().IsAlpha()==1) 
    {cout << "gskymap::FillData: resolution is not well defined " << endl;exit(1);}
  double sms = tsms->GetString().Atof();
  in.getline(iline,1024);
  tok = TString(iline).Tokenize(TString(' '));
  TObjString* ttheta_1  = (TObjString*)tok->At(1);
  if(ttheta_1->GetString().IsAlpha()==1) 
    {cout << "gskymap::FillData: theta_1 is not well defined " << endl;exit(1);}
  double theta_1 = ttheta_1->GetString().Atof();
  in.getline(iline,1024);
  tok = TString(iline).Tokenize(TString(' '));
  TObjString* ttheta_2  = (TObjString*)tok->At(1);
  if(ttheta_2->GetString().IsAlpha()==1) 
    {cout << "gskymap::FillData: theta_2 is not well defined " << endl;exit(1);}
  double theta_2 = ttheta_2->GetString().Atof();
  in.getline(iline,1024);
  tok = TString(iline).Tokenize(TString(' '));
  TObjString* tphi_1  = (TObjString*)tok->At(1);
  if(tphi_1->GetString().IsAlpha()==1) 
    {cout << "gskymap::FillData: phi_1 is not well defined " << endl;exit(1);}
  double phi_1 = tphi_1->GetString().Atof();
  in.getline(iline,1024);
  tok = TString(iline).Tokenize(TString(' '));
  TObjString* tphi_2  = (TObjString*)tok->At(1);
  if(tphi_2->GetString().IsAlpha()==1) 
    {cout << "gskymap::FillData: phi_2 is not well defined " << endl;exit(1);}
  double phi_2 = tphi_2->GetString().Atof();
  in.getline(iline,1024);
  tok = TString(iline).Tokenize(TString(' '));
  TObjString* tlenght  = (TObjString*)tok->At(1);
  if(tlenght->GetString().IsAlpha()==1) 
    {cout << "gskymap::FillData: lenght is not well defined " << endl;exit(1);}
  int lenght = tlenght->GetString().Atoi();

  //sms = double(int(1000*sms))/1000.;

  //cout << sms << " " << theta_1 << " " << theta_2 << " " << phi_1 << " " << phi_2 << endl;

  int L = skymap::size();
  if(L!=lenght) 
    {cout << "gskymap::FillData: lenght is not consistent" << endl;exit(1);}

  while (1) {

    in.getline(iline,1024);
    if (!in.good()) break;
    //cout << ++cnt << " " << iline << endl;
    tok = TString(iline).Tokenize(TString(' '));

    if (tok->GetEntries()==1) continue;

    TObjString* tindex = (TObjString*)tok->At(0);
    TObjString* tvalue = (TObjString*)tok->At(1);

    if (tindex->GetString().IsAlpha()==1) 
      {cout << "gskymap::FillData: index " << tindex->GetString().Data() << " is not well defined " << endl;exit(1);}
    //if (tvalue->GetString().IsAlpha()==1) continue;
    //  {cout << "gskymap::FillData: value " << tvalue->GetString().Data() << " is not well defined " << endl;exit(1);}

    int index = tindex->GetString().Atoi();
    double value = tvalue->GetString().Atof();

    if(index>=0&&index<L)
      skymap::set(index,value);
    else 
      {cout << "gskymap::FillData: index " << index << " not allowed - max " << L-1 << endl;}

  }
  in.close();

  FillData();

  return;
}

//______________________________________________________________________________
void
gskymap::FillData() {
//
// Fill 2D histo with the internal skymap data
//

  int size=int(2*180*resolution*2*360*resolution);

  double* phi = new double[size]; 
  double* theta = new double[size]; 
  double* binc = new double[size]; 

  int k=0;
  for(int i=0;i<2*180*resolution;i++) {
    for(int j=0;j<2*360*resolution;j++) {
      phi[k]=(double)j/(2*resolution);
      theta[k]=(double)i/(2*resolution);
      binc[k]=skymap::get(theta[k],phi[k]);
      k++;
    }
  }

  if (coordinate.CompareTo("GEOGRAPHIC")==0) {
    for (int i=0;i<k;i++) CwbToGeographic(phi[i],theta[i],phi[i],theta[i]);
  } 
  if (coordinate.CompareTo("CELESTIAL")==0) {
    for (int i=0;i<k;i++) CwbToCelestial(phi[i],theta[i],phi[i],theta[i]);
  } 

  FillData(k, phi, theta, binc);

  delete [] phi;
  delete [] theta;
  delete [] binc;
}

//______________________________________________________________________________
void
gskymap::FillData(int size, double* phi, double* theta, double* binc) {
//
// Fill 2D histo using input arrays
//
// Input: size  - size of the input arrays
//        phi   - arrays of phi values
//        theta - arrays of theta values
//        binc  - arrays of amplitudes 
//

  for (int i=0;i<size;i++) {
    double x,y;
    if (coordinate.CompareTo("CWB")==0) {
      if(phi[i]<0) phi[i]+=360.;
      if(phi[i]>360) phi[i]-=360.;
      if(theta[i]<0) theta[i]+=180.;
      if(theta[i]>180) theta[i]-=180.;
    }
    if (coordinate.CompareTo("GEOGRAPHIC")==0) {
      if(phi[i]<-180) phi[i]+=360.;
      if(phi[i]>180) phi[i]-=360.;
      if(theta[i]<-90) theta[i]+=180.;
      if(theta[i]>90) theta[i]-=180.;
      phi[i]+=180.;
      theta[i]+=90.;
    }
    if (coordinate.CompareTo("CELESTIAL")==0) {
      if(phi[i]<0) phi[i]+=360.;
      if(phi[i]>360) phi[i]-=360.;
      if(theta[i]<-90) theta[i]+=180.;
      if(theta[i]>90) theta[i]-=180.;
      theta[i]+=90.;
    }
    if (projection.CompareTo("hammer")==0) {
      ProjectHammer((phi[i]-180), (theta[i]-90), x, y);
      x+=180+1;
      y+=90+1;
    } else
    if (projection.CompareTo("sinusoidal")==0) {
      ProjectSinusoidal((phi[i]-180), (theta[i]-90), x, y);
      x+=180+1;
      y+=90+1;
    } else
    if (projection.CompareTo("parabolic")==0) {
      ProjectParabolic((phi[i]-180), (theta[i]-90), x, y);
      x+=180+1;
      y+=90+1;
    } else {
      x=phi[i];
      y=theta[i];
    }
    h2->SetBinContent(int(x*resolution)+1,int(y*resolution)+1,binc[i]);

    if (coordinate.CompareTo("GEOGRAPHIC")==0) {
      CwbToGeographic(phi[i],theta[i],phi[i],theta[i]);
      //phi[i]-=180.;
      //theta[i]+=90.;
    }
    if (coordinate.CompareTo("CELESTIAL")==0) {
      CwbToCelestial(phi[i],theta[i],phi[i],theta[i]);
    }
  }
}

//______________________________________________________________________________
void
gskymap::Draw(int dpaletteId, Option_t* goption) {
//
// draw skymap 
//
// Input: dpaletteId - color palette used for the skymap
//        goption    - graphical options 
//

  if(goff) return;

  CreateCanvas();

  TGaxis::SetMaxDigits(3);

  if(changed) {FillData();changed=false;}

  if (dpaletteId==DUMMY_PALETTE_ID) {
    if (paletteId!=0) {
      SetPlotStyle(paletteId);
    } else {
      gStyle->SetPalette(kBird,0);
    }
  } else {
    if (dpaletteId!=0) {
      SetPlotStyle(dpaletteId);
    } else {
      gStyle->SetPalette(kBird,0);
    }
  }

  canvas->cd();
  canvas->SetLogz(isLogz);
  h2->Draw(goption);
  // change palette's width
  canvas->Update();
  TPaletteAxis *palette = (TPaletteAxis*)h2->GetListOfFunctions()->FindObject("palette");
  if(palette) {
    palette->SetX1NDC(0.91);
    palette->SetX2NDC(0.933);
    palette->SetTitleOffset(0.92);
    palette->GetAxis()->SetTickSize(0.01);
    canvas->Modified();
  }

  Double_t rlonait[360*4],rlatait[360*4];
  Double_t kfix;
  Int_t kk,k;
  double mdistance=8.;	// min distance (degrees) between polyline segments (avoid solid line when saved) 

// draw lines of constant longitude
  if (isGridy) {
    ClearGridy();
    for (kk=0;kk<=360;kk+=15){
      kfix=kk;
      int klen=0;
      for (k=0;k<180*4;k++){
        double x,y;
        double phi = kfix;
        double theta = k/4.;
        if (projection.CompareTo("hammer")==0) {
          ProjectHammer((phi-180), (theta-90), x, y);
        } else
        if (projection.CompareTo("sinusoidal")==0) {
          ProjectSinusoidal((phi-180), (theta-90), x, y);
        } else
        if (projection.CompareTo("parabolic")==0) {
          ProjectParabolic((phi-180), (theta-90), x, y);
        } else {
          x=phi-180;
          y=theta-90;
        }
        if (coordinate.CompareTo("CWB")==0) {x+=180;y+=90;}
        //if (coordinate.CompareTo("CELESTIAL")==0) {x+=180;}
        if (coordinate.CompareTo("CELESTIAL")==0) {x+=180;x=360-x;}  //TO BE CHECKED
        rlonait[klen] = x;
        rlatait[klen] = y;
        double distance = klen>0 ? sqrt(pow(x-rlonait[klen-1],2)+pow(y-rlatait[klen-1],2)) : mdistance; 
        if(distance>=mdistance || k==(180*4-1)) klen++;
      }
      TPolyLine *pL = new TPolyLine(klen,rlonait,rlatait);
      pL->SetLineColor(colorGridx);
      if(kk==0 || kk==360) pL->SetLineStyle(1); else pL->SetLineStyle(3);
      pL->Draw();
      gridyL.push_back(pL);
    }
  }

/// draw lines of constant latitude
  if (isGridx) {
    ClearGridx();
    //for (kk=10;kk<180;kk+=10){
    for (kk=15;kk<180;kk+=15){
      kfix=kk;
      int klen=0;
      for (k=0;k<360*4;k++){
        double x,y;
        double phi = k/4.;
        double theta = kfix;
        if (projection.CompareTo("hammer")==0) {
          ProjectHammer((phi-180), (theta-90), x, y);
        } else
        if (projection.CompareTo("sinusoidal")==0) {
          ProjectSinusoidal((phi-180), (theta-90), x, y);
        } else
        if (projection.CompareTo("parabolic")==0) {
          ProjectParabolic((phi-180), (theta-90), x, y);
        } else {
          x=phi-180;
          y=theta-90;
        }
        if (coordinate.CompareTo("CWB")==0) {x+=180;y+=90;}
        //if (coordinate.CompareTo("CELESTIAL")==0) {x+=180;}
        if (coordinate.CompareTo("CELESTIAL")==0) {x+=180;x=360-x;}  // TO BE CHECKED
        rlonait[klen] = x;
        rlatait[klen] = y;
        double distance = klen>0 ? sqrt(pow(x-rlonait[klen-1],2)+pow(y-rlatait[klen-1],2)) : mdistance; 
        if(distance>=mdistance || k==(360*4-1)) klen++;
      }
      TPolyLine *pL = new TPolyLine(klen,rlonait,rlatait);
      pL->SetLineColor(colorGridy);
      if(kk==0 || kk==180) pL->SetLineStyle(1); else pL->SetLineStyle(3);
      pL->Draw();
      gridxL.push_back(pL);
    }
  }

/// draw labels of constant latitude
/*
  if ((projection.CompareTo("")!=0)&&(coordinate.CompareTo("GEOGRAPHIC")==0)) {
    ClearAxisLabel();
    double phi[100],theta[100];
    int k=0;
    //for (kk=0;kk<=360;kk+=40) {phi[k]=kk-180;theta[k]=0;k++;}
    //for (kk=0;kk<=180;kk+=10) {phi[k]=-180;theta[k]=kk-90;k++;}
    for (kk=0;kk<=360;kk+=15) {phi[k]=kk-180;theta[k]=0;k++;}
    for (kk=0;kk<=180;kk+=15) {phi[k]=-180;theta[k]=kk-90;k++;}
    char label[16];
    for(kk=0;kk<k;kk++) {
      if(theta[kk]!=0) 
        sprintf(label,"%d",(int)theta[kk]); 
      else
        sprintf(label,"%d",(int)phi[kk]+180); 
      double x,y;
      if (projection.CompareTo("hammer")==0) {
        ProjectHammer(phi[kk], theta[kk], x, y);
      } else
      if (projection.CompareTo("sinusoidal")==0) {
        ProjectSinusoidal(phi[kk], theta[kk], x, y);
      } else
      if (projection.CompareTo("parabolic")==0) {
        ProjectParabolic(phi[kk], theta[kk], x, y);
      } else {
        x=phi[kk];y=theta[kk];
      }
      TText *tP = new TText(x,y,label);
      tP->SetTextSize(0.04);
      tP->SetTextColor(kBlack);
      tP->SetTextFont(42);
      if(theta[kk]>0)  tP->SetTextAlign(31);
      if(theta[kk]==0) tP->SetTextAlign(32);
      if(theta[kk]<0)  tP->SetTextAlign(33);
      tP->Draw();
      axisT.push_back(tP);
    }
  }
*/

/// draw line of galactic disk
  if (gpsGalacticDisk>=0.) {
    ClearGalacticDisk();

    double gphi = gpsGalacticDisk>0. ? skymap::RA2phi(0., gpsGalacticDisk) : 0.;

    int L=4*360;
    wavearray<float> th(L);
    wavearray<float> ph(L);
    for (int l=0;l<L;l++) {
      double phi,theta;
      GalacticToEquatorial(l/4.,0.,phi,theta);
      ph.data[l]=phi+gphi;
      th.data[l]=theta;
      ph.data[l]=fmod(ph.data[l],360);
    }

    wavearray<double> x(L);
    wavearray<double> y(L);
    int P=0;
    for (int l=0;l<L;l++) {
      if((fabs(ph.data[l]-ph.data[l==0?L-1:l-1])<180)&&(l<L-1)) {
        x.data[P]=ph.data[l]; y.data[P]=th.data[l];

        if (projection.CompareTo("hammer")==0) {
          ProjectHammer((x.data[P]-180), (-y.data[P]), x.data[P], y.data[P]);
          x.data[P]=x.data[P]+180; y.data[P]=-y.data[P];
        } else if (projection.CompareTo("sinusoidal")==0) {
          ProjectSinusoidal((x.data[P]-180), (-y.data[P]), x.data[P], y.data[P]);
          x.data[P]=x.data[P]+180; y.data[P]=-y.data[P];
        } else if (projection.CompareTo("parabolic")==0) {
          ProjectParabolic((x.data[P]-180), (-y.data[P]), x.data[P], y.data[P]);
          x.data[P]=x.data[P]+180; y.data[P]=-y.data[P];
        } else {
        }
        if (coordinate.CompareTo("GEOGRAPHIC")==0)
          {x.data[P]=x.data[P]-180;}
        if (coordinate.CompareTo("CWB")==0)
          {y.data[P]=90+y.data[P];}
        if (coordinate.CompareTo("CELESTIAL")==0)
          {x.data[P]=360-x.data[P];y.data[P]=90+y.data[P];}
        P++;
      } else {
        TPolyLine *pL = new TPolyLine(P,x.data,y.data);
        pL->SetLineColor(colorGalacticDisk);
        pL->Draw();
        gdL.push_back(pL);
        P=0;
      }
    }
  }

  if (drawWorldMap && (coordinate.CompareTo("CELESTIAL")!=0)) {
    ClearWorldMap();
    if (wm_size==0) {
      // geographic coordinates
      wm_size = ReadWorlMapCoastLine(wm_lon, wm_lat);
    }
    for (int n=0;n<wm_size;n++) {
      double x,y;
      double phi = wm_lon[n];
      double theta = wm_lat[n];
      if (coordinate.CompareTo("CWB")==0) {phi+=180;if(phi>180) phi=phi-360;}
      if (coordinate.CompareTo("CELESTIAL")==0) {phi+=180;if(phi>180) phi=phi-360;phi=360-phi;}  //TOBE CHECK
      if (projection.CompareTo("hammer")==0) {
        ProjectHammer(phi, theta, x, y);
      } else
      if (projection.CompareTo("sinusoidal")==0) {
        ProjectSinusoidal(phi, theta, x, y);
      } else
      if (projection.CompareTo("parabolic")==0) {
        ProjectParabolic(phi, theta, x, y);
      } else {
        x=phi;
        y=theta;
      }
      if (coordinate.CompareTo("CWB")==0) {GeographicToCwb(x+180,y,x,y);}
      if (coordinate.CompareTo("CELESTIAL")==0) {CelestialToCwb(x+180,y,x,y);}  // TO BE CHECKED
      TMarker *mP = new TMarker(x,y, 1); // WARNING : style=20,size=0.1 has issue when printed
      mP->SetMarkerColor(kGray+1);
      mP->Draw();
      wmM.push_back(mP);
    }
  }

  if (coordinate.CompareTo("CELESTIAL")==0) {
    ReverseXAxis(h2);
  }

  TGaxis::SetMaxDigits();
}

//______________________________________________________________________________
void 
gskymap::DrawMarker(double ra, double dec, double gps, int marker, Size_t msize, Color_t mcolor) {
//
// Draw a marker in the ra,dec coordinates (see ROOT TMarker input parameters)
//
// Input: ra,dec     - celestial coordinates
//        gps        - gps time
//        mcolor     - color of marker
//        msize      - size of marker
//        mstyle     - style of marker
//

  double phi,theta;
  GeographicToCwb(ra,dec,phi,theta);
  phi = gps>0 ? skymap::RA2phi(phi, gps) : phi;
  if(coordinate.CompareTo("GEOGRAPHIC")==0) CwbToGeographic(phi,theta,phi,theta);
  if(coordinate.CompareTo("CELESTIAL")==0) CwbToCelestial(phi,theta,phi,theta);
  DrawMarker(phi, theta, marker, msize, mcolor);
}

//______________________________________________________________________________
void 
gskymap::DrawMarker(double phi, double theta, int marker, Size_t msize, Color_t mcolor) {
//
// Draw a marker in the phi,theta coordinates (see ROOT TMarker input parameters)
//
// Input: phi,theta  - skymap coordinates
//        mcolor  - color of marker
//        msize   - size of marker
//        mstyle  - style of marker
//

  if(goff) return;

  double x,y;
  if (projection.CompareTo("hammer")==0) {
    //ProjectHammer((phi-180), -(theta-90), x, y);
    ProjectHammer(phi, theta, x, y);
  } else
  if (projection.CompareTo("sinusoidal")==0) {
    //ProjectSinusoidal((phi-180), -(theta-90), x, y);
    ProjectSinusoidal(phi, theta, x, y);
  } else
  if (projection.CompareTo("parabolic")==0) {
    //ProjectParabolic((phi-180), -(theta-90), x, y);
    ProjectParabolic(phi, theta, x, y);
  } else {
    x=phi;
    y=theta;
  }

  TMarker *mP = new TMarker(x,y, marker);
  mP->SetMarkerSize(msize);
  mP->SetMarkerColor(mcolor);
  mP->Draw();
}

//______________________________________________________________________________
void 
gskymap::DrawText(double ra, double dec, double gps, TString text, double tsize, Color_t tcolor) {
//
// Draw a text to ra,dec coordinates
//
// Input: ra,dec     - celestial coordinates
//        text       - text   
//        tsize      - size of text
//        tcolor     - color of text
//

  double phi,theta;
  GeographicToCwb(ra,dec,phi,theta);
  phi = gps>0 ? skymap::RA2phi(phi, gps) : phi;
  if(coordinate.CompareTo("GEOGRAPHIC")==0) CwbToGeographic(phi,theta,phi,theta);
  if(coordinate.CompareTo("CELESTIAL")==0) CwbToCelestial(phi,theta,phi,theta);
  DrawText(phi, theta, text, tsize, tcolor);
}

//______________________________________________________________________________
void 
gskymap::DrawText(double phi, double theta, TString text, double tsize, Color_t tcolor) {
//
// Draw a text to phi,theta coordinates
//
// Input: phi,theta  - skymap coordinates
//        text       - text   
//        tsize      - size of text
//        tcolor     - color of text
//

  if(goff) return;

  double x,y;
  if (projection.CompareTo("hammer")==0) {
    //ProjectHammer((phi-180), -(theta-90), x, y);
    ProjectHammer(phi, theta, x, y);
  } else
  if (projection.CompareTo("sinusoidal")==0) {
    //ProjectSinusoidal((phi-180), -(theta-90), x, y);
    ProjectSinusoidal(phi, theta, x, y);
  } else
  if (projection.CompareTo("parabolic")==0) {
    //ProjectParabolic((phi-180), -(theta-90), x, y);
    ProjectParabolic(phi, theta, x, y);
  } else {
    x=phi;
    y=theta;
  }

  TText *tP = new TText(x,y,text);
  tP->SetTextSize(tsize);
  tP->SetTextColor(tcolor);
  tP->Draw();
}

//______________________________________________________________________________
void 
gskymap::ProjectHammer(Double_t l, Double_t b, Double_t &Al, Double_t &Ab)
{
//
// Static function
// Convert Right Ascension, Declination to X,Y using an AITOFF projection.
// This procedure can be used to create an all-sky map in Galactic
// coordinates with an equal-area Aitoff projection.  Output map
// coordinates are zero longitude centered.
// Also called Hammer-Aitoff projection (first presented by Ernst von Hammer in 1892)
// source: GMT
// code from  Ernst-Jan Buis
// http://en.wikipedia.org/wiki/Hammer_projection
//
// INPUT : Geographic Coordinates
// LONGITUDE (W:-Pi,E:+Pi)
// LATITUDE  (N:+Pi/2,S:-Pi/2)
//

   Double_t x, y;

   Double_t alpha2 = (l/2)*TMath::DegToRad();
   Double_t delta  = b*TMath::DegToRad();
   Double_t r2     = TMath::Sqrt(2.);
   Double_t f      = 2*r2/TMath::Pi();
   Double_t cdec   = TMath::Cos(delta);
   Double_t denom  = TMath::Sqrt(1. + cdec*TMath::Cos(alpha2));
   x      = cdec*TMath::Sin(alpha2)*2.*r2/denom;
   y      = TMath::Sin(delta)*r2/denom;
   x     *= TMath::RadToDeg()/f;
   y     *= TMath::RadToDeg()/f;
   //  x *= -1.; // for a skymap swap left<->right
   Al = x;
   Ab = y;

   return;
}

//______________________________________________________________________________
void 
gskymap::ProjectSinusoidal(Double_t l, Double_t b, Double_t &Al, Double_t &Ab)
{
//
// Static function
// code from  Ernst-Jan Buis
//
// INPUT : Geographic Coordinates
// LONGITUDE (W:-Pi,E:+Pi)
// LATITUDE  (N:+Pi/2,S:-Pi/2)
//

   Al = l*cos(b*TMath::DegToRad());
   Ab = b;
   return;
}

//______________________________________________________________________________
void 
gskymap::ProjectParabolic(Double_t l, Double_t b, Double_t &Al, Double_t &Ab)
{
//
// Static function
// code from  Ernst-Jan Buis
//
// INPUT : Geographic Coordinates
// LONGITUDE (W:-Pi,E:+Pi)
// LATITUDE  (N:+Pi/2,S:-Pi/2)
//

   Al = l*(2.*TMath::Cos(2*b*TMath::DegToRad()/3) - 1);
   Ab = 180*TMath::Sin(b*TMath::DegToRad()/3);
   return;
}

//______________________________________________________________________________
void
gskymap::SetPlotStyle(int paletteId) {
//
// set the color palette used for the skymap plot
//
// Input: paletteId  - palette number (1,2,3,4,5)
//
// -----------------------------------------------------------------
// http://ultrahigh.org/2007/08/20/making-pretty-root-color-palettes/
// -----------------------------------------------------------------
//

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

//______________________________________________________________________________
int  
gskymap::ReadWorlMapCoastLine(double*& wm_lon, double*& wm_lat) {
//
// Read the World Coastline map 
//
// return in longitudes/latitudes arrays the world cordinates
//
// -----------------------------------------------------------------
// Coastline map (shapefile) download from :
// http://www.naturalearthdata.com/downloads/10m-physical-vectors/10m-coastline/
// hape2Text transforms ESRI shapefiles into Ascii  files
// http://www.zonums.com/shp2text.html
// -----------------------------------------------------------------
//

  char ifileName[256];
  sprintf(ifileName,"%s/%s",worldMapPath.Data(),WorlMapCoastLine);

  ifstream in;
  in.open(ifileName,ios::in);
  if (!in.good()) {cout << "gskymap::ReadWorlMapCoastLine: Error Opening File : " << ifileName << endl;exit(1);}

  int cnt=0;
  char iline[1024];
  int size = WM_ENTRIES;
  wm_lon = new double[size];
  wm_lat = new double[size];
  while (1) {

    in.getline(iline,1024);
    if (!in.good()) break;
    //cout << ++cnt << " " << iline << endl;
    TObjArray* tok = TString(iline).Tokenize(TString(' '));

    if (tok->GetEntries()==1) continue;

    TObjString* tra   = (TObjString*)tok->At(0);
    TObjString* tdec  = (TObjString*)tok->At(1);

    if (tdec->GetString().IsAlpha()==1) continue;
    //cout << tra->GetString().Data() << " " << tdec->GetString().Data() << endl;

    double ra   = tra->GetString().Atof();
    double dec  = tdec->GetString().Atof();

    if (fabs(dec)<0.001) continue;

    wm_lon[cnt] = ra;
    wm_lat[cnt] = dec;
    cnt++;

  }
  cout << cnt << endl;

  in.close();

  drawWorldMap = false;

  return size;
}

//______________________________________________________________________________
void 
gskymap::ClearAxisLabel() {
//
// clear the axis labels
//

  for(int i=0; i<(int)axisT.size(); i++) {if(axisT[i]) delete axisT[i];}
  axisT.clear();
  return;
}

//______________________________________________________________________________
void 
gskymap::ClearGalacticDisk() {
//
// clear galactic disk
//

  for(int i=0; i<(int)gdL.size(); i++) {if(gdL[i]) delete gdL[i];}
  gdL.clear();
  return;
}

//______________________________________________________________________________
void 
gskymap::ClearWorldMap() {
//
// clear world map
//

  for(int i=0; i<(int)wmM.size(); i++) {if(wmM[i]) delete wmM[i];}
  wmM.clear();
  return;
}

//______________________________________________________________________________
void 
gskymap::ClearGridx() {
//
// clear grids along the x axis 
//

  for(int i=0; i<(int)gridxL.size(); i++) {if(gridxL[i]) delete gridxL[i];}
  gridxL.clear();
}

//______________________________________________________________________________
void 
gskymap::ClearGridy() {
//
// clear grids along the y axis 
//

  for(int i=0; i<(int)gridyL.size(); i++) {if(gridyL[i]) delete gridyL[i];}
  gridyL.clear();
  return;
}

//______________________________________________________________________________
void 
gskymap::Print(TString pname) {
//
// save plot to file
//
// Input: pname  - output file name
//

  if(canvas==NULL) {cout << "gskymap::Print: Error - canvas not initialized " << endl;exit(1);}

  TGaxis::SetMaxDigits(3);

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
    pname.ReplaceAll(".PNG",".png");
    canvas->Print(pname);
  }
*/

  pname.ReplaceAll(".PNG",".png");
  canvas->Print(pname);

  TGaxis::SetMaxDigits();
}

//______________________________________________________________________________
void
gskymap::DumpSkyMap(char* fname) {
//
// dump skymap to ascii file
//
// Input: fname  - output file name
//
// format :  
//
// header :         sky_res  -> area of a pixel (degrees) 
//                  theta_1  -> theta min
//                  theta_2  -> theta max
//                  phi_1    -> phi min
//                  phi_2    -> phi max
//                  lenght   -> number of entries 
//
// list of data:    sky_index value
// 

  int L = skymap::size();

  double x;
  FILE* fp;             // dump file

  if((fp = fopen(fname, "w")) == NULL) {
    cout << "netevent::DumpSkyMap() error: cannot open file " << fname <<". \n";
    return;
  };
  fprintf(fp,"wat %s\n",watversion('f'));
  fprintf(fp,"class skymap\n");
  x = cos(this->theta_1*PI/180.)-cos(this->theta_2*PI/180.);
  x*= (this->phi_2-this->phi_1)*180/PI/L;
  x = double(int(1000*x))/1000.;
  fprintf(fp,"sky_res %f\n",sqrt(fabs(x)));
  fprintf(fp,"theta_1 %f\n",this->theta_1);
  fprintf(fp,"theta_2 %f\n",this->theta_2);
  fprintf(fp,"phi_1 %f\n",this->phi_1);
  fprintf(fp,"phi_2 %f\n",this->phi_2);
  fprintf(fp,"lenght %d\n",L);

  for(int l=0;l<L;l++) fprintf(fp,"%d %e\n",l,this->get(l));

  fclose(fp);
  return;
}

//______________________________________________________________________________
void
gskymap::ReverseXAxis(TH2D *h2) {
//
// reverse the x axis
//

  char xlabel[16];
  for (int i=1;i<=h2->GetNbinsX();i++) {
    //cout << i << " " << h2->GetXaxis()->GetBinCenter(i) << endl;
    double xvalue = h2->GetXaxis()->GetBinCenter(i)+h2->GetXaxis()->GetBinWidth(i)/2.;
    int nbin=int(h2->GetNbinsX()/4.); 
    if(i%nbin==0) {
      sprintf(xlabel,"%d",int(360-xvalue));
      //cout << i << " " << xvalue << endl;
      h2->GetXaxis()->SetBinLabel(i,xlabel);
      h2->GetXaxis()->SetLabelSize(0.06);
      h2->GetXaxis()->LabelsOption("h");
      h2->GetXaxis()->SetNdivisions(604);
      h2->GetXaxis()->SetLabelFont(42);
      h2->GetXaxis()->SetLabelOffset(0.010);
      h2->GetXaxis()->SetTitleOffset(1.5);
    }
  }
}

//______________________________________________________________________________
void gskymap::DumpObject(const char* file, const char* name)
{
//
// dump skymap into root file
//

  TFile* rfile = new TFile(const_cast<char*>(file), "RECREATE");
  this->Write(const_cast<char*>(name)); // write this object
  rfile->Close();
}

//______________________________________________________________________________
void gskymap::LoadObject(const char* file, const char* name)
{
//
// load gskymap object from root file
//

  TFile* rfile = new TFile(const_cast<char*>(file));
  gskymap* gsm = (gskymap*)rfile->Get(const_cast<char*>(name));
  if(gsm==NULL) {
    cout << "gskymap::LoadObject : Error - input file don't contains object gskymap" << endl;
    exit(1);
  }
  *this=*gsm;

  gridxL.clear();
  gridyL.clear();
  gdL.clear();
  wmM.clear();
  axisT.clear();

  FillData();
  rfile->Close();
}

//______________________________________________________________________________
void gskymap::Plot()
{
//
// Used to draw skymap from TBrowser
//

  gridxL.clear();
  gridyL.clear();
  gdL.clear();
  wmM.clear();
  axisT.clear();

  FillData();
  Draw();
}

