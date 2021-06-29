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


// draw the network antenna pattern (read config from the user_parameters.C file) : used with the cwb_draw_antpat 
{

  #define RESOLUTION  2 
  #define COORDINATES "Geographic"
  #define WORLD_MAP_DIR "$CWB_GWAT/data/"

  // polarization=0 -> |Fx|                     DPF
  // polarization=1 -> |F+|                     DPF
  // polarization=2 -> |Fx|/|F+|                DPF
  // polarization=3 -> sqrt((|F+|^2+|Fx|^2)/nIFO)      DPF
  // polarization=4 -> |Fx|^2                   DPF
  // polarization=5 -> |F+|^2                   DPF
  // polarization=6 -> Fx                       only with 1 detector
  // polarization=7 -> F+                       only with 1 detector
  // polarization=8 -> F1x/F2x                  only with 2 detectors
  // polarization=9 -> F1+/F2+                  only with 2 detectors
  // polarization=10 -> sqrt(|F1+|^2+|F1x|^2)/sqrt(|F2+|^2+|F2x|^2)     only with 2 detectors
  // polarization=11 -> The same as (10) but averaged over psi          only with 2 detectors


  bool save_plot=false;
  int polarization=3;
  if(gSystem->Getenv("CWB_ANTPAT_POLARIZATION")!=NULL) {
    TString cwb_antpat_polarization=TString(gSystem->Getenv("CWB_ANTPAT_POLARIZATION"));
    if(cwb_antpat_polarization.CompareTo("")!=0) {
      if(!cwb_antpat_polarization.IsFloat()) {cout<< "Error : CWB_ANTPAT_POLARIZATION is not a number" << endl;exit(1);}
      polarization=cwb_antpat_polarization.Atoi();
    }
  }
  if(gSystem->Getenv("CWB_ANTPAT_SAVE_PLOT")!=NULL) {
    TString cwb_save_plot=TString(gSystem->Getenv("CWB_ANTPAT_SAVE_PLOT"));
    if(cwb_save_plot.CompareTo("")!=0) {
      if(!cwb_save_plot.IsFloat()) {cout<< "Error : CWB_ANTPAT_SAVE_PLOT is not a number" << endl;exit(1);}
      if((cwb_save_plot.Atoi()==0)||(cwb_save_plot.Atoi()==1)) save_plot=cwb_save_plot.Atoi();
    }
  }
  TString cwb_antpat_projection="hammer";
  if(gSystem->Getenv("CWB_ANTPAT_PROJECTION")!=NULL) {
    cwb_antpat_projection=TString(gSystem->Getenv("CWB_ANTPAT_PROJECTION"));
    if(cwb_antpat_projection.Contains("rect")) cwb_antpat_projection="";
  }

  double rad2deg = 180./TMath::Pi();
  double deg2rad = TMath::Pi()/180.;

  gnetwork* gNET = new gnetwork;

  gskymap* gSM = gNET->GetGskymap();
  gSM->SetOptions(cwb_antpat_projection,COORDINATES,RESOLUTION/2);

  TString ifos[60];
  for(int n=0; n<nIFO; n++) ifos[n]=ifo[n];
  char ifostr[64]="";
  for(int n=0; n<nIFO; n++) sprintf(ifostr,"%s %s",ifostr,ifo[n]);
  cout << "Network : " << ifostr << endl;

  TString world_map = gSystem->ExpandPathName(WORLD_MAP_DIR);
  gSM->SetWorldMapPath(world_map.Data());
  gSM->SetWorldMap();

  TH2D* h2 = (TH2D*)gSM->GetHistogram();
  h2->GetXaxis()->SetTitleSize(0.05);
  h2->GetXaxis()->SetLabelSize(0.05);
  h2->GetYaxis()->SetTitleSize(0.05);
  h2->GetYaxis()->SetLabelSize(0.05);
  h2->GetYaxis()->SetLabelFont(42);
  h2->GetYaxis()->SetLabelFont(42);
  h2->GetXaxis()->SetTitleFont(42);
  h2->GetYaxis()->SetTitleFont(42);

  if(nIFO>1) h2->GetZaxis()->SetRangeUser(0,1.0); 
  else h2->GetZaxis()->SetRangeUser(0,1.0);

  if(polarization==2) h2->GetZaxis()->SetRangeUser(0,1.0);

  detector* pD[3];
  for(int i=0; i<nIFO; i++) pD[i] = new detector(ifo[i]); // built in detector
  for(int i=0; i<nIFO; i++) gNET->add(pD[i]);

  int palette = 0;
  bool btitle = true;
  gNET->DrawAntennaPattern(polarization,palette,btitle);
  gNET->DrawSites(kBlack,2.0);
  gNET->DrawSitesArms(1000000,kWhite,3.0);
  gNET->DrawSitesShortLabel(kBlack);

  if(save_plot) {
    char ofName[256];
    sprintf(ofName,"%s/AntPat_%s",dump_dir,data_label);

    if(polarization==0) sprintf(ofName,"%s%s.png",ofName,"_Fc");
    if(polarization==1) sprintf(ofName,"%s%s.png",ofName,"_Fp");
    if(polarization==2) sprintf(ofName,"%s%s.png",ofName,"_Fc_over_Fp");
    if(polarization==3) sprintf(ofName,"%s%s.png",ofName,"_Sqrt_Fp2_plus_Fc2");

    cout << "Write : " << ofName << endl;
    gSM->Print(ofName);
    exit(0);
  }
}

