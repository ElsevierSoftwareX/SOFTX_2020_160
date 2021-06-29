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


#define MAX_GW	100


void DrawSensitivityPE(TString type, TString ifname_ff, TString ifname_of, TString ifname_re, TString odir, TString label) {

  bool check = (type!="phase")&&(type!="time")&&(type!="amp");
       check = check&&(type!="dphase50")&&(type!="dtime50")&&(type!="damp50");
       check = check&&(type!="dphase90")&&(type!="dtime90")&&(type!="damp90");
  if(check) {cout << "DrawSensitivityPE.C - Input type option error" << endl;exit(1);}

  // init blind colors 
  Color_t color[4];
  color[0] = CWB::Toolbox::getTableau10BlindColor("DarkOrange1");
  color[1] = CWB::Toolbox::getTableau10BlindColor("DeepSkyBlue4");
  color[2] = CWB::Toolbox::getTableau10BlindColor("DarkGray");
  color[3] = CWB::Toolbox::getTableau10BlindColor("SandyBrown");

  int N=0;
  double event[MAX_GW];
  double zero[MAX_GW]; for(int i=0;i<MAX_GW;i++) zero[i]=0;
  char dummy[256];
  char gw_name[MAX_GW][256];

  ifstream in_ff;
  in_ff.open(ifname_ff.Data(),ios::in);
  if (!in_ff.good()) {cout << "DrawSensitivityPE.C - Error Opening File : " << ifname_ff.Data() << endl;exit(1);}
  double dff_phase[MAX_GW],dff_time[MAX_GW],dff_amp[MAX_GW];
  double ff_dphase50[MAX_GW],ff_dtime50[MAX_GW],ff_damp50[MAX_GW];
  double ff_dphase90[MAX_GW],ff_dtime90[MAX_GW],ff_damp90[MAX_GW];
  N=0;
  while(1) {
    in_ff >> gw_name[N] >> dummy >> dummy 
	  >> dummy >> dff_phase[N] >> dummy >> dff_time[N] >> dummy >> dff_amp[N]
	  >> dummy >> ff_dphase50[N] >> dummy >> ff_dtime50[N] >> dummy >> ff_damp50[N]
	  >> dummy >> ff_dphase90[N] >> dummy >> ff_dtime90[N] >> dummy >> ff_damp90[N];
    if(!in_ff.good()) break;
    //cout <<" "<< gw_name[N] <<" "<< dff_phase[N] <<" "<< dff_time[N] <<" "<< dff_amp[N] << endl;
    event[N] = N+1;
    dff_time[N]/=1000.;		// sec -> msec
    ff_dtime50[N]*=1000.;	// sec -> msec
    ff_dtime90[N]*=1000.;	// sec -> msec
    N++;
  }
  int NFF=N;
  in_ff.close();

  ifstream in_of;
  in_of.open(ifname_of.Data(),ios::in);
  if (!in_of.good()) {cout << "DrawSensitivityPE.C - Error Opening File : " << ifname_of.Data() << endl;exit(1);}
  double dof_phase[MAX_GW],dof_time[MAX_GW],dof_amp[MAX_GW];
  double of_dphase50[MAX_GW],of_dtime50[MAX_GW],of_damp50[MAX_GW];
  double of_dphase90[MAX_GW],of_dtime90[MAX_GW],of_damp90[MAX_GW];
  N=0;
  while(1) {
    in_of >> gw_name[N] >> dummy >> dummy 
          >> dummy >> dof_phase[N] >> dummy >> dof_time[N] >> dummy >> dof_amp[N]
          >> dummy >> of_dphase50[N] >> dummy >> of_dtime50[N] >> dummy >> of_damp50[N]
          >> dummy >> of_dphase90[N] >> dummy >> of_dtime90[N] >> dummy >> of_damp90[N]; 
    if(!in_of.good()) break;
    //cout <<" "<< gw_name[N] <<" "<< dof_phase[N] <<" "<< dof_time[N] <<" "<< dof_amp[N] << endl;
    dof_time[N]/=1000.;		// sec -> msec
    of_dtime50[N]*=1000.;	// sec -> msec
    of_dtime90[N]*=1000.;	// sec -> msec
    N++;
  }
  int NOF=N;
  in_of.close();

  ifstream in_re;
  in_re.open(ifname_re.Data(),ios::in);
  if (!in_re.good()) {cout << "DrawSensitivityPE.C - Error Opening File : " << ifname_re.Data() << endl;exit(1);}
  double dre_phase[MAX_GW],dre_time[MAX_GW],dre_amp[MAX_GW];
  double re_dphase50[MAX_GW],re_dtime50[MAX_GW],re_damp50[MAX_GW];
  double re_dphase90[MAX_GW],re_dtime90[MAX_GW],re_damp90[MAX_GW];
  N=0;
  while(1) {
    in_re >> gw_name[N] >> dummy >> dummy
          >> dummy >> dre_phase[N] >> dummy >> dre_time[N] >> dummy >> dre_amp[N]
          >> dummy >> re_dphase50[N] >> dummy >> re_dtime50[N] >> dummy >> re_damp50[N]
          >> dummy >> re_dphase90[N] >> dummy >> re_dtime90[N] >> dummy >> re_damp90[N]; 
    if(!in_re.good()) break;
    //cout <<" "<< gw_name[N] <<" "<< dre_phase[N] <<" "<< dre_time[N] <<" "<< dre_amp[N] << endl;
    dre_time[N]/=1000.;		// sec -> msec
    re_dtime50[N]*=1000.;	// sec -> msec
    re_dtime90[N]*=1000.;	// sec -> msec
    N++;
  }
  int NRE=N;
  in_re.close();

  TCanvas* canvas = new TCanvas("sensitivity", "sensitivity", 300,40, 1000, 600);
  canvas->SetGrid();

  gStyle->SetGridStyle(3);
  gStyle->SetGridWidth(1);
  gStyle->SetFrameBorderMode(0);
  gStyle->SetTitleTextColor(kBlack);
  gStyle->SetTitleFont(12,"D");

  double xmax=0,ymax=0,ymin=1e20;
  for(int i=0;i<NFF;i++) {
    if(type=="phase") {
      if(dff_phase[i]<ymin) ymin=dff_phase[i];
      if(dff_phase[i]>ymax) ymax=dff_phase[i];
      if(dof_phase[i]<ymin) ymin=dof_phase[i];
      if(dof_phase[i]>ymax) ymax=dof_phase[i];
      if(dre_phase[i]<ymin) ymin=dre_phase[i];
      if(dre_phase[i]>ymax) ymax=dre_phase[i];
    }
    if(type=="time") {
      if(dff_time[i]<ymin) ymin=dff_time[i];
      if(dff_time[i]>ymax) ymax=dff_time[i];
      if(dof_time[i]<ymin) ymin=dof_time[i];
      if(dof_time[i]>ymax) ymax=dof_time[i];
      if(dre_time[i]<ymin) ymin=dre_time[i];
      if(dre_time[i]>ymax) ymax=dre_time[i];
    }
    if(type=="amp") {
      if(dff_amp[i]<ymin) ymin=dff_amp[i];
      if(dff_amp[i]>ymax) ymax=dff_amp[i];
      if(dof_amp[i]<ymin) ymin=dof_amp[i];
      if(dof_amp[i]>ymax) ymax=dof_amp[i];
      if(dre_amp[i]<ymin) ymin=dre_amp[i];
      if(dre_amp[i]>ymax) ymax=dre_amp[i];
    }
    if(type=="dphase50") {
      if(ff_dphase50[i]<ymin) ymin=ff_dphase50[i];
      if(ff_dphase50[i]>ymax) ymax=ff_dphase50[i];
      if(of_dphase50[i]<ymin) ymin=of_dphase50[i];
      if(of_dphase50[i]>ymax) ymax=of_dphase50[i];
      if(re_dphase50[i]<ymin) ymin=re_dphase50[i];
      if(re_dphase50[i]>ymax) ymax=re_dphase50[i];
    }
    if(type=="dtime50") {
      if(ff_dtime50[i]<ymin) ymin=ff_dtime50[i];
      if(ff_dtime50[i]>ymax) ymax=ff_dtime50[i];
      if(of_dtime50[i]<ymin) ymin=of_dtime50[i];
      if(of_dtime50[i]>ymax) ymax=of_dtime50[i];
      if(re_dtime50[i]<ymin) ymin=re_dtime50[i];
      if(re_dtime50[i]>ymax) ymax=re_dtime50[i];
    }
    if(type=="damp50") {
      if(ff_damp50[i]<ymin) ymin=ff_damp50[i];
      if(ff_damp50[i]>ymax) ymax=ff_damp50[i];
      if(of_damp50[i]<ymin) ymin=of_damp50[i];
      if(of_damp50[i]>ymax) ymax=of_damp50[i];
      if(re_damp50[i]<ymin) ymin=re_damp50[i];
      if(re_damp50[i]>ymax) ymax=re_damp50[i];
    }
    if(type=="dphase90") {
      if(ff_dphase90[i]<ymin) ymin=ff_dphase90[i];
      if(ff_dphase90[i]>ymax) ymax=ff_dphase90[i];
      if(of_dphase90[i]<ymin) ymin=of_dphase90[i];
      if(of_dphase90[i]>ymax) ymax=of_dphase90[i];
      if(re_dphase90[i]<ymin) ymin=re_dphase90[i];
      if(re_dphase90[i]>ymax) ymax=re_dphase90[i];
    }
    if(type=="dtime90") {
      if(ff_dtime90[i]<ymin) ymin=ff_dtime90[i];
      if(ff_dtime90[i]>ymax) ymax=ff_dtime90[i];
      if(of_dtime90[i]<ymin) ymin=of_dtime90[i];
      if(of_dtime90[i]>ymax) ymax=of_dtime90[i];
      if(re_dtime90[i]<ymin) ymin=re_dtime90[i];
      if(re_dtime90[i]>ymax) ymax=re_dtime90[i];
    }
    if(type=="damp90") {
      if(ff_damp90[i]<ymin) ymin=ff_damp90[i];
      if(ff_damp90[i]>ymax) ymax=ff_damp90[i];
      if(of_damp90[i]<ymin) ymin=of_damp90[i];
      if(of_damp90[i]>ymax) ymax=of_damp90[i];
      if(re_damp90[i]<ymin) ymin=re_damp90[i];
      if(re_damp90[i]>ymax) ymax=re_damp90[i];
    }
  }

  TH2F *frame = new TH2F("frame","",1000,0.0,NFF+1,1000,0.9*ymin,1.1*ymax);
  if(type=="phase")    frame->SetTitle("phase sensitivity");
  if(type=="time")     frame->SetTitle("time sensitivity");
  if(type=="amp")      frame->SetTitle("amp sensitivity");
  if(type=="dphase50") frame->SetTitle("phase-50cr sensitivity");
  if(type=="dtime50")  frame->SetTitle("time-50cr sensitivity");
  if(type=="damp50")   frame->SetTitle("amp-50cr sensitivity");
  if(type=="dphase90") frame->SetTitle("phase-90cr sensitivity");
  if(type=="dtime90")  frame->SetTitle("time-90cr sensitivity");
  if(type=="damp90")   frame->SetTitle("amp-90cr sensitivity");
  frame->SetStats(0);
  if(type=="phase")    frame->GetYaxis()->SetTitle("phase sensitivity per degree");
  if(type=="time")     frame->GetYaxis()->SetTitle("time sensitivity per msec");
  if(type=="amp")      frame->GetYaxis()->SetTitle("amp sensitivity per amp factor");
  if(type=="dphase50") frame->GetYaxis()->SetTitle("phase-50cr sensitivity (degree)");
  if(type=="dtime50")  frame->GetYaxis()->SetTitle("time-50cr sensitivity (msec)");
  if(type=="damp50")   frame->GetYaxis()->SetTitle("amp-50cr sensitivity (amp factor)");
  if(type=="dphase90") frame->GetYaxis()->SetTitle("phase-90cr sensitivity (degree)");
  if(type=="dtime90")  frame->GetYaxis()->SetTitle("time-90cr sensitivity (msec)");
  if(type=="damp90")   frame->GetYaxis()->SetTitle("amp-90cr sensitivity (amp factor)");
  frame->GetXaxis()->CenterTitle(kTRUE);
  frame->GetXaxis()->SetRangeUser(0,NFF+1);
  frame->GetYaxis()->SetTitleOffset(1.10);
  frame->GetXaxis()->CenterTitle(kTRUE);
  frame->GetYaxis()->CenterTitle(kTRUE);
  frame->GetXaxis()->SetTitleFont(132);
  frame->GetXaxis()->SetLabelFont(132);
  frame->GetYaxis()->SetTitleFont(132);
  frame->GetYaxis()->SetLabelFont(132);
  frame->GetXaxis()->SetTitleSize(0.045);
  frame->GetXaxis()->SetLabelSize(0.045);
  frame->GetYaxis()->SetTitleSize(0.045);
  frame->GetYaxis()->SetLabelSize(0.045);
  frame->GetYaxis()->SetLabelOffset(0.01);
  frame->LabelsOption("x");
  if(N<10) {
    frame->GetXaxis()->SetMoreLogLabels();
    frame->GetXaxis()->SetNoExponent();
  }
  frame->Draw();

  canvas->SetBottomMargin(0.15);
  //frame->GetXaxis()->SetNdivisions(-221);
  frame->GetXaxis()->SetNdivisions(120);
  frame->GetXaxis()->SetTitle("");
  frame->GetXaxis()->SetLabelOffset(.05);
  frame->GetXaxis()->SetLabelSize(0.03);
  frame->GetXaxis()->CenterLabels();
  for(int n=0;n<NFF;n++) frame->GetXaxis()->ChangeLabel(n+1,45,-1,-1,-1,-1,gw_name[n]);
  frame->GetXaxis()->ChangeLabel(N+1,45,-1,-1,-1,-1," ");

  double dff[MAX_GW];
  if(type=="phase")    for(int i=0;i<NFF;i++) dff[i]=dff_phase[i];
  if(type=="time")     for(int i=0;i<NFF;i++) dff[i]=dff_time[i];
  if(type=="amp")      for(int i=0;i<NFF;i++) dff[i]=dff_amp[i];
  if(type=="dphase50") for(int i=0;i<NFF;i++) dff[i]=ff_dphase50[i];
  if(type=="dtime50")  for(int i=0;i<NFF;i++) dff[i]=ff_dtime50[i];
  if(type=="damp50")   for(int i=0;i<NFF;i++) dff[i]=ff_damp50[i];
  if(type=="dphase90") for(int i=0;i<NFF;i++) dff[i]=ff_dphase90[i];
  if(type=="dtime90")  for(int i=0;i<NFF;i++) dff[i]=ff_dtime90[i];
  if(type=="damp90")   for(int i=0;i<NFF;i++) dff[i]=ff_damp90[i];

  TGraph* gdff = new TGraph(NFF,event,dff);
  gdff->SetLineWidth(2);
  gdff->SetLineColor(color[2]);
  gdff->SetMarkerColor(color[2]);
  gdff->SetMarkerStyle(20);
  gdff->Draw("LPSAME");

  double dof[MAX_GW];
  if(type=="phase")    for(int i=0;i<NFF;i++) dof[i]=dof_phase[i];
  if(type=="time")     for(int i=0;i<NFF;i++) dof[i]=dof_time[i];
  if(type=="amp")      for(int i=0;i<NFF;i++) dof[i]=dof_amp[i];
  if(type=="dphase50") for(int i=0;i<NFF;i++) dof[i]=of_dphase50[i];
  if(type=="dtime50")  for(int i=0;i<NFF;i++) dof[i]=of_dtime50[i];
  if(type=="damp50")   for(int i=0;i<NFF;i++) dof[i]=of_damp50[i];
  if(type=="dphase90") for(int i=0;i<NFF;i++) dof[i]=of_dphase90[i];
  if(type=="dtime90")  for(int i=0;i<NFF;i++) dof[i]=of_dtime90[i];
  if(type=="damp90")   for(int i=0;i<NFF;i++) dof[i]=of_damp90[i];

  TGraph* gdof = new TGraph(NFF,event,dof);
  gdof->SetLineWidth(2);
  gdof->SetLineColor(color[1]);
  gdof->SetMarkerColor(color[1]);
  gdof->SetMarkerStyle(20);
  gdof->Draw("LPSAME");

  double dre[MAX_GW];
  if(type=="phase")    for(int i=0;i<NFF;i++) dre[i]=dre_phase[i];
  if(type=="time")     for(int i=0;i<NFF;i++) dre[i]=dre_time[i];
  if(type=="amp")      for(int i=0;i<NFF;i++) dre[i]=dre_amp[i];
  if(type=="dphase50") for(int i=0;i<NFF;i++) dre[i]=re_dphase50[i];
  if(type=="dtime50")  for(int i=0;i<NFF;i++) dre[i]=re_dtime50[i];
  if(type=="damp50")   for(int i=0;i<NFF;i++) dre[i]=re_damp50[i];
  if(type=="dphase90") for(int i=0;i<NFF;i++) dre[i]=re_dphase90[i];
  if(type=="dtime90")  for(int i=0;i<NFF;i++) dre[i]=re_dtime90[i];
  if(type=="damp90")   for(int i=0;i<NFF;i++) dre[i]=re_damp90[i];

  TGraph* gdre = new TGraph(NFF,event,dre);
  gdre->SetLineWidth(2);
  gdre->SetLineColor(color[0]);
  gdre->SetMarkerColor(color[0]);
  gdre->SetMarkerStyle(20);
  gdre->Draw("LPSAME");

  TLegend *leg = new TLegend(0.1102204,0.7456446,0.3226453,0.8867596,NULL,"brNDC");
  leg->SetBorderSize(1);
  leg->SetTextAlign(22);
  leg->SetTextFont(12);
  leg->SetLineColor(1);
  leg->SetLineStyle(1);
  leg->SetLineWidth(1);
  leg->SetFillColor(0);
  leg->SetFillStyle(1001);
  leg->SetTextSize(0.04);
  leg->SetLineColor(kBlack);
  leg->SetFillColor(kWhite);
  leg->AddEntry(gdff,"Fitting Factor","lp");
  leg->AddEntry(gdof,"Overlap Factor","lp");
  leg->AddEntry(gdre,"Residual Energy","lp");
  leg->Draw();

  // dump pvalue plot 
  if(odir!="") odir=odir+"/";
  TString ofname;
  if(type=="phase")    ofname=odir+"/"+"A_phase_sensitivity.png";
  if(type=="dphase50") ofname=odir+"/"+"B_dphase_50cr_sensitivity.png";
  if(type=="dphase90") ofname=odir+"/"+"C_dphase_90cr_sensitivity.png";
  if(type=="time")     ofname=odir+"/"+"D_time_sensitivity.png";
  if(type=="dtime50")  ofname=odir+"/"+"E_dtime_50cr_sensitivity.png";
  if(type=="dtime90")  ofname=odir+"/"+"F_dtime_90cr_sensitivity.png";
  if(type=="amp")      ofname=odir+"/"+"G_amp_sensitivity.png";
  if(type=="damp50")   ofname=odir+"/"+"H_damp_50cr_sensitivity.png";
  if(type=="damp90")   ofname=odir+"/"+"I_damp_90cr_sensitivity.png";
  cout << ofname << endl;
  canvas->Print(ofname);
}

