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


#define MAX_EVT  	100		// max number of events

#define NLOOP	 	100000		// used to compute the multiple trial p-value

//#define SHOW_OFFSOURCE_DISTRIBUTION

//#define CHANGE_LABELS

#define DRAW_SIGMA_MEDIAN

#define DRAW_NAMES	4		// max number of displayed names

#define YMIN		0.001		// minimum y-axis p-value 

#define SIGMA1    	0.682689492137
#define SIGMA2    	0.954499736104
#define SIGMA3    	0.997300203937

double GetProbability(int N, int k, double p);
double InverseProbability(int N, int k, double p);

double InverseProbabilityBH(int N, int k, double p);
TGraph* DrawProbabilityBH(int N, double p, Color_t color, int style=9);
TGraph* DrawStepGraph(int N, double* x, double* y, Color_t color=kBlack);
TGraphAsymmErrors* DrawBelt(int N, double pinf, double psup, Color_t color);


void DrawPvaluesPE(TString ifname, TString odir, TString side, TString label="", bool title=true, double FDR=0) {

  if(FDR>0) {
    if(FDR>1)             {cout << endl << "DrawPvaluesPE Error, FDR must be [0,1]" << endl << endl;exit(1);}
    if(fmod(10*FDR,1)!=0) {cout << endl << "DrawPvaluesPE Error, FDR must be a multiple of 0.1" << endl << endl;exit(1);}
  }

  gRandom->SetSeed(150914);

  // init blind colors 
  Color_t color[4];
  color[0] = CWB::Toolbox::getTableau10BlindColor("DarkOrange1");
  color[1] = CWB::Toolbox::getTableau10BlindColor("DeepSkyBlue4");
  color[2] = CWB::Toolbox::getTableau10BlindColor("DarkGray");
  color[3] = CWB::Toolbox::getTableau10BlindColor("SandyBrown");

// -----------------------------> READ P-VALUES

  // read values from ifname file
  ifstream in;
  in.open(ifname.Data(),ios::in);
  if (!in.good()) {cout << "DrawPvaluesPE Error Opening File : " << ifname.Data() << endl;exit(1);}
  int NEVT=0;
  char   gw_name[MAX_EVT][256];
  double median[MAX_EVT],l99[MAX_EVT],l90[MAX_EVT],l50[MAX_EVT],u50[MAX_EVT],u90[MAX_EVT],u99[MAX_EVT];
  double match[MAX_EVT],pvalue[MAX_EVT];
  char   dummy[256];
  cout << endl;
  cout << "READ EVENTS" << endl;
  cout << endl;
  while(1) {
    in >> gw_name[NEVT] >> dummy >> dummy >> match[NEVT] >> dummy >> pvalue[NEVT] >> dummy >> median[NEVT]
       >> dummy >> l99[NEVT] >> dummy >> l90[NEVT] >> dummy >> l50[NEVT] >> dummy >> u50[NEVT] >> dummy >> u90[NEVT] >> dummy >> u99[NEVT];
    if(!in.good()) break;
    cout << NEVT <<" "<< gw_name[NEVT] <<" "<< match[NEVT] <<" "<< pvalue[NEVT] << endl;
    NEVT++;
  }
  cout << endl;
  in.close();

// -----------------------------> SORT P-VALUES

  // sort pvalues
  int index[MAX_EVT];
  TMath::Sort(NEVT,pvalue,index,false);

  double pval_min=1e20;
  TString sside[MAX_EVT];
  double sindex[MAX_EVT];
  double spvalue[MAX_EVT];
  TString sgw_name[MAX_EVT];
  cout << endl;
  for(int n=0;n<NEVT;n++) {
    int k=index[n];
    sindex[n]=n+1;
    spvalue[n]=pvalue[k];
    sgw_name[n]=gw_name[k];
  }
  cout << endl;
  cout << "SORT EVENTS" << endl;
  cout << endl;
  for(int n=0;n<NEVT;n++) cout << n << "\t" << sindex[n] << "\t" << sgw_name[n] << "\t" << spvalue[n] << endl;;
  cout << endl;

// -----------------------------> DRAW P-VALUES

  // create plots
  gStyle->SetFrameBorderMode(0);     // remove the red box around canvas
  gROOT->ForceStyle();
  gStyle->SetTitleFont(12,"D");
  TCanvas *canvas = new TCanvas("pvalue", "pvalue", 300,40, 1000, 600);
  canvas->Clear();
  canvas->SetLogx();
  canvas->SetLogy();
  canvas->SetGridx();
  canvas->SetGridy();
  canvas->SetFillColor(kWhite);

  // draw pvalues
  TGraph* gr = new TGraph(NEVT,sindex,spvalue);
  if(title) gr->SetTitle("p-values distribution ("+label+")"); else gr->SetTitle("");
  //gr->GetXaxis()->SetRangeUser(1,double(NEVT));
  //double xMAX = (1+NEVT/int(pow(10,int(log10(NEVT)))))*int(pow(10,int(log10(NEVT))));
  double xMAX = NEVT;
  gr->GetHistogram()->GetXaxis()->SetLimits(1,xMAX);
  gr->GetYaxis()->SetRangeUser(spvalue[0],1.0);
  if(spvalue[0]<YMIN) gr->GetYaxis()->SetRangeUser(YMIN,1.0);
  gr->GetXaxis()->SetTitle("Cumulative number of events");
  gr->GetYaxis()->SetTitle("p-value");
  gr->GetYaxis()->SetRangeUser(0,1);
  gr->GetXaxis()->SetTitleOffset(0.80);
  gr->GetYaxis()->SetTitleOffset(0.80);
  gr->GetXaxis()->CenterTitle(kTRUE);
  gr->GetYaxis()->CenterTitle(kTRUE);
  gr->GetXaxis()->SetTitleFont(132);
  gr->GetXaxis()->SetLabelFont(132);
  gr->GetYaxis()->SetTitleFont(132);
  gr->GetYaxis()->SetLabelFont(132);
  gr->GetXaxis()->SetTitleSize(0.045);
  gr->GetXaxis()->SetLabelSize(0.045);
  gr->GetYaxis()->SetTitleSize(0.045);
  gr->GetYaxis()->SetLabelSize(0.045);
  gr->GetYaxis()->SetLabelOffset(0.01);
  gr->GetYaxis()->SetNdivisions(3);

  if(NEVT<10) {
    gr->GetXaxis()->SetMoreLogLabels();
    gr->GetXaxis()->SetNoExponent();
  }
  gr->SetLineColor(kGray+1);
  gr->SetLineWidth(1);
  gr->SetMarkerSize(1.3);
  gr->SetMarkerColor(kGreen+2);
//  gr->SetMarkerColor(color[0]);
  gr->SetMarkerStyle(20);
  gr->Draw("APL");

// -----------------------------> SET GRID COLOR

#ifdef SET_GRID_COLOR
  // set grid color
  gr->GetXaxis()->SetAxisColor(16);
  gr->GetYaxis()->SetAxisColor(16);

  // fix frame lines color
  double xmin = gr->GetHistogram()->GetXaxis()->GetXmin();
  double xmax = gr->GetHistogram()->GetXaxis()->GetXmax();
  //double ymin = pow(10,gPad->GetUymin());
  double ymin = 0.0;
  double ymax = gPad->GetUymax();
  cout << "xmin=" << xmin << " xmax=" << xmax << endl;
  cout << "ymin=" << ymin << " ymax=" << ymax << endl;

  TLine* frameLine[4];
  frameLine[0] = new TLine(xmin,ymin,xmax,ymin);
  frameLine[1] = new TLine(xmin,ymin,xmin,ymax);
  frameLine[2] = new TLine(xmax,ymin,xmax,ymax);
  frameLine[3] = new TLine(0.,ymax,xmax,ymax);
  for(int i=0;i<2;i++) frameLine[i]->Draw();
#endif

// -----------------------------> CHANGE X LABELS

#ifdef CHANGE_LABELS
  canvas->SetBottomMargin(0.15);
  //gr->GetHistogram()->GetXaxis()->SetNdivisions(-221);
  gr->GetHistogram()->SetMinimum(YMIN);
  gr->GetXaxis()->SetNdivisions(120);
  gr->GetXaxis()->SetTitle("");
  gr->GetXaxis()->SetLabelOffset(.05);
  gr->GetXaxis()->SetLabelSize(0.03);
  gr->GetXaxis()->CenterLabels();
  for(int n=0;n<NEVT;n++) gr->GetXaxis()->ChangeLabel(n+1,45,-1,-1,-1,-1,sgw_name[n]);
  gr->GetXaxis()->ChangeLabel(NEVT+1,45,-1,-1,-1,-1," ");
#endif

// -----------------------------> DRAW MEDIAN AND CONFIDENCE REGION

#ifdef DRAW_SIGMA_MEDIAN
  TGraph* gmean    = DrawProbabilityBH(NEVT, 0.5, kRed, 9);
//  TGraph* gmean    = DrawProbabilityBH(NEVT, 0.5, color[0], 9);
//  TGraph* glsigma1 = DrawProbabilityBH(NEVT, (1-SIGMA1)/2, 14);
//  TGraph* gusigma1 = DrawProbabilityBH(NEVT, (1+SIGMA1)/2, 14);
//  TGraph* glsigma2 = DrawProbabilityBH(NEVT, (1-SIGMA2)/2, 16);
//  TGraph* gusigma2 = DrawProbabilityBH(NEVT, (1+SIGMA2)/2, 16);
//  TGraph* glsigma3 = DrawProbabilityBH(NEVT, (1-SIGMA3)/2, 18);
//  TGraph* gusigma3 = DrawProbabilityBH(NEVT, (1+SIGMA3)/2, 18);

  // draw 90% confidence 
  TGraph* gcl90inf = DrawProbabilityBH(NEVT, 0.05, kBlue);
  TGraph* gcl90sup = DrawProbabilityBH(NEVT, 0.95, kBlue);

/*
  TGraphAsymmErrors* gbelt3 = DrawBelt(NEVT, (1-SIGMA3)/2, (1+SIGMA3)/2, 18);
  TGraphAsymmErrors* gbelt2 = DrawBelt(NEVT, (1-SIGMA2)/2, (1+SIGMA2)/2, 17);
  TGraphAsymmErrors* gbelt1 = DrawBelt(NEVT, (1-SIGMA1)/1, (1+SIGMA1)/2, 16);
*/

  // draw again events marker, must be in first level
  gr->Draw("P");
#endif

// -----------------------------> FDR

  TF1 *fdr[30];
  int fdr_id=-1;
  if(FDR>0) {
    gr->GetYaxis()->SetNdivisions(10);
    for(int n=1;n<=30;n++) {
      double alpha;
      if(n>1  && n<=10) alpha=n*YMIN;
      if(n>10 && n<=20) alpha=(n-10)*0.01;
      if(n>20 && n<=30) alpha=(n-20)*0.1;

      fdr[n-1] = new TF1("func", "[0]*x", 1,MAX_EVT);
      fdr[n-1]->SetNpx(10000);
      fdr[n-1]->SetParameter(0,alpha/NEVT);
      fdr[n-1]->SetLineWidth(1);
      fdr[n-1]->SetLineColor(kGray);
      fdr[n-1]->SetLineStyle(1);
      //cout << n << " ALPHA " << alpha << endl;
      //if((alpha==0.01)||(alpha==0.1)) fdr[n-1]->Draw("SAME");	// draw main FRD lines
      if(fabs(alpha-FDR)<1.e-6) {				// draw main FRD threshold
        //fdr->SetLineColor(kBlack);
        fdr[n-1]->SetLineStyle(1);
        fdr[n-1]->Draw("SAME");
        fdr_id=n-1;
      }
    }

    // find events below the FDR threshold
    int nmin;
    int nfdr=0;
    double prob_min=1e20;
    double cl1[MAX_EVT];
    double prob[MAX_EVT];
    cout << endl;
    cout << "FDR EVENTS" << endl;
    cout << endl;
    for(int n=0;n<NEVT;n++) {
      if(spvalue[n]>=fdr[fdr_id]->Eval(n+1)) continue;
      prob[n] = GetProbability(NEVT,n+1,spvalue[n]);
      cl1[n] = InverseProbability(NEVT,n+1,0.01);
      cout << n << "\t" << sgw_name[n] << "\tpvalue = " << spvalue[n] << "\tprob = " << prob[n] << "\tcl1 = " << cl1[n] << endl;
      if(prob[n]<prob_min) {prob_min=prob[n];nmin=n;}
      nfdr++;
    }

    if(nfdr>0) {

      cout << endl;
      cout << nmin << " " << sgw_name[nmin] << "\tpvalue = " << spvalue[nmin] 
           << "\tprob = " << prob_min << "\tlog10(prob_min) " << log10(prob_min) << endl;
      cout << endl;

      double pmin[NLOOP];
      for(int l=0;l<NLOOP;l++) {	// repeat the generation of NEVT

        // generate NEVT pvalues 
        for(int i=0;i<NEVT;i++) pvalue[i]=gRandom->Uniform(0,1);

        // sort pvalues
        TMath::Sort(NEVT,pvalue,index,false);

        pmin[l]=1e20;
        for(int n=0;n<NEVT;n++) {
          int k=index[n];
          if(pvalue[k]>=fdr[fdr_id]->Eval(n+1)) continue;
          prob[n] = GetProbability(NEVT,n+1,pvalue[k]);
          if(prob[n]<pmin[l]) pmin[l]=prob[n];
        }
      }

      // compute ntrials p-value
      double PVALUE=0;
      for(int l=0;l<NLOOP;l++) {	// repeat the generation of NtRIALS
         if(pmin[l]<prob_min) PVALUE+=1;
      }
      PVALUE/=NLOOP;
      cout << endl;
      cout << "PVALUE " << PVALUE << endl;
      cout << endl;

      if(title) gr->SetTitle(TString(gr->GetTitle())+TString::Format(" - significance (#alpha=%.1f): %.4f",FDR,PVALUE)); else gr->SetTitle("");

#ifdef SHOW_OFFSOURCE_DISTRIBUTION
      TH1F* hist = new TH1F("hist","hist",100,1e-5,0);
      for(int l=0;l<NLOOP;l++) {	// repeat the generation of NtRIALS
         hist->Fill(log10(pmin[l]));
      }
      hist->Draw();

      double xline[2]={log10(prob_min),log10(prob_min)};
      double yline[2]={0,NLOOP};
      TGraph* line = new TGraph(2,xline,yline);
      line->SetLineColor(kRed);
      line->Draw("same");

      gStyle->SetLineColor(kBlack);
      canvas->SetLogx(false);
      return;
#endif
    }

    // draw alpha lbels
    TF1 *f3 = new TF1("f3","log10(x)",YMIN,0.999);
    TGaxis *A3 = new TGaxis(double(NEVT),0.0,double(NEVT),1.0,"f3",505,"=SG");
    A3->SetTitle("#alpha");
    A3->CenterTitle(kTRUE);
    A3->SetTitleFont(132);
    A3->SetLabelFont(132);
    A3->SetLabelSize(0.045);
    A3->SetTickSize(0.02);
    A3->SetLabelOffset(0.03);
    A3->SetTitleSize(0.045);
    A3->SetTitleOffset(0.4);
    A3->Draw();
    canvas->Update();
  }

// -----------------------------> DRAW GW NAMES

#ifdef DRAW_NAMES
  float xrange = NEVT;
  float yrange = spvalue[0]==0 ? YMIN : spvalue[0];

  TBox* box;
  TLatex* latex;
  TArrow* arrow;
  cout << endl;
  for(int n=0;n<DRAW_NAMES;n++) {
  
    float xfactor = 1 + (0.05 * 10./xrange);
    float yfactor =1 - (0.30 * log10(yrange)/-3.0);
    double xpos = sindex[n]*xfactor; 
    double ypos = spvalue[n]*yfactor; 
    if(ypos<YMIN) ypos=1.1*(2-yfactor)*YMIN;
    //cout << n << " " << sindex[n] << " " << sgw_name[n] << "\t" << spvalue[n] << "\t" << xpos << "\t" << ypos << endl;

    // draw white box
    box = new TBox(xpos, ypos*(2-yfactor), xpos*1.33*xfactor, ypos*(2-yfactor)*yfactor);
    box->Draw();
    box->SetFillColor(kWhite);

    // draw label
    latex = new TLatex(xpos, ypos, sgw_name[n].Data());
    latex->SetTextFont(132);
    latex->SetTextSize(0.035);
    latex->SetTextColor(kBlack);
    latex->Draw();

    // draw arrow for events with pvalue<YMIN
    ypos = spvalue[n]*yfactor; 
    if(ypos<YMIN) {
      xpos = xpos*1.14*xfactor;
      ypos = YMIN;
      arrow = new TArrow(xpos,ypos*(2-yfactor),xpos,ypos,0.05,"|>");
      arrow->SetAngle(60);
      arrow->SetLineWidth(2);
      arrow->SetLineColor(kBlue);
      arrow->SetFillColor(kBlue);
      arrow->SetArrowSize(0.012);
      arrow->Draw();
    }
  }
  cout << endl;
#endif

// -----------------------------> LEGEND

  // draw legend
  TLegend *leg = new TLegend(0.501002,0.1184669,0.8897796,0.3501742,NULL,"brNDC");
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
  leg->AddEntry(gr,"observed","lp");
  leg->AddEntry(gmean,"null hypothesis","lp");
//  leg->AddEntry(glsigma1,"1 sigma","lp");
//  leg->AddEntry(glsigma2,"2 sigma","lp");
  //leg->AddEntry(gcl,"1% confidence","lp");
  leg->AddEntry(gcl90inf,"90% confidence region","lp");
  if(fdr_id>=0) leg->AddEntry(fdr[fdr_id],TString::Format("FDR (#alpha = %.1f)",FDR),"lp");
  leg->Draw();

  // dump pvalue plot 
  if(odir!="") odir=odir+"/";
  TString ofname=odir+"/"+gSystem->BaseName(ifname);
  TString tag = TString::Format("_pvalue.png");
  ofname.ReplaceAll(".txt",tag);
  cout << ofname << endl;
  
  canvas->Print(ofname);

//  exit(0);
}

double GetProbability(int N, int k, double p) {

  double prob=1;
  for(int i=0;i<k;i++) {
    prob-=TMath::Binomial(N,i)*pow(p,i)*pow(1-p,N-i);
  }
  return prob;
}

double InverseProbability(int N, int k, double p) {

#define LMAX		100
#define PRECISION       1.e-6

  int L=100;

  double x;
  double xmin=0;
  double xmax=1;

  double a=0.;
  for(int i=0;i<LMAX;i++) {
    x = gRandom->Uniform(xmin,xmax);
    //cout << i << "\txmin " << xmin << "\txmax " << xmax << "\tx " << x << endl; 
    a=GetProbability(N, k, x);
    if(a<p) xmin=x; else xmax=x;
    if(fabs(a-p)<PRECISION) break;
  }
  return x;
}

double InverseProbabilityBH(int N, int k, double p) {

#define LMAX		100
#define PRECISION       1.e-6

  int L=100;

  double x;
  double xmin=0;
  double xmax=1;

  double a=0.;
  for(int i=0;i<LMAX;i++) {
    x = gRandom->Uniform(xmin,xmax);
    //cout << i << "\txmin " << xmin << "\txmax " << xmax << "\tx " << x << endl; 
    a=ROOT::Math::inc_beta(x, k, N-k+1);
    if(a<p) xmin=x; else xmax=x;
    if(fabs(a-p)<PRECISION) break;
  }
  return x;
}

TGraph* DrawProbabilityBH(int N, double p, Color_t color, int style=9) {
// BENJAMINI-HOCHBERG

  double* x = new double[N];
  for(int i=0;i<N;i++) x[i]=i+1;

  double* prob = new double[N];
  for(int k=1;k<=N;k++) prob[k-1]=InverseProbabilityBH(N, k, p);

  TGraph* gprob = new TGraph(N,x,prob);
  gprob->SetLineColor(color);
//  gprob->SetLineWidth(2);
  gprob->SetLineStyle(style);
  gprob->Draw("same");
  gprob->GetXaxis()->SetRangeUser(1,N);
  gprob->GetYaxis()->SetRangeUser(0,1);

  return gprob;
}

TGraph* DrawStepGraph(int N, double* x, double* y, Color_t color) {

  int M=N-N%2;

  double* X = new double[2*M+N%2];
  double* Y = new double[2*M+N%2];

  for(int i=0;i<M;i++) {
    X[2*i+0]=x[i+0]; Y[2*i+0]=y[i];
    X[2*i+1]=x[i+1]; Y[2*i+1]=y[i];
  }
  if(N%2) X[2*M]=x[N-1]; Y[2*M]=y[N-1];

  TGraph* gprob = new TGraph(2*M+N%2,X,Y);
  gprob->SetLineColor(color);
  gprob->SetLineWidth(2);
  gprob->Draw("same");

  return gprob;
}

TGraphAsymmErrors* DrawBelt(int N, double pinf, double psup, Color_t color) {

  double* x = new double[N];
  for(int i=0;i<N;i++) x[i]=i+1;
  double* ex = new double[N];
  for(int i=0;i<N;i++) ex[i]=0;

  double* y = new double[N];
  for(int k=1;k<=N;k++) y[k-1]=InverseProbabilityBH(N, k, 0.5);

  double* yinf = new double[N];
  for(int k=1;k<=N;k++) yinf[k-1]=y[k-1]-InverseProbabilityBH(N, k, pinf);

  double* ysup = new double[N];
  for(int k=1;k<=N;k++) ysup[k-1]=InverseProbabilityBH(N, k, psup)-y[k-1];

  TGraphAsymmErrors* egr;
  egr = new TGraphAsymmErrors(N,x,y,ex,ex,yinf,ysup);
  egr->SetMarkerStyle(20);
  egr->SetMarkerSize(0);
  egr->SetLineColor(15);
  egr->SetFillColorAlpha(color, 0.1);
//  egr->SetFillStyle(1001);
//  egr->SetFillColor(color);    
  egr->Draw("3SAME");

  return egr;
}

