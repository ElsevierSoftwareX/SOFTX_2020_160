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


//
// Draw median50/90 sky localization error region angle vs SNR
// Note : this macro is used to generate the PRC report (cwb_report merge_label prc)
// Author : Gabriele Vedovato
//
// input lst file names format :
// output root merge file_name              mdc_name           mdc_type   line_color   line_style
// wave_ADV_SIM_BRST_L1H1V1_run1.M1.root    WNB250_100_0d100   1       	  6            0
// note : mdc_type = type[1]+1, where type[1] is the parameter declared in the waveburst tree 
//
// command : root -l 'root 'macro/DrawMedianPRCvsSNR.C("file.lst","plot title","MEDIAN50","output_plot.png")'
//

#define APPLY_FIT

#define MAX_FILES 20		// num max of file in the input file list

//------------------------------------------------------
// Thresholds Settings                                  
//------------------------------------------------------

#define TREE_NAME "waveburst"

#define MIN_SNR_NET 		5	// SNRnet inf
#define MAX_SNR_NET 		100	// SNRnet sup
#define NSNR    		20	// number of snr bins

#define MIN_NSNR 		20	// minimum number of entries in the SNR bin
#define MAX_MEDIAN50		50	// max median50 to be included in the fit
#define MAX_MEDIAN90		100	// max median90 to be included in the fit

#define LINX				// set Logx to false

Double_t 
medianfit(Double_t *x, Double_t *par);		// fit function
int 
GetMedian(TString ifName, int type, TString pType, wavearray<double>& snr, 
          wavearray<double>& median, wavearray<double>& entries,
          float T_win, int pp_inetcc, float T_cor, int pp_irho, float T_cut, float T_vED, float T_pen, float T_ifar);

void 
DrawMedianPRCvsSNR(TString ilist_file, TString pTitle, TString pType, TString ofName,
                   float T_win=0, int pp_inetcc=0, float T_cor=0, int pp_irho=0, 
                   float T_cut=0, float T_vED=0, float T_pen=0, float T_ifar=0) {

  if((ofName!="")&&(!ofName.Contains(".png"))) {
    cout << "Error Output File : " << ilist_file.Data() << endl;
    cout << "Must have .png extension" << endl;
    exit(1);
  } else {
    ofName.ReplaceAll(".png",".gif");
  }

  TString ptitle="Sky Localization : median error area vs SNR";

  if((pType!="MEDIAN50")&&(pType!="MEDIAN90")&&(pType!="")) {
    ptitle=pType;
  } else {
    if(pType=="") pType="MEDIAN50";
    if(pType=="MEDIAN50")   ptitle=ptitle+"50% error region";
    if(pType=="MEDIAN90")   ptitle=ptitle+"90% error region";
  }

  TString ifile[MAX_FILES];
  TString name[MAX_FILES];
  int type[MAX_FILES];
  int colors[MAX_FILES];
  int lstyle[MAX_FILES];

  // Open list
  ifstream in;
  in.open(ilist_file.Data(),ios::in);
  if (!in.good()) {cout << "Error Opening File : " << ilist_file.Data() << endl;exit(1);}

  int size=0;
  char str[1024];
  int fpos=0;
  while(true) {
    in.getline(str,1024);
    if (!in.good()) break;
    if(str[0] != '#') size++;
  }
  cout << "size " << size << endl;
  in.clear(ios::goodbit);
  in.seekg(0, ios::beg);
  if (size==0) {cout << "Error : File " << ilist_file.Data() << " is empty" << endl;exit(1);}

  char sfile[256];
  char sname[256];
  int stype;
  int scolors;
  int sltyle;
  int NFILE=0;
  while(true) {
    in >> sfile >> sname >> stype >> scolors >> sltyle;
    if(!in.good()) break;
    if(sfile[0]=='#') continue;
    ifile[NFILE]=sfile;
    name[NFILE]=sname;
    type[NFILE]=stype;
    colors[NFILE]=scolors;
    lstyle[NFILE]=sltyle;
    cout << NFILE+1 << " " << ifile[NFILE] << " " << name[NFILE] << " " 
         << type[NFILE] << " " << colors[NFILE] << " " << lstyle[NFILE] << endl;
    NFILE++;
  }
  in.close();

#ifdef APPLY_FIT
  TF1* f1 = NULL;
  f1=new TF1("medianfit",medianfit,0.5,4,3);
  f1->SetParNames("A","B","C");
  f1->SetLineWidth(2);
  f1->SetNpx(10000);
#endif

  // get medians
  double max_median=0;
  int nevents[MAX_FILES];
  wavearray<double>* snr = new wavearray<double>[MAX_FILES];
  wavearray<double>* median = new wavearray<double>[MAX_FILES];
  wavearray<double>* entries = new wavearray<double>[MAX_FILES];
  for(int n=0;n<NFILE;n++) {
    cout << "Process File : " << ifile[n].Data() << endl;
    nevents[n]=GetMedian(ifile[n], type[n], pType, snr[n], median[n], entries[n],
                         T_win, pp_inetcc, T_cor, pp_irho, T_cut, T_vED, T_pen, T_ifar);
    // find max median
    for(int i=0;i<snr[n].size();i++) {
      if(median[n][i]>max_median) max_median=median[n][i];
      //cout << n << " " << i << "\t snr " << snr[n][i] << "\t median " << median[n][i] << endl;
    }
    //cout << n << " max_median " << max_median << endl;
  }

  // create plots
  gStyle->SetFrameBorderMode(0);     // remove the red box around canvas
  gROOT->ForceStyle();               

  gStyle->SetTitleFont(72);
  gStyle->SetMarkerColor(50);
  gStyle->SetLineColor(kWhite);
  gStyle->SetTitleW(0.98);     
  gStyle->SetTitleH(0.05);     
  gStyle->SetTitleY(0.98);     
  gStyle->SetFillColor(kWhite);
  gStyle->SetLineColor(kWhite);
  gStyle->SetTitleFont(12,"D");

  TCanvas *canvas = new TCanvas("median", "median", 300,40, 600, 600);
  canvas->Clear();                                                       
  canvas->ToggleEventStatus();                                           
#ifdef LINX
  canvas->SetLogx(false);                                                     
#else
  canvas->SetLogx(true);         
#endif                                            
  //canvas->SetLogy();                                                     
  canvas->SetGridx();                                                    
  canvas->SetGridy();                                                    
  canvas->SetFillColor(kWhite);                                          

  double A[MAX_FILES];
  double B[MAX_FILES];
  double C[MAX_FILES];
  double chi2[MAX_FILES];
  //TGraph* gr[MAX_FILES];
  //for(int n=0;n<NFILE;n++) gr[n] = new TGraph(snr[n].size(),snr[n].data,median[n].data);
  TGraphErrors* gr[MAX_FILES];
  for(int n=0;n<NFILE;n++) {
    wavearray<double> esnr=snr[n];esnr=0; 
    wavearray<double> emedian=median[n];emedian=0;
    //for(int i=0;i<emedian.size();i++) emedian[i]=1./sqrt(entries[n][i]); 
    gr[n] = new TGraphErrors(snr[n].size(),snr[n].data,median[n].data,esnr.data,emedian.data);
  }
  for(int n=0;n<NFILE;n++) {                                            
    gr[n]->SetLineWidth(0);                                            
    gr[n]->SetLineColor(colors[n]);                                    
    gr[n]->SetLineStyle(1);                                    
    gr[n]->SetMarkerColor(colors[n]);                                  
    gr[n]->SetMarkerStyle(20);                                  
    gr[n]->SetMarkerSize(1);                                  

#ifdef APPLY_FIT
    f1->SetLineColor(colors[n]);
    f1->SetLineStyle(lstyle[n]);
//    f1->SetLineWidth(1);
    f1->SetNpx(100);
    gr[n]->Fit("medianfit");
    chi2[n]=f1->GetChisquare();
    A[n]=f1->GetParameter(0);
    B[n]=f1->GetParameter(1);
    C[n]=f1->GetParameter(2);
#endif
  }                                                                    

  TMultiGraph* mg = new TMultiGraph();
  TString gTitle = "Sky Localization "+pType+" vs SNR";
  if(pTitle!="") gTitle = gTitle+" : "+pTitle; 
  mg->SetTitle(gTitle.Data()); 
  for(int n=0;n<NFILE;n++) mg->Add(gr[n]);  
  mg->Paint("AP");                        

  mg->GetHistogram()->GetXaxis()->SetRangeUser(MIN_SNR_NET,MAX_SNR_NET);
//  mg->GetHistogram()->GetYaxis()->SetRangeUser(0,max_median);

  mg->GetHistogram()->GetXaxis()->SetTitle("Injected SNRnet");
  if(pType=="MEDIAN50") mg->GetHistogram()->GetYaxis()->SetTitle("Median50 error angle (degrees)");
  else                  mg->GetHistogram()->GetYaxis()->SetTitle("Median90 error angle (degrees)");

  mg->GetHistogram()->GetXaxis()->SetLabelSize(0.04);
  mg->GetHistogram()->GetYaxis()->SetLabelSize(0.04);
  mg->GetHistogram()->GetXaxis()->SetTitleSize(0.04);
  mg->GetHistogram()->GetYaxis()->SetTitleSize(0.04);
  mg->GetHistogram()->GetXaxis()->SetTitleOffset(1.2);
  mg->GetHistogram()->GetYaxis()->SetTitleOffset(1.2);
  mg->GetHistogram()->GetXaxis()->CenterTitle(true);
  mg->GetHistogram()->GetYaxis()->CenterTitle(true);
  mg->GetHistogram()->GetXaxis()->SetLabelFont(42);
  mg->GetHistogram()->GetXaxis()->SetTitleFont(42);
  mg->GetHistogram()->GetYaxis()->SetLabelFont(42);
  mg->GetHistogram()->GetYaxis()->SetTitleFont(42);
  mg->Draw("AP");

  // draw the legend

  TLegend* leg;                
  double hleg = 0.8-NFILE*0.05; 
  leg = new TLegend(0.5793173,hleg,0.9915385,0.8721805,NULL,"brNDC");

  leg->SetBorderSize(1);
  leg->SetTextAlign(22);
  leg->SetTextFont(12); 
  leg->SetLineColor(1); 
  leg->SetLineStyle(1); 
  leg->SetLineWidth(1); 
  leg->SetFillColor(0); 
  leg->SetFillStyle(1001);
  leg->SetTextSize(0.05); 
  leg->SetLineColor(kBlack);
  leg->SetFillColor(kWhite);

  for(int n=0;n<NFILE;n++) {
    double median_snr15 = A[n]+pow(fabs(B[n])/15.,1)+pow(C[n]/15.,2);
    char legLabel[256];
    //sprintf(legLabel,"%s : %2.1f (%g)",name[n].Data(),median_snr10,nevents[n]);
    sprintf(legLabel,"SNR(15) = %2.1f",median_snr15);
    leg->AddEntry(gr[n],legLabel,"lp");
  }
  leg->Draw();

  // save plot

  if(ofName!="") {
    TString gfileName=ofName;
    canvas->Print(gfileName);
    TString pfileName=gfileName;
    pfileName.ReplaceAll(".gif",".png");
    char cmd[1024];
    sprintf(cmd,"convert %s %s",gfileName.Data(),pfileName.Data());
    cout << cmd << endl;
    gSystem->Exec(cmd);
    sprintf(cmd,"rm %s",gfileName.Data());
    cout << cmd << endl;
    gSystem->Exec(cmd);
    //gSystem->Exit(0);
  }
}

Double_t 
medianfit(Double_t *x, Double_t *par) {
  // constrain par[1] to be >=0
  double y = par[0]+pow(fabs(par[1])/x[0],1)+pow(par[2]/x[0],2);
  return y;
}

int 
GetMedian(TString ifName, int type, TString pType, wavearray<double>& snr, 
          wavearray<double>& median, wavearray<double>& entries,
          float T_win, int pp_inetcc, float T_cor, int pp_irho, float T_cut, float T_vED, float T_pen, float T_ifar) {

  // open input tree 
  TFile *file0 = TFile::Open(ifName);
  cout << ifName.Data() << endl;     
  if (file0==NULL) {cout << "Error Opening file" << endl;exit(1);}
  TTree* itree = (TTree *) gROOT->FindObject(TREE_NAME);          
  if (itree==NULL) {cout << "Error Opening tree" << endl;exit(1);}

  // get detector list
  TList* list = itree->GetUserInfo();
  int nIFO=list->GetSize();          
  if (nIFO==0) {cout << "Error : no ifo present in the tree" << endl;exit(1);}
  for (int n=0;n<list->GetSize();n++) {                                       
    detector* pDetector;                                                      
    pDetector = (detector*)list->At(n);                                       
    detectorParams dParams = pDetector->getDetectorParams();                  
    //pDetector->print();                                                       
  }                                                                           

  // define selection cuts
  char cut[1024];         
  char tmp[1024];         
  sprintf(cut,"abs(time[0]-time[%d])<%f && netcc[%d]>%f && rho[%d]>%f",
          nIFO,T_win,pp_inetcc,T_cor,pp_irho,T_cut);
  if(T_vED>0)  {strcpy(tmp,cut);sprintf(cut,"%s && neted[0]/ecor<%f",tmp,T_vED);}
  if(T_pen>0)  {strcpy(tmp,cut);sprintf(cut,"%s && penalty>%f",tmp,T_pen);}
  if(type>0)   {strcpy(tmp,cut);sprintf(cut,"%s && type[1]==%d",tmp,type);}
  if(T_ifar>0) {strcpy(tmp,cut);sprintf(cut,"%s && ifar>(24.*3600.*365.)*%f",tmp,T_ifar);}

#ifdef LINX
  double fsnr = (MAX_SNR_NET-MIN_SNR_NET)/NSNR;
#else
  double fsnr = TMath::Exp(TMath::Log(MAX_SNR_NET-MIN_SNR_NET)/NSNR);
#endif
  entries.resize(NSNR);    
  median.resize(NSNR);    
  snr.resize(NSNR);    
  int nsnr=0;
  double max_median = (pType=="MEDIAN50") ? MAX_MEDIAN50 : MAX_MEDIAN90;
  double infsnr=MIN_SNR_NET;
#ifdef LINX
  double supsnr=fsnr+infsnr;
#else
  double supsnr=fsnr*infsnr;
#endif
  char SNRnet[256]="";
  strcpy(tmp,"iSNR[0]");
  for(int i=1;i<nIFO;i++) {sprintf(SNRnet,"%s+iSNR[%d]",tmp,i);strcpy(tmp,SNRnet);}
  sprintf(SNRnet,"sqrt(%s)",tmp);
  int nevents=0;		
  for(int n=0;n<NSNR;n++) {            

    char icut[512];
    sprintf(icut,"%s && %s>%f && %s<=%f",cut,SNRnet,infsnr,SNRnet,supsnr);
    //sprintf(icut,"%s && sqrt(likelihood/%d)>%f && sqrt(likelihood/%d)<=%f",cut,nIFO,infsnr,nIFO,supsnr);
    //cout << icut << endl;
    itree->Draw("erA[0]",icut,"goff");

    int size = itree->GetSelectedRows();
    nevents+=size;
    double* erA = itree->GetV1();

    if(size>MIN_NSNR) { 
      int* index = new int[size];
      TMath::Sort(size,erA,index,false);

      entries[nsnr] = size;
      median[nsnr]  = (pType=="MEDIAN50") ? erA[index[int(size*0.5)]] : erA[index[int(size*0.9)]];	// median
      snr[nsnr]     = (supsnr+infsnr)/2.;		// SNRnet
      if(median[nsnr]<max_median) {
        cout << "Selected entries : " << entries[nsnr] << "\t SNRnet : " 
             << snr[nsnr] << "\t median : " << median[nsnr] << endl;
        nsnr++;
      }

      delete [] index;
    }
    
#ifdef LINX
    infsnr += fsnr;
    supsnr += fsnr;
#else
    infsnr *= fsnr;
    supsnr *= fsnr;
#endif
  }

  entries.resize(nsnr);    
  median.resize(nsnr);    
  snr.resize(nsnr);    

  return nevents;		// return total number of selected events
}

