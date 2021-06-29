//
// Draw "Monster Event Display" vs "NN Frame Display"
// Author : Gabriele Vedovato
//

#define nIFO 3

#define WATPLOT                          // monster event display
#define PLOT_FRAME                       // frame display

TH2F* hist2D=NULL;
TCanvas* canvas;  
watplot* WTS;

void PlotFrame(std::vector<double>* nn_frame, int cid, bool save);

void DisplayFrame(TString ifName, int ID, bool save=false) {		// ID = entry to be displayed 

#ifdef WATPLOT
    WTS = new watplot(const_cast<char*>("wts"));
#endif

   TChain sigTree("waveburst");
   sigTree.Add(ifName);
   netevent signal(&sigTree,nIFO);
   int sig_entries = signal.GetEntries();
   cout << "sig entries : " << sig_entries << endl;
   std::vector<double>* sig_frame = new vector<double>;
   signal.fChain->SetBranchAddress("nnframe", &sig_frame);
   netcluster* pwc = new netcluster();
   signal.fChain->SetBranchAddress("netcluster",&pwc);

   float sig_rho;
   float sig_netcc;
   float sig_likelihood;

   for(int i=0; i<sig_entries; i++) {

     if(i!=ID) continue;

     signal.GetEntry(i);
     //signal.Show(i);

     sig_rho = signal.rho[0];
     sig_netcc = signal.netcc[1];                      
     sig_likelihood = signal.likelihood;
     cout << sig_frame->size() << endl; 
     //for(int j=0;j<sig_frame->size();j++) cout << j << " " << (*sig_frame)[j] << endl;
     cout << sig_rho << " " << sig_netcc << " " << sig_likelihood << endl; 

#ifdef PLOT_FRAME                       // frame display
     PlotFrame(sig_frame, i, save);
#endif

#ifdef WATPLOT                          // monster event display
     WTS->canvas->cd();
     WTS->plot(pwc, 1, nIFO, 'L', 0, const_cast<char*>("COLZ")); // use ID=1
     if(save) {
       char fname[1024];
       sprintf(fname, "monster_event_display_%d.png",i);
       cout << "write " << fname << endl;
       WTS->canvas->Print(fname);
     }
     WTS->clear();
#endif

   }
}


void 
PlotFrame(std::vector<double>* nn_frame, int cid, bool save) {

  // get dimension of the frame 
  int nframe = nn_frame->size();
  if(nframe != pow(int(sqrt(nframe)),2)) {
    cout << "PlotFrame - Error : size is not a square number " << endl;
    gSystem->Exit(1);                                                  
  }                                                                    
  nframe = sqrt(nframe);                                               

  if(hist2D) { delete hist2D; hist2D=NULL; }
  hist2D=new TH2F("frame", "frame", nframe, 0, nframe, nframe, 0, nframe);

  hist2D->SetStats(kFALSE);
  hist2D->SetTitleFont(12);
  hist2D->SetFillColor(kWhite);

  hist2D->GetXaxis()->SetNdivisions(506);
  hist2D->GetXaxis()->SetLabelFont(42);  
  hist2D->GetXaxis()->SetLabelOffset(0.014);
  hist2D->GetXaxis()->SetTitleOffset(1.4);  
  hist2D->GetYaxis()->SetTitleOffset(1.2);  
  hist2D->GetYaxis()->SetNdivisions(506);   
  hist2D->GetYaxis()->SetLabelFont(42);     
  hist2D->GetYaxis()->SetLabelOffset(0.01); 
  hist2D->GetZaxis()->SetLabelFont(42);     
  hist2D->GetZaxis()->SetNoExponent(false); 
  hist2D->GetZaxis()->SetNdivisions(506);   

  hist2D->GetXaxis()->SetTitleFont(42);
  hist2D->GetXaxis()->SetTitle("Time");
  hist2D->GetXaxis()->CenterTitle(true);     
  hist2D->GetYaxis()->SetTitleFont(42);      
  hist2D->GetYaxis()->SetTitle("Frequency"); 
  hist2D->GetYaxis()->CenterTitle(true);         

  hist2D->GetZaxis()->SetTitleOffset(0.6);
  hist2D->GetZaxis()->SetTitleFont(42);   
  hist2D->GetZaxis()->CenterTitle(true);  

  hist2D->GetXaxis()->SetLabelSize(0.03);
  hist2D->GetYaxis()->SetLabelSize(0.03);
  hist2D->GetZaxis()->SetLabelSize(0.03);

  for(int i=0;i<nframe;i++) {
    for(int j=0;j<nframe;j++) {
      //cout << i << " " << j << " " << (*nn_frame)[i*nframe+j] << endl;
      hist2D->SetBinContent(i+1,j+1,(*nn_frame)[i*nframe+j]);           
    }                                                                   
  }                                                                     

  if(canvas) { delete canvas; canvas=NULL; }
  canvas= new TCanvas("frame", "frame", 200, 20, 600, 600);
  canvas->Clear();                                     
  canvas->ToggleEventStatus();                         
  canvas->SetGridx();                                  
  canvas->SetGridy();                                  
  canvas->SetFillColor(kWhite);                        
  canvas->SetRightMargin(0.10);                        
  canvas->SetLeftMargin(0.10);                         
  canvas->SetBottomMargin(0.13);                       
  canvas->SetBorderMode(0);                            

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
  gStyle->SetCanvasColor(kWhite);
  gStyle->SetStatBorderSize(1);

  hist2D->Draw("COLZ");

  // change palette's width
  canvas->Update();
  TPaletteAxis *palette = (TPaletteAxis*)hist2D->GetListOfFunctions()->FindObject("palette");
  palette->SetX1NDC(0.91);
  palette->SetX2NDC(0.933);
  palette->SetTitleOffset(0.92);
  palette->GetAxis()->SetTickSize(0.01);
  canvas->Modified();

  if(save) {
    char fname[1024];
    sprintf(fname, "nn_frame_display_%d.png",cid);
    cout << "write " << fname << endl;
    canvas->Print(fname);
  }
}


