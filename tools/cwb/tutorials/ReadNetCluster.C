// this macro shows how to read the netcluster from output ROOT file produced with CWB_Plugin_WF.C
// netcluster must be enabled in the config/user_parameters.C with the following option:
//
//   options += "wf_output_enable=cluster "; // enable save event crate/cluster into the output root file
//
// Author : G.Vedovato

// ------------------------------------------------------------
//  structure
// ------------------------------------------------------------

struct XPARMS {                         // structure for output for estimated parameters

  std::vector<TString> ifo; 		// network ifo names

  float isnr;				// injected network SNR
  float osnr;				// reconstructed network SNR

  double seggps;                        // segment start time

  netcluster nc;      			// netcluster
};

// ------------------------------------------------------------
// functions
// ------------------------------------------------------------

int ReadDataFromROOT(TString ifile, XPARMS* xparms);
void PlotNetCluster(watplot* WTS, netcluster* pwc, double seggps, TString odir);
void DumpNetCluster(netcluster* pwc, int cid, int nifo, TString odir);


// ------------------------------------------------------------
// global variables
// ------------------------------------------------------------

int  gENTRY;

watplot* WTS;
int nIFO;

TCanvas* canvas;

void ReadNetCluster(TString ifroot, TString odir=".", int entry=0) {
//
// ifroot:	input root file name (event parameters + nectcluster)
// odir:	output directory where to dump the plots likelihood + null
// entry:	event entry in the input root file to be investigated 
//

  // init
  canvas=NULL;
  WTS=NULL;

  gROOT->SetBatch(true);	// macro is executed in batch mode

  gENTRY = entry;

  CWB::Toolbox::checkFile(ifroot);

  // read xparms from input root file
  XPARMS xparms;
  int err = ReadDataFromROOT(ifroot, &xparms);
  if(err>0) exit(1);

  PlotNetCluster(WTS, &xparms.nc, xparms.seggps, odir);

  DumpNetCluster(&xparms.nc, 1, nIFO, odir);

  exit(0);
}

int ReadDataFromROOT(TString ifile, XPARMS* xparms) { 
//
// ----------------------------------------------------
// Read netcluster
// ----------------------------------------------------

  TFile* froot = new TFile(ifile,"READ");
  if(froot==NULL) {cout << "ReadDataFromROOT - Error opening file " << ifile << endl;return 1;}
  TTree* itree = (TTree*)froot->Get("waveburst");
  if(itree==NULL) {cout << "ReadDataFromROOT - Error : no waveburst present in the file" << endl;return 1;}

  // get detector list
  TList* list = itree->GetUserInfo();
  nIFO=list->GetSize();
  if (nIFO==0) {cout << "ReadDataFromROOT - Error : no ifo present in the tree" << endl;return 1;}
  for (int n=0;n<list->GetSize();n++) {
    detector* pDetector;
    pDetector = (detector*)list->At(n);
    xparms->ifo.push_back(pDetector->Name);
    detectorParams dParams = pDetector->getDetectorParams();
    //pDetector->print();                                                       
  }

  itree->SetBranchAddress("ndim",&nIFO);
  itree->GetEntry(gENTRY);

  double seggps=0;			   // segment start time
  float crate;				   // netcluster::rate
  float iSNR[NIFO_MAX];
  float oSNR[NIFO_MAX];	
  std::vector<netpixel>* cluster;

  cluster = new std::vector<netpixel>;

  itree->SetBranchAddress("iSNR",iSNR);
  itree->SetBranchAddress("oSNR",oSNR);

  itree->SetBranchAddress("seggps",&seggps);
  itree->SetBranchAddress("crate",&crate);
  itree->SetBranchAddress("cluster",&cluster);

  itree->GetEntry(gENTRY);

  xparms->isnr = 0;
  for(int n=0;n<nIFO;n++) xparms->isnr += iSNR[n];
  xparms->isnr = sqrt(xparms->isnr);
  xparms->osnr = 0;
  for(int n=0;n<nIFO;n++) xparms->osnr += oSNR[n];
  xparms->osnr = sqrt(xparms->osnr);

  xparms->seggps=seggps;

  cout << "segment start time " << int(seggps) << endl;
  cout << "cluster rate       " << crate << endl;
  cout << "cluster size       " << cluster->size() << endl;
 
  if(cluster->size()==0) {cout << "ReadDataFromROOT - Error : no clusters in the root file !!!" << endl;return 1;}

  xparms->nc.rate=crate;
  std::vector<int> pList;
  for(int i=0;i<cluster->size();i++) {
    xparms->nc.pList.push_back((*cluster)[i]);
    pList.push_back(i);
  }
  xparms->nc.cList.push_back(pList);
  clusterdata cd;  // dummy cluster data
  xparms->nc.cData.push_back(cd);

  if(cluster) delete cluster;
  if(itree)   delete itree;
  if(froot)   delete froot;

  return 0;
}

void PlotNetCluster(watplot* WTS, netcluster* pwc, double seggps, TString odir) {
//
// ----------------------------------------------------
// Plot netcluster
// ----------------------------------------------------

  if(WTS==NULL) WTS = new watplot(const_cast<char*>("wtswrc"));

  WTS->canvas->cd();

  char fname[1024];

  TString xtitle = TString::Format("Time (sec) : GPS OFFSET = %.3f",seggps);

  // dump likelihood plot
  sprintf(fname, "%s/l_tfmap_scalogram.%s", odir.Data(), "png");
  WTS->plot(pwc, 1, nIFO, 'L', 0, const_cast<char*>("COLZ"),256,true);
  WTS->hist2D->GetXaxis()->SetTitle(xtitle);
  WTS->canvas->Print(fname);
  WTS->clear();

  // dump null plot
  sprintf(fname, "%s/n_tfmap_scalogram.%s", odir.Data(), "png");
  WTS->plot(pwc, 1, nIFO, 'N', 0, const_cast<char*>("COLZ"),256,true);
  WTS->hist2D->GetXaxis()->SetTitle(xtitle);
  WTS->canvas->Print(fname);
  WTS->clear();

  return;
}

void DumpNetCluster(netcluster* pwc, int cid, int nifo, TString odir) {
//
// pwc    : pointer to netcluster object
// cid    : cluster id
// nifo   : number of detectors
//

  char ofname[1024];
  sprintf(ofname, "%s/tfmap_dump.%s", odir.Data(), "txt");

  FILE *fp;
  if((fp = fopen(ofname, "w")) == NULL ) {
     cout << "DumpNetCluster error: cannot open file " << ofname <<". \n";
     return;
  };
  cout << endl << "ofile name = " << ofname << endl << endl;

  double RATE = pwc->rate;                      	// original rate

  std::vector<int>* vint = &(pwc->cList[cid-1]); 	// pixel list

  int V = vint->size();                         	// cluster size
  if(!V) return;                                               

  netpixel* pix = pwc->getPixel(cid,0);

  for(int n=0; n<V; n++) {
    netpixel* pix = pwc->getPixel(cid,n);
    if(!pix->core) continue;            

    double like = pix->likelihood>0. ? pix->likelihood : 0.;
    double null = pix->null>0. ? pix->null : 0.;
  
    double time = int(pix->time/pix->layers)/double(pix->rate); 	// central pixel time
    double dt = 1./pix->rate;
    time -= dt/2.;							// begin pixel time

    double freq = pix->frequency*pix->rate/2.;
    double df = 1./(2.*dt);

    char sout[1024];
    sprintf(sout,"pixel %*d\ttime = %*.3f (sec) \tfreq = %*.3f (Hz) \tnull = %*.3f \tlikelihood = %*.3f \tRes = %*.3f x %*.3f (sec,Hz)\n",7,n,7,time,7,freq,7,null,7,like,7,dt,7,df);
    cout << sout;
    fprintf(fp,"%s\n",sout);
  }                                                                   

  fclose(fp);


  return;
}

