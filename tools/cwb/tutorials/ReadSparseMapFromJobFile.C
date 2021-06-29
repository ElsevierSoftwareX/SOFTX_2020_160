//
// Read & Plot Sparse Map From Job Files
// Author : Gabriele Vedovato

#define JFILE	"data/supercluster_931158223_75_ADV_SIM_EOBNRv2_L1H1V1_2G_OVERLAP_job1.root"

//#define SAVE_PLOT		// save plot to disk

watplot* WTS;

void 
ReadSparseMapFromJobFile(int ifoId=0, int resId=0, bool build_core=true, bool mask_core=false) {

  TString fName = JFILE;

  detector* pD[NIFO_MAX];                       //! pointers to detectors
  WSeries<double>* pTF[NIFO_MAX];               //! pointers to WSeries
  network NET;                                  //! network
  CWB::config cfg;				//! configuration

  TFile* jfile = new TFile(fName);
  if(jfile==NULL) {cout << "Error opening root file : " << fName.Data() << endl;exit(1);}
  // read network object
  if(jfile->Get("network")!=NULL) {
    // read network object
    NET = *(network*)jfile->Get("network");
  } else {
    cout << "Error : net is not contained in root file " << fName.Data() << endl;
    exit(1);
  }
  // read config object
  CWB::config* pcfg = (CWB::config*)jfile->Get("config");
  cfg = *pcfg;

  int nIFO=cfg.nIFO;			 // number of detectors 
  int nRES = cfg.l_high-cfg.l_low+1;     // number of frequency resolution levels

  if((ifoId<0)||(ifoId>=nIFO)) {
    cout << "ifoId non available - must be [0:" << nIFO-1 << "]" << endl;
    cout << "root 'ReadSparseMapFromJobFile(ifoId, resId, build_core=true/false, mask_core=false/true)'" << endl;
    exit(1);
  }
  if((resId<0)||(resId>=nRES)) {
    cout << "resId non available - must be [0:" << nRES-1 << "]" << endl;
    cout << "root 'ReadSparseMapFromJobFile(ifoId, resId, build_core=true/false, mask_core=false/true)'" << endl;
    exit(1);
  }

  for(int n=0; n<nIFO; n++) pD[n] = NET.getifo(n);	// get detectors
  //for(int n=0; n<nIFO; n++) pD[n]->print();

  int ifactor = 0;

  // read sparse map from job file
  cout << "Loading sparse TF map ... " << endl;
  for(int n=0; n<nIFO; n++) {
    pD[n]->sclear();   // clear vector with sparse maps
    for(int i=0; i<nRES; i++) {
      char swname[32];
      if(cfg.simulation) sprintf(swname,"sparse/%s-level:%d:%d",cfg.ifo[n],ifactor,i+cfg.l_low);
      else               sprintf(swname,"sparse/%s-level:%d",cfg.ifo[n],i+cfg.l_low);
      SSeries<double>* psw = (SSeries<double>*)jfile->Get(swname);
      if(psw==NULL) {
        cout << "sparse map " << swname << " not exist in job file" << endl;exit(1);
      }
      SSeries<double> SS = *psw;
      pD[n]->vSS.push_back(SS);
      delete psw;
    }
    cout<<endl;
  }

  // plot sparse maps
  WTS = new watplot(const_cast<char*>("WTS"));
  for(int n=0; n<nIFO; n++) {
    if(n!=ifoId) continue;
    for(int i=0; i<nRES; i++) {
      if(i!=resId) continue;
      // bild_core=true  -> rebuild TF map using only core pixels
      // bild_core=false -> rebuild TF map using all pixels
      pD[n]->vSS[i].Expand(build_core);	
      SSeries<double>* vSS = &pD[n]->vSS[i];
      cout << "Num Slices : " << vSS->GetSlices() << endl;
      cout << "Num Layers : " << vSS->GetLayers() << endl;
      cout << "Halo Slice : " << vSS->GetHaloSlice() << endl;
      cout << "Halo Layer : " << vSS->GetHaloLayer() << endl;

      // Plot WDM Scalogram
      double start = vSS->start();
      double stop  = vSS->start()+vSS->size()/vSS->rate();
      double flow  = cfg.fLow;
      double fhigh = cfg.fHigh;
      WTS->plot(vSS, 2, start, stop,const_cast<char*>("COLZ"));
      WTS->hist2D->GetYaxis()->SetRangeUser(flow, fhigh);

      char title[128];
      sprintf(title,"IFO = %s  -  Level = %d  -  dT = %g (sec)  -  dF = %g (Hz)",
              cfg.ifo[n],i+cfg.l_low,vSS->GetTimeResolution(),vSS->GetFreqResolution());
      WTS->hist2D->SetTitle(title);

      if(mask_core) {
        // set core pixels to white
        int xsize=WTS->hist2D->GetNbinsX();
        int ysize=WTS->hist2D->GetNbinsY();
        wavearray<int> index = vSS->GetSparseIndex();
        cout << "index size : " << index.size() << endl;
        for(int m=0;m<index.size();m++) {
          int k = vSS->GetSlice(index[m]);
          int j = vSS->GetLayer(index[m]);
          // write pixel 2 times because in watplot is duplicated
          // hist2D start from 1
          WTS->hist2D->SetBinContent(2*k+0,j+1,0);
          WTS->hist2D->SetBinContent(2*k+1,j+1,0);
        }
      }

#ifdef SAVE_PLOT
      // dump spectrum to disk
      char h2name[32];sprintf(h2name,"h%s-level:%d.png",cfg.ifo[n],cfg.l_high-i);
      WTS->canvas->Print(h2name);
      exit(0); 
#endif

      if(i==resId) return;
    }
  }

  jfile->Close();

}
