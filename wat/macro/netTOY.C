// 2G coherent WaveBurst production script
// S.Klimenko, University of Florida, June 15, 2012
{
   int  i,j,n,m,k;

  gSystem->Exec("date");

  if(search=='E' || search=='E') cout<<"\n un-modeled search (Energy): "<<search<<endl;
  if(search=='b' || search=='B') cout<<"\n un-modeled search (single stream): "<<search<<endl;
  if(search=='r' || search=='R') cout<<"\n un-modeled search (dual stream): "<<search<<endl;
  if(search=='i' || search=='I') cout<<"\n elliptical polarisation: "<<search<<endl;
  if(search=='g' || search=='G') cout<<"\n circular polarisation: "<<search<<endl;
  if(search=='s' || search=='S') cout<<"\n linear polarisation: "<<search<<endl;

  cout<<" network of ";
  for(i=0; i<nIFO; i++) cout<<ifo[i]<<" ";
  cout<<" detectors\n\n";

  // network initialization

  detector*           pD[nIFO];  // poiners to detectors
  WSeries<double>*   pTF[nIFO];  // pointer to WSeries
  wavearray<double>* hot[nIFO];  // temporary time series
  network            NET;        // network

  for(i=0; i<nIFO; i++)  pD[i] = new detector(ifo[i]);
  for(i=0; i<nIFO; i++) NET.add(pD[i]); 

  if(healpix) NET.setSkyMaps(4);
  else        NET.setSkyMaps(angle,Theta1,Theta2,Phi1,Phi2);

  NET.setAntenna(); 
  NET.constraint(delta,gamma);
  NET.setDelay(refIFO);
  NET.Edge = waveoffset;
  NET.netCC = netCC;
  NET.netRHO = netRHO;
  NET.EFEC = EFEC;
  NET.optim = optim;
  NET.precision = precision;

  TString wat_dir = gSystem->Getenv("HOME_ALUNO");  
  if(l_high==9)     NET.setMRAcatalog(wat_dir+"/OverlapCatalog5.bin");
  else if(l_low==3) NET.setMRAcatalog(wat_dir+"/OverlapCatalog8-1024.bin");
  else              NET.setMRAcatalog(wat_dir+"/OverlapCatalog16-1024.bin");
  NET.setAcore(Acore);         

  double mTau=NET.getDelay("MAX");  // maximum time delay 
  double dTau=NET.getDelay();       // time delay difference 
  cout<<"maximum time delay between detectors: "<<mTau<<endl;
  cout<<"       maximum time delay difference: "<<dTau<<endl;
  cout<<"                  skymap search mode: "<<mode<<endl;
  cout<<"           skymap angular resolution: "<<angle<<endl;
  cout<<"          skymap size in polar angle: "<<Theta1<<", "<<Theta2<<endl;
  cout<<"      skymap size in azimuthal angle: "<<Phi1<<", "<<Phi2<<endl;
  cout<<"                    netRHO and netCC: "<<netRHO<<", "<<netCC<<endl;
  cout<<"                       acor and e2or: "<<NET.acor<<", "<<NET.e2or<<endl;
  cout<<"              regulator delta, gamma: "<<delta<<" "<<gamma<<endl<<endl;

  injection mdc(nIFO);
  livetime live;
  netevent netburst(nIFO);
  variability wavevar;
  wavenoise noiserms;

  // time_frequency initialization

  int nRES = l_high-l_low+1;     // number of frequency resolution levels
  Meyer<double> B(1024);         // set wavelet for resampling
  WDM<double>* pwdm[nRES];       // wavelet pointers: wdm[0] - l_high, wdm[nRES-1] l_low 
  WDM<double> WD2(1024,1024,6,12);
  WDM<double> WD32(32,32,6,12);
  WDM<double> WD4(1024/2,1024/2,6,12);

  for(i=l_high; i>=l_low; i--) {
     pwdm[l_high-i] = new WDM<double>(1<<i,1<<i,6,10);
     NET.add(pwdm[l_high-i]);
  }

  gSystem->Exec("date");
  gSystem->Exec("hostname");
  gRandom->SetSeed(0);

  runID = 1;
  cout<<"job ID: "<<runID<<endl;
  //cout<<"Input : "<<input_dir<<"\n label: "<<input_label<<endl;
  //cout<<"Output: "<<output_dir<<"\n label: "<<output_label<<endl;
  NET.setRunID(runID);

  // read input segment and injection lists

  char file[512], tdf00[512], tdf90[512], buFFer[1024];
  int rnID = int(gRandom->Rndm(13)*1.e9);   // random name ID
  rnID = 4;

  double Tb=931159200.;                     // hardcoded test setup for bursts
  //double Tb=993107016.;                     // hardcoded test setup for NSNS
  double dT=616.;                           // hardcoded test setup
  double Te=Tb+dT;

// read and dump data on local disk (nodedir)

  wavearray<double> x;           // temporary time series
  wavearray<double> y;           // temporary time series

  x.start(Tb);
  x.stop(Te);
  x.edge(waveoffset);
  x.resize(1024*16*616);
  x.rate(1024*16);

  wavearray<double>* px;
  double gfactor=8.*sqrt(nIFO);           // SNR of simulated injections
  TString jname = wat_dir+"/test.root";   // read data from test root file
  //TString jname = wat_dir+"/NSNS_993107016.root";

  // open input job file
  TFile* jfile = new TFile(jname);
  if(jfile==NULL||!jfile->IsOpen())
    {cout << "Error : file " << jname << " not found" <<  endl;exit(1);}

  for(i=0; i<nIFO; i++) {

    sprintf(file,"strain/%s",ifo[i]);
    px = (wavearray<double>*)jfile->Get(file);     // read ifo strain from temporary job file
    x = *px; delete px;
    sprintf(file,"mdc/%s",ifo[i]);
    px = (wavearray<double>*)jfile->Get(file);     // read injection strain from temporary job file
    (*px)*=gfactor;
    x.add(*px);                                    // add injections
    delete px;

    sprintf(file,"%s/%s_%d_%s_%d_%d.dat",nodedir,ifo[i],int(Tb),output_label,runID,rnID);

    pTF[i] = pD[i]->getTFmap();
    hot[i] = pD[i]->getHoT();
    //pTF[i]->Forward(x,B,levelR);              
    //pTF[i]->getLayer(y,0);
    //y->DumpBinary(file);
    x->DumpBinary(file);

    fprintf(stdout,"start=%f duration=%f rate=%f\n",
	    x.start(),x.size()/x.rate(),x.rate());
  }

  size_t lags = 0;
  double factor = 1.0;          // strain factor
  double R = x.rate();          // data rate
  x.resize(0);

  char tmpFile[1024];
  char outFile[1024];
  char endFile[1024];
  char outDump[1024];
  char endDump[1024];
  char out_CED[1024];
  char end_CED[1024];
  char command[1024];
  FileStat_t fstemp;

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// if simulation==1, loop on the injection strain factors
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  if(!simulation) nfactor = 1;

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// data conditioning
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  lagSize=1;

  wavearray<double> htmp;
  WSeries<double> wtmp;
  WSeries<double> utmp;
  WSeries<double> ntmp;

  for(i=0; i<nIFO; i++) {
     sprintf(file,"%s/%s_%d_%s_%d_%d.dat",nodedir,ifo[i],int(Tb),output_label,runID,rnID);
     hot[i]->ReadBinary(file);
     hot[i]->rate(R);
     cout<<hot[i]->size()<<" "<<hot[i]->rate()<<endl;
     
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++
// REGRESSION stage should be included here
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++

     pTF[i]->Forward(*hot[i],WD4);
     regression rr(*pTF[i],"target",16,fHigh);
     rr.add(*hot[i],"target");   
     rr.setFilter(8);
     rr.setMatrix(NET.Edge,0.95);
     rr.solve(0.,10,'h');
     rr.apply(0.8);
     *hot[i] = rr.getClean();

     pTF[i]->Forward(*hot[i],WD2);
     pTF[i]->setlow(fLow);
     pD[i]->white(60.,0,NET.Edge,20.);
     pTF[i]->white(pD[i]->nRMS,1);
     pTF[i]->white(pD[i]->nRMS,-1);
     pTF[i]->sethigh(fHigh);
     pD[i]->bandPass(16.,fHigh);              // band pass filtering     
     wtmp = *pTF[i];
     pTF[i]->Inverse();
     wtmp.Inverse(-2);
     *hot[i] = *pTF[i];
     *hot[i] += wtmp;
     *hot[i] *= 0.5;
     
      cout<<"After "<<ifo[i]<<" data conditioning"<<endl; 
      gSystem->Exec("date"); // GetProcInfo();
    }

    if(!simulation) {                     // setup lags
      lags = NET.setTimeShifts(lagSize,lagStep,lagOff,lagMax,lagFile,lagMode,lagSite);
      cout<<"lag step: "<<lagStep<<endl;
      cout<<"number of time lags: "<<lags<<endl;
    }
    else if(!lags) lags = NET.setTimeShifts();

    double TL = NET.setVeto(gap);
    cout<<"live time in zero lag: "<<TL<<endl<<endl;  // set veto array 
    if(TL <= 0.) exit(1);                             // exit if live time is zero

    lagSize=1;

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// NETPIXEL stage: start of the coherent search
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    int npix;
    double Eo;
    bool append = false;              // file append for savemode
    wavearray<int> nP(NET.nLag);      // pixel counter
    wavearray<int> nC(NET.nLag);      // cluster counter
    TH1F* hic[nRES];
    for(j=0; j<nRES; j++) hic[j] = NULL;
    hic[3] = new TH1F(" c "," c ",1000,0.,50.);


    cout<<"Start coherent search: \n"; gSystem->Exec("date");

    for(j=1; j<nRES; j++) {           // loop over TF resolutions
	  
       for(i=0; i<nIFO; i++) {        // produce TF maps with max over the sky energy
          NET.getifo(i)->getTFmap()->maxEnergy(*hot[i],*pwdm[j],mTau,4);
       }

       gSystem->Exec("date");           // GetProcInfo();
       Eo = NET.THRESHOLD(bpp);         // threshold on pixel energy
       cout<<"thresholds in units of noise variance: Eo="<<Eo<<" Em="<<Eo*2<<endl;

       for(k=0; k<1; k++) {

          npix = NET.getNetworkPixels(k,Eo,Eo*2,hic[j]);
          cout<<"lag="<<k<<" pixels: "<<npix<<"  ";
	  n=1; m=1;	    
          cout<<"clusters: "<<NET.cluster(1,1)<<" ";

	  sprintf(file,"%s/pix_%d_%s_%d_%d_%d.lag",
		  nodedir,int(Tb),output_label,int(k),n,rnID);

	  if(!append) nP.data[k]  = NET.getwc(k)->write(file,0);
	  else        nP.data[k] += NET.getwc(k)->write(file,1);
          nC.data[k]+=NET.getwc(k)->csize();
          NET.getwc(k)->clear();
	  cout<<"stored: "<<nC.data[k]<<"|"<<nP.data[k]<<endl;
          
       }
       append = true;                
       gSystem->Exec("date"); // GetProcInfo();
    }	  	  

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// SUPERCLUSTER stage
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   gSystem->Exec("date"); // GetProcInfo();

   netcluster* pwc;
   netcluster wc;
   FILE* fp;
   size_t count; 

   int nnn = 0;
   int mmm = 0;

   TH2F* hist[2];
   hist[0] = new TH2F(" 0 "," 0 ",101,0.,1.01,100,0.,1);           // subNetCut performance histogram
   hist[1] = new TH2F(" 1 "," 1 ",400,-0.,1.1,400,-0.,1.1);           // subNetCut performance histogram
   TH1F* his = new TH1F(" 3 "," 3 ",1000,0.,3.);           // subNetCut performance histogram

// set time-delay filters

   lagSize=1;
   for(k=0; k<lagSize; k++){
      gSystem->Exec("date"); // GetProcInfo();
      sprintf(file,"%s/pix_%d_%s_%d_%d_%d.lag",
              nodedir,int(Tb),output_label,int(k),1,rnID);
      
      wc.read(file);
      cout<<file<<endl;
      cout<<"process lag "<<k<<" "<<NET.e2or<<endl;
      cout<<"loaded clusters|pixels: "<<wc.csize()<<"|"<<wc.size()<<endl;
      
      m = wc.supercluster('L', NET.e2or, 6.5, false);
      cout<<"super  clusters|pixels: "<<wc.esize(0)<<"|"<<wc.psize(0)<<endl;
      
      pwc = NET.getwc(k);
      pwc->cpf(wc, false);            

      gSystem->Exec("date"); // GetProcInfo();

      for(i=0; i<nIFO; i++) pD[i]->sclear();

      for(j=0; j<nRES; j++) { // loop over TF resolutions
         pwdm[j]->setTDFilter(12, 1);   // set low-rate TD filters
         for(i=0; i<nIFO; i++) {
            pTF[i]->Forward(*hot[i],*pwdm[j]);
            pD[i]->addSTFmap(pwc);
         }
      } 

      NET.setDelayIndex(hot[0]->rate()); 
      pwc->setcore(false);               // release all pixels
      

      n = 0;
      while(1){
         count = pwc->loadTDampSSE(NET, 'a', 10000, 300);            
         n += NET.subNetCut(k,subnet,0.33,hist[0]);         
         cout<<"selected pixels: "<<n<<", fraction: "<<n/double(pwc->psize(1)+pwc->psize(-1))<< endl;
         if(count<10000) break;
      }
      
      nnn += pwc->psize(-1);
      mmm += pwc->psize(1)+pwc->psize(-1);

      pwc->defragment(Tgap,Fgap);

      cout<<"events in the buffer: "<<NET.events()<<"|"<<nnn<<"|"<<nnn/double(mmm)<<"\n";

      gSystem->Exec("date"); // GetProcInfo();

      pwc->setcore(false);
      pwc->setcuts();

      for(j=0; j<nRES; j++) pwdm[j]->setTDFilter(12, 4);   // set low-rate TD filters
      NET.setDelayIndex(hot[0]->rate()*4); 
      NET.pOUT=1;
 
      gSystem->Exec("date"); // GetProcInfo();

       while(1){
         count = pwc->loadTDampSSE(NET, 'a', 10000,9999);
         NET.likelihood2G(search,k,930,hist[1]);            
         cout<<"processed pixels "<<count<<endl;
         if(count<10000) break;
      }
      gSystem->Exec("date"); // GetProcInfo();
      return;  
   }
}







