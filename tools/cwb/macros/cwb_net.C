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


// Coherent WaveBurst production script for GW network
{
  CWB::Toolbox TB;

  TB.checkFile(gSystem->Getenv("CWB_ROOTLOGON_FILE"));
  TB.checkFile(gSystem->Getenv("CWB_PARAMETERS_FILE"));
  TB.checkFile(gSystem->Getenv("CWB_UPARAMETERS_FILE"));


  int  i,j,n,m;

  if(search=='E' || search=='E') cout<<"\n un-modeled search (Energy): "<<SEARCH()<<endl;
  if(search=='b' || search=='B') cout<<"\n un-modeled search (single stream): "<<SEARCH()<<endl;
  if(search=='r' || search=='R') cout<<"\n un-modeled search (dual stream): "<<SEARCH()<<endl;
  if(search=='i' || search=='I') cout<<"\n elliptical polarisation: "<<SEARCH()<<endl;
  if(search=='g' || search=='G') cout<<"\n circular polarisation: "<<SEARCH()<<endl;
  if(search=='s' || search=='S') cout<<"\n linear polarisation: "<<SEARCH()<<endl;

  // Check single detector mode 
  // if nIFO=1 the analysis is done as a network of 2 equal detectors
  bool singleDetector=false;
  if(nIFO==1) {
    cout << "------> cWB in Sigle Detector Mode !!!" << endl;
    singleDetector=true;
    sprintf(ifo[1],"%s",ifo[0]);
    strcpy(refIFO,ifo[0]);
    strcpy(channelNamesRaw[1],channelNamesRaw[0]);
    strcpy(channelNamesMDC[1],channelNamesMDC[0]);
    strcpy(frFiles[2],frFiles[1]);  // mdc frFiles
    strcpy(frFiles[1],frFiles[0]);

    for(int i=0;i<nDQF;i++) DQF[i+nDQF]=DQF[i];
    nDQF*=2;

    eDisbalance = false;   // disable energy disbalance
    lagSize = 1;
    lagOff  = 0;
    mode    = 1;           // 1 - exclude duplicate delay configurations
    delta   = 0;           // weak regulator
    GAMMA() = 0;           // force regulator not to be hard

    nIFO++;
  }


  cout<<" network of ";
  for(i=0; i<nIFO; i++) cout<<ifo[i]<<" ";
  cout<<" detectors\n\n";

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// declarations
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  Meyer<double> B(1024);       // set wavelet for resampling
  Meyer<double> S(1024,2);     // set wavelet for production

  netcluster wc;           
  netcluster* pwc;

  wavearray<double> x;         // temporary time series
  WSeries<double> wM;          // mdc WSeries
  WSeries<float>   v[nIFO];    // noise variability
  detector*       pD[nIFO];    // pointers to detectors
  WSeries<double>* pTF[nIFO];  // pointers to WSeries

  for(i=0; i<nIFO; i++) pD[i] = new detector(ifo[i]);

  network         NET;         // network

  for(i=0; i<nIFO; i++) NET.add(pD[i]); 
  NET.setSkyMaps(angle,Theta1,Theta2,Phi1,Phi2);
  NET.setAntenna(); 
  NET.constraint(delta,GAMMA());
  NET.setDelay(refIFO);
  NET.Edge = segEdge;
  NET.netCC = netCC;
  NET.netRHO = netRHO;
  NET.EFEC = EFEC;
  NET.precision = precision;
  NET.nSky = nSky;  
  NET.eDisbalance = eDisbalance;

  double mTau=NET.getDelay("MAX");  // maximum time delay 
  double dTau=NET.getDelay();       // time delay difference 
  cout<<"maximum time delay between detectors: "<<mTau<<endl;
  cout<<"       maximum time delay difference: "<<dTau<<endl;
  cout<<"                  skymap search mode: "<<mode<<endl;
  cout<<"           skymap angular resolution: "<<angle<<endl;
  cout<<"          skymap size in polar angle: "<<Theta1<<", "<<Theta2<<endl;
  cout<<"      skymap size in azimuthal angle: "<<Phi1<<", "<<Phi2<<endl;
  cout<<"                    netRHO and netCC: "<<netRHO<<", "<<netCC<<endl;
  cout<<"              regulator delta, local: "<<delta<<" "<<NET.local<<endl<<endl;

  injection mdc(nIFO);
  livetime live;
  netevent netburst(nIFO,Psave);
  variability wavevar;
  wavenoise noiserms;

  gSystem->Exec("/bin/date");
  gSystem->Exec("/bin/hostname");
  gRandom->SetSeed(0);

  // parse input

  int dump_infos_and_exit = TString(gSystem->Getenv("CWB_DUMP_INFOS_AND_EXIT")).Atoi();
  int dump_sensitivity_and_exit = TString(gSystem->Getenv("CWB_DUMP_SENSITIVITY_AND_EXIT")).Atoi();

  TString srunID = TString(gSystem->Getenv("CWB_JOBID"));
  runID = srunID.Atoi();

  cout<<"job ID: "<<runID<<endl;
  cout<<"Output: "<<output_dir<<"\n  label: "<<data_label<<endl;

  NET.setRunID(runID);

  // ---------------------------------------------------------------------
  // CWB HISTORY INIT
  // ---------------------------------------------------------------------

  CWB::Toolbox histTB;
  char job_stage[256];
  if(simulation==1) sprintf(job_stage,"SIMULATION");
  else              sprintf(job_stage,"PRODUCTION");

  char* cwbBuffer = histTB.getEnvCWB();
  if(cwbBuffer!=NULL) {
    history.AddHistory(job_stage, "CWB_ENV", cwbBuffer);
    delete [] cwbBuffer;
  }

  history.AddHistory(job_stage, "WATVERSION", watversion('s'));
  history.AddHistory(job_stage, "WORKDIR", work_dir);
  history.AddHistory(job_stage, "DATALABEL", data_label);

  char cmd_line[512]="";
  for(int i=0;i<gApplication->Argc();i++) sprintf(cmd_line,"%s %s",cmd_line,gApplication->Argv(i));
  history.AddHistory(job_stage, "CMDLINE", cmd_line);

  char* rootlogonBuffer = histTB.readFile("rootlogon.C");
  if(rootlogonBuffer!=NULL) {
    history.AddHistory(job_stage, "ROOTLOGON", rootlogonBuffer);
    delete [] rootlogonBuffer;
  }

  // save configuration files
  TString cwb_parameters_name = TString(gSystem->Getenv("CWB_PARAMETERS_FILE"));
  for(int i=0;i<gApplication->Argc()-1;i++) {  // skip last argument (net.C) 
    if(TString(gApplication->Argv(i)).Contains(".C")) {
      char* parametersBuffer = histTB.readFile(gApplication->Argv(i));
      if(parametersBuffer!=NULL) {
        if(TString(gApplication->Argv(i))==cwb_parameters_name) {
          history.AddHistory(job_stage, "PARAMETERS", parametersBuffer);
        } else {
          history.AddLog(job_stage, parametersBuffer);
        }
      }
      delete [] parametersBuffer;
    }
  }

  history.AddLog(job_stage, "START JOB");
  int job_start_time=CWB::Time("now").GetGPS();

  // ---------------------------------------------------------------------
  // JOB INIT
  // ---------------------------------------------------------------------

  #include <vector>

  // cwb toolbox
  CWB::Toolbox slagTB;
  CWB::Toolbox frTB[nIFO+1];   // the first nIFO positions contains the ifo frFiles, nIFO+1 contains the sim frFile
  int nfrFiles[nIFO+1];   
  vector<waveSegment> detSegs;
  vector<waveSegment> cat1List,cat2List;
  int jobID;int slagID;int segID[20];
  double mdcShift=0.;
  vector<TString> ifos(nIFO);
  for(int n=0;n<nIFO;n++) ifos[n]=ifo[n];

  if((simulation)&&((slagSize>1)||(lagSize!=1))) {
    cout << "Error : slagSize<=1 & lagSize==1 in simulation mode !!!" << endl;exit(1);
  }

  if(slagSize>0) {   // SLAG Segments

    cout << endl << "START SLAG Init ..." << endl << endl;

    // get zero lag merged dq cat 1 list
    // used only to get slagSegs && slagJobList
    cat1List=slagTB.readSegList(nDQF, DQF, CWB_CAT1);

    // get number/list of possible super lag jobs 
    vector<waveSegment> slagJobList=slagTB.getSlagJobList(cat1List, segLen);  
    int slagSegs=slagJobList.size();  

    // get super lag list
    vector<slag> slagList=slagTB.getSlagList(nIFO, slagSize, slagSegs, slagOff, slagMin, slagMax, slagSite);

    // init slag structures
    slag SLAG=slagTB.getSlag(slagList,runID);
    if(SLAG.jobId!=runID) {cout << "jobID " << runID << " not found in the slag list !!!" << endl;exit(1);}
    slagID = SLAG.slagId[0];
    jobID  = SLAG.jobId;
    for(n=0; n<nIFO; n++) segID[n]=SLAG.segId[n];
    cout << "SuperLag=" << slagID << " jobID=" << jobID;
    for(n=0; n<nIFO; n++) cout << " segID[" << ifo[n] << "]=" << segID[n];cout << endl;

    // set slag shifts into the DQF structures
    for(int n=0;n<nIFO;n++) ifos[n]=ifo[n];
    slagTB.setSlagShifts(SLAG, ifos, segLen, nDQF, DQF);  

    // get shifted merged dq cat 1 list
    cat1List=slagTB.readSegList(nDQF, DQF, CWB_CAT1);

    // extract detector's slag segments range from dq cat 1
    detSegs=slagTB.getSegList(SLAG, slagJobList, segLen, segMLS, segEdge, cat1List);
    cout << "Slag Segments" << endl;

    cout << endl << "END SLAG Init ..." << endl << endl;

  } else {     // Standard Segments      

    jobID  = runID;
    slagID = 0;
    for(int n=0;n<nIFO;n++) segID[n]=0;

    // get shifted merged dq cat 1 list
    cat1List=slagTB.readSegList(nDQF, DQF, CWB_CAT1);

    // extract detector's segments range from dq cat 1
    detSegs=slagTB.getSegList(jobID, nIFO, segLen, segMLS, segEdge, cat1List);
    cout << "Standard Segments" << endl;
  }

  // store slags in NET class
  float slagShift[20];
  for(n=0; n<nIFO; n++) slagShift[n]=(segID[n]-segID[0])*segLen;
  slagShift[nIFO]=slagID;
  netburst.setSLags(slagShift);

  // print detector's segments for this job
  cout.precision(14);
  for(int i=0;i<nIFO;i++) { 
    cout << "detSegs_dq1[" << ifo[i] << "] Range : " << detSegs[i].start << "-" << detSegs[i].stop << endl;
  }
  if(detSegs.size()==0) {cout << "no segments found for this job, job terminated !!!" << endl;exit(1);}

  // set & get frame file list
  frfile FRF[nIFO+1];
  for(n=0; n<nIFO; n++) {
    nfrFiles[n]=frTB[n].frl2FrTree(frFiles[n]);
    cout << ifo[n] << " -> nfrFiles : " << nfrFiles[n] << endl;
    FRF[n] = frTB[n].getFrList(detSegs[n].start-dataShift[n], detSegs[n].stop-dataShift[n], segEdge); // dataShift
  }
  if(simulation) {
    // set frame file list
    nfrFiles[nIFO]=frTB[nIFO].frl2FrTree(frFiles[nIFO]);
    cout << "MDC " << " -> nfrFiles : " << nfrFiles[nIFO] << endl;
    if(mdc_shift.startMDC<0) {  
      // read mdc range from mdc frl files
      waveSegment mdc_range  = frTB[nIFO].getFrRange();
      cout << "mdc_range : " << mdc_range.start << " " << mdc_range.stop << endl;
      mdc_shift.startMDC=mdc_range.start;
      mdc_shift.stopMDC=mdc_range.stop;
    }
    mdcShift  = frTB[nIFO].getMDCShift(mdc_shift, detSegs[0].start);
    cout << "mdcShift : " << mdcShift << endl;
    FRF[nIFO] = frTB[nIFO].getFrList(detSegs[0].start-mdcShift, detSegs[0].stop-mdcShift, segEdge);
  }

  // get shifted merged dq cat 2 list
  cat2List=slagTB.readSegList(nDQF, DQF, CWB_CAT2);
  // store cat2List into network class
  NET.segList=cat2List;		

  // check if seg+cat2 data length is not zero
  vector<waveSegment> detSegs_dq2;
  detSegs_dq2.push_back(detSegs[0]);
  detSegs_dq2 = slagTB.mergeSegLists(detSegs_dq2,cat2List);
  for(int i=0;i<detSegs_dq2.size();i++) {
    cout << "detSegs_dq2[" << i << "] Range : " << detSegs_dq2[i].start << "-" << detSegs_dq2[i].stop << endl;
  }
  double detSegs_ctime = slagTB.getTimeSegList(detSegs_dq2);
  cout << "live time after cat 2 : " << detSegs_ctime << endl;
  if(detSegs_ctime<segTHR) {cout << "job segment live time after cat2 < " << segTHR << " sec, job terminated !!!" << endl;exit(1);}

  double Tb=detSegs[0].start;
  double Te=detSegs[0].stop;
  double dT = Te-Tb;                        // WB segment duration

  // read input framelist file

  char file[512], tdf00[512], tdf90[512], buFFer[1024];
  int rnID = int(gRandom->Rndm(13)*1.e9);   // random name ID

  if(simulation) {                          // reag MDC log file
    i=NET.readMDClog(injectionList,double(long(Tb))-mdcShift);  
    printf("GPS: %16.6f saved,  injections: %d\n",double(long(Tb)),i);
    frTB[nIFO].shiftBurstMDCLog(NET.mdcList, ifos, mdcShift);  
    for(int i=0;i<NET.mdcTime.size();i++) NET.mdcTime[i]+=mdcShift; 

    // check if seg+cat2+inj data length is not zero
    vector<waveSegment> mdcSegs(NET.mdcTime.size());
    for(int k=0;k<NET.mdcTime.size();k++) {mdcSegs[k].start=NET.mdcTime[k]-gap;mdcSegs[k].stop=NET.mdcTime[k]+gap;}
    vector<waveSegment> mdcSegs_dq2 = slagTB.mergeSegLists(detSegs_dq2,mdcSegs);
    double mdcSegs_ctime = slagTB.getTimeSegList(mdcSegs_dq2);
    cout << "live time in zero lag after cat2+inj : " << mdcSegs_ctime << endl;
    if(mdcSegs_ctime==0) {cout << "job segment with zero cat2+inj live time in zero lag, job terminated !!!" << endl;exit(1);}
  } 

  if(dump_infos_and_exit) exit(0);

  if(mask>0.) NET.setSkyMask(mask,skyMaskFile);

  for(i=0; i<nIFO; i++) {
    frTB[i].readFrames(FRF[i],channelNamesRaw[i],x);  
    x.start(x.start()+dataShift[i]);                // dataShift
    x.start(x.start()-segLen*(segID[i]-segID[0]));  // SLAG
    if(singleDetector) TB.resampleToPowerOfTwo(x);
    sprintf(file,"%s/%s_%d_%s_%d_%d.dat",
	    nodedir,ifo[i],int(Tb),data_label,runID,rnID);
    if(dump_sensitivity_and_exit) {     // used only to save sensitivity
      sprintf(file,"%s/sensitivity_%s_%d_%s_job%d.txt",dump_dir,ifo[i],int(Tb),data_label,runID);
      cout << endl << "Dump Sensitivity : " << file << endl << endl;
      TB.makeSpectrum(file, x);
      continue;
    }
    if(dcCal[i]>0.) x*=dcCal[i];               // DC correction
    if(fResample>0) {
      x.FFT(1); x.resize(fResample/x.rate()*x.size()); x.FFT(-1); x.rate(fResample); //RESAMPLING
    }
    pTF[i] = pD[i]->getTFmap();
    pTF[i]->Forward(x,B,levelR);
    pTF[i]->getLayer(x,0);
    pTF[i]->Forward(x,S,levelD);
    pTF[i]->DumpBinary(file);
    n = pTF[i]->size();

    fprintf(stdout,"start=%f duration=%f rate=%f\n",
	    x.start(),x.size()/x.rate(),x.rate());
    if(i>0 && pTF[0]->start() != x.start()) {
      cout << "net.C - Error : data not synchronized" << endl;
      cout << ifo[i] << " " << x.start() << " != " << ifo[0] << " " << pTF[0]->start() << endl;
      exit(1);
    }
    if(i>0 && pTF[0]->rate()  != x.rate()) {
      cout << "net.C - Error : data have different rates" << endl;
      cout << ifo[i] << " " << x.rate() << " != " << ifo[0] << " " << pTF[0]->rate() << endl;
      exit(1);
    }
 
    if(simulation) {
      frTB[nIFO].readFrames(FRF[nIFO],channelNamesMDC[i],x);  // SLAG
      x.start(x.start()+mdcShift);  // SLAG
      sprintf(file,"%s/mdc%s_%d_%s_%d_%d.dat",
	      nodedir,ifo[i],int(Tb),data_label,runID,rnID);
      cout<<file<<endl;
      if(fResample>0) {
        x.FFT(1); x.resize(fResample/x.rate()*x.size()); x.FFT(-1); x.rate(fResample); //RESAMPLING
      }
      wM.Forward(x,B,levelR); 
      wM.getLayer(x,0); 
      wM.Forward(x,S,levelD);
      wM.DumpBinary(file);
      m = wM.size();

      fprintf(stdout,"start=%f duration=%f rate=%f\n",
	      x.start(),x.size()/x.rate(),x.rate());
      if(i>0 && wM.start() != x.start()) exit(1);
      if(i>0 && wM.rate() != x.rate())   exit(1);
      if(n != m)   exit(1);
    }
    if(singleDetector) break;
  }
  if(dump_sensitivity_and_exit) exit(0);

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

  for(int iii=0; iii<nfactor; iii++) {
    
    if(simulation) {
      factor=factors[iii];
      cout<<"factor="<<factor<<endl;
      char sim_label[512];
      sprintf(sim_label,"%d_%d_%s_%g_job%d",int(Tb),int(dT),data_label,factor,runID);

      sprintf(outFile,"%s/wave_%s.root",nodedir,sim_label);
      sprintf(endFile,"%s/wave_%s.root",output_dir,sim_label);
      sprintf(tmpFile,"%s/wave_%s.root.tmp",nodedir,sim_label);
      sprintf(outDump,"%s/wave_%s.txt",nodedir,sim_label);
      sprintf(endDump,"%s/wave_%s.txt",output_dir,sim_label);
      sprintf(out_CED,"%s/ced_%s",nodedir,sim_label);
      sprintf(end_CED,"%s/ced_%s",output_dir,sim_label);

      if(!gSystem->GetPathInfo(endFile,fstemp)) {
	printf("The file %s already exists - skip\n",endFile);
	fflush(stdout);
	TFile rf(endFile); 
	if(!rf.IsZombie()) continue;
      }

    }
    else {
      char prod_label[512];
      sprintf(prod_label,"%d_%d_%s_slag%d_lag%d_%d_job%d",
              int(Tb),int(dT),data_label,slagID,lagOff,lagSize,runID);

      sprintf(outFile,"%s/wave_%s.root",nodedir,prod_label);
      sprintf(endFile,"%s/wave_%s.root",output_dir,prod_label);
      sprintf(tmpFile,"%s/wave_%s.root.tmp",nodedir,prod_label);
      sprintf(outDump,"%s/wave_%s.txt",nodedir,prod_label);
      sprintf(endDump,"%s/wave_%s.txt",output_dir,prod_label);
      sprintf(out_CED,"%s/ced_%s",nodedir,prod_label);
      sprintf(end_CED,"%s/ced_%s",output_dir,prod_label);
    }

    cout<<"output file on the node: "<<outFile<<endl;
    cout<<"final output file name : "<<endFile<<endl;
    cout<<"temporary output file  : "<<tmpFile<<endl;

    gSystem->Exec("/bin/date"); GetProcInfo();

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// data conditioning
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    for(i=0; i<nIFO; i++) {
      sprintf(file,"%s/%s_%d_%s_%d_%d.dat",
	      nodedir,ifo[i],int(Tb),data_label,runID,rnID);
      pTF[i]->setLevel(levelD);
      pTF[i]->ReadBinary(file);
      
      if(simulation) {
	sprintf(file,"%s/mdc%s_%d_%s_%d_%d.dat",
		nodedir,ifo[i],int(Tb),data_label,runID,rnID);
	wM.setLevel(levelD);
	wM.ReadBinary(file); wM*=factor;
	pTF[i]->add(wM);
	wM*=1./factor;
      }

      pTF[i]->lprFilter(2,0,Tlpr,4.);
      pTF[i]->setlow(fLow);
      pD[i]->white(60.,1,8.,20.);
      if(simulation) pD[i]->setsim(wM,NET.getmdcTime(),10.,8.,true);   
      pTF[i]->Inverse(levelD-levelF);
      pTF[i]->lprFilter(2,0,Tlpr,4.);
      pTF[i]->Forward(levelD-levelF);
      pTF[i]->sethigh(fHigh);
      v[i] = pTF[i]->variability();
      pD[i]->bandPass();                // band pass filtering
            
      cout<<"After "<<ifo[i]<<" data conditioning"<<endl; 
      gSystem->Exec("/bin/date"); GetProcInfo();

      if(singleDetector) {*pD[1]=*pD[0]; break;}
    }
/* MLAG
    if(!simulation) {                     // setup lags
      lags = NET.setTimeShifts(lagSize,lagStep,lagOff,lagMax,lagFile,lagMode,lagSite);
      cout<<"lag step: "<<lagStep<<endl;
      cout<<"number of time lags: "<<lags<<endl;
    }
    else if(!lags) lags = NET.setTimeShifts();

    double TL = NET.setVeto(gap);
    cout<<"live time in zero lag: "<<TL<<endl<<endl;  // set veto array 
    if(TL <= 0.) exit(1);                             // exit if live time is zero
*/
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// initialization of the output files 
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    TFile *froot = new TFile(tmpFile, "RECREATE");
    TTree* net_tree = netburst.setTree();
    TTree* live_tree= live.setTree();

    if(simulation) {
      TTree* mdc_tree = mdc.setTree();
    }
    else {
      TTree* var_tree = wavevar.setTree();
      TTree* noise_tree = noiserms.setTree();
    }

    //if(dump) netburst.dopen(outDump,"w");  

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// start of the coherent search
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    Long_t xid,xsize,xflags,xmt;
    int xestat = gSystem->GetPathInfo(outDump,&xid,&xsize,&xflags,&xmt);
    if (xestat==0) {
      sprintf(command,"/bin/rm %s",outDump);
      gSystem->Exec(command);
    }

    int mlagSize=lagOff+lagSize;
    int mlagOff=lagOff;

    if(mlagStep==0) mlagStep=lagSize; // if mlagStep=0 -> standard lag analysis

    for(int mlag=mlagOff;mlag<mlagSize;mlag+=mlagStep) {

      lagOff  = mlag;
      lagSize = lagOff+mlagStep<=mlagSize ? mlagStep : mlagSize-lagOff;
      if(lagSize==0) continue;

      cout << "lagSize : " << lagSize << " lagOff : " << lagOff << endl;


      if(!simulation) {                     // setup lags
        lags = NET.setTimeShifts(lagSize,lagStep,lagOff,lagMax,lagFile,lagMode,lagSite);
        cout<<"lag step: "<<lagStep<<endl;
        cout<<"number of time lags: "<<lags<<endl;
      }
      else if(!lags) lags = NET.setTimeShifts();
  
      double TL = NET.setVeto(gap);
      cout<<"live time in zero lag: "<<TL<<endl<<endl;  // set veto array 
      if(TL <= 0.) exit(1);                             // exit if live time is zero
  
  
      double Ao;
      wavearray<int> np(NET.nLag);      // pixel counter
      bool append = false;              // file append for savemode
      int ceddir = 0;                   // flag if ced directory exists
  
      cout<<"Start coherent search: "; gSystem->Exec("/bin/date");
  
      for(i=levelD; i>=l_low; i--) {  // loop over TF resolutions
  	  
        if(i<=l_high) {
  	    	    	    
       	  sprintf(tdf00,"%s/Meyer1024wat482_00_L%1d.dat",filter_dir,i);
  	  NET.setDelayFilters(tdf00);
  	  if(i==l_high) {
  	    NET.setDelayIndex();
  	    NET.setIndexMode(1);
    	  }
  
  	  Ao = NET.threshold(bpp,dTau);
  	  NET.set2or(x2or*Ao*Ao);
  	  cout<<"pixel threshold in units of noise rms: "<<Ao<<endl;
  	  cout<<"2 OR  threshold in units of noise var: "<<x2or*Ao*Ao<<endl;
  	
  	  cout<<"total    pixels: "<<NET.coherence(Ao)<<"  ";
  	    
  	  n = size_t(2.*Tgap*pD[0]->getTFmap()->resolution(0)+0.1);
  	  m = size_t(Fgap/pD[0]->getTFmap()->resolution(0)+0.0001);
  	    
  	  cout<<"clusters: "<<NET.cluster(n,m)<<"  ";
  	  cout<<"selected pixels: "<<NET.likelihood('E',Acore)<<"\n";
  
  	  for(j=0; j<NET.nLag; j++) {
  	    sprintf(file,"%s/pix_%d_%s_%d_%d_%d.lag",
  		    nodedir,int(Tb),data_label,int(j),iii,rnID);
  	    wc = *(NET.getwc(j));
  	    if(!append) np.data[j]  = wc.write(file,0);
  	    else        np.data[j] += wc.write(file,1);
  	    cout<<wc.csize()<<"|"<<wc.size()<<"|"<<np.data[j]<<" "; 
  	  }
  	  cout<<endl;
  	  append = true;                
        }
	  
        if(i>l_low) NET.Inverse(1); 
        gSystem->Exec("/bin/date"); GetProcInfo();
	  
      }

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// supercluster analysis
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      for(j=0; j<lags; j++){
        sprintf(file,"%s/pix_%d_%s_%d_%d_%d.lag",
  	        nodedir,int(Tb),data_label,int(j),iii,rnID);
        wc.read(file);
        m = wc.supercluster('L',NET.e2or,true);
        pwc = NET.getwc(j); pwc->cpf(wc,true);
        cout<<m<<"|"<<pwc->size()<<" ";
        sprintf(command,"/bin/rm %s",file);
        gSystem->Exec(command);
      }
      cout<<endl;
      cout<<"events in the buffer: "<<NET.events()<<"\n";
      gSystem->Exec("/bin/date"); GetProcInfo();


//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// likelihood
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      for(i=l_low; i<=l_high; i++) {
        sprintf(tdf00,"%s/Meyer1024wat482_00%s_L%1d.dat",filter_dir,filter,i);
        sprintf(tdf90,"%s/Meyer1024wat482_90%s_L%1d.dat",filter_dir,filter,i);
        NET.setDelayFilters(tdf00,tdf90);

        if(i==l_low) {
	  NET.setDelayIndex();
	  NET.setIndexMode(mode);
        }

        cout<<"selected core pixels: "<<NET.likelihood(SEARCH(),Acore)<<" for level "<<i<<"\n";
        cout<<"rejected weak pixels: "<<NET.netcut(netRHO,'r',0,1)<<"\n";  // remove weak glitches
        cout<<"rejected loud pixels: "<<NET.netcut(netCC,'c',0,1)<<"\n";   // remove loud glitches
        cout<<"events in the buffer: "<<NET.events()<<"\n";

        if(cedDump) {
	  cout<<"dump ced into "<<out_CED<<"\n";
	  //if(CED(&NET,out_CED,cedRHO)) ceddir = 1; 
          //if(netburst.ced(&NET,out_CED,cedRHO,factor)) ceddir = 1; 
          CWB::ced ced(&NET,&netburst,out_CED);
          ced.SetOptions(cedRHO);
          if(singleDetector) ced.SetChannelName(channelNamesRaw[0]);
          bool fullCED = singleDetector ? false : true; 
          if(ced.Write(factor,fullCED)) ceddir = 1;
        }

        if(i<l_high) NET.Forward(1);
      }

      gSystem->Exec("/bin/date");
      gSystem->Exec("/bin/date"); GetProcInfo();
      cout<<"\nSearch done\n";
      cout<<"reconstructed events: "<<NET.events()<<"\n";

      if(simulation) NET.printwc(0);

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// save data in root file
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      if(dump) netburst.dopen(outDump,"w");

      //live.output(live_tree,&NET);
      live.output(live_tree,&NET,slagShift);   // FIX LIVETIME

      if(simulation) {
        netburst.output(net_tree,&NET,factor);
        mdc.output(mdc_tree,&NET,factor);
      }
      else {
        netburst.output(net_tree,&NET);
        for(i=0; i<nIFO; i++) {
//	  wavevar.output(var_tree,&v[i],i+1,segEdge);
	  noiserms.output(noise_tree,&pD[i].nRMS,i+1,R/2);
        }
      }

      history.AddLog(job_stage, "STOP JOB"); 
      history.Write();

      froot->Write();
      if(dump) netburst.dclose();

    } // end mlag loop

    froot->Close();

    sprintf(command,"/bin/mv %s %s", tmpFile, outFile);
    if(!cedDump) gSystem->Exec(command);
    sprintf(command,"/bin/mv %s %s",outFile,endFile);
    if(!cedDump) gSystem->Exec(command);
    sprintf(command,"/bin/mv %s %s",outDump,endDump);
    if(!cedDump && dump) gSystem->Exec(command);
    xestat = gSystem->GetPathInfo(end_CED,&xid,&xsize,&xflags,&xmt);
    if (xestat==0) {
      sprintf(command,"/bin/mv %s/* %s/.",out_CED,end_CED);
    } else {
      sprintf(command,"/bin/mv %s %s",out_CED,end_CED);
    }
    if(cedDump && ceddir) gSystem->Exec(command);

  }

// clean-up temporary data files

  for(i=0; i<nIFO; i++) {
    sprintf(file,"%s/%s_%d_%s_%d_%d.dat",nodedir,ifo[i],int(Tb),data_label,runID,rnID);
    sprintf(command,"/bin/rm %s",file);
    gSystem->Exec(command);

    if(simulation) {
      sprintf(file,"%s/mdc%s_%d_%s_%d_%d.dat",nodedir,ifo[i],int(Tb),data_label,runID,rnID);
      sprintf(command,"/bin/rm %s",file);
      gSystem->Exec(command);
    }
    if(singleDetector) break;  
  }
  cout<<"Stopping the job "<<runID<<endl;
  gSystem->Exec("/bin/date");
  int job_stop_time=CWB::Time("now").GetGPS();
 
  int job_elapsed_time  = (job_stop_time-job_start_time);
  int job_elapsed_hour  = int(job_elapsed_time/3600);
  int job_elapsed_min   = int((job_elapsed_time-3600*job_elapsed_hour)/60);
  int job_elapsed_sec   = int(job_elapsed_time-3600*job_elapsed_hour-60*job_elapsed_min);
  int job_data_size_sec = int(detSegs[0].stop-detSegs[0].start);
  double job_speed_factor = double(job_data_size_sec)/double(job_elapsed_time); 
  cout << endl;
  printf("Job Elapsed Time - %02d:%02d:%02d (hh:mm:ss)\n",job_elapsed_hour,job_elapsed_min,job_elapsed_sec);
  printf("Job Speed Factor - %2.1fX\n",job_speed_factor);
  cout << endl;
 
//  return 0;
  exit(0);
}
