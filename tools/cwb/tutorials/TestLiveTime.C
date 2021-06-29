// This macro compare the livetime computed by the CWB pipeline and the approximate one
// computed CWB::Toolbox::getLiveTime
// usage : root -n -l -b ${CWB_PARMS_FILES} 'TestLiveTime.C("live.root",3)'
// WARNING : this macro works only for standard lags (not for slag) !!!
// Author : Gabriele Vedovato


void TestLiveTime(TString liveName, int nIFO=3, TString fName="live_bech.txt") {
// liveName : live root file name
// nIFO     : number of detectors
// fName    : output benchmark file name 

  CWB::Toolbox TB;

  // check if liveName exist
  TB.checkFile(liveName);

  // check if parameter files exist
  TB.checkFile(gSystem->Getenv("CWB_ROOTLOGON_FILE"));
  TB.checkFile(gSystem->Getenv("CWB_PARAMETERS_FILE"));
  TB.checkFile(gSystem->Getenv("CWB_UPARAMETERS_FILE"));

  // compute cat2List & jobList
  vector<waveSegment> cat1List=TB.readSegList(nDQF, DQF, CWB_CAT1);
  vector<waveSegment> cat2List=TB.readSegList(nDQF, DQF, CWB_CAT2);
  vector<waveSegment> jobList=TB.getJobList(cat1List, cat2List, segLen, segMLS, segTHR, segEdge);
  cout<<"Number of standard jobs : " << jobList.size() <<endl<<endl;

  // open livetime root file
  TChain live("liveTime");
  live.Add(liveName);

  wavearray<double> Trun(500000); Trun = 0.;
  wavearray<double>* Wlag = new wavearray<double>[nIFO+1];
  wavearray<double>* Wslag = new wavearray<double>[nIFO+1];
  wavearray<double> Tdlag;
  wavearray<double> Tlag;
  wavearray<double> Rlag;
  wavearray<double> Elag;
  wavearray<double> Olag;

  int lag_number=-1;			// get all lag except zero lag
  int slag_number=0;
  // get total livetime
  double liveTOT=TB.getLiveTime(nIFO,live,Trun,Wlag,Wslag,Tlag,Tdlag,lag_number,slag_number);
  cout << "total livetime from cWB : " << (int)liveTOT << endl;

  // Wlag don't contains lag zero, a new lag array is made which contains the zero lag
  vector<wavearray<double> > lag(nIFO+1);
  for(int i=0;i<nIFO+1;i++) {
    lag[i].resize(Wlag[i].size()+1);
    lag[i][0]=0;			// add lag zero
    for(int j=0;j<Wlag[i].size();j++) lag[i][j+1]=Wlag[i][j];
  }

  // build a sub jobList which contains only the processed set
  vector<waveSegment> sjobList;
  for(int i=0;i<jobList.size();i++) {
    int jobID = jobList[i].index; 
    if(Trun[jobID-1]==0) continue;	// skip jobID with 0 livetime
    sjobList.push_back(jobList[i]);
  }

  FILE *fP;
  if((fP = fopen(fName.Data(), "w")) == NULL) {
    cout << "cannot open output file " << fName.Data() <<". \n";
    gSystem->Exit(1);
  };
  cout << "Write output file : " << fName.Data() << endl;
  fprintf(fP,"# lagID\t liveCWB\t liveTB\t\t 100*(TB-CWB)/CWB\n");

  // get cWB livetime for each lag and compare with the approximate value from CWB::Toolbox::getLiveTime
  vector<double> shiftList(nIFO);
  for(int i=0;i<lag[0].size();i++) {

    lag_number = lag[nIFO][i];

    cout << endl;
    cout << i << " lagID " << lag_number <<  " -> ";
    for(int j=0;j<nIFO;j++) shiftList[j]=lag[j][i];
    for(int j=0;j<nIFO;j++) cout << shiftList[j] << " ";
    cout << endl;

    double liveCWB=TB.getLiveTime(nIFO,live,Trun,Wlag,Wslag,Tlag,Tdlag,lag_number,slag_number);
    cout << "livetime from CWB : " << (int)liveCWB << endl;

    double liveTB = TB.getLiveTime(sjobList, cat2List, shiftList);
    cout << "livetime from TB  : " << (int)liveTB << endl;

    double liveDif = 100*(liveTB-liveCWB)/liveCWB;
    cout << "livetime diff 100*(TB-CWB)/CWB " << liveDif << endl;

    fprintf(fP,"%d\t %f\t %f\t %g\t \n",lag_number,liveCWB,liveTB,liveDif);
  }
  if(fP!=NULL) fclose(fP);

  delete [] Wlag;
  delete [] Wslag;

  exit(0);
}
