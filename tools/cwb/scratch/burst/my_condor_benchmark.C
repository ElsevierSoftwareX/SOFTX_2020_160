#include <vector>

#define NSTAGE 7

void my_condor_benchmark(int start = 6) {

  TString STAGE[NSTAGE] = {"FULL","INIT","STRAIN","CSTRAIN","COHERENCE","SUPERCLUSTER","LIKELIHOOD"};

  CWB::Toolbox TB;

 // TB.checkFile(gSystem->Getenv("CWB_ROOTLOGON_FILE"));
//  TB.checkFile(gSystem->Getenv("CWB_PARAMETERS_FILE"));
 // TB.checkFile(gSystem->Getenv("CWB_UPARAMETERS_FILE"));
  //TB.checkFile(gSystem->Getenv("OUTPUT_DIR"));
  TB.checkFile(output_dir);
  //TString OUTPUT_DIR=TString(gSystem->Getenv("OUTPUT_DIR")); 
  TString OUTPUT_DIR=TString(output_dir); 
  // char job_label[512];sprintf(job_label,"%s%s",data_label,sfactor);

  cout << "Starting reading output directory ..." << endl;
  vector<TString> fileList = TB.getFileListFromDir(OUTPUT_DIR,".root","wave_","",true);
  int nfile = fileList.size();
  float* jobId = new float[nfile];
  float* istat = new float[nfile];
  float* istatc = new float[nfile];
  float* sstat = new float[nfile];
  float* sstatc = new float[nfile];
  float* tstat = new float[nfile];
  float* prev = new float[nfile];
  for(int n=0;n<nfile;n++) {jobId[n]=n;istat[n]=0.0;istatc[n]=0.0;sstat[n]=0.0;sstatc[n]=0.0;tstat[n]=0.0;prev[n]=0.0;}
 /* for(int n=0;n<nfile;n++) {
    //cout << n << " " << fileList[n].Data()<< endl;

    if (n%10==0) cout << "cwb_condor benchmark - " << n << "/" << fileList.size() << " files" << endl;

    vector<float> bench = CWB::Toolbox::getJobBenchmark(fileList[n],jstage,
                                        cwb_bench_res,cwb_bench_factor,cwb_bench_name);
    jobId[n]=bench[0];
    istat[n]=bench[1];
  }*/
 // float b[nfile];
 //  TChain live("liveTime");
  double liveTot = 0.;
  double liveZero = 0.;

  wavearray<double> Trun(500000); Trun = 0.;
  wavearray<double>* Wlag = new wavearray<double>[nIFO+1];
  wavearray<double>* Wslag = new wavearray<double>[nIFO+1];
  wavearray<double> Tdlag;
  wavearray<double> Tlag;

  for(int j=start;j<8;j++){
          for(int n=0;n<nfile;n++) {
				 TChain live("liveTime");  
                live.Add(fileList[n]);    
                liveTot = TB.getLiveTime(nIFO,live,Trun,Wlag,Wslag,Tlag,Tdlag,-1,-1);
                liveTot += TB.getLiveTime(nIFO,live,Trun,Wlag,Wslag,Tlag,Tdlag,0,0);
                vector<float> bench = CWB::Toolbox::getJobBenchmark(fileList[n],j,"JET");
                jobId[n]=bench[0];
                istat[n]=bench[1]-prev[n];
                istatc[n]=bench[1];
                tstat[n]=liveTot;
                sstat[n]=liveTot/(bench[1]-prev[n]);
                sstatc[n]=liveTot/bench[1];
          //    b[n]=bench[1];
                if ((n<200)||(n%1000==0)) {cout<<n<<" : Jobid: "<<bench[0] <<" Stageid: "<<j<< " JET = "<<istat[n]/3600. <<" (Cumulative: "<<bench[1]/3600.<<") Live: "<<liveTot/3600. << "  Speed: "<<sstat[n]<<" (Cumulative: "<<sstatc[n]<<")"<< endl;}
                liveTot = 0.0;
				 Trun = 0.; 
				prev[n] = bench[1];
				live.Clear(); 			  
         }

  cout<<"Njobs: "<< n <<" Average JET @stageid"<<j<< " = "<<TMath::Mean(nfile,istat)/3600.<<" +/- "<<TMath::RMS(nfile,istat)/3600.<< " [hr]  Average Livetime = "<<TMath::Mean(nfile,tstat)/3600<<" +/- "<<TMath::RMS(nfile,tstat)/3600<<" [hr] Average Throughput = "<<TMath::Mean(nfile,sstat)<<" +/- "<<TMath::RMS(nfile,sstat)<<endl;
  }
 cout<<"Njobs: "<< n <<" Average Cumulative JET @stageid"<<j<< " = "<<TMath::Mean(nfile,istatc)/3600.<<" +/- "<<TMath::RMS(nfile,istatc)/3600.<< " [hr]  Average Livetime = "<<TMath::Mean(nfile,tstat)/3600<<" +/- "<<TMath::RMS(nfile,tstat)/3600<<" [hr] Average Cumulative Throughput = "<<TMath::Mean(nfile,sstatc)<<" +/- "<<TMath::RMS(nfile,sstatc)<<endl;


  
}
