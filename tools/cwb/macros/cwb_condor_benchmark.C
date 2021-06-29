/*
# Copyright (C) 2019 Gabriele Vedovato, Marco Drago
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


// draw benchmark plots of computational load : used by the cwb_condor command

#include <vector>

#define NSTAGE 7

TCanvas* canvas;
TH1F* hplot;
TGraph* gplot;
TH2F* h2plot;

void cwb_condor_benchmark() {

  TString STAGE[NSTAGE] = {"FULL","INIT","STRAIN","CSTRAIN","COHERENCE","SUPERCLUSTER","LIKELIHOOD"};

  CWB::Toolbox TB;

  TB.checkFile(gSystem->Getenv("CWB_ROOTLOGON_FILE"));
  TB.checkFile(gSystem->Getenv("CWB_PARAMETERS_FILE"));
  TB.checkFile(gSystem->Getenv("CWB_UPARAMETERS_FILE"));

  // get benchmark options
  TString cwb_jstage_name = "FULL";
  TString cwb_fstage_name = "FULL";
  TString cwb_plot_type   = "HIST";
  TString cwb_plot_save   = "";
  TString cwb_bench_name  = "JET"; 		// JET/SET/MEM/SET
  double  cwb_bench_min   = -1;			// min bench value
  double  cwb_bench_max   = -1;			// max bench value
  int     cwb_bench_res   = -1;			// res bench value (disabled)
  int     cwb_bench_factor= -1;			// factor bench value (disabled)
  TString cwb_bench_opts=TString(gSystem->Getenv("CWB_BENCH_OPTS"));
  if(cwb_bench_opts!="") {
    int bstage=false;
    TString option="";
    // get fstage
    option = TB.getParameter(cwb_bench_opts,"--fstage");
    if(option!="") {
      cwb_fstage_name = option;
      cwb_fstage_name.ToUpper();
      bstage=false;
      for(int i=0;i<NSTAGE;i++) if(cwb_fstage_name==STAGE[i]) bstage=true;
      if(!bstage) {
        cout << "cwb_condor_benchmark.C : fstage " 
        << cwb_fstage_name.Data() << " not correct" << endl;
        cout << "select : " << endl;
        for(int i=0;i<NSTAGE;i++) cout << STAGE[i] << endl;
        gSystem->Exit(1); 
      }
    } 
    // get jstage
    option = TB.getParameter(cwb_bench_opts,"--jstage");
    if(option!="") {
      cwb_jstage_name = option;
      cwb_jstage_name.ToUpper();
      bstage=false;
      for(int i=0;i<NSTAGE;i++) if(cwb_jstage_name==STAGE[i]) bstage=true;
      if(!bstage) {
        cout << "cwb_condor_benchmark.C : jstage " 
        << cwb_jstage_name.Data() << " not correct" << endl;
        cout << "select : " << endl;
        for(int i=0;i<NSTAGE;i++) cout << STAGE[i] << endl;
        gSystem->Exit(1); 
      }
    } 
    // get res
    option = TB.getParameter(cwb_bench_opts,"--res");
    if(option!="") cwb_bench_res = option.Atoi();
    // get factor
    option = TB.getParameter(cwb_bench_opts,"--factor");
    if(option!="") cwb_bench_factor = option.Atoi();
    // get bench
    option = TB.getParameter(cwb_bench_opts,"--bench");
    if(option!="") {
      cwb_bench_name=option;
      cwb_bench_name.ToUpper();
      if((cwb_bench_name!="MEM")   && (cwb_bench_name!="JET")   &&
         (cwb_bench_name!="JFS")   && (cwb_bench_name!="SET")   &&
         (cwb_bench_name!="PSIZE") && (cwb_bench_name!="CSIZE") &&
         (cwb_bench_name!="THR")   && (cwb_bench_name!="GT")) {
        cout << "cwb_condor_benchmark.C : Error - "
             << "currently only JET/MEM/JFS/SET/THR/PSIZE/CSIZE/GT is implemented" << endl;
        gSystem->Exit(1); 
      }
    }
    if((cwb_bench_res>=0) && 
       (cwb_bench_name!="THR")&&(cwb_bench_name!="PSIZE")&&(cwb_bench_name!="CSIZE")) {
      cout << "cwb_condor_benchmark.C : Error - " 
           << "'--res' parm must declared only with bench=THR/PSIZE/CSIZE" << endl;
      gSystem->Exit(1); 
    }
    if(((cwb_bench_name=="THR")  &&(cwb_bench_res<0)) ||
       ((cwb_bench_name=="PSIZE")&&(cwb_bench_res<0)) ||
       ((cwb_bench_name=="CSIZE")&&(cwb_bench_res<0))) {
      cout << "cwb_condor_benchmark.C : Error - "
           << "when bench=THR/PSIZE/CSIZE the '--res' parm must be > 0" << endl;
      gSystem->Exit(1); 
    }
    if(((cwb_bench_name=="THR")  &&(cwb_jstage_name!="COHERENCE")) ||
       ((cwb_bench_name=="PSIZE")&&(cwb_jstage_name!="COHERENCE")) ||
       ((cwb_bench_name=="GT")   &&(cwb_jstage_name!="COHERENCE")) ||
       ((cwb_bench_name=="CSIZE")&&(cwb_jstage_name!="COHERENCE"))) {
      cout << "cwb_condor_benchmark.C : Error - "
           << "bench=THR/PSIZE/CSIZE/GT is available only in COHERENCE stage" << endl;
      gSystem->Exit(1); 
    }
    // get plot type
    option = TB.getParameter(cwb_bench_opts,"--plot");
    if(option!="") {
      cwb_plot_type=option;
      cwb_plot_type.ToUpper();
      if((cwb_plot_type!="HIST")&&(cwb_plot_type!="GRAPH")&&(cwb_plot_type!="HIST2")) {
        cout << "cwb_condor_benchmark.C : currently only HIST/GRAPH/HIST2 is implemented" << endl;
        gSystem->Exit(1); 
      }
    }
    // get plot save
    option = TB.getParameter(cwb_bench_opts,"--save");
    if(option!="") {
      cwb_plot_save=option;
      if(option.EndsWith(".png")) {	    // if ends with .png then it is a file
        // get directory
        if(option.Contains("/")) {
          option.Remove(option.Last('/'));  // strip file name
        } else option=".";
      }
      CWB::Toolbox::checkFile(option); 	    // check if dir exist
    }

    // get min
    option = TB.getParameter(cwb_bench_opts,"--min");
    if(option!="") cwb_bench_min = option.Atof();
    // get max
    option = TB.getParameter(cwb_bench_opts,"--max");
    if(option!="") cwb_bench_max = option.Atof();
    if(cwb_bench_min==0) cwb_bench_min=-1; 
    if(cwb_bench_max==0) cwb_bench_max=-1; 
    if(cwb_bench_min>0 && cwb_bench_max>0 && cwb_bench_min>=cwb_bench_max) {
      cout << "cwb_condor_benchmark.C : Error min must be < max " << endl;
      gSystem->Exit(1); 
    }
  }

  // if bench = JET/SET/JFS/MEM the values of the current stage are reported in the next stage 
  bool uns = false;		      // uns : use next stage
  if(cwb_bench_name=="JET") uns = true;
  if(cwb_bench_name=="SET") uns = true;
  if(cwb_bench_name=="JFS") uns = true;
  if(cwb_bench_name=="MEM") uns = true;

  CWB_STAGE jstage=CWB_STAGE_FINISH;                             
  if(cwb_jstage_name=="INIT")         jstage = uns ? CWB_STAGE_STRAIN : CWB_STAGE_INIT;
  if(cwb_jstage_name=="STRAIN")       jstage = uns ? CWB_STAGE_CSTRAIN : CWB_STAGE_STRAIN;
  if(cwb_jstage_name=="CSTRAIN")      jstage = uns ? CWB_STAGE_COHERENCE : CWB_STAGE_CSTRAIN;
  if(cwb_jstage_name=="COHERENCE")    jstage = uns ? CWB_STAGE_SUPERCLUSTER: CWB_STAGE_COHERENCE;
  if(cwb_jstage_name=="SUPERCLUSTER") jstage = uns ? CWB_STAGE_LIKELIHOOD : CWB_STAGE_SUPERCLUSTER;
  if(cwb_jstage_name=="LIKELIHOOD")   jstage = uns ? CWB_STAGE_SAVE : CWB_STAGE_LIKELIHOOD;
  if(cwb_jstage_name=="FULL") {       jstage = CWB_STAGE_FINISH;
                                      if(cwb_bench_name=="SET") cwb_bench_name = "JET";}

  TString fstage="wave_";                             
  if(cwb_fstage_name=="INIT")         fstage = "init_";
  if(cwb_fstage_name=="STRAIN")       fstage = "strain_";
  if(cwb_fstage_name=="CSTRAIN")      fstage = "cstrain_";
  if(cwb_fstage_name=="COHERENCE")    fstage = "coherence_";
  if(cwb_fstage_name=="SUPERCLUSTER") fstage = "supercluster_";
  if(cwb_fstage_name=="LIKELIHOOD")   fstage = "wave_";
  if(cwb_fstage_name=="FULL")         fstage = "wave_";

  if(cwb_jstage_name==cwb_fstage_name) jstage = CWB_STAGE_FINISH;

  if(nfactor<=0) nfactor=1;     // fix nfactor when nfactor is not defined
  // get the number of job submit by condor
  char full_condor_dir[1024];
  sprintf(full_condor_dir,"%s/%s",work_dir,condor_dir);
  char condor_dag_file[1024];
  sprintf(condor_dag_file,"%s/%s%s.dag",full_condor_dir,data_label,"");
  Long_t id,size=0,flags,mt;
  int estat = gSystem->GetPathInfo(condor_dag_file,&id,&size,&flags,&mt);
  vector<int> jobList;
  if (estat==0) jobList=TB.getCondorJobList(full_condor_dir, data_label);
  else {
    cout << endl << "cwb_condor_benchmark.C : error - dag file not exist !!!"
         << endl << condor_dag_file << endl << endl;
    if (!online) gSystem->Exit(1);
  }
  int ncondor_jobs = jobList.size();
  int max_jobid=0;
  for(int i=0;i<ncondor_jobs;i++) if(jobList[i]>max_jobid) max_jobid=jobList[i];

  // factor label
  char sfactor[32]="";
  if(simulation) {
    double factor=factors[nfactor-1]; 
    if(simulation==3) {
      if(factor<0)  sprintf(sfactor,"_n%g_",fabs(factor));
      if(factor==0) sprintf(sfactor,"_z%g_",factor);
      if(factor>0)  sprintf(sfactor,"_p%g_",factor);
    } else if(simulation==4) {
      int ioffset = int(factors[0])<=0 ? 1 : int(factors[0]);       
      //ioffset+=nfactor-1; 
      sprintf(sfactor,"_%i_",ioffset);
    } else          sprintf(sfactor,"_%g_",factor);
  }
  char job_label[512];
  if(fstage=="wave_") sprintf(job_label,"%s%s",data_label,sfactor);
  else                sprintf(job_label,"%s",data_label);

  cout << "Starting reading output directory ..." << endl;
  vector<TString> fileList = TB.getFileListFromDir(output_dir,".root",fstage,job_label,true);
  if (online) fileList = TB.getFileListFromDir(output_dir,".root",fstage,"",true);
  int nfile = fileList.size();
  float* jobId = new float[nfile];
  float* istat = new float[nfile];
  for(int n=0;n<nfile;n++) {jobId[n]=n;istat[n]=0;}
  for(int n=0;n<nfile;n++) {
    //cout << n << " " << fileList[n].Data()<< endl;

    if (n%100==0) cout << "cwb_condor benchmark - " << n << "/" << fileList.size() << " files" << endl;

    vector<float> bench = CWB::Toolbox::getJobBenchmark(fileList[n],jstage,
                                        cwb_bench_res,cwb_bench_factor,cwb_bench_name);
    jobId[n]=bench[0];
    istat[n]=bench[1];
  }

  //for(int n=0;n<nfile;n++) {if(istat[n]>0) cout << jobId[n] << " " << istat[n] << endl;}

  if(cwb_plot_save!="") gROOT->SetBatch(true);

  canvas = new TCanvas("cwb_condor_benchmark", "LVC experiment", 300,40, 800, 600);
  canvas->Clear();
  canvas->ToggleEventStatus();
  canvas->SetGridx(true);
  canvas->SetGridy(true);
  canvas->SetFillColor(kWhite);

  gStyle->SetTitleH(0.050);
  gStyle->SetTitleW(0.98);
  gStyle->SetTitleY(0.98);
  gStyle->SetTitleFont(72);
//  gStyle->SetLineColor(kWhite);
  gStyle->SetTitleFillColor(kWhite);
  gStyle->SetNumberContours(256);
  gStyle->SetMarkerStyle(7);
  gStyle->SetMarkerSize(2);
  gStyle->SetPalette(1,0);

  if((cwb_bench_name=="JET")||(cwb_bench_name=="SET")) {
    for(int n=0;n<nfile;n++) istat[n]/=3600.;	// sec -> hour
  }
  if(cwb_bench_name=="JFS") {
    for(int n=0;n<nfile;n++) istat[n]/=(1024.*1024.);	// bytes -> MB
  }
  if((cwb_bench_name=="THR")&&(cwb_bench_max<=0)) {
    for(int n=0;n<nfile;n++) if(istat[n]>50) istat[n]=50;
  }

  // compute istat min/max
  double istat_min=100000;
  double istat_max=0; 
  for(int n=0;n<nfile;n++) {
    if(istat[n]>0 && istat[n]<istat_min) istat_min=istat[n];
    if(istat[n]>0 && istat[n]>istat_max) istat_max=istat[n];
  }
  if(cwb_bench_max>0) istat_max=cwb_bench_max;

  //for(int n=0;n<nfile;n++) {if(istat[n]>0) cout << jobId[n] << " " << istat[n] << endl;}

  TString ptitle;
  if(cwb_bench_name=="SET")   ptitle="SET - Elapsed Time in Stage "+cwb_jstage_name;
  if(cwb_bench_name=="JET")   ptitle="JET - Job Elapsed Time at the end of Stage "+cwb_jstage_name;
  if(cwb_bench_name=="JFS")   ptitle="JFS - Job File Size at the end of Stage "+cwb_jstage_name;
  if(cwb_bench_name=="MEM")   ptitle="MEM - Virtual Memory used at the end of Stage "+cwb_jstage_name;
  if(cwb_bench_name=="THR")   {ptitle=cwb_jstage_name+" - Threshold @ ResolutionID : ";
                               ptitle+=cwb_bench_res;ptitle+=" & FactorID : ";
                               ptitle+=cwb_bench_factor;}
  if(cwb_bench_name=="GT")    {ptitle=cwb_jstage_name+" - Gating Time @ FactorID : ";
                               ptitle+=cwb_bench_factor;}
  if(cwb_bench_name=="PSIZE") {ptitle=cwb_jstage_name+" - Number of pixels per lag @ ResolutionID : ";
                               ptitle+=cwb_bench_res;ptitle+=" & FactorID : ";
                               ptitle+=cwb_bench_factor;}
  if(cwb_bench_name=="CSIZE") {ptitle=cwb_jstage_name+" - Number of clusters per lag @ ResolutionID : ";
                               ptitle+=cwb_bench_res;ptitle+=" & FactorID : ";
                               ptitle+=cwb_bench_factor;}
  TString ytitle;
  if(cwb_bench_name=="SET")   ytitle="hour";
  if(cwb_bench_name=="JET")   ytitle="hour";
  if(cwb_bench_name=="JFS")   ytitle="MB";
  if(cwb_bench_name=="MEM")   ytitle="MB";
  if(cwb_bench_name=="THR")   ytitle="energy";
  if(cwb_bench_name=="GT")    ytitle="sec";
  if(cwb_bench_name=="PSIZE") ytitle="count";
  if(cwb_bench_name=="CSIZE") ytitle="count";

  char ofile_name[1024]="";
  if(cwb_plot_save!="") {
    if(cwb_plot_save.EndsWith(".png")) {	        // user defined file name
      sprintf(ofile_name,"%s",cwb_plot_save.Data());    
    } else {						// user define directory
      TString extraLabel="";
      if(cwb_bench_res>=0) {extraLabel+="_RESID_";extraLabel+=cwb_bench_res;}
      if(cwb_bench_factor>=0) {extraLabel+="_FACTORID_";extraLabel+=cwb_bench_factor;}
      // file name is provided by macro
      sprintf(ofile_name,"%s/benchmark_%s_%s_%s_%s%s.png",cwb_plot_save.Data(),data_label,
              cwb_bench_name.Data(),cwb_jstage_name.Data(),cwb_plot_type.Data(),extraLabel.Data());
    } 
  }

  if(cwb_plot_type=="GRAPH") {
    // sort data
    float* x = new float[max_jobid+1];
    float* y = new float[max_jobid+1];
    for(int n=0;n<=max_jobid;n++) {x[n]=n;y[n]=0;}
    for(int n=0;n<nfile;n++) {
      if(istat[n]>0) y[int(jobId[n])]=istat[n];
    }
    int N=0;
    for(int n=0;n<=max_jobid;n++) {
      if(cwb_bench_name=="GT") {x[N]=x[n];y[N]=y[n];N++;}
      else                     if(y[n]>0) {x[N]=x[n];y[N]=y[n];N++;}
    } 

    gplot = new TGraph(N,x,y);
    gplot->SetName("benchmark");
    gplot->SetTitle("graph");
    gplot->SetLineColor(kRed);
    gplot->SetMarkerColor(kRed);
    gplot->Draw("APL");
    gplot->GetHistogram()->SetXTitle("Job id");
    gplot->GetHistogram()->SetYTitle(ytitle);
    gplot->SetTitle(ptitle);
    if(cwb_plot_save!="") {
      cout << "cwb_condor_benchmark.C : dump -> " << ofile_name << endl;
      canvas->Print(ofile_name);
      exit(0);
    } 
  }

  if(cwb_plot_type=="HIST") {
    if((cwb_bench_name=="THR")||(cwb_bench_name=="PSIZE")||(cwb_bench_name=="CSIZE")||(cwb_bench_name=="GT")) {
      canvas->SetLogy(true);
      istat_max+=1;
    }
    hplot = new TH1F("benchmark","hist",100,istat_min,istat_max);
    for(int n=0;n<nfile;n++) {
      if(cwb_bench_name=="GT") hplot->Fill(istat[n]);	
      else                     if(istat[n]>0) hplot->Fill(istat[n]);	
    }
    hplot->Draw("HIST");
    hplot->SetXTitle(ytitle);
    hplot->SetYTitle("#counts");
    hplot->SetTitle(ptitle);
    hplot->SetLineColor(kRed);
    hplot->SetFillColor(kRed);
    hplot->GetYaxis()->SetTitleOffset(1.3);
    if(cwb_plot_save!="") {
      cout << "cwb_condor_benchmark.C : dump -> " << ofile_name << endl;
      canvas->Print(ofile_name);
      exit(0);
    }
  }

  if(cwb_plot_type=="HIST2") {

    int NX = sqrt(max_jobid); NX=NX-NX%10; if(NX==0) NX=max_jobid;
    int NY = max_jobid/NX; if(max_jobid%NX) NY+=1;

    h2plot = new TH2F("benchmark","hist2",NX,0,NX,NY,0,NX*NY);
    for(int n=0;n<nfile;n++) {
      int job = jobId[n]-1;
      int nx  = job%NX;
      int ny  = job;
      double z = istat[n];
      if(cwb_bench_max>0) if(istat[n]>cwb_bench_max) z=cwb_bench_max; 
      //cout << job << " " << nx << " " << ny << " " << istat[n] << " " << endl;
      if(istat[n]) h2plot->Fill(nx,ny,z);	
    }
    canvas->SetLogz(false);
    h2plot->Draw("colz");
    char h2title[256];sprintf(h2title," : condor jobs = %d/%d",nfile,ncondor_jobs);
    h2plot->SetTitle(ptitle+h2title);
    h2plot->SetStats(kFALSE);
    h2plot->SetXTitle("");
    h2plot->SetYTitle("job#");
    h2plot->SetZTitle(ytitle);
    h2plot->SetLineColor(kRed);
    h2plot->SetFillColor(kRed);
    h2plot->GetYaxis()->SetTitleOffset(1.3);
    h2plot->GetZaxis()->SetTitleOffset(0.5);
    h2plot->GetZaxis()->SetRangeUser(0,int(istat_max)+1);
    //h2plot->GetZaxis()->SetNdivisions(-309);
    h2plot->GetYaxis()->CenterTitle(true);
    h2plot->GetZaxis()->CenterTitle(true);
    if(cwb_plot_save!="") {
      cout << "cwb_condor_benchmark.C : dump -> " << ofile_name << endl;
      canvas->Print(ofile_name);
      exit(0);
    }
  }
}
