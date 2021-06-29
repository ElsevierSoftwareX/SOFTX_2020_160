// Create Range Plot  
// Note : this macro is used to generate the CBC report
// Author : Francesco Salemi

#include "TROOT.h"                               
#include "TSystem.h"                             
#include "TFile.h"                             
#include "TTree.h" 
#include "TTreeIndex.h"                             
#include <fstream>                               
#include <iostream>                              
#include "TGraphErrors.h"                              
#include "TMultiGraph.h"                         
#include "TCanvas.h"                             
#include "TH1F.h"                                
#include "TMath.h"                               
#include "TH1F.h"                                
#include "TH2F.h"                                
#include <TComplex.h>                            
#include <TStyle.h>                              
#include <TRandom.h>                             
#include "TVector3.h"                            
#include "TRotation.h"                           
#include "Math/Vector3Dfwd.h"                    
#include "Math/Rotation3D.h"                     

#define MINMAXIFAR 4278.
#define VT 0.40017445


TGraphErrors* CreateGraphRadiusIFAR(char* sim_file_name, char* mdc_file_name, TString SEL, float shell_volume, Color_t color=kBlue, TString opt="default",double liveTot=1e6, float T_ifar=0.0, float T_win=0.2, int TRIALS = 1, int nIFO = 2){

  float myifar, netcc[3];
  float rho[2];
  double mytime[6];
  float factor, distance, mchirp;
  float mass[2];
  float spin[6];
 // float chi[3];

  float chirp[6];
  float range[2];
 
  CWB::CBCTool cbcTool;
 
  TFile* filein = new TFile(sim_file_name);
  TTree* sim_org = nullptr; 
  filein->GetObject("waveburst",sim_org);
  gROOT->cd();
  if(!sim_org->GetListOfBranches()->FindObject("chip")) {  
		cout << "Adding Chi_p branch to wave tree"<<endl;
		cbcTool.AddChip(sim_file_name, "waveburst");	
  }
  TTree* sim = sim_org->CopyTree(SEL);
  sim->SetBranchAddress("mass", mass);
  sim->SetBranchAddress("factor", &factor);
  sim->SetBranchAddress("range", range);
  sim->SetBranchAddress("chirp", chirp);
  sim->SetBranchAddress("rho", rho);
  sim->SetBranchAddress("netcc", netcc);
  sim->SetBranchAddress("ifar", &myifar);
  sim->SetBranchAddress("time", mytime);
  sim->SetBranchAddress("spin", spin);

 
  TFile* filein2 = new TFile(mdc_file_name);
  TTree* mdc_org = nullptr; 
  filein2->GetObject("mdc",mdc_org);
  gROOT->cd();
  SEL.ReplaceAll("chirp[0]","mchirp");
  if(!mdc_org->GetListOfBranches()->FindObject("chip")) {
        cout << "Adding Chi_p branch to mdc tree" << endl;
        cbcTool.AddChip(mdc_file_name, "mdc");
  }
  TTree* mdc = mdc_org->CopyTree(SEL);
  mdc->SetBranchAddress("time",mytime);
  mdc->SetBranchAddress("mass", mass);
  mdc->SetBranchAddress("factor", &factor);
  mdc->SetBranchAddress("distance", &distance);
  mdc->SetBranchAddress("mchirp", &mchirp);
  mdc->SetBranchAddress("spin", spin);
 
  double dV = 0.0;
  std::vector<double> vdv;
  std::vector<double> vifar;
  std::vector<double> vx;
  std::vector<double> vV;
  std::vector<double> veV;
  std::vector<double> vR;
  std::vector<double> veR;
  std::vector<double> vsifar;
  std::vector<double> vseifar;
  std::vector<double> vVr;
  std::vector<double> veVr;
  std::vector<double> vRr;
  std::vector<double> veRr;
  std::vector<double> vsrho;
  TGraphErrors* co_gr = nullptr;
  
  int nevts = (int)mdc->GetEntries();
  if (nevts == 0) {
	cout <<"No events left after cut!"<<endl;
	vR.push_back(0.0);
	vsifar.push_back(0.0);
	co_gr = new TGraphErrors(vR.size(), &vsifar[0], &vR[0], 0, 0);	
	return co_gr;
  } 
  float shell_volume_per_injection = shell_volume/nevts;
  int ifactor;
  float maxIFAR = 0.0;
  float CYS = 86400.*365.25;
  if (sim->GetListOfBranches()->FindObject("ifar")) {
        maxIFAR =
            TMath::CeilNint(sim->GetMaximum("ifar") / CYS );
       // cout << "Maximum empirically estimated IFAR : " << maxIFAR << " [years]" << endl;
  } else {
        cout << "Missing ifar branch: either use cbc_plots or add it "
                "to wave tree."
             << endl;
        exit(1);
  }
  if(opt.Contains("rho")) {
	sim->BuildIndex("0","rho[1]*10000");   //BEWARE rho[1] could be search depandent!
  } else if (opt.Contains("time")) { 
  
    sim->BuildIndex("time[0]","rho[1]*1000");   //BEWARE rho[1] could be search depandent!        
  
  }else {
	sim->BuildIndex("ifar","rho[1]*1000");   //BEWARE rho[1] could be search depandent!
  }
  TTreeIndex *I1=(TTreeIndex*)sim->GetTreeIndex(); // get the tree index
  Long64_t* index1=I1->GetIndex(); //create an array of entries in sorted order
  
  //cout << nevts << " injected signals " << sim->GetEntries("ifar>0") << " recovered signals" <<endl; 
  int countv = 0;
  int countvifar = 0;
  bool Error = false;
  double mpcTscale = 1.0;  
  for (int g = 0; g < (int)sim->GetEntries(); g++) {
                sim->GetEntry(index1[g]);
                ifactor = (int)factor - 1;
		//cout << "g=" << g << " trueindex=" <<index1[g] <<" IFAR=" << myifar/CYS << " RHO1=" << rho[1] <<endl;
                if (myifar <= T_ifar * CYS) {
			countvifar++;				
			//cout << g << " " <<index1[g] <<" " << myifar/CYS << " ";	
                        continue;
		
                }

		if ((mytime[0] - mytime[nIFO]) < -T_win ||
                    (mytime[0] - mytime[nIFO]) > 2 * T_win) {
                        countv++;
                        continue;
                } // NOT checking for detector 1 and 2: very small bias...

		if (opt.Contains("DDistrVolume")){
			dV = shell_volume_per_injection;
		}
		else if (opt.Contains("DDistrUniform")){
			dV = pow(range[1], 2) * shell_volume_per_injection;
		}
		else if (opt.Contains("DDistrChirpMass")){
			dV = pow(range[1], 2) * shell_volume_per_injection * pow(chirp[0] / 1.22,5./6.);
		}
		else if (opt.Contains("Redshift")){
				dV = VT;
				mpcTscale = Tscale*1e-9;   //Conversion from Gpc^3 to Mpc^3
		}
		else {
			if(!Error){
			    cout << "No defined distance distribution? "
			    "WARNING: Assuming uniform in volume"
			    << endl;
			    Error=1;
			 }   
			 dV = shell_volume_per_injection;
                }
		//vdv.push_back(dV+internal_volume);
		vdv.push_back(dV);
  		vifar.push_back(myifar);
		opt.Contains("time") ? vx.push_back(time[0]) : vx.push_back(rho[1]);
			 	
  }
  
  //cout << endl;
  //cout << countvifar << " events vetoed by T_ifar : " << T_ifar << endl;
  //cout << countv << " events vetoed by T_win" << endl;
  //cout << vdv.size() << " events selected" << endl;
       
	//cout << "DEB2" << endl;
	vV.push_back(vdv[vdv.size()-1]);
	veV.push_back(pow(vdv[vdv.size()-1], 2));
	if(vifar[vdv.size()-1]< MINMAXIFAR*CYS ){vsifar.push_back(vifar[vdv.size()-1]);}
	else {vsifar.push_back(MINMAXIFAR*CYS);}
	vVr.push_back(vdv[vdv.size()-1]);
	veVr.push_back(pow(vdv[vdv.size()-1], 2));
	vsrho.push_back(vx[vdv.size()-1]);
	//cout << "vifar[vdv.size()-1] = " << vifar[vdv.size()-1] << endl; 
	//cout << "vsifar[0] = " << vsifar[0] << endl;  
	// break;
	int mcount_ifar = 0;
	int mcount_rho = 0;
	for (int i = vdv.size()-1; i >= 0; i--) {
		      if (vifar[i] == 0) {
		              continue;
		      }
		      if (vifar[i] > vsifar[0]) {
			      vV[0] += vdv[i];
                              veV[0] += pow(vdv[i], 2);
	
		      }	else if (vifar[i] == vsifar[mcount_ifar]) {
		              vV[mcount_ifar] += vdv[i];
		              veV[mcount_ifar] += pow(vdv[i], 2);
		      } else {
		              vsifar.push_back(vifar[i]);
		              vseifar.push_back(TMath::Sqrt(TMath::Nint(liveTot * vifar[i])));	
		              vV.push_back(vV[mcount_ifar] + vdv[i]);
		              veV.push_back(veV[mcount_ifar] + pow(vdv[i], 2));
		              mcount_ifar++;
		      }
		     // cout << i << " " << vsifar[mcount_ifar] << " " << vV[mcount_ifar] << endl;	
		      if (vx[i] == vsrho[mcount_rho]) {
		              vVr[mcount_rho] += vdv[i];
		              veVr[mcount_rho] += pow(vdv[i], 2);
		      } else {
		              vsrho.push_back(vx[i]);
		              vVr.push_back(vVr[mcount_rho] + vdv[i]);
		              veVr.push_back(veVr[mcount_rho] + pow(vdv[i], 2));
		              mcount_rho++;
		      }
	}
	//cout << "Length of ifar/volume vector: " << vV.size() << endl;
	//cout << "Length of rho/volume vector: " << vVr.size() << endl;
//	float 
	for (int i = 0; i < (int) vV.size(); i++) {
		veV[i] = TMath::Sqrt(veV[i]);
		vsifar[i] /= (TRIALS * CYS);
		vR.push_back(
		          pow(3. / (4. * TMath::Pi()* mpcTscale) * vV[i], 1. / 3.));
		      veR.push_back(pow(3. / (4 * TMath::Pi() * mpcTscale), 1. / 3.) *
		                    1. / 3 * pow(vV[i], -2. / 3.) * veV[i]);
	//	cout << i << " " <<  vsifar[i] << " " << vR[i] << " " <<vV[i] <<endl; 	
	}

	for (int i = 0; i < (int) vVr.size(); i++) {
		      veVr[i] = TMath::Sqrt(veVr[i]);
					vRr.push_back(
		          pow(3. / (4. * TMath::Pi()* mpcTscale) * vVr[i], 1. / 3.));
		      veRr.push_back(pow(3. / (4 * TMath::Pi() * mpcTscale), 1. / 3.) *
		                    1. / 3 * pow(vVr[i], -2. / 3.) * veVr[i]);
		//cout <<  vsrho[i] << " ";
		 	
	}
  //cout << endl;

	//std::vector<TGraphErrors> v_gr(2);
  //TGraphErrors *co_gr = new TGraphErrors[2];
  //TGraphErrors co_gr;
  //TGraphErrors* co_gr = new TGraphErrors(vR.size(), &vsifar[0], &vR[0], 0, &veR[0]);
  if(opt.Contains("rho")) {co_gr = new TGraphErrors(vRr.size(), &vsrho[0], &vRr[0], 0, &veRr[0]);}
  else {co_gr = new TGraphErrors(vR.size(), &vsifar[0], &vR[0], 0, &veR[0]);}
  //co_gr[0] = TGraphErrors(vR.size(), &vsifar[0], &vR[0], 0, &veR[0]);
  //co_gr->GetYaxis()->SetTitle("Sensitive Distance [Mpc]");
  //co_gr->GetXaxis()->SetTitle("Inverse False Alarm Rate [yr]");
 // co_gr->GetXaxis()->SetLimits(MINIFAR,MAXIFAR);
 // co_gr->GetYaxis()->SetRangeUser(MINRADIUS,MAXRADIUS);
  co_gr->GetXaxis()->SetTitleOffset(1.3);
  co_gr->GetYaxis()->SetTitleOffset(1.25);
  co_gr->GetXaxis()->SetTickLength(0.01);
  co_gr->GetYaxis()->SetTickLength(0.01);
  co_gr->GetXaxis()->CenterTitle(kTRUE);
  co_gr->GetYaxis()->CenterTitle(kTRUE);
  co_gr->GetXaxis()->SetTitleFont(42);
  co_gr->GetXaxis()->SetLabelFont(42);
  co_gr->GetYaxis()->SetTitleFont(42);
  co_gr->GetYaxis()->SetLabelFont(42);
  co_gr->SetMarkerStyle(20);
  co_gr->SetMarkerSize(1.0);
  co_gr->SetMarkerColor(1);
  co_gr->SetLineColor(color);
  co_gr->SetLineWidth(3);
  co_gr->SetTitle("");
  co_gr->SetFillColor(color);
  co_gr->SetFillStyle(3001);
  //co_gr->Draw("aple3");

  filein->Close();
  filein2->Close();
  delete filein, filein2;
  delete mdc, mdc_org, sim, sim_org;
  vifar.clear(); vdv.clear(); vx.clear(); vsifar.clear(); vseifar.clear(); vsrho.clear(); 
  vR.clear(); veR.clear(); veRr.clear(); vV.clear(); veV.clear(); vRr.clear(); vVr.clear(); 
  veVr.clear();
  return co_gr;

    //exit(0);

}

