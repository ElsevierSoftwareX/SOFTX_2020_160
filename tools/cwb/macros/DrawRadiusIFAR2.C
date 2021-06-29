/*
# Copyright (C) 2019 Francesco Salemi
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
// Draw Radius vs IFAR
// Note : this macro is used to generate the CBC report
// Author : Francesco Salemi
/*
#include "TROOT->h"                               
#include "TSystem->h"                             
#include "TFile->h"                             
#include "TTree->h"                             
#include <fstream>                               
#include <iostream>                              
#include "TGraph->h"                              
#include "TMultiGraph->h"                         
#include "TCanvas->h"                             
#include "TH1F->h"                                
#include "TMath->h"                               
#include "TH1F->h"                                
#include "TH2F->h"                                
#include <TComplex->h>                            
#include <TStyle->h>                              
#include <TRandom->h>                             
#include "TVector3->h"                            
#include "TRotation->h"                           
#include "Math/Vector3Dfwd->h"                    
#include "Math/Rotation3D->h"                     
*/

#define MINIFAR 0.01
#define MAXIFAR 1000.
#define MINRADIUS 0
#define MAXRADIUS 1000
//#define SEL ""
//#define SEL "mass[0]+mass[1]>75.75"
//#define SEL "mass[0]+mass[1]<27.25"

//void DrawRadiusIFAR2(char* sim_file_name, char* mdc_file_name, float shell_volume, char* netdir, double liveTot, float T_ifar=0.0, float T_win=0.2, float internal_volume=0.0, int TRIALS = 1){
TGraphErrors* DrawRadiusIFAR2(char* sim_file_name, char* mdc_file_name, TString SEL, float shell_volume){

  TCanvas* co_canvas;
  co_canvas = new TCanvas("sd2", "SD2", 3, 47, 1000, 802);
  co_canvas->SetGridx();
  co_canvas->SetGridy();
  co_canvas->SetLogx();


  float myifar, ecor, m1, m2, netcc[3], neted, penalty;
  float rho[2];
  double mytime[6];
  float factor, distance, mchirp;
  float mass[2];
  float spin[6];
  float chi[3];

  float chirp[6];
  float range[2];
 
  TFile* filein = new TFile(sim_file_name);
  TTree* sim_org = nullptr; 
  filein->GetObject("waveburst",sim_org);
  gROOT->cd();
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

  
 /* 
  TChain sim("waveburst");
  sim->Add(sim_file_name);
  sim->SetBranchAddress("mass", mass);
  sim->SetBranchAddress("factor", &factor);
  sim->SetBranchAddress("range", range);
  sim->SetBranchAddress("chirp", chirp);
  sim->SetBranchAddress("rho", rho);
  sim->SetBranchAddress("netcc", netcc);
  sim->SetBranchAddress("ifar", &myifar);
  sim->SetBranchAddress("time", mytime);
  sim->SetBranchAddress("spin", spin);
*/

  TFile* filein2 = new TFile(mdc_file_name);
  TTree* mdc_org = nullptr; 
  filein2->GetObject("mdc",mdc_org);
  gROOT->cd();
  TTree* mdc = mdc_org->CopyTree(SEL);
  mdc->SetBranchAddress("time",mytime);
  mdc->SetBranchAddress("mass", mass);
  mdc->SetBranchAddress("factor", &factor);
  mdc->SetBranchAddress("distance", &distance);
  mdc->SetBranchAddress("mchirp", &mchirp);
  mdc->SetBranchAddress("spin", spin);
 
  
  int nevts = (int)mdc->GetEntries(); 
  float shell_volume_per_injection = shell_volume/nevts;

  float maxIFAR = 0.0;
  float CYS = 86400.*365.25;
  if (sim->GetListOfBranches()->FindObject("ifar")) {
        maxIFAR =
            TMath::CeilNint(sim->GetMaximum("ifar") / CYS );
        cout << "Maximum empirically estimated IFAR : " << maxIFAR << " [years]" << endl;
  } else {
        cout << "Missing ifar branch: either use cbc_plots or add it "
                "to wave tree->->->"
             << endl;
        exit(1);
  }
  double dV = 0.0;
  std::vector<double> vdv;
  std::vector<double> vifar;
  std::vector<double> vMtot;

  sim->BuildIndex("ifar","rho[1]");   //BEWARE rho[1] should be search depandent!
  TTreeIndex *I1=(TTreeIndex*)sim->GetTreeIndex(); // get the tree index
  Long64_t* index1=I1->GetIndex(); //create an array of entries in sorted order

  cout << nevts << " injected signals " << sim->GetEntries("ifar>0") << " recovered signals" <<endl; 
  int countv = 0;
  int countvifar = 0; 
  for (int g = 0; g < (int)sim->GetEntries(); g++) {
                sim->GetEntry(index1[g]);
		//cout << g << " " <<index1[g] <<" " << myifar/CYS << " ";
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

		dV = pow(range[1], 2) * shell_volume_per_injection * pow(chirp[0] / 1.22,5./6.);
		vdv.push_back(dV+internal_volume);
                vifar.push_back(myifar);
		vMtot.push_back(mass[0]+mass[1]);
			 	
  }
  cout << endl;
  cout << countvifar << " events vetoed by T_ifar : " << T_ifar << endl;
  cout << countv << " events vetoed by T_win" << endl;
  cout << vdv.size() << " events selected" << endl;
  //      Int_t *mindex = new Int_t[vdv->size()];
  //      TMath::Sort((int)vdv->size(), &vifar[0], mindex, true);
        std::vector<double> vV;
        std::vector<double> veV;
	std::vector<double> vR;
        std::vector<double> veR;
        std::vector<double> vsifar;
	std::vector<double> vseifar;
       
//cout << "DEB2" << endl;
        vV.push_back(vdv[vdv.size()-1]);
        veV.push_back(pow(vdv[vdv.size()-1], 2));
        vsifar.push_back(vifar[vdv.size()-1]);
    //cout << "DEB3" << endl;   
        // break;
        int mcount = 0;
        for (int i = vdv.size()-1; i >= 0; i--) {
                if (vifar[i] == 0) {
                        continue;
                }
                if (vifar[i] == vsifar[mcount]) {
                        vV[mcount] += vdv[i];
                        veV[mcount] += pow(vdv[i], 2);
                } else {
                        vsifar.push_back(vifar[i]);
                       	vseifar.push_back(TMath::Sqrt(TMath::Nint(liveTot * vifar[i])));	
                        vV.push_back(vV[mcount] + vdv[i]);
                        veV.push_back(veV[mcount] + pow(vdv[i], 2));
                        mcount++;
                }
        }
        cout << "Length of ifar/volume vector: " << vV.size() << endl;
	float Tscale = 1.0;
        for (int i = 0; i < vV.size(); i++) {
                veV[i] = TMath::Sqrt(veV[i]);
                vfar[i] *= TRIALS * CYS;
                vefar[i] *= TRIALS * CYS;
		vsifar[i] /= (TRIALS * CYS);
		vR.push_back(
                    pow(3. / (4. * TMath::Pi()* Tscale) * vV[i], 1. / 3.));
                veR.push_back(pow(3. / (4 * TMath::Pi() * Tscale), 1. / 3.) *
                              1. / 3 * pow(vV[i], -2. / 3.) * veV[i]);
		//cout <<  vsifar[i] << " " << 	
        }





  TGraphErrors *co_gr = new TGraphErrors(vR.size(), &vsifar[0], &vR[0], 0, &veR[0]);
  co_gr->GetYaxis()->SetTitle("Sensitive Distance [Mpc]");
  co_gr->GetXaxis()->SetTitle("Inverse False Alarm Rate [yr]");
  co_gr->GetXaxis()->SetLimits(MINIFAR,MAXIFAR);
  co_gr->GetYaxis()->SetRangeUser(MINRADIUS,MAXRADIUS);
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
  co_gr->SetLineColor(kBlue);
  co_gr->SetTitle("");
  co_gr->SetFillColor(kBlue);
  co_gr->SetFillStyle(3001);
 // h->Draw("e2same"); 
  co_gr->Draw("aple3");
 

 

 /* char fname[1024];  
  sprintf(fname, "%s/ROC_IFAR2.png", netdir);
  co_canvas->Update();
  co_canvas->SaveAs(fname);
*/
 return co_gr;
    //exit(0);

}

