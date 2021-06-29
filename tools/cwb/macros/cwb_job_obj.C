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


// Coherent WaveBurst production script for GW network (obsolete)

{

  #include <vector>

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// declarations
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  TFile* jfile;
  TDirectory* cdsim;
  TDirectory* cdprod;

  char tdf00[1024];
  char tdf90[1024];
  char end_CED[1024];

  netcluster wc;           
  netcluster* pwc;
  detector* pD[NIFO_MAX];    
  WSeries<double>* pTF[NIFO_MAX];  
  WSeries<double> wM;          // mdc WSeries

  jfile = gROOT->GetFile(); 
  jfile->ReOpen("UPDATE");
  if(jfile==NULL) {cout << "Error opening input root file " << endl;exit(1);}
  jfile->ls();

  // read network object
  network NET = *(network*)jfile->Get("net");
  // read config object
  CWB::config config = *(CWB::config*)jfile->Get("config");
  // restore config parameters in CINT
  config.Export();  

  if(simulation) cdsim  = (TDirectory*)jfile->Get("sim");
  else           cdprod = (TDirectory*)jfile->Get("prod"); 

  nSky=1000;
  NET.nSky=nSky;

  cedDump=true;
  TObjArray* token = TString(jfile->GetName()).Tokenize(TString("/"));
  TObjString* path_tok = (TObjString*)token->At(token->GetEntries()-1);
  TString jname = path_tok->GetString().Data();
  cout << jname.Data() << endl;

  // print ifo infos
  cout << "nIFO : " << nIFO << endl;
  for(int n=0; n<nIFO; n++) pD[n] = NET.getifo(n);
  for(int n=0; n<nIFO; n++) pD[n]->print();

  // restore skymaps 
  NET.setSkyMaps(angle);
  NET.setAntenna();
  NET.setDelay(pD[0]->Name);

/* not needed with fixed network copy constructor  !!!
  // restore network NDM & ifoName lists
  vectorD v; v.clear();
  for(int n=0; n<nIFO; n++) {
    NET.NDM.push_back(v);
    for(int m=0; m<nIFO; m++) NET.NDM[n].push_back(0);
    NET.ifoName.push_back(pD[n]->Name);
  }
  cout << "NDM size " << NET.NDM.size() << endl;
*/
  // there is a issue in network class (no TClass for char*)
  NET.ifoName.clear();
  for(int n=0; n<nIFO; n++) NET.ifoName.push_back(pD[n]->Name);

  // declare netevent
  netevent netburst(nIFO,Psave);

  // print network infos
  double mTau=NET.getDelay("MAX");  // maximum time delay 
  double dTau=NET.getDelay();       // time delay difference 
  cout<<"maximum time delay between detectors: "<<mTau<<endl;
  cout<<"       maximum time delay difference: "<<dTau<<endl;
  cout<<"                    netRHO and netCC: "<<NET.netRHO<<", "<<NET.netCC<<endl;
  cout<<"              regulator delta, gamma: "<<NET.delta <<" " <<NET.gamma<<endl<<endl;

  gSystem->Exec("/bin/date");
  gSystem->Exec("/bin/hostname");

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// data conditioning
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  if(!simulation) nfactor = 1;
  for(int iii=0; iii<nfactor; iii++) {
    double factor=factors[iii];
    cout<<"factor="<<factor<<endl;
    char sfactor[32];
    sprintf(sfactor,"factor_%g",factor);
    // build end_CED file name
    sprintf(end_CED,"%s/ced_%s_%g",ced_dir,jname.ReplaceAll(".root","").Data(),factor);

    for(int i=0; i<nIFO; i++) {
      pTF[i] = (WSeries<double>*)jfile->Get(TString("proc/")+pD[i]->Name);
      if(pTF[i]==NULL || simulation) {  // process raw data
        pTF[i] = (WSeries<double>*)jfile->Get(TString("raw/")+pD[i]->Name);
        if(pTF[i]==NULL) {cout << "Error : data not present, job terminated!!!" << endl;return;}
        pTF[i] = pD[i]->getTFmap();
        *pTF[i] = *(WSeries<double>*)jfile->Get(TString("raw/")+pD[i]->Name);
        if(simulation) {
          wM = *(WSeries<double>*)jfile->Get(TString("raw/")+ifo[i]+TString(":mdc"));
          wM*=factor;
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
        pD[i]->bandPass();                
        cout<<"After "<<pD[i]->Name<<" data conditioning"<<endl; 
        gSystem->Exec("/bin/date"); GetProcInfo();
      } else { // process proc data
        pTF[i] = pD[i]->getTFmap();
        *pTF[i] = *(WSeries<double>*)jfile->Get(TString("proc/")+pD[i]->Name);
      }        
    }

    // set the decomposition level
    for(int i=0; i<nIFO; i++) pTF[i]->Inverse(levelD-l_low);
    NET.setTimeShifts();

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// supercluster analysis
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    int lags=1;
    for(int j=0; j<lags; j++){
      wc.clear();
      TString wcdir = simulation ? TString("sim/")+sfactor : TString("prod/lag_")+j;
      TIter nextkey(((TDirectory*)jfile->Get(wcdir))->GetListOfKeys());
      TKey *key;while(key=(TKey*)nextkey()) wc.append(*(netcluster*)key->ReadObj());
      int m = wc.supercluster('L',NET.e2or,true);
      pwc = NET.getwc(j); pwc->cpf(wc,true);
      cout<<m<<"|"<<pwc->size()<<" ";
    }

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// likelihood
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    TDirectory* cdced = NULL;
    if(cedDump) {cdced = simulation ? cdsim->mkdir("ced") : cdprod->mkdir("ced"); }
  
    for(int i=l_low; i<=l_high; i++) {
      sprintf(tdf00,"%s/Meyer1024wat482_00%s_L%1d.dat",filter_dir,filter,i);
      sprintf(tdf90,"%s/Meyer1024wat482_90%s_L%1d.dat",filter_dir,filter,i);
      NET.setDelayFilters(tdf00,tdf90);

      if(i==l_low) {
        //NET.pOUT=true;
        NET.setDelayIndex();
        NET.setIndexMode(mode);
      }

      cout<<endl;
      cout<<"selected core pixels: "<<NET.likelihood(SEARCH(),Acore)<<" for level "<<i<<"\n";
      cout<<"rejected weak pixels: "<<NET.netcut(NET.netRHO,'r',0,1)<<"\n";  // remove weak glitches
      cout<<"rejected loud pixels: "<<NET.netcut(NET.netCC,'c',0,1)<<"\n";   // remove loud glitches
      cout<<"events in the buffer: "<<NET.events()<<"\n";

      if(cedDump) {
        if(cdced==NULL) {
          cout<<"dump ced into "<<end_CED<<"\n";
          CWB::ced ced(&NET,&netburst,end_CED);
        } else {
          //gBenchmark->Start("ced");
          CWB::ced ced(&NET,&netburst,cdced);
          //gBenchmark->Show("ced");
          //gBenchmark->Reset();
        }
        ced.SetOptions(cedRHO);
        ced.Write(factor);
      }

      if(i<l_high) NET.Forward(1);
    }

    gSystem->Exec("/bin/date"); GetProcInfo();
    cout<<"\nSearch done\n";
    cout<<"reconstructed events: "<<NET.events()<<"\n";

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// save data in root file
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    if(simulation) cdsim->cd(); else cdprod->cd();

    TTree* net_tree = netburst.setTree();

    // build output data txt name  (in the /dev/shm directory)
    char ofile[256];
    gRandom->SetSeed(0);
    int rnID = int(gRandom->Rndm(13)*1.e9);
    UserGroup_t* uinfo = gSystem->GetUserInfo();
    TString uname = uinfo->fUser;
    gSystem->Exec(TString("mkdir -p /dev/shm/")+uname);
    sprintf(ofile,"/dev/shm/%s/waveburst-%d.txt",uname.Data(),rnID);

    dump=true;

    if(dump) netburst.dopen(ofile,"w");
    netburst.output(net_tree,&NET);
    jfile->Write();
    if(dump) netburst.dclose();

    TMacro macro(ofile);
    macro.Write("output.txt"); 
    gSystem->Exec(TString("rm ")+ofile);
  }

  jfile->ReOpen("READ");
  gROOT->RefreshBrowsers(); 
  gROOT->ForceStyle(0);
}


