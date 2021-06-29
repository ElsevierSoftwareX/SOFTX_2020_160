//
// Make BRST HF [512:5120] Hz Set for O1 (proposal)
// How to use it :
// root -b Make_BRST_O1_HF.C
//
// 1) TIME & FFT plots are generated for each waveform and saved under the ODIR/plots directory
// 2) A texi file html_index.texi is created under ODIR
// 3) texi file is converted into ODIR/html_index/index.html
//
// Author : Gabriele Vedovato
//

#define ODIR	"brst_o1_hf"	// output plot directory

CWB::mdc* MDC;

void Make_BRST_O1_HF() {

  gROOT->SetBatch(true);

  MDC = new CWB::mdc();

  vector<mdcpar> par;

  // -------------------------------------------
  // SGQ3
  // -------------------------------------------

  par.resize(2);
  double F_Q3[5] = {849,1615,2000,2477,3067};
  for(int n=0;n<5;n++) {
    par[0].name="frequency"; par[0].value=F_Q3[n];
    par[1].name="Q"; par[1].value=3.;
    MDC->AddWaveform(MDC_SGE, par);
  }

  // -------------------------------------------
  // SGQ9
  // -------------------------------------------

  par.resize(2);
  double F_Q9[15] = {849,1053,1304,1415,1615,1797,2000,2226,2477,2756,3067,3413,3799,4225,5000};
  double Q_Q9[15] = {8.9,9,9,9,9,9,9,9,9,9,9,9,9,9,9};
  for(int n=0;n<15;n++) {
    par[0].name="frequency"; par[0].value=F_Q9[n];
    par[1].name="Q"; par[1].value=Q_Q9[n];
    if((F_Q9[n]==1053)||(F_Q9[n]==2000)||(F_Q9[n]==3067)||(F_Q9[n]==5000)) {
      MDC->AddWaveform(MDC_SGL, par);	// linear polarized for back compatibility
    } else {
      MDC->AddWaveform(MDC_SGE, par);
    }
  }

  // -------------------------------------------
  // SGQ100
  // -------------------------------------------

  par.resize(2);
  double F_Q100[6] = {849,1615,2000,2477,3067,5000};
  for(int n=0;n<6;n++) {
    par[0].name="frequency"; par[0].value=F_Q100[n];
    par[1].name="Q"; par[1].value=100.;
    MDC->AddWaveform(MDC_SGE, par);
  }

  // -------------------------------------------
  // WNB
  // -------------------------------------------

  par.resize(6);
  double F_WNB[7] = {1000,2000,3500,1000,1000,2000,3500};
  double B_WNB[7] = {10,100,100,1000,1000,1000,1000};
  double D_WNB[7] = {0.1,0.1,0.1,0.1,0.01,0.01,0.01};
  for(int n=0;n<7;n++) {
    par[0].name="frequency"; par[0].value=F_WNB[n];
    par[1].name="bandwidth"; par[1].value=B_WNB[n];
    par[2].name="duration";  par[2].value=D_WNB[n];
    for(int m=0;m<30;m++) {
      par[3].name="pseed";     par[3].value=100000+n*100+m;
      par[4].name="xseed";     par[4].value=100001+n*100+m;
      par[5].name="mode";      par[5].value=0;  // asymmetric
      MDC->AddWaveform(MDC_WNB, par);
    }
  }

  // -------------------------------------------
  // RDET200
  // -------------------------------------------

  par.resize(2);
  double F_T200[3] = {1590,2090,2590};
  for(int n=0;n<3;n++) {
    par[0].name="frequency"; par[0].value=F_T200[n];
    par[1].name="tau"; par[1].value=0.2;
    MDC->AddWaveform(MDC_RDE, par);
  }

  // -------------------------------------------
  // RDEQ9
  // -------------------------------------------

  char wf_name[128];
  par.resize(2);
  double F2_Q9[3] = {2000,3067,5000};
  for(int n=0;n<3;n++) {
    par[0].name="frequency"; par[0].value=F2_Q9[n];
    par[1].name="tau"; par[1].value=9./(TMath::TwoPi()*F2_Q9[n]);
    sprintf(wf_name,"RDE%.0fQ9",F2_Q9[n]);
    MDC->AddWaveform(MDC_RDE, par, wf_name);
  }

  MDC->Print();

  vector<waveform> wfList = MDC->GetWaveformList();

  // loop over the waveform list
  for(int n=0;n<wfList.size();n++) {

    // Get the first waveform hp,hx components 
    waveform wf = MDC->GetWaveform(wfList[n].name,0);

    if(wf.status==false) {
      cout << "Error : Waveform " << wf.name << " not exist in the MDC pool !!!" << endl << endl;
      gSystem->Exit(1);
    }

    cout << n << " " << wf.name << " size : " << wf.hp.size() << " rate : " << wf.hp.rate() 
         << " start : " << (int)wf.hp.start() << endl;

    wf.hp.start(0);	// set start to 0 (needed by draw Method)
    wf.hx.start(0);

    watplot* plot = NULL;

    gSystem->mkdir(TString::Format("%s",ODIR));
    gSystem->mkdir(TString::Format("%s/plots",ODIR));

    // draw hp,hx in time domain
    plot=MDC->Draw(wf.hp,MDC_TIME,"ALP ZOOM");
    MDC->Draw(wf.hx,MDC_TIME,"same",kRed);
    if(plot) plot->graph[0]->SetTitle(TString::Format("%s : h+(black), hx(red))",wf.name.Data()));
    if(plot) plot->canvas->SaveAs(TString::Format("%s/plots/%s_TIME.png",ODIR,wf.name.Data()));

    // draw hp,hx in frequency domain
    plot = MDC->Draw(wf.hp,MDC_FFT,"ALP NOLOGX");	// draw hp in frequency domain
    MDC->Draw(wf.hx,MDC_FFT,"same NOLOGX",kRed);
    plot->graph[0]->GetHistogram()->GetXaxis()->SetRangeUser(512,5120);
    if(plot) plot->graph[0]->SetTitle(TString::Format("%s : h+(black), hx(red))",wf.name.Data()));
    if(plot) plot->canvas->SaveAs(TString::Format("%s/plots/%s_FFT.png",ODIR,wf.name.Data()));
  }

  TString texiName=TString::Format("%s/html_index.texi",ODIR);
  bool overwrite=CWB::Toolbox::checkFile(texiName,true);
  if(!overwrite) gSystem->Exit(0);
  ofstream out;
  out.open(texiName,ios::out);
  out << "@c include texi macros" << endl;
  out << "@include macros.texi" << endl;
  out << "@include mathjax.texi" << endl << endl;
  for(int i=0;i<wfList.size();i++) {
    out << "@center @txtfont{"<<wfList[i].name<<", red, h1}" << endl;
    out << "@multitable @columnfractions .5 .5" << endl;
    out << "@item @center @displayimage{../plots,"<<TString::Format("%s_TIME",wfList[i].name.Data())<<",480}"<<endl;
    out << "@tab  @center @displayimage{../plots,"<<TString::Format("%s_FFT",wfList[i].name.Data())<<",480}"<<endl;
    out << "@end multitable" << endl;
  }
  out.close();

  // convert texi into html
  TString cwb_scripts = TString(gSystem->Getenv("CWB_SCRIPTS"));
  //TString exec_cmd = TString::Format("%s/cwb_mkhtml.csh %s wheader;rm %s",
  //                   cwb_scripts.Data(),texiName.Data(),texiName.Data());
  TString exec_cmd = TString::Format("%s/cwb_mkhtml.csh %s wheader",
                     cwb_scripts.Data(),texiName.Data());
  int ret=gSystem->Exec(exec_cmd);
  if(ret) {
    cout << "Make_BRST_O1_HF.C : Error while executing cwb_mkhtml html_index.texi !!!" << endl;
    gSystem->Exit(1);
  }
  cout << endl << "New html file created : " << ODIR << "/html_index/index.html" << endl << endl;

  gSystem->Exit(0);
}
