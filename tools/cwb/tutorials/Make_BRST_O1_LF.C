//
// Make BRST LF [32:1024] Hz Set for O1 (proposal)
// How to use it :
// root -b Make_BRST_O1_LF.C
//
// 1) TIME & FFT plots are generated for each waveform and saved under the ODIR/plots directory
// 2) A texi file html_index.texi is created under ODIR
// 3) texi file is converted into ODIR/html_index/index.html
//
// Author : Gabriele Vedovato
//

#define ODIR	"brst_o1_lf"	// output plot directory

CWB::mdc* MDC;

void Make_BRST_O1_LF() {

  gROOT->SetBatch(true);

  MDC = new CWB::mdc();

  vector<mdcpar> par;

  // -------------------------------------------
  // GA
  // -------------------------------------------

  par.resize(1);
  double F_GA[4]   = {0.0001,0.001,0.0025,0.004};
  for(int n=0;n<4;n++) {
    par[0].name="duration"; par[0].value=F_GA[n];
    MDC->AddWaveform(MDC_GA, par);
  }

  // -------------------------------------------
  // SGQ3
  // -------------------------------------------

  par.resize(2);
  double F_Q3[3] = {70,235,849};
  for(int n=0;n<3;n++) {
    par[0].name="frequency"; par[0].value=F_Q3[n];
    par[1].name="Q"; par[1].value=3.;
    MDC->AddWaveform(MDC_SGE, par);
  }

  // -------------------------------------------
  // SGQ9
  // -------------------------------------------

  par.resize(2);
  double F_Q8d9[8] = {70,100,153,235,361,554,849,1053};
  double Q[8] = {8.9,8.9,8.9,8.9,8.9,8.9,8.9,9};
  for(int n=0;n<8;n++) {
    par[0].name="frequency"; par[0].value=F_Q8d9[n];
    par[1].name="Q"; par[1].value=Q[n];
    if((F_Q8d9[n]==70)||(F_Q8d9[n]==100)||(F_Q8d9[n]==235)||(F_Q8d9[n]==361)||(F_Q8d9[n]==1053)) {
      MDC->AddWaveform(MDC_SGL, par);	// linear polarized for back compatibility
    } else {
      MDC->AddWaveform(MDC_SGE, par);
    }
  }

  // -------------------------------------------
  // SGQ100
  // -------------------------------------------

  par.resize(2);
  double F_Q100[3] = {70,235,849};
  for(int n=0;n<3;n++) {
    par[0].name="frequency"; par[0].value=F_Q100[n];
    par[1].name="Q"; par[1].value=100.;
    MDC->AddWaveform(MDC_SGE, par);
  }

  // -------------------------------------------
  // WNB
  // -------------------------------------------

  par.resize(6);
  double F_WNB[2] = {100,250};
  double B_WNB[2] = {100,100};
  double D_WNB[2] = {0.1,0.1};
  for(int n=0;n<2;n++) {
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
    plot = MDC->Draw(wf.hp,MDC_FFT,"ALP NOLOGX");	
    if(wf.name=="GA0d1") MDC->Draw(wf.hx,MDC_FFT,"same NOLOGX NOLOGY",kRed);
    else                 MDC->Draw(wf.hx,MDC_FFT,"same NOLOGX",kRed);
    plot->graph[0]->GetHistogram()->GetXaxis()->SetRangeUser(32,1024);
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
    cout << "Make_BRST_O1_LF.C : Error while executing cwb_mkhtml html_index.texi !!!" << endl;
    gSystem->Exit(1);
  }
  cout << endl << "New html file created : " << ODIR << "/html_index/index.html" << endl << endl;

  gSystem->Exit(0);
}
