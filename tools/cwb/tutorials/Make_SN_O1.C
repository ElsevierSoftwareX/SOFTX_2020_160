//
// Make BRST LF [32:1024] Hz Set for O1 (proposal)
// How to use it :
// root -b Make_SN_O1.C
//
// 1) TIME & FFT plots are generated for each waveform and saved under the ODIR/plots directory
// 2) A texi file html_index.texi is created under ODIR
// 3) texi file is converted into ODIR/html_index/index.html
//
// Author : Gabriele Vedovato
//

// ------------------------------------------------------------------------------
// Scheidegger2010 : R4E1FC
// https://wiki.ligo.org/viewauth/Bursts/O1Waveforms
// https://svn.ligo.caltech.edu/svn/snsearch/Raw_Waveforms/Scheidegger2010/
// ------------------------------------------------------------------------------
/*
#define SNNAME  "R4E1FC_L"
#define SNDIR   "/home/waveburst/WAVEFORMS/SN/Scheidegger2010/WaveformFiles"
#define ODIR	"sn_o1/R4E1FC_L"	// output plot directory
#define NTH	14
#define NPH	23
*/

// ------------------------------------------------------------------------------
// Ott2013 : s27fheat1p05
// https://wiki.ligo.org/viewauth/Bursts/O1Waveforms
// https://svn.ligo.caltech.edu/svn/snsearch/Raw_Waveforms/Ott2013/WaveformFiles/
// ------------------------------------------------------------------------------

#define SNNAME  "s27fheat1p05"
#define SNDIR   "/home/waveburst/WAVEFORMS/SN/Ott2013/WaveformFiles"
#define ODIR	"sn_o1/s27fheat1p05"	// output plot directory
#define NTH	18
#define NPH	27


#define HEADER	// add cwb header to html 

CWB::mdc* MDC;

void Plot_SN_O1(vector<waveform> wfList, int ID, int id); 
void WriteHtml_SN_O1(ofstream& out, vector<waveform> wfList, int ID, int id, int index);

void Make_SN_O1() {

  gROOT->SetBatch(true);

  MDC = new CWB::mdc();

  // read SN waveform file list
  vector<TString> snList = CWB::Toolbox::getFileListFromDir(SNDIR, "plus.txt",SNNAME,"",true);
  if(snList.size()==0) {
    cout << "No SN files are present on the directory : " << SNDIR << endl;
    return;
  }
  for(int j=0;j<snList.size();j++) {
    //cout << j << " " << snList[j] << endl;
    TString sn_name = SNNAME;
    TString hp_name = snList[j];
    TString hc_name = snList[j]; hc_name.ReplaceAll("plus.txt","cross.txt");
    TString name = gSystem->BaseName(hp_name.Data()); name.ReplaceAll("-plus.txt","");
    TString sth  = name;    sth.Remove(0,NTH);sth.Remove(5,sth.Sizeof()); 
    TString sph  = name;    sph.Remove(0,NPH);sph.Remove(5,sph.Sizeof()); 
    vector<mdcpar> sn_par(4);
    sn_par[0].name="name";  sn_par[0].svalue=name;
    sn_par[1].name="theta"; sn_par[1].value=sth.Atof();
    sn_par[2].name="phi";   sn_par[2].value=sph.Atof();
    sn_par[3].name="hrss";  sn_par[3].value=-1;		// do not normalize hrss
    MDC->AddWaveform(sn_name,hp_name,hc_name,16384,sn_par);
  }

  MDC->Print(1);

  vector<waveform> wfList = MDC->GetWaveformList();

  gSystem->mkdir(TString::Format("%s/plots",ODIR),true);

  // loop over the waveform list
  for(int n=0;n<wfList.size();n++) {
    Plot_SN_O1(wfList,n,0);
    for(int j=0;j<(int)wfList[n].list.size();j++) {
      Plot_SN_O1(wfList,n,j);
    }
  }

  TString texiName=TString::Format("%s/html_index.texi",ODIR);
  bool overwrite=CWB::Toolbox::checkFile(texiName,true);
  if(!overwrite) gSystem->Exit(0);
  int index=1;
  ofstream out;
  out.open(texiName,ios::out);
  out << "@c include texi macros" << endl;
  out << "@include macros.texi" << endl;
  out << "@include mathjax.texi" << endl << endl;
  for(int n=0;n<wfList.size();n++) {
    WriteHtml_SN_O1(out,wfList,n,0,index++);
    for(int j=0;j<(int)wfList[n].list.size();j++) {
      WriteHtml_SN_O1(out,wfList,n,j,index++);
    }
  }
  out.close();

  // convert texi into html
  TString cwb_scripts = TString(gSystem->Getenv("CWB_SCRIPTS"));
  //TString exec_cmd = TString::Format("%s/cwb_mkhtml.csh %s wheader;rm %s",
  //                   cwb_scripts.Data(),texiName.Data(),texiName.Data());
#ifdef HEADER
  TString exec_cmd = TString::Format("%s/cwb_mkhtml.csh %s wheader",
                     cwb_scripts.Data(),texiName.Data());
#else
  TString exec_cmd = TString::Format("%s/cwb_mkhtml.csh %s",
                     cwb_scripts.Data(),texiName.Data());
#endif
  int ret=gSystem->Exec(exec_cmd);
  if(ret) {
    cout << "Make_SN_O1.C : Error while executing cwb_mkhtml html_index.texi !!!" << endl;
    gSystem->Exit(1);
  }
  cout << endl << "New html file created : " << ODIR << "/html_index/index.html" << endl << endl;

  gSystem->Exit(0);
}

void Plot_SN_O1(vector<waveform> wfList, int ID, int id) {

  waveform wf = MDC->GetWaveform(wfList[ID].name,id);
  if(wf.status==false) {
    cout << "Error : Waveform " << wf.name << " not exist in the MDC pool !!!" << endl << endl;
    gSystem->Exit(1);
  }

  cout << ID << " " << wf.name << " size : " << wf.hp.size() << " rate : " << wf.hp.rate() 
       << " start : " << (int)wf.hp.start() << endl;

  wf.hp.start(0);	// set start to 0 (needed by draw Method)
  wf.hx.start(0);

  watplot* plot = NULL;
  TString name = wf.par[0].svalue;

  plot=MDC->Draw(wf.hp,MDC_TIME,"ALP");
  MDC->Draw(wf.hx,MDC_TIME,"same",kRed);
  if(plot) plot->graph[0]->SetTitle(TString::Format("%s : h+(black), hx(red))",name.Data()));
  if(plot) plot->canvas->SaveAs(TString::Format("%s/plots/%s_TIME.png",ODIR,name.Data()));

  plot = MDC->Draw(wf.hp,MDC_FFT,"ALP");	// draw hp in frequency domain
  MDC->Draw(wf.hx,MDC_FFT,"same",kRed);
  plot->graph[0]->GetHistogram()->GetXaxis()->SetRangeUser(0,8*1024);
  if(plot) plot->graph[0]->SetTitle(TString::Format("%s : h+(black), hx(red))",name.Data()));
  if(plot) plot->canvas->SaveAs(TString::Format("%s/plots/%s_FFT.png",ODIR,name.Data()));

}

void WriteHtml_SN_O1(ofstream& out, vector<waveform> wfList, int ID, int id, int index) {

  waveform wf = MDC->GetWaveform(wfList[ID].name,id);
  TString parms = TString::Format("hrss @@ 10kpc : %s = %0.3g - %s = %0.3g", 
                                  wf.par[4].name.Data(),wf.par[4].value,
                                  wf.par[5].name.Data(),wf.par[5].value);
  TString name = wf.par[0].svalue;
  out << "@center @txtfont{"<<TString::Format("%d - %s",index,name.Data()).Data() <<", red, h1}" << endl;
  out << "@center @txtfont{"<<parms<<", blue, h2}" << endl;
  out << "@multitable @columnfractions .5 .5" << endl;
  out << "@item @center @displayimage{../plots,"
      <<TString::Format("%s_TIME",name.Data())<<",480}"<<endl;
  out << "@tab  @center @displayimage{../plots,"
      <<TString::Format("%s_FFT",name.Data())<<",480}"<<endl;
  out << "@end multitable" << endl;
}
