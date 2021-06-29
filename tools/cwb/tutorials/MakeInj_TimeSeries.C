CWB::mdc *MDC; 

#define APPROX		"EXP02"
#define XML_FILE        "../EXP02/freffix/posterior_samples_S190521g_EXP02_map.xml"
#define GEO_CENTER_TIME	1242442967.473284864

#define BEG_TIME	1242442957
#define END_TIME	1242442977

void MakeInj_TimeSeries() {

  #define N_IFOS 3

  TString ifo[N_IFOS]={"L1","H1","V1"};

  // set inspiral parms

  MDC = new CWB::mdc(N_IFOS,ifo);
  TString inspOptions="";
  inspOptions+="--xml "+TString(XML_FILE)+ " ";
  inspOptions+= "--dir ./ ";
  MDC->SetInspiral("PE",inspOptions);

  // get mdc buffer draw waveform
  wavearray<double> w[N_IFOS];
  for(int i=0;i<N_IFOS;i++) {
    w[i].resize(10*16384);
    w[i].rate(16384);
    w[i].start(GEO_CENTER_TIME-5);

    MDC->Get(w[i],ifo[i]);     // get waveform

    //w[i].start(0);
    //MDC->DrawTime(w[i],"ALP ZOOM");
    //MDC->DrawFFT(w[i],"ALP ZOOM");
    //return

    TString ofName = TString::Format("cWB_S190521g_%s_map_%.3f_%s.txt",APPROX,w[i].start(),ifo[i].Data());
    ofstream out;
    out.open(ofName.Data(),ios::out);
    if (!out.good()) {cout << "Error Opening Output File : " << ofName.Data() << endl;exit(1);}
    cout << "Create Output File : " << ofName.Data() << endl;
    out.precision(19);
    for(int j=0;j<6*w[i].rate();j+=8) {
      double dt = 1./w[i].rate();
      double t = w[i].start()+j*dt;
      out << t << "\t" << w[i][j] << endl;
    }
    out.close();
  }
  exit(0); 
}
