//
// Print Antenna Pattern Components
// Author : Gabriele Vedovato


//#define L1_ENABLED
//#define H1_ENABLED
//#define V1_ENABLED
//#define H2_ENABLED
//#define G1_ENABLED
//#define T1_ENABLED
//#define A1_ENABLED
//#define A2_ENABLED
#define O1_ENABLED
#define N1_ENABLED
//#define E1_ENABLED
//#define J1_ENABLED

void TestPrintAntennaPattern(double inj_theta, double inj_phi, double inj_psi=0., bool rad=false) {

  int nIFO=0;
  TString ifo[12];
  TString ifoName[12];
#ifdef H1_ENABLED
   ifo[nIFO]="H1";   // LHO1
   ifoName[nIFO++]="LHO1";   // LHO1
#endif
#ifdef L1_ENABLED
   ifo[nIFO]="L1";   // LLO
   ifoName[nIFO++]="LLO";   // LLO
#endif
#ifdef G1_ENABLED
   ifo[nIFO]="G1";  // GEO
   ifoName[nIFO++]="GEO";  // GEO
#endif
#ifdef V1_ENABLED
   ifo[nIFO]="V1";   // VIRGO
   ifoName[nIFO++]="VIRGO";   // VIRGO
#endif
#ifdef T1_ENABLED
   ifo[nIFO]="T1";  // TAMA
   ifoName[nIFO++]="TAMA";  // TAMA
#endif
#ifdef H2_ENABLED
   ifo[nIFO]="H2";   // LHO2
   ifoName[nIFO++]="LHO2";   // LHO2
#endif
#ifdef J1_ENABLED
   ifo[nIFO]="J1";   // LCGT
   ifoName[nIFO++]="LCGT";   // LCGT
#endif
#ifdef A1_ENABLED
   ifo[nIFO]="A1";   // AIGO
   ifoName[nIFO++]="AIGO1";   // AIGO
#endif
#ifdef A2_ENABLED
   ifo[nIFO]="A2";   // AIGO
   ifoName[nIFO++]="AIGO2";   // AIGO
#endif
#ifdef O1_ENABLED
   ifo[nIFO]="O1";   // AURIGA
   ifoName[nIFO++]="AURIGA";   // AURIGA
#endif
#ifdef N1_ENABLED
   ifo[nIFO]="N1";   // NAUTILUS
   ifoName[nIFO++]="NAUTILUS";   // NAUTILUS
#endif
#ifdef E1_ENABLED
   ifo[nIFO]="E1";   // EXPLORER
   ifoName[nIFO++]="EXPLORER";   // EXPLORER
#endif

  char ifostr[32]="";
  for(int n=0; n<nIFO; n++) {
    sprintf(ifostr,"%s %s",ifostr,ifo[n].Data());
  }

  char ifoNamestr[32]="";
  sprintf(ifoNamestr,"%s",ifoName[0].Data());
  for(int n=1; n<nIFO; n++) {
    sprintf(ifoNamestr,"%s %s",ifoNamestr,ifoName[n].Data());
  }

  cout << "NETWORK : " << ifoNamestr << endl;

  gnetwork gNET(nIFO,ifo);
  //for(int n=0;n<nIFO;n++) gNET.GetSite(ifo[n]);

  double Pi=TMath::Pi();
  for(int n=0;n<nIFO;n++) {

    double theta = inj_theta;
    double phi = inj_phi;
    double psi = inj_psi;

    if(rad) {
      theta = acos(inj_theta);
      theta*= 180/Pi;

      phi = inj_phi > 0 ? inj_phi : 2*Pi+inj_phi;
      phi*= 180/Pi;

      psi = inj_psi;
      psi*= 180/Pi;
    }

    cout << "phi : " << phi << " theta : " << theta << " psi : " << psi << endl;

    double Fp = gNET.GetAntennaPattern(ifo[n],phi,theta,psi,true);
    double Fc = gNET.GetAntennaPattern(ifo[n],phi,theta,psi,false);

    cout << "IFO : " << ifo[n] << " Fp " << Fp << " Fc : " << Fc << endl;

    cout << "Delay : " << ifo[n] << "-" << ifo[0] << " " 
                       << gNET.GetDelay(ifo[n],ifo[0],phi,theta) << endl; 
  }

  double Fp_dpf = gNET.GetAntennaPattern(inj_phi,inj_theta,inj_psi,true);
  double Fc_dpf = gNET.GetAntennaPattern(inj_phi,inj_theta,inj_psi,false);

  cout << "Network : " << " Fp_dpf " << Fp_dpf << " Fc_dpf : " << Fc_dpf << " Fp_dpf/Fc_dpf : " << Fp_dpf/Fc_dpf << endl;

  exit(0);
}

