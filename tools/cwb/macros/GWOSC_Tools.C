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

// this macro is a collection of tools used manage the CWOSC data used for the cwb_gwosc command (see $CWB_SCRIPTS/Makefile.gwosc)
// these functions work for the GWOSC GWTC-1 events

#include <lal/LALSimInspiral.h>

#define MAX_IFOS	3

void ConvertSamples(TString ifname, TString ofname, double gps);
void Posterior2XML(TString gwname, int seed=-1);
void ConvertPSD(TString ifname);

void GWOSC_Tools(TString tool, TString spar1="", TString spar2="", double dpar1=0) {

  if(tool=="ConvertSamples") 	ConvertSamples(spar1, spar2, dpar1);
  if(tool=="Posterior2XML") 	Posterior2XML(spar1,int(dpar1));
  if(tool=="ConvertPSD") 	ConvertPSD(spar1);

}

void ConvertSamples(TString ifname, TString ofname, double gps) {
//
// this macro is used to convert public samples (IMRPhenomPv2) format to the format used by cWB
//
// input:	ifname: name of input posterior sample file, obtained from the public hdf5 file and converted with h5dump command
//			eg: h5dump -o GW150914_GWTC-1.txt -y -w 400 GW150914_GWTC-1.hdf5
//		gps:	gps time used to generate the xml file (the gps time is not included in the public hdf5 file)
// output:	ofname:	name of the converted output file
//

  std::map<TString, double> psample;
  std::vector<TString> pname;

  TString name="";

  name="dummy"; 	  psample[name] = 0; 	pname.push_back(name);		// defined to avoid a bug in mdc::Posterior2XML

  name="time"; 		  psample[name] = gps; 	pname.push_back(name);		// coa time

  name="costheta_jn"; 	  psample[name] = 0; 	pname.push_back(name);		// simTable->inclination
  name="distance"; 	  psample[name] = 0; 	pname.push_back(name);		// luminosity distance
  name="ra"; 		  psample[name] = 0; 	pname.push_back(name);	
  name="dec"; 		  psample[name] = 0; 	pname.push_back(name);
  name="m1"; 		  psample[name] = 0; 	pname.push_back(name);		// mass1 in the detector frame
  name="m2"; 		  psample[name] = 0; 	pname.push_back(name);		// mass2 in the detector frame
  name="a1"; 		  psample[name] = 0; 	pname.push_back(name);		// spin1
  name="a2"; 		  psample[name] = 0; 	pname.push_back(name);		// spin2
  name="costilt1"; 	  psample[name] = 0; 	pname.push_back(name);		// used to compute spin components s1x,s1y,s1z
  name="costilt2"; 	  psample[name] = 0; 	pname.push_back(name);		// used to compute spin components s2x,s2y,s2z

  if(ifname.Contains("GW170817")) {
    name="lambda1"; 	  psample[name] = 0; 	pname.push_back(name);		// lambda1
    name="lambda2"; 	  psample[name] = 0; 	pname.push_back(name);		// lambda2
  }

  name="q"; 		  psample[name] = 0; 	pname.push_back(name);		// mass ratio
  name="mc"; 		  psample[name] = 0; 	pname.push_back(name);		// chirp mass

  name="lal_pnorder"; 	  psample[name] = 8; 	pname.push_back(name);		// pseudoFourPN
  name="lal_amporder"; 	  psample[name] = 0; 	pname.push_back(name);
  name="lal_approximant"; psample[name] = IMRPhenomPv2;	pname.push_back(name);	// IMRPhenomPv2
  name="flow"; 		  psample[name] = 20;	pname.push_back(name);		// simTable->f_lower
  name="f_ref"; 	  psample[name] = 20;	pname.push_back(name);		// simTable->f_final

  name="phase"; 	  psample[name] = 0; 	pname.push_back(name);		// simTable->coa_phase
  name="psi"; 		  psample[name] = 0; 	pname.push_back(name);		// polarization
  name="phi12"; 	  psample[name] = 0; 	pname.push_back(name);		// used to compute spin components s1x,s1y,s1z
  name="phi_jl"; 	  psample[name] = 0; 	pname.push_back(name);		// used to compute spin components s2x,s2y,s2z


  ifstream in;
  in.open(ifname.Data(),ios::in);
  if(!in.good()) {cout << "Error Opening File : " << ifname.Data() << endl;exit(1);}

  ofstream out;
  out.open(ofname.Data(),ios::out);
  if(!out.good()) {cout << "Error Opening File : " << ofname.Data() << endl;exit(1);}
  out.precision(14);

  // write header
  for(int i = 0; i < pname.size(); i++) out << pname[i] << "\t";
  out << endl;

  int  entries=0;
  char line[1024];
  bool data=false;
  while (1) {
    in >> line; if(!in.good()) break;
    if(TString(line).Contains("{")) {
      in >> line; if(!in.good()) break;
      data=true;
    }
    if(data) {
      psample["costheta_jn"]	= TString(line).Atof(); in >> line; if(!in.good()) break;
      psample["distance"]	= TString(line).Atof(); in >> line; if(!in.good()) break;
      psample["ra"]		= TString(line).Atof(); in >> line; if(!in.good()) break;
      psample["dec"]		= TString(line).Atof(); in >> line; if(!in.good()) break;
      psample["m1"]		= TString(line).Atof(); in >> line; if(!in.good()) break;
      psample["m2"]		= TString(line).Atof(); in >> line; if(!in.good()) break;
      if(ifname.Contains("GW170817")) {
        psample["lambda1"]	= TString(line).Atof(); in >> line; if(!in.good()) break;
        psample["lambda2"]	= TString(line).Atof(); in >> line; if(!in.good()) break;
      }
      psample["a1"]		= TString(line).Atof(); in >> line; if(!in.good()) break;
      psample["a2"]		= TString(line).Atof(); in >> line; if(!in.good()) break;
      psample["costilt1"]	= TString(line).Atof(); in >> line; if(!in.good()) break;
      psample["costilt2"]	= TString(line).Atof(); in >> line; if(!in.good()) break;

      double m1 = psample["m1"];
      double m2 = psample["m2"];
      psample["mc"] = pow(m1*m2,3./5.)/pow(m1+m2,1./5.);
      psample["q"]  = m2<m1 ? m2/m1 : m1/m2;

      // write data
      for(int i = 0; i < pname.size(); i++) out << psample[pname[i]] << "\t";
      out << endl;

      entries++; 
/*
      cout << "----------------------------------------------------------" << endl;
      cout << entries << endl;
      cout << "----------------------------------------------------------" << endl;

      for (int i = 0; i < pname.size(); i++) {
        cout << pname[i] << "\t" << psample[pname[i]] << endl;
      }
*/
    }
    if(TString(line).Contains("}")) data=false;

//    if(entries>3) break;
  }

  cout << endl << "ConvertSamples.C - converted entries: " << entries << endl;

  cout << endl << "ConvertSamples.C - created file: " << ofname << endl << endl;

  out.close();
  in.close();

  exit(0);
}

void Posterior2XML(TString gwname, int seed) {
//
// extract one sample from the posterior samples file and convert into the xml format
// the public files are at: https://dcc.ligo.org/LIGO-P1800370/public
// how to download, eg: wget https://dcc.ligo.org/public/0157/P1800370/005/GW150914_GWTC-1.hdf5 -q -O GW150914_GWTC-1.hdf5
//
// input:	gwname: GW name
//		seed:	seed used to randomly select a sample from the gwname+"_GWTC-1.dat" posterior sample file
//			if <0 then the first sample is selected
//

  gRandom->SetSeed(seed);

  TString options = "";

  options += "--ninjections 1 ";
  options += TString::Format("--seed %d ",seed);

  options += "--source "+gwname+" ";
  options += "--waveform IMRPhenomPv2pseudoFourPN ";
  options += "--f_ref 20 ";
  options += "--f_lower 20 ";

  TString posteriorFile  = gwname+"_GWTC-1.dat";

  TString xmlFile        = gwname+"_GWTC-1.xml";
  options += "--clb_file "+gwname+"_GWTC-1.clb ";

  CWB::mdc::Posterior2XML(posteriorFile, xmlFile, options);

  exit(0);
}

void ConvertPSD(TString ifname) {
//
// convert the GWOSC PSD curves to the format used by the the CWB_Plugin_SimNoise.C plugin
// the public files are at: https://dcc.ligo.org/LIGO-P1900011/public
// how to download, eg: wget https://dcc.ligo.org/public/0158/P1900011/001/GWTC1_GW150914_PSDs.dat -q -O GWTC1_GW150914_PSDs.dat 
// 
// input:	ifname:	input public PSD file
//

  ifstream in;
  in.open(ifname.Data(), ios::in);
  if (!in.good()) {cout << "Error Opening File : " << ifname << endl;exit(1);}

  int size=0;
  char str[1024];
  int fpos=0;
  while(true) {
    in.getline(str,1024);
    if (!in.good()) break;
    if(str[0] != '#') size++;
  }
  cout << "size " << size << endl;
  in.clear(ios::goodbit);
  in.seekg(0, ios::beg);

  in.getline(str,1024);	// get header
  TString header = str;

  int iH1 = header.Index("LIGO_Hanford_PSD"); 
  int iL1 = header.Index("LIGO_Livingston_PSD"); 
  int iV1 = header.Index("Virgo_PSD"); 

  double  iIFO[MAX_IFOS];
  TString sIFO[MAX_IFOS];

  int nIFO=0;
  if(iH1>=0) {iIFO[nIFO]=iH1;sIFO[nIFO]="H1";nIFO++;}
  if(iL1>=0) {iIFO[nIFO]=iL1;sIFO[nIFO]="L1";nIFO++;}
  if(iV1>=0) {iIFO[nIFO]=iV1;sIFO[nIFO]="V1";nIFO++;}

  if(nIFO==0 || nIFO>MAX_IFOS) {cout << "Error: no IFOs in PSDs file or not allowed Format" << endl; exit(1);}

  double* freq = new double[size];
  double* psd[MAX_IFOS]; for(int n=0;n<MAX_IFOS;n++) psd[n] = new double[size];

  int lines=0;
  while(1) {
    if(nIFO==1) in >> freq[lines] >> psd[0][lines];
    if(nIFO==2) in >> freq[lines] >> psd[0][lines] >> psd[1][lines];
    if(nIFO==3) in >> freq[lines] >> psd[0][lines] >> psd[1][lines] >> psd[2][lines];
    if (!in.good()) break;
    size=lines;
    lines++;
  }
  in.close();

  int idx[3];
  TMath::Sort(nIFO,iIFO,idx,false);

  for(int n=0;n<nIFO;n++) {
    int m=idx[n];
    TString ofname = ifname;
    ofname.ReplaceAll("_PSDs.dat","_ASD_"+sIFO[m]+".dat");
    cout << ofname << endl;
    ofstream out;
    out.open(ofname.Data(),ios::out);
    if(!out.good()) {cout << "Error Opening Output File : " << ofname << endl;exit(1);}
    out << 0 << "\t" << 10*sqrt(psd[m][0]) << endl;		// blind the detector at low frequencies
    for(int i=0;i<size;i++) out << freq[i] << "\t" << sqrt(psd[m][i]) << endl;
    out << 8192 << "\t" << 10*sqrt(psd[m][size-1]) << endl;	// blind the detector at high frequencies	
    out.close();
  }

  exit(0);
}

