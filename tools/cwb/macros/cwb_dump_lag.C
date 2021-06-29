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


// list lag used in production : used by the cwb_dump command

{
  CWB::Toolbox TB;

  TB.checkFile(gSystem->Getenv("CWB_ROOTLOGON_FILE"));
  TB.checkFile(gSystem->Getenv("CWB_PARAMETERS_FILE"));
  TB.checkFile(gSystem->Getenv("CWB_UPARAMETERS_FILE"));

  if(simulation) {
    cout << "cwb_dump lag : simulation is !=0 -> only zero lag is defined!!!" << endl << endl;
    exit(0);
  }

  // check input user configuration
  CWB::config cfg;
  cfg.Import();
  cfg.Check();

  // write lag list to file
  if(lagMode[0]!='r') {
    lagFile = new char[256];
    sprintf(lagFile,"%s/%s.lag",dump_dir,data_label);
    cout << "Write lag list : " << lagFile << endl;
  }

  // define network
  detector* pD[NIFO_MAX];                   	// pointers to detectors
  for(int i=0; i<nIFO; i++) {
    if(strlen(ifo[i])>0) pD[i] = new detector(ifo[i]);        // built in detector
    else                 pD[i] = new detector(detParms[i]);   // user define detector
  }

  network NET;                          	// network
  for(int i=0; i<nIFO; i++) NET.add(pD[i]);	// add deetctors to network object

  wavearray<double> x(segLen*16384);        
  pD[0]->TFmap=x;
 
  // setup lags
  int lags = NET.setTimeShifts(lagSize,lagStep,lagOff,lagMax,lagFile,lagMode,lagSite);
  cout<<"lag step: "<<lagStep<<endl;
  cout<<"number of time lags: "<<lags<<endl;

/*
  // print selected lags
  printf("%8s ","lag");
  for(int n=0; n<nIFO; n++) printf("%12.12s%2s","ifo",NET.getifo(n)->Name);
  printf("\n");
  for(int m=0; m<NET.getifo(0)->lagShift.size(); m++) {
    printf("%8d ",m);
    for(int n=0; n<nIFO; n++) printf("%14.5f",NET.getifo(n)->lagShift.data[m]);
    printf("\n");
  }
*/
  exit(0);
}
