#include <vector>

#define HEN_LIST "config/HEN_list_2009_2010.txt"

#define S6A_START 931035615   // S6A start: Jul 07 2009 21:00:00 UTC
#define S6A_STOP  935798415   // S6A last : Sep 01 2009 00:00:00 UTC

#define ONSOURCE_TIME  1000   // sec


void CreateSegmentHEN_S6A_L1H1V1() {


  CWB::Toolbox TB;

  TB.checkFile(gSystem->Getenv("CWB_PARAMETERS_FILE"));
  gROOT->Macro(gSystem->Getenv("CWB_PARAMETERS_FILE"));  
  TB.checkFile(gSystem->Getenv("CWB_UPARAMETERS_FILE"));
  gROOT->Macro(gSystem->Getenv("CWB_UPARAMETERS_FILE"));  

  nDQF=12;	// exclude the HEN DQ

  double PGPS=0;
  double GPS=0;
  float RA;
  float DEC;
  float RADIUS;

  ifstream in;
  in.open(HEN_LIST,ios::in);
  if (!in.good()) {
    cout << "CreateSegmentHEN_S6A_L1H1V1.C - Error Opening File : " << HEN_LIST << endl;
    gSystem->Exit(1);
  }

  ofstream out;
  char ofile[512];sprintf(ofile,"%s/HEN_S6A_OnSource_Segments.txt",input_dir);
  out.open(ofile,ios::out);
  if (!out.good()) {cout << "CreateSegmentHEN_S6A_L1H1V1.C - Error : Error Opening File : " << ofile << endl;exit(1);}

  char str[1024];
  while(true) {
    in.getline(str,1024);
    if (!in.good()) break;
    if(str[0] == '#') continue;
    TObjArray* token = TString(str).Tokenize(TString(" "));
    if(token->GetEntries()!=16) {
      cout << "CreateSegmentHEN_S6A_L1H1V1.C - bad line format : " << str << endl;
      gSystem->Exit(1);
    }
    TObjString* otoken;
    TString stoken;

    otoken = (TObjString*)token->At(2);
    stoken = otoken->GetString();
    GPS = stoken.Atof();

    otoken = (TObjString*)token->At(3);
    stoken = otoken->GetString();
    RA = stoken.Atof();

    otoken = (TObjString*)token->At(4);
    stoken = otoken->GetString();
    DEC = stoken.Atof();

    otoken = (TObjString*)token->At(5);
    stoken = otoken->GetString();
    RADIUS = stoken.Atof();

    if(GPS<S6A_START || GPS>S6A_STOP) continue;
    if(PGPS && fabs(GPS-PGPS)<=ONSOURCE_TIME) {
      cout << "Trigger discarted because fabs(GPS-PGPS)<=" << ONSOURCE_TIME << endl;
      continue;
    } 
    PGPS=GPS;  

    cout << "CreateSegmentHEN_S6A_L1H1V1.C - trigger found at gps : " 
         << int(GPS) << " ra : " << RA << " dec : " << DEC << " radius : " << RADIUS << endl;

    out << int(GPS-ONSOURCE_TIME/2) << " " << int(GPS+ONSOURCE_TIME/2) << endl;

    delete token;
  }
  in.close();
  out.close();

  // add to DQF the on source trigger list segments
  strcpy(DQF[nDQF].ifo, ifo[0]);
  sprintf(DQF[nDQF].file, "%s",ofile);
  DQF[nDQF].cat    = CWB_CAT0;
  DQF[nDQF].shift  = 0.;
  DQF[nDQF].invert = false;
  DQF[nDQF].c4     = false;
  nDQF++;

  // get the on source trigger list segments after dq cat1
  vector<waveSegment> cat2ListHEN=TB.readSegList(nDQF, DQF, CWB_CAT1);

  // get number/list of available super lag segments
  vector<waveSegment> slagJobList=TB.getSlagJobList(cat2ListHEN, int(segLen));

  // compute cat1ListHEN
  vector<waveSegment> cat1ListHEN;

  int key;
  int J=0;
  for(int i=0;i<cat2ListHEN.size();i++) {
    // discart segment with len < 2*segMLS
    int len = int(cat2ListHEN[i].stop)-int(cat2ListHEN[i].start);
    if(len<2*segMLS) continue;
    //cout << i << " " << int(cat2ListHEN[i].start) << " " << int(cat2ListHEN[i].stop) << endl;
    //cout << int(cat2ListHEN[i].start) << " " << int(cat2ListHEN[i].stop) << endl;
    key = int(cat2ListHEN[i].start);
    for(int j=J;j<slagJobList.size();j++) {
      if(key>slagJobList[j].start && key<slagJobList[j].stop) {
        //cout << j << " " << key << " " << int(slagJobList[j].start) << " " << int(slagJobList[j].stop) << endl;
        cat1ListHEN.push_back(slagJobList[j]);
        J=j+1;
        continue;
      }
    }
    key = int(cat2ListHEN[i].stop);
    for(int j=J;j<slagJobList.size();j++) {
      if(key>slagJobList[j].start && key<slagJobList[j].stop) {
        //cout << j << " " << key << " " << int(slagJobList[j].start) << " " << int(slagJobList[j].stop) << endl;
        cat1ListHEN.push_back(slagJobList[j]);
        J=j+1;
        continue;
      }
    }
  }

  for(int i=0;i<cat2ListHEN.size();i++) {
    //cout << i << " " << int(cat1ListHEN[2*i].start) << " " << int(cat2ListHEN[i].start) << " " 
    //                 << int(cat2ListHEN[i].stop)  << " " << int(cat1ListHEN[2*i+1].stop) << endl;
    //cout << int(cat1ListHEN[i].start) << " " << int(cat1ListHEN[i].stop) << endl;
  }

  // write cat1ListHEN
  ofstream out1;
  char ofile1[512];sprintf(ofile1,"%s/HEN_S6A_L1H1V1_OnSource_Segments_Cat1.txt",input_dir);
  cout << "write file : " << ofile1 << endl;
  out1.open(ofile1,ios::out);
  if (!out1.good()) {cout << "CreateSegmentHEN_S6A_L1H1V1.C - Error : Error Opening File : " << ofile1 << endl;exit(1);}
  for(int i=0;i<cat1ListHEN.size();i++) out1 << int(cat1ListHEN[i].start) << " " << int(cat1ListHEN[i].stop) << endl;
  out1.close();

  // write cat2ListHEN
  ofstream out2;
  char ofile2[512];sprintf(ofile2,"%s/HEN_S6A_L1H1V1_OnSource_Segments_Cat2.txt",input_dir);
  cout << "write file : " << ofile2 << endl;
  out2.open(ofile2,ios::out);
  if (!out2.good()) {cout << "CreateSegmentHEN_S6A_L1H1V1.C - Error : Error Opening File : " << ofile2 << endl;exit(1);}
  for(int i=0;i<cat2ListHEN.size();i++) {
    // discart segment with len < 2*segMLS
    int len = int(cat2ListHEN[i].stop)-int(cat2ListHEN[i].start);
    if(len<2*segMLS) continue;
    out2 << int(cat2ListHEN[i].start) << " " << int(cat2ListHEN[i].stop) << endl;
  }
  out2.close();

  gSystem->Exit(0);
}

