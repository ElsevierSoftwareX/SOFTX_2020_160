//
// Merge Segment Lists
// Author : Gabriele Vedovato


#include "/home/vedovato/Y2/coherent/waveburst/wat/wat-5.4.0/network.hh"
#include "TString.h"
#include <fstream>
#include <stdlib.h>
#include <string.h>
#include <vector>
#include "TMath.h"

using namespace std;

struct dqfile {
  int ifo;
  TString file;
  int cat;
  double shift;
  bool invert;
  bool c4;
};

double mergeLists(std::vector<waveSegment>& ilist1, std::vector<waveSegment>& ilist2) {

  int i=0;
  int j=0;
  int n1=ilist1.size();
  int n2=ilist2.size();
  double ctime=0;
  double start=0;
  double stop=0;

  std::vector<waveSegment> olist;
  waveSegment SEG;
  SEG.index = 0;

  while (i<n1 && j<n2) {
    if (ilist2[j].stop<=ilist1[i].start) j++;
    else if (ilist1[i].stop <= ilist2[j].start) i++;
    else {
      if (ilist2[j].start < ilist1[i].start) start=ilist1[i].start;
      else start = ilist2[j].start;
      if (ilist2[j].stop > ilist1[i].stop) stop=ilist1[i].stop;
      else stop=ilist2[j].stop;
      if (ilist2[j].stop >= ilist1[i].stop) i++;
      else j++;
      ctime+=stop-start;

      SEG.index++;
      SEG.start = start;
      SEG.stop  = stop;
      olist.push_back(SEG);
    }
  }

  ilist2=olist;
  olist.clear();

  return ctime;
}

double readSEGlist(dqfile DQF, std::vector<waveSegment>& olist) {

// Open list 

  ifstream in;
  in.open(DQF.file,ios::in);
  if (!in.good()) {cout << "Error Opening File : " << DQF.file << endl;exit(1);}

  int size=0;
  char str[1024];
  int fpos=0;
  while(true) {
    in.getline(str,1024);
    if (!in.good()) break;
    if(str[0] != '#') size++;
  }
  //cout << "size " << size << endl;
  in.clear(ios::goodbit);
  in.seekg(0, ios::beg);

  int N=0;
  int dummy;
  double* start = new double[size];
  double* stop  = new double[size];
  while (1) {
    fpos=in.tellg();
    in.getline(str,1024);
    if(str[0] == '#') continue;
    in.seekg(fpos, ios::beg);

    if(DQF.c4) in >> dummy >> start[N] >> stop[N] >> dummy;
    else    in >> start[N] >> stop[N];
    if (!in.good()) break;
    start[N]+=DQF.shift;
    stop[N]+=DQF.shift;
    if(stop[N]<=start[N]) {
      cout << "Error Ranges : " << start[N] << " " << stop[N] << endl;
      exit(1);
    }
    N++;
    fpos=in.tellg();
    in.seekg(fpos+1, ios::beg);
  }
  in.close();

// sort list

  Int_t *id = new Int_t[N];
  TMath::Sort(N,start,id,false);
  for(int i=1;i<N;i++) {
    bool flag=true;
    if(start[id[i]]<=0) flag=false;
    if(start[id[i]]>stop[id[i]]) flag=false;
    if(start[id[i]]<=stop[id[i-1]]) flag=false;
    if(!flag) {
      cout << "Error in segment list file : " << DQF.file.Data() << endl;
      cout << i << " " << start[id[i]] << " " << stop[id[i]] << " flag " << flag << endl;
      exit(1);
    }
  }

// write list to output vector

  double ctime=0.;
  waveSegment SEG;
  SEG.index = 0;
  olist.clear();
  if(DQF.invert) { // invert bad with good periods
    SEG.index++;
    SEG.start = 0;
    SEG.stop  = start[id[0]];
    olist.push_back(SEG);
    for(int i=0;i<N-1;i++) {
      SEG.index++;
      SEG.start = stop[id[i]];
      SEG.stop  = start[id[i+1]];
      ctime+=SEG.stop-SEG.start;
      olist.push_back(SEG);
    }
    SEG.index++;
    SEG.start = stop[id[N-1]];
    SEG.stop  = 4000000000;
    olist.push_back(SEG);
  } else {
    for(int i=0;i<N;i++) {
      SEG.index++;
      SEG.start = start[id[i]];
      SEG.stop  = stop[id[i]];
      ctime+=stop[id[i]]-start[id[i]];
      olist.push_back(SEG);
    }
  }

  delete [] start;
  delete [] stop;

  return ctime;
}

double readSEGlist(int nDQF, dqfile* DQF, int cat, std::vector<waveSegment>& olist) {

  double ctime=0.;

  int ndqf=0;
  dqfile dqf[nDQF];
  for(int n=0;n<nDQF;n++) {
    if(DQF[n].cat<=cat) dqf[ndqf++]=DQF[n];  
  }

  std::vector<waveSegment>* list = new std::vector<waveSegment>[ndqf];

  for(int n=0;n<ndqf;n++) readSEGlist(dqf[n],list[n]);

  olist=list[0];
  for(int n=1;n<ndqf;n++) ctime = mergeLists(list[n],olist);

  for(int n=0;n<ndqf;n++) list[n].clear();

  return ctime;
}

int dumpList(std::vector<waveSegment> list, TString fName, bool c4=false) {

  FILE* fP;
  if((fP = fopen(fName.Data(), "w")) == NULL) {
    cout << "cannot open output file " << fName.Data() <<". \n";
    exit(1);
  };
  cout << "Write output file : " << fName.Data() << endl;

  if(c4) 
    fprintf(fP,"# seg   start           stop            duration\n");
  for(int i=0;i<list.size();i++) {
    if(c4) {
      fprintf(fP,"%-d\t%-d\t%-d\t%-d\n",
                  int(list[i].index)-1, 
                  int(list[i].start), 
                  int(list[i].stop), 
                  int(list[i].stop-list[i].start));
    } else {
      //fprintf(fP,"%-d\t%-d\n",
      fprintf(fP,"%d %d\n",
                  int(list[i].start), 
                  int(list[i].stop));
    }
  }
  if(fP!=NULL) fclose(fP);
  return 0;
}


#define WRITE_OFILE

void MergeSegments14() {

  int nDQF=12;
  dqfile DQF[12]={
                  {0 ,"S6D_R9_segments/S6D_OFFLINE_L1SCIENCE.txt"        , 0, 0., false, false},
                  {0 ,"S6D_R9_segments/S6D_OFFLINE_L1_DQCAT1SEGMENTS.txt", 1, 0., true , false},
                  {0 ,"S6D_R9_segments/S6D_OFFLINE_L1_DQCAT2SEGMENTS.txt", 2, 0., true , false},
                  {0 ,"S6D_R9_segments/S6D_OFFLINE_L1_DQCAT4SEGMENTS.txt", 1, 0., true , false},
                  {1 ,"S6D_R9_segments/S6D_OFFLINE_H1SCIENCE.txt"        , 0, 0., false, false},
                  {1 ,"S6D_R9_segments/S6D_OFFLINE_H1_DQCAT1SEGMENTS.txt", 1, 0., true , false},
                  {1 ,"S6D_R9_segments/S6D_OFFLINE_H1_DQCAT2SEGMENTS.txt", 2, 0., true , false},
                  {1 ,"S6D_R9_segments/S6D_OFFLINE_H1_DQCAT4SEGMENTS.txt", 1, 0., true , false},
                  {2 ,"S6D_R9_segments/S6D_OFFLINE_V1SCIENCE.txt"        , 0, 0., false, false},
                  {2 ,"S6D_R9_segments/S6D_OFFLINE_V1_DQCAT1SEGMENTS.txt", 1, 0., true , false},
                  {2 ,"S6D_R9_segments/S6D_OFFLINE_V1_DQCAT2SEGMENTS.txt", 2, 0., true , false},
                  {2 ,"S6D_R9_segments/S6D_OFFLINE_V1_DQCAT4SEGMENTS.txt", 1, 0., true , false},
                 };


  int cat = 1;
  std::vector<waveSegment> olist;
  double ctime=readSEGlist(nDQF, DQF, cat, olist);
  cout << "ctime : " << int(ctime) << " sec " << ctime/3600. << " h " << ctime/86400. << " day" << endl;

#ifdef WRITE_OFILE
  dumpList(olist,"olist.txt", false);
#endif

  olist.clear();

  exit(0);
}
