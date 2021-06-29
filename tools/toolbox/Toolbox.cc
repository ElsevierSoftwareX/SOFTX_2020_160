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

#include "Toolbox.hh"
#include "watversion.hh"
#include "TTreeFormula.h"
#include "GToolbox.hh"
#include "TThread.h"
#include <Riostream.h>

#define CCAST(PAR) const_cast<char*>(PAR)                                                        \

// Without this macro the THtml doc for CWB::Toolbox can not be generated
ClassImp(CWB::Toolbox)

int compareSegments(waveSegment a, waveSegment b) {return a.start < b.start;}

// definitions used by CWB::Toolbox::mergeCWBTrees with threads

#define MAX_THREADS 16

#define MAX_TREE_SIZE	100000000000LL

struct MergeParms {
  int threadID;
  vector<TString> fileList;
  bool simulation;
  TString odir;
  TString label;
  bool brms; 
  bool bvar; 
  bool bpsm;
};

void *MergeHandle(void *ptr) {

   MergeParms mp = *((MergeParms*) ptr);

   // dump thread file list
   char ofName[1024];
   sprintf(ofName,"%s/mergeList.T%d.txt",mp.odir.Data(),mp.threadID);

   CWB::Toolbox::dumpFileList(mp.fileList, ofName);

   // exec macro merge
   char cmd[1024];
   sprintf(cmd,"root -l -b '%s/cwb_merge_thread.C(\"%s\",%d,\"%s\",\"%s\",%d,%d,%d)'",
           gSystem->ExpandPathName("$CWB_MACROS"), ofName, mp.simulation, 
           mp.odir.Data(), mp.label.Data(), mp.brms, mp.bvar, mp.bpsm);
   gSystem->Exec(cmd);

   // remove thread file list
   gSystem->Exec(TString::Format("rm %s",ofName));

   return 0;
}

//______________________________________________________________________________
vector<waveSegment> 
CWB::Toolbox::readSegments(TString ifile) {
//
// read segment list from file  
//  
// Input:  ifile - input file name
//
// Output: return the waveSegment list   
//

  // Open file segment list                                        
  ifstream in;                                                     
  in.open(ifile.Data(),ios::in);                                   
  if(!in.good()) {
    cout << "CWB::Toolbox::readSegments - Error Opening File : " << ifile << endl;
    gSystem->Exit(1);
  }

  char str[1024];
  int fpos=0;
  int index=0;
  double start;
  double stop;
  waveSegment seg;
  vector<waveSegment> iseg;
  while (1) {
    fpos=in.tellg();
    in.getline(str,1024);
    if(str[0] == '#') continue;
    in.seekg(fpos, ios::beg);

    in >> start >> stop;
    if (!in.good()) break;
    fpos=in.tellg();
    in.seekg(fpos+1, ios::beg);

    seg.index=index++; seg.start=start; seg.stop=stop; iseg.push_back(seg);
  }

  in.close();

  return iseg;
}

//______________________________________________________________________________
vector<waveSegment> 
CWB::Toolbox::unionSegments(vector<waveSegment>& ilist) {
//
// Join & sort a waveSegment list 
//  
// Input:  ilist - waveSegment list
//
// Output: return the joined waveSegment list   
//
// ilist
// xxxxxxxx         xxxxx
//      xxxxxxxx  xxx      xxx
//
// output list
// xxxxxxxxxxxxx  xxxxxxx  xxx
//

  vector<waveSegment> vsegs;
  int n = ilist.size();
  if(n==0) return vsegs;

  waveSegment* isegs = new waveSegment[n];
  waveSegment* osegs = new waveSegment[n+1];

  for(int i=0;i<n+1;i++) {osegs[i].start=0;osegs[i].stop=0;}
  for(int i=0;i<n;i++) {
    if(ilist[i].start>ilist[i].stop) {
      cout << "CWB::Toolbox::unionSegments - Error : start must be <= stop - start = " 
           << ilist[i].start << " stop = " << ilist[i].stop << endl;
      exit(1);
    }
    if(ilist[i].start<0) {
      cout << "CWB::Toolbox::unionSegments - Error : start must be positive - start = " 
           << ilist[i].start << endl;
      exit(1);
    }
    isegs[i] = ilist[i];
  }

  std::sort(isegs, isegs + n, compareSegments);
  double right_most = -1;
  int cnt = 0;
  for (int i = 0 ; i < n ; i++) {
    if (isegs[i].start > right_most) {
      right_most = isegs[i].stop;
      ++cnt;
      osegs[cnt].start = isegs[i].start;
      osegs[cnt].stop  = isegs[i].stop;
    }
    if (isegs[i].stop > right_most) {
      right_most      = isegs[i].stop;
      osegs[cnt].stop = isegs[i].stop;
    }
  }

  int vcnt=0;
  waveSegment seg;
  for(int i=0;i<cnt+1;i++) if(osegs[i].stop>0) {seg=osegs[i];seg.index=vcnt++;vsegs.push_back(seg);}

  delete [] isegs;
  delete [] osegs;  

  return vsegs;
}

//______________________________________________________________________________
vector<waveSegment> 
CWB::Toolbox::sortSegments(vector<waveSegment>& ilist) {
//
// sort a waveSegment list 
//  
// Input:  ilist - waveSegment list
//
// Output: return the sorted waveSegment list   
//

  int size = ilist.size();

  double* start = new double[size];
  double* stop  = new double[size];
  for(int i=0;i<size;i++) {
    start[i] = ilist[i].start;
    stop[i]  = ilist[i].stop;
  }

// sort list

  Int_t *id = new Int_t[size];
  TMath::Sort(size,start,id,false);
  for(int i=1;i<size;i++) {
    bool flag=true;
    if(start[id[i]]<=0) flag=false;
    if(start[id[i]]>stop[id[i]]) flag=false;
    if(start[id[i]]<stop[id[i-1]]) flag=false;
    if(!flag) {
      cout.precision(14);
      cout << "CWB::Toolbox::invertSegments - Error in segment list (duplicated veto) : " << endl;
      cout << id[i-1]+1 << " " << start[id[i-1]] << " " << stop[id[i-1]] << " flag " << flag << endl;
      cout << id[i]+1 << " " << start[id[i]] << " " << stop[id[i]] << " flag " << flag << endl;
      gSystem->Exit(1);
    }
  }

// write list to output vector

  vector<waveSegment> olist;
  double ctime=0.;
  waveSegment SEG;
  SEG.index = 0;
  olist.clear();

  // sort list
  for(int i=0;i<size;i++) {
    SEG.index++;
    SEG.start = start[id[i]];
    SEG.stop  = stop[id[i]];
    ctime+=stop[id[i]]-start[id[i]];
    olist.push_back(SEG);
  }

  delete [] id;
  delete [] start;
  delete [] stop;

  return olist;
}

//______________________________________________________________________________
vector<waveSegment> 
CWB::Toolbox::invertSegments(vector<waveSegment>& ilist) {
//
// invert a waveSegment list 
//  
// Input:  ilist - waveSegment list
//
// Output: return the inverted waveSegment list   
//
// ilist
// xxxxxxxx         xxxxx
//
// output list
//         xxxxxxxxx     xxxxxxxxxxxxxxxxxxxxxxx
//

  int size = ilist.size();

  double* start = new double[size];
  double* stop  = new double[size];
  for(int i=0;i<size;i++) {
    start[i] = ilist[i].start;
    stop[i]  = ilist[i].stop;
  }

// sort list

  Int_t *id = new Int_t[size];
  TMath::Sort(size,start,id,false);
  for(int i=1;i<size;i++) {
    bool flag=true;
    if(start[id[i]]<=0) flag=false;
    if(start[id[i]]>stop[id[i]]) flag=false;
    if(start[id[i]]<stop[id[i-1]]) flag=false;
    if(!flag) {
      cout.precision(14);
      cout << "CWB::Toolbox::invertSegments - Error in segment list (duplicated veto) : " << endl;
      cout << id[i-1]+1 << " " << start[id[i-1]] << " " << stop[id[i-1]] << " flag " << flag << endl;
      cout << id[i]+1 << " " << start[id[i]] << " " << stop[id[i]] << " flag " << flag << endl;
      gSystem->Exit(1);
    }
  }

// write list to output vector

  vector<waveSegment> olist;
  double ctime=0.;
  waveSegment SEG;
  SEG.index = 0;
  olist.clear();

  // invert bad with good periods
  SEG.index++;
  SEG.start = 0;
  SEG.stop  = start[id[0]];
  olist.push_back(SEG);
  for(int i=0;i<size-1;i++) {
    SEG.index++;
    SEG.start = stop[id[i]];
    SEG.stop  = start[id[i+1]];
    ctime+=SEG.stop-SEG.start;
    olist.push_back(SEG);
  }
  SEG.index++;
  SEG.start = stop[id[size-1]];
  SEG.stop  = 4000000000;
  olist.push_back(SEG);

  delete [] id;
  delete [] start;
  delete [] stop;

  return olist;
}

//______________________________________________________________________________
void 
CWB::Toolbox::blackmanharris (double* window, int n) {
//
// blackmanharris window
//
// Input/Output - window      - array of double intialized (must be initialized)
//                              return the window values
// Input        - n           - window length
//

  for (int i = 0; i < n; i++)
  {
    double f = 0.0;
    f = ((double) i) / ((double) (n - 1));
    window[i] = 0.35875 -
      0.48829 * cos(2.0 * TMath::Pi() * f) +
      0.14128 * cos(4.0 * TMath::Pi() * f) -
      0.01168 * cos(6.0 * TMath::Pi() * f);
  }

  double norm = 0;
  for (int i=0;i<n;i++) norm += pow(window[i],2);
  norm /= n;
  for (int i=0;i<n;i++) window[i] /= sqrt(norm);
}

//______________________________________________________________________________
vector<waveSegment> 
CWB::Toolbox::mergeSegLists(vector<waveSegment>& ilist1, vector<waveSegment>& ilist2) {
//
// Merge 2 segment lists 
// NOTE : the input lists must be sorted !!!
//
// Input:  ilist1 - first waveSegment list
//         ilist2 - second waveSegment list
//
// Output: return the merged waveSegment list   
//
// ilist1
// xxxxxxxx    xxxxxxxxxxxx  
//
// ilist2
//      xxxx  xxxxxx           xxx
//
// olist
//      xxx    xxxxx           
//

  int i=0;
  int j=0;
  int n1=ilist1.size();
  int n2=ilist2.size();
  double ctime=0;
  double start=0;
  double stop=0;

  vector<waveSegment> olist;
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

  //ilist2=olist;
  //olist.clear();

  return olist;
}

//______________________________________________________________________________
vector<waveSegment> 
CWB::Toolbox::readSegList(dqfile DQF) {
//
// Read DQ file 
// Format A : #entry  start stop stop-start
// Format B : start stop
// list is sorted, check if start<stop, check if entries are overlapped  
//
// Input:  DQF   - DQ structure
//
// Output: return the list of time ranges 
//

  vector<waveSegment> olist;

// Open list 

  ifstream in;
  in.open(DQF.file,ios::in);
  if (!in.good()) {cout << "CWB::Toolbox::readSegList - Error Opening File : " << DQF.file << endl;gSystem->Exit(1);}

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
  if (size==0) {cout << "CWB::Toolbox::readSegList - Error : File " << DQF.file << " is empty" << endl;gSystem->Exit(1);}

  int N=0;
  float dummy;
  double* start = new double[size];
  double* stop  = new double[size];
  while (1) {
    fpos=in.tellg();
    in.getline(str,1024);
    if(str[0] == '#') continue;
    in.seekg(fpos, ios::beg);

    if(DQF.c4) in >> dummy >> start[N] >> stop[N] >> dummy;
    else       in >> start[N] >> stop[N];
    if (!in.good()) break;
    fpos=in.tellg();
    in.seekg(fpos+1, ios::beg);
    start[N]+=DQF.shift;
    stop[N]+=DQF.shift;
    if(stop[N]<=start[N]) {
      cout.precision(14);
      cout << "CWB::Toolbox::readSegList - Error Ranges : " << start[N] << " " << stop[N] << endl;
      gSystem->Exit(1);
    }
    N++;
    if(N>size) {
      cout << "CWB::Toolbox::readSegList - Error Max Range : " << N << " " << size << endl;
      gSystem->Exit(1);
    }
  }
  in.close();

// sort list

  Int_t *id = new Int_t[N];
  TMath::Sort(N,start,id,false);
  for(int i=1;i<N;i++) {
    bool flag=true;
    if(start[id[i]]<=0) flag=false;
    if(start[id[i]]>stop[id[i]]) flag=false;
    if(start[id[i]]<stop[id[i-1]]) flag=false;
    if(!flag) {
      cout.precision(14);
      cout << "CWB::Toolbox::readSegList - Error in segment list file (duplicated veto) : " << DQF.file << endl;
      cout << id[i-1]+1 << " " << start[id[i-1]] << " " << stop[id[i-1]] << " flag " << flag << endl;
      cout << id[i]+1 << " " << start[id[i]] << " " << stop[id[i]] << " flag " << flag << endl;
      gSystem->Exit(1);
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

  delete [] id;
  delete [] start;
  delete [] stop;

  return olist;
}

//______________________________________________________________________________
vector<waveSegment> 
CWB::Toolbox::readSegList(int nDQF, dqfile* DQF, CWB_CAT dqcat) {
//
// Read & Merge DQ files with DQF.cat<=dqcat
//
//
// Input:  nDQF   - number of DQ files
//         dqfile - DQ structure array
//         dqcat  - max DQ cat 
//
// Output: return the list of time ranges which pass the merged max DQ cat conditions
// 

  vector<waveSegment> olist;
 
  int ndqf=0;
  dqfile dqf[nDQF];
  for(int n=0;n<nDQF;n++) {					// Read DQ files with DQF.cat<=dqcat
    if(DQF[n].cat<=dqcat) dqf[ndqf++]=DQF[n];  
  }
  if(ndqf==0) {
    cout << "CWB::Toolbox::readSegList - no CWB_CAT=" << dqcat << " files in the list" << endl;
    gSystem->Exit(1);
  };

  vector<waveSegment>* list = new vector<waveSegment>[ndqf];

  for(int n=0;n<ndqf;n++) list[n]=readSegList(dqf[n]);

  olist=list[0];
  for(int n=1;n<ndqf;n++) olist = mergeSegLists(list[n],olist); // Merge DQ files 

  for(int n=0;n<ndqf;n++) list[n].clear();

  delete [] list;

  return olist;							// return the list of time ranges
}

//______________________________________________________________________________
int 
CWB::Toolbox::dumpSegList(vector<waveSegment> list, TString fName, bool c4) {
//
// Dump to file a waveSegment list 
//
// Input:  list  - waveSegment list
//         fName - output file name
//         c4    - output format
//                 false : start stop
//                 true  : index start stop (stop-start)
//

  FILE* fP;
  if((fP = fopen(fName.Data(), "w")) == NULL) {
    cout << "cannot open output file " << fName.Data() <<". \n";
    gSystem->Exit(1);
  };
  cout << "Write output file : " << fName.Data() << endl;

  if(c4) 
    fprintf(fP,"# seg   start           stop            duration\n");
  for(int i=0;i<(int)list.size();i++) {
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

//______________________________________________________________________________
double 
CWB::Toolbox::getTimeSegList(vector<waveSegment> list) {
//
// Input: list - segment list
//
// Return the length in sec of the sum of the segments
//

  double segtime=0;
  for(int i=0;i<(int)list.size();i++) segtime+=list[i].stop-list[i].start;
  return segtime;
}

//______________________________________________________________________________
void
CWB::Toolbox::dumpJobList(vector<waveSegment> ilist, TString fName, double segLen, double segMLS, double segEdge) {
//
// dump to file the job segment list
//
// Input: ilist      - waveSegment list
//        fName      - output file name
//        segLen     - Segment length [sec]
//        segMLS     - Minimum Segment Length after DQ_CAT1 [sec]
//        segEdge    - wavelet boundary offset [sec]
//

  vector<waveSegment> olist;
  olist=getJobList(ilist, segLen, segMLS, segEdge);
  dumpSegList(olist, fName, false);
  olist.clear();

  return;
}

//______________________________________________________________________________
vector<waveSegment>
CWB::Toolbox::getJobList(vector<waveSegment> ilist, double segLen, double segMLS, double segEdge) {
//
// build the job segment list
//
// the job segments are builded starting from the input ilist
// each segment must have a minimum length of segMLS+2segEdge and a maximum length of segLen+2*segEdge
// in order to maximize the input live time each segment with lenght<2*(segLen+segEdge) is 
// divided in 2 segments with length<segLen+2*segEdge
// segEdge     : xxx 
// segMLS      : -------
// segLen      : ---------------
// input seg   : ----------------------------
// output segA : xxx---------xxx
// output segB :             xxx----------xxx 
//
// Input: ilist      - number of detectors
//        segLen     - Segment length [sec]
//        segMLS     - Minimum Segment Length after DQ_CAT1 [sec]
//        segEdge    - wavelet boundary offset [sec]
//
// Output: Return the job segment list
//

  if(segMLS>segLen) {
    cout << endl << "CWB::Toolbox::getJobList - Error : segMLS must be <= segLen" << endl;
    cout << "segMLS = " << segMLS << " - segLen = " << segLen << endl << endl;
    exit(1);
  }

  int start,stop,len;
  int remainder,half;
  int n;
  int size = ilist.size();
  int lostlivetime=0;
  vector<waveSegment> olist;

  waveSegment SEG;
  SEG.index=0;
  olist.clear();

  for(int i=0;i<size;i++) {
    start = fmod(ilist[i].start,1)!=0 ? int(ilist[i].start+1) : ilist[i].start; // 1.x -> 2.0
    stop  = fmod(ilist[i].stop ,1)!=0 ? int(ilist[i].stop)    : ilist[i].stop;  // 1.x -> 1.0
    start += segEdge;
    stop  -= segEdge;
    len=stop-start;
    if(len<=0) continue;

    n=len/segLen;
    if(n==0) {
      if(len<segMLS) {
        //printf("too small a segment = %d %d %d %d\n",i,start-8,stop+8,len);
        lostlivetime+=len;
        continue;
      }
      SEG.index++;
      SEG.start=start;
      SEG.stop=stop;
      olist.push_back(SEG);
      continue;
    }
    if(n==1) {
      if(len>segLen) {
        remainder=len;
        half=int(remainder/2);
        if(half>=segMLS) {
          SEG.index++;
          SEG.start=start;
          SEG.stop=SEG.start+half;
          olist.push_back(SEG);
          SEG.index++;
          SEG.start=SEG.stop;
          SEG.stop=stop;
          olist.push_back(SEG);
        } else {
          SEG.index++;
          SEG.start=start;
          SEG.stop=SEG.start+segLen;
          olist.push_back(SEG);
        }
      } else {
        SEG.index++;
        SEG.start=start;
        SEG.stop=stop;
        olist.push_back(SEG);
      }
      continue;
    }

    for(int j=0;j<n-1;j++) {
      SEG.index++;
      SEG.start=segLen*j+start;
      SEG.stop=SEG.start+segLen;
      olist.push_back(SEG);
    }
    remainder=stop-SEG.stop;
    half=int(remainder/2);
    if(half>=segMLS) {
      SEG.index++;
      SEG.start=SEG.stop;
      SEG.stop=SEG.start+half;
      olist.push_back(SEG);
      SEG.index++;
      SEG.start=SEG.stop;
      SEG.stop=stop;
      olist.push_back(SEG);
    } else {
      SEG.index++;
      SEG.start=SEG.stop;
      SEG.stop=SEG.start+segLen;
      olist.push_back(SEG);
    }
  }

  printf("Toolbox::getJobList : lost livetime after building of the standard job list = %d sec\n",lostlivetime);
  return olist;
}

//______________________________________________________________________________
vector<waveSegment>
CWB::Toolbox::getJobList(vector<waveSegment> cat1List, vector<waveSegment> cat2List, 
                         double segLen, double segMLS, double segTHR, double segEdge) {
//
// build the job segment list
//
// 1) build job list starting from cat1List
// 2) select jobs which have lenght < segTHR after cat2List
//
// Input: ilist      - number of detectors
//        segLen     - Segment length [sec]
//        segMLS     - Minimum Segment Length after DQ_CAT1 [sec]
//        segTHR     - Minimum Segment Length after DQ_CAT2 [sec]
//        segEdge    - wavelet boundary offset [sec]
//
// Output: Return the job segment list
//

  // build job list starting from cat1
  vector<waveSegment> jobList1=getJobList(cat1List, segLen, segMLS, segEdge);

  // select jobs which have after cat2 lenght < segTHR
  vector<waveSegment> jobList2;
  for(int i=0;i<jobList1.size();i++) {
    // check if job segment after cat2 is < segTHR
    vector<waveSegment> detSegs_dq2;
    detSegs_dq2.push_back(jobList1[i]);
    detSegs_dq2 = mergeSegLists(detSegs_dq2,cat2List);
    double detSegs_ctime = getTimeSegList(detSegs_dq2);
    if(detSegs_ctime<segTHR) {
      cout << "CWB::Toolbox::getJobList : job segment=" << i+1 
           << " live time after cat2 : " << detSegs_ctime << endl;
      cout << "                           segTHR=" << segTHR 
           << " sec -> job removed !!!" << endl;
    } else {
      jobList2.push_back(jobList1[i]);
    }
  }

  return jobList2;
}

//______________________________________________________________________________
vector<slag> 
CWB::Toolbox::getSlagList(size_t  nIFO,   size_t slagSize, int slagSegs, int slagOff, 
                          size_t slagMin, size_t slagMax, size_t* slagSite, char* slagFile) {
//
// Get SuperLags list 
//
//
// Input: nIFO     - number of ifos
//        slagSize - number of super lags (=1 if simulation>0) - if slagSize=0 -> Standard Segments
//        slagSegs - number of segments
//        slagOff  - first slag id (slagOff=0 - include zero slag )
//        slagMin  - select the minimum available slag distance (see example) : slagMin must be <= slagMax
//        slagMax  - select the maximum available slag distance (see example)
//        slagSite - site index starting with 0 (same definition as lagSite)
//        slagFile - user slag file list (default = NULL)
//                   format : slagID   slagID_ifo1   slagID_ifo2   ...
//                   lines starting with # are skipped
//                   Example with 3 ifos :
//
//                   # SLAG         ifo[0]        ifo[1]        ifo[2]
//                        0             0             0             0
//                        1             0             1            -1
//                        2             0            -1             1
//
// Return the slag list
//
// Built-in slags
// The built-in algorithm generates the slag list according the increasing value of the slag distance
// The entries with the same distance are randomly reshuffled
// The slag distance is defined as the sum over the ifos of the absolute values of the shifts 
//  
// Example with 3 ifos :
//
//          SLAG         ifo[0]        ifo[1]        ifo[2]
//             0             0             0             0
//             1             0             1            -1
//             2             0            -1             1
//             3             0             2             1
//             4             0             2            -1
//             5             0            -1             2
//             6             0            -2            -1
//             7             0             1             2
//             8             0             1            -2
//             9             0            -1            -2
//
// the distance is :
// 0 for SLAG 0
// 2 for SLAG 1,2
// 3 for SLAG 3,4,5,6,7,8,9
//
// slagMin,slagMax select the minimum/maximum available distance
// ex: if slagMin=3,slagMax=3 the available slags are 3,4,5,6,7,8,9
// slagOff select the slag offset of the available slags
// ex: if slagMin=3,slagMax=3,slagOff=2 the available slags are 5,6,7,8,9
// slagSize define the number of selected slags
// ex: if slagMin=3,slagMax=3,slagOff=2,slagSize=2 then the selected slags are 5,6
//  

  if(slagSize<1) slagSize=1;

  if((int)slagMax>slagSegs) {
    cout << "CWB::Toolbox::makeSlagList : Error -  slagMax must be < slagSegs" << endl;
    gSystem->Exit(1);
  }

  if(slagSegs<=0.) {
    cout << "CWB::Toolbox::makeSlagList : Error -  slagSegs must be positive" << endl;
    gSystem->Exit(1);
  }

  if((slagMax>0)&&(slagMin>slagMax)) {
    cout << "CWB::Toolbox::getSlagList : Error - slagMin must be < slagMax" << endl;
    gSystem->Exit(-1);
  }

  if(slagSite!=NULL) for(int n=0; n<(int)nIFO; n++) {
    if(slagSite[n] >= nIFO) {
      cout << "CWB::Toolbox::getSlagList : Error slagSite - value out of range " << endl;
      gSystem->Exit(-1);
    }
  }

  vector<slag> slagList;
  vector<int> id;

  if(slagFile) { 	// read list of slags from file

    cout << endl << "CWB::Toolbox::getSlagList : read slags from file -> " << slagFile << endl;

    ifstream in;
    in.open(slagFile, ios::in);
    if(!in.good()) {
      cout << "CWB::Toolbox::getSlagList : Error Opening File : " << slagFile << endl;
      gSystem->Exit(1);
    }

    char str[1024];
    int fpos=0;
    int id;
    int size = slagOff+slagSize;
    while(true) {
      fpos=in.tellg();
      in.getline(str,1024);
      if(str[0] == '#') continue;
      in.seekg(fpos, ios::beg);
      slag SLAG;
      in >> id; 
      if (!in.good()) break;
      // check duplicated slag id
      for(int j=0; j<(int)slagList.size(); j++) {
        if(id==slagList[j].jobId) {
          cout << "CWB::Toolbox::getSlagList : duplicated slag id " 
               << " or wrong format in slag file -> " << slagFile << endl;
          gSystem->Exit(1);
        }
      }
      SLAG.jobId=id;
      for(int n=0; n<(int)nIFO; n++) {
        in >> id; SLAG.slagId.push_back(id);
      } 
      slagList.push_back(SLAG);
      if (slagList.size()==size) break;
      if (!in.good()) break;
      fpos=in.tellg();
      in.seekg(fpos+1, ios::beg);
    }

    in.close();

    if(size>slagList.size()) {
      cout << "CWB::Toolbox::getSlagList : Error - slagOff+slagSize " << size << " is > entries in slagFile " << slagList.size() << endl;
      gSystem->Exit(1);
    }

  } else {	// built-in slag generation

    cout << endl << "CWB::Toolbox::getSlagList : built-in slags" << endl;

    // find number of non co-located detectors
    vector<int> ifo;
    for(int n=0;n<(int)nIFO;n++) {
      if(slagSite!=NULL) {
        bool unique=true;
        for(int m=0;m<(int)ifo.size();m++) if((int)slagSite[n]==ifo[m]) unique=false;
        if(unique) ifo.push_back(slagSite[n]);
      } else {
        ifo.push_back(n);
      }
    }
    //for(int m=0;m<(int)ifo.size();m++) cout << m << " " << ifo[m] << endl;

    // add slag 0
    int jobId=0;
    if(slagMin==0) {
      slag SLAG;
      SLAG.jobId=jobId++;
      for(int n=0; n<(int)nIFO; n++) SLAG.slagId.push_back(0);
      slagList.push_back(SLAG);
    }
    for(int m=(int)slagMin;m<=(int)slagMax;m++) {
      //int size = slagSize ? slagOff+slagSize : 0;
      int size = 0;				// get all slags contained in the range
      int nifo = ifo.size();			// number of non co-located detectors
      vector<slag> slags;
      getSlagList(slags,size,m,nifo,id);	// get slag list located at range m
      // add slags to slagList in random mode
      gRandom->SetSeed(m+1);     		// random seed is initialized with range
      while(slags.size()>0) {   
        int k = int(gRandom->Uniform(0,int(slags.size())));
        slags[k].jobId=jobId++;
        slagList.push_back(slags[k]);
        slags.erase(slags.begin()+k);		// remove slag from slags list
      } 
    }
    // set slags according to slagSite setup 
    if(slagSite!=NULL) {
      for(int j=0; j<(int)slagList.size(); j++) {
        slag SLAG = slagList[j]; 
        slagList[j].slagId.resize(nIFO);
        for(int n=0; n<(int)nIFO; n++) {
          for(int m=0;m<(int)ifo.size();m++) 
            if((int)slagSite[n]==ifo[m]) slagList[j].slagId[n] = SLAG.slagId[m];
        }
      }
    }
  }

  if((slagSize>slagList.size()-slagOff)) {
    cout << "CWB::Toolbox::getSlagList : Error - the number of available slags = " 
         << slagList.size() << " are less than the requested slag (slagSize) = " << slagSize << endl;
    gSystem->Exit(-1);
  }

  // print slag list
  cout << endl;
  printf("%14s ","SLAG");
  for(int n=0; n<nIFO; n++) printf("%8sifo[%d]","",n);
  printf("\n");
  for(int i=(int)slagOff; i<int(slagOff+slagSize); i++) {
    printf("%14d", slagList[i].jobId);
    for (int n=0; n<nIFO; n++) printf("%14d",slagList[i].slagId[n]);
    printf("\n");
  }
  cout << endl;

  // fill slag list
  int jobID=0;
  slag SLAG;
  vector<slag> slist;
  for(int j=(int)slagOff; j<int(slagOff+slagSize); j++){
    int nseg=0;
    for(int k=1; k<=slagSegs; k++){
      bool check=true;
      for (int n=0; n<(int)nIFO; n++) if((k+slagList[j].slagId[n])>slagSegs) check=false;
      for (int n=0; n<(int)nIFO; n++) if((k+slagList[j].slagId[n])<0) check=false;
      if(check){
        SLAG.jobId=++jobID;
        SLAG.slagId.clear();
        SLAG.segId.clear();
        SLAG.slagId.push_back(slagList[j].jobId);
        SLAG.slagId.push_back(++nseg);  // this number identify the segment number in a slag
        for(int n=0; n<(int)nIFO; n++) {
          SLAG.slagId.push_back(slagList[j].slagId[n]);
          SLAG.segId.push_back(slagList[j].slagId[n]+k);
        }
        slist.push_back(SLAG);
      }
    }
  }
  slagList.clear();

  return slist;
}

//______________________________________________________________________________
void
CWB::Toolbox::getSlagList(vector<slag>& slagList, int slagSize, int slagRank, int nifo, vector<int>& id) {
//
// Get SuperLags list associate to the slagRank value
//
// Input: slagSize - number of super lags 
//        slagRank - slag distance 
//        nifo     - number of ifos
//        id       - temporary auxiliary vector id
//
//  Return the slag list for the slagRank value
//

  if(slagSize && ((int)slagList.size()==slagSize)) return;
  for(int j=-slagRank;j<=slagRank;j++) {             
    if(j==0) continue;                               
    bool unique=true;
    for(int n=0; n<(int)id.size(); n++) if(j==id[n]) unique=false;
    if(!unique) continue;
    if(nifo==2) {
      id.push_back(j);
      int m=0;
      for(int n=0; n<(int)id.size(); n++) m+=abs(id[n]);
      if(m==slagRank) {
        slag SLAG;
        SLAG.jobId=slagList.size();
        SLAG.slagId.push_back(0);  // set first detector with slag shift = 0
        for(int n=0; n<(int)id.size(); n++) SLAG.slagId.push_back(id[n]);
        slagList.push_back(SLAG);
        if(slagSize && ((int)slagList.size()==slagSize)) return;
      }
      id.resize(id.size()-1);
      continue;
    } else {
      id.push_back(j);
      getSlagList(slagList, slagSize, slagRank, nifo-1, id);
      id.resize(id.size()-1);
    }
  }
  return;
}

//______________________________________________________________________________
void 
CWB::Toolbox::dumpSlagList(vector<slag> slagList, TString slagFile, bool slagOnly) {
//
// dump slag list to file 
//
// Input: slagList   - slag list
//        slagFile   - output file name
//        slagOnly   - true:  write only slagID
//                            Example with 3 ifos :
//
//                            # SLAG         ifo[0]        ifo[1]        ifo[2]
//                                 0             0             0             0
//                                 1             0             1            -1
//                                 2             0            -1             1
//
//                   - false: write slagID + segID
//

  if(((int)slagList.size()>0)&&(slagFile!="")) {

    FILE *fP=NULL;
    if((fP = fopen(slagFile.Data(), "w")) == NULL) {
      cout << "CWB::Toolbox::makeSlagList : Error - cannot open file " << slagFile.Data() << endl;
      gSystem->Exit(1);
    }

    int nIFO = slagList[0].segId.size();

    // write header

    fprintf(fP,"#");for (int n=0;n<=nIFO+1;n++) fprintf(fP,"--------------");fprintf(fP,"\n");
    fprintf(fP,"# Super Lags List - %d jobs\n",(int)slagList.size());
    fprintf(fP,"#");for (int n=0;n<=nIFO+1;n++) fprintf(fP,"--------------");fprintf(fP,"\n");
    fprintf(fP,"#         nIFO %13d\n",int(nIFO));
    fprintf(fP,"#");for (int n=0;n<=nIFO+1;n++) fprintf(fP,"--------------");fprintf(fP,"\n");
    if(!slagOnly) fprintf(fP,"#%13s","jobId"); else fprintf(fP,"#");
    fprintf(fP,"%14s","slagId");
    for(int n=0; n<nIFO; n++) fprintf(fP,"%11s[%1d]","segID",int(n));
    fprintf(fP,"\n");
    fprintf(fP,"#");for (int n=0;n<=nIFO+1;n++) fprintf(fP,"--------------");fprintf(fP,"\n");

    // write jobs
    for(int j=0; j<(int)slagList.size(); j++){
      if(slagOnly) {
        if(slagList[j].slagId[1]==1) {
          fprintf(fP,"%14d", slagList[j].slagId[0]);
          for (int n=0; n<nIFO; n++) fprintf(fP,"%14d",slagList[j].slagId[n+2]);
          fprintf(fP,"\n");
        }
      } else {
        fprintf(fP,"%14d", slagList[j].jobId);                               // jobID
        fprintf(fP,"%14d", slagList[j].slagId[0]);                           // slagID
        for (int n=0; n<nIFO; n++) fprintf(fP,"%14d",slagList[j].segId[n]);  // segID
        fprintf(fP,"\n");
      }
    }

    if(fP!=NULL) fclose(fP);
  }

  return;
}

//______________________________________________________________________________
int 
CWB::Toolbox::createDagFile(vector<int> jobList, TString condor_dir, TString label, 
                            vector<TString> jobFiles, TString stage, int jobmin, int jobmax) {
// 
// produce the condor dag file for stages
//
// Input: jobList     - list of job ID
//        condor_dir  - dag,sub condor directory 
//        label       - label used for dag file = 'condor_dir'/'label'.dag
//        jobFiles    - name of job files created by the previous stage
//        stage       - analysis stage (Ex : SUPERCLUSTER, LIKELIHOOD)
//        jobmin      - beg jobList index array
//        jobmax      - end jobList index array
//   

  vector<waveSegment> _jobList(jobList.size());
  for(int i=0;i<(int)jobList.size();i++) _jobList[i].index=jobList[i];
  return createDagFile(_jobList, condor_dir, label, jobFiles, stage, jobmin, jobmax);
}

//______________________________________________________________________________
int 
CWB::Toolbox::createDagFile(vector<waveSegment> jobList, TString condor_dir, TString label, 
                            int jobmin, int jobmax) {
// 
// produce the condor dag file for CWB_STAGE_FULL stage
//
// Input: jobList     - list of job : only waveSegment::index is used
//        condor_dir  - dag,sub condor directory 
//        label       - label used for dag file = 'condor_dir'/'label'.dag
//        jobmin      - beg jobList index array
//        jobmax      - end jobList index array
//

  vector<TString> jobFiles;
  return createDagFile(jobList, condor_dir, label, jobFiles, "CWB_STAGE_FULL", jobmin, jobmax);
}

//______________________________________________________________________________
int 
CWB::Toolbox::createDagFile(vector<waveSegment> jobList, TString condor_dir, TString label, 
                            vector<TString> jobFiles, TString stage, int jobmin, int jobmax) {
// 
// produce the condor dag file for stages
//
// Input: jobList     - list of job : only waveSegment::index is used
//        condor_dir  - dag,sub condor directory 
//        label       - label used for dag file = 'condor_dir'/'label'.dag
//        jobFiles    - name of job files created by the previous stage
//        stage       - analysis stage (Ex : SUPERCLUSTER, LIKELIHOOD)
//        jobmin      - beg jobList index array
//        jobmax      - end jobList index array
//   

  if(jobFiles.size()==0) {          // fill jobFile with CWB_UPARAMETERS_FILE env
    TString cwb_uparameters_file;
    if(gSystem->Getenv("CWB_UPARAMETERS_FILE")!=NULL) {
      cwb_uparameters_file=TString(gSystem->Getenv("CWB_UPARAMETERS_FILE"));
      if(cwb_uparameters_file!="") checkFile(cwb_uparameters_file);
      jobFiles.resize(jobList.size());
      for(int n=0;n<(int)jobFiles.size();n++) jobFiles[n]=cwb_uparameters_file;
    } else {
      cout << "CWB::Toolbox::createDagFile : Error - CWB_UPARAMETERS_FILE env is not defined" << endl;
      gSystem->Exit(1);
    }
  }
  if(jobFiles.size()!=jobList.size()) {
    cout << "CWB::Toolbox::createDagFile : Error - jobFiles size " << jobFiles.size() << 
            "is != jobList size " << jobList.size() << endl;
    gSystem->Exit(1);
  }

  int start = jobmin<1 ? 1 : jobmin;
  int end   = jobmax<1||jobmax>(int)jobList.size() ? jobList.size() : jobmax;

  char ofile[1024];
  sprintf(ofile,"%s/%s.dag",condor_dir.Data(),label.Data());
  ofstream out;
  out.open(ofile,ios::out);

  // read jobList   
  int jID=0;
  for(int n=end;n>=start;n--) {
    if(jobFiles[n-1]=="") continue;	
    jID = jobList[n-1].index; 
    char ostring[1024];
    sprintf(ostring,"JOB A%i %s/%s.sub",jID,condor_dir.Data(),label.Data());
    out << ostring << endl;
    sprintf(ostring,"VARS A%i PID=\"%i\" CWB_UFILE=\"%s\" CWB_STAGE=\"%s\"",
                    jID,jID,jobFiles[n-1].Data(),stage.Data());
    out << ostring << endl;
    if(gSystem->Getenv("_USE_OSG")!=NULL) {
      sprintf(ostring,"SCRIPT POST A%i %s/cwb_net_osg_post.sh %i",jID,gSystem->Getenv("CWB_SCRIPTS"),jID);
      out << ostring << endl;
    } else {
      sprintf(ostring,"RETRY A%i 3000",jID);
      out << ostring << endl;
    }
  }
  out.close();
  return 0;
}

//______________________________________________________________________________
int 
CWB::Toolbox::createDagFile(vector<slag> slagList, TString condor_dir, TString label, 
                            int jobmin, int jobmax) {
// 
// produce the condor dag file for CWB_STAGE_FULL stage
//
// Input: slagList    - list of jobId : only slag::index is used
//        condor_dir  - dag,sub condor directory 
//        label       - label used for dag file = 'condor_dir'/'label'.dag
//        jobmin      - beg jobList index array
//        jobmax      - end jobList index array
//   

  vector<TString> jobFiles;
  return createDagFile(slagList, condor_dir, label, jobFiles, "CWB_STAGE_FULL", jobmin, jobmax);
}

//______________________________________________________________________________
int 
CWB::Toolbox::createDagFile(vector<slag> slagList, TString condor_dir, TString label, 
                            vector<TString> jobFiles, TString stage, int jobmin, int jobmax) {
// 
// produce the condor dag file for stages
//
// Input: slagList    - list of jobId : only slag::index is used
//        condor_dir  - dag,sub condor directory 
//        label       - label used for dag file = 'condor_dir'/'label'.dag
//        jobFiles    - name of job files created by the previous stage
//        stage       - analysis stage (Ex : SUPERCLUSTER, LIKELIHOOD)
//        jobmin      - beg jobList index array
//        jobmax      - end jobList index array
//   

  if(jobFiles.size()==0) {          // fill jobFile with CWB_UPARAMETERS_FILE env
    TString cwb_uparameters_file;
    if(gSystem->Getenv("CWB_UPARAMETERS_FILE")!=NULL) {
      cwb_uparameters_file=TString(gSystem->Getenv("CWB_UPARAMETERS_FILE"));
      if(cwb_uparameters_file!="") checkFile(cwb_uparameters_file);
      jobFiles.resize(slagList.size()); 
      for(int n=0;n<(int)jobFiles.size();n++) jobFiles[n]=cwb_uparameters_file; 
    } else {
      cout << "CWB::Toolbox::createDagFile : Error - CWB_UPARAMETERS_FILE env is not defined" << endl;
      gSystem->Exit(1);
    }
  }
  if(jobFiles.size()!=slagList.size()) {
    cout << "CWB::Toolbox::createDagFile : Error - jobFiles size " << jobFiles.size() << 
            "is != slagList size " << slagList.size() << endl;
    gSystem->Exit(1);
  }

  int start = jobmin<1 ? 1 : jobmin;
  int end   = jobmax<1||jobmax>(int)slagList.size() ? slagList.size() : jobmax;

  char ofile[1024];
  sprintf(ofile,"%s/%s.dag",condor_dir.Data(),label.Data());
  ofstream out;
  out.open(ofile,ios::out);

  // read slagList   // SLAG
  int jID=0;
  for(int n=end;n>=start;n--) {
    if(jobFiles[n-1]=="") continue;	
    jID = slagList[n-1].jobId; 
    char ostring[1024];
    sprintf(ostring,"JOB A%i %s/%s.sub",jID,condor_dir.Data(),label.Data());
    out << ostring << endl;
    sprintf(ostring,"VARS A%i PID=\"%i\" CWB_UFILE=\"%s\" CWB_STAGE=\"%s\"",
                    jID,jID,jobFiles[n-1].Data(),stage.Data());
    out << ostring << endl;
    sprintf(ostring,"RETRY A%i 3000",jID);
    out << ostring << endl;
    if(gSystem->Getenv("_USE_OSG")!=NULL) {
      sprintf(ostring,"SCRIPT POST A%i %s/cwb_net_osg_post.sh",jID,gSystem->Getenv("CWB_SCRIPTS"));
      out << ostring << endl;
    }
  }
  out.close();
  return 0;
}

//______________________________________________________________________________
int 
CWB::Toolbox::createSubFile(TString label, TString condor_dir, TString out_dir, TString err_dir, 
                            TString log_dir, TString ext, TString condor_tag) {
// 
// produce the condor sub file 
//
// Input: label       - label used for dag file = 'condor_dir'/'label'.dag
//        condor_dir  - dag,sub condor directory 
//        out_dir     - out condor files directory
//        err_dir     - err condor files directory
//        ext         - condor directory
//        condor_tag  - Define a Unique Tag for Condor Jobs 
//   

  char ofile[1024];
  if(ext=="")
    sprintf(ofile,"%s/%s.sub",condor_dir.Data(),label.Data());
  else
    sprintf(ofile,"%s/%s.sub.%s",condor_dir.Data(),label.Data(),ext.Data());

  FILE *fP=NULL;
  if((fP = fopen(ofile, "w")) == NULL) {
    cout << "CWB::Toolbox::createSubFile : Error - cannot open file " << ofile << endl;
    gSystem->Exit(1);
  }

  fprintf(fP,"universe = vanilla\n");
  fprintf(fP,"getenv = true\n");
  fprintf(fP,"priority = $(PRI)\n");
  fprintf(fP,"on_exit_hold = ( ExitCode != 0 )\n");
  fprintf(fP,"request_memory = 3000\n");
  fprintf(fP,"executable = cwb.sh\n");
  fprintf(fP,"job_machine_attrs = Machine\n");
  fprintf(fP,"job_machine_attrs_history_length = 5\n");
  fprintf(fP,"requirements = target.machine =!= MachineAttrMachine1 && target.machine =!= MachineAttrMachine2 && target.machine =!= MachineAttrMachine3 && target.machine =!= MachineAttrMachine4 && target.machine =!= MachineAttrMachine5\n");
  fprintf(fP,"environment = CWB_JOBID=$(PID);CWB_UFILE=$(CWB_UFILE);CWB_STAGE=$(CWB_STAGE)\n");
  if(condor_tag!="") fprintf(fP,"accounting_group = %s\n",condor_tag.Data());
  fprintf(fP,"output = %s/$(PID)_%s_$(CWB_STAGE).out\n",out_dir.Data(),label.Data());
  fprintf(fP,"error = %s/$(PID)_%s_$(CWB_STAGE).err\n",err_dir.Data(),label.Data());
  if(gSystem->Getenv("_USE_OSG")!=NULL) {
    TString home_dir= TString(gSystem->Getenv("HOME"));
    //cout<<home_dir.Data()<<endl;
    fprintf(fP,"log = %s/condor/%s.log\n",home_dir.Data(),label.Data());
    fprintf(fP,"should_transfer_files = YES\n");
    fprintf(fP,"when_to_transfer_output = ON_EXIT\n");
    fprintf(fP,"transfer_input_files = %s/.rootrc,  %s/%s.tgz\n",home_dir.Data(), condor_dir.Data(),label.Data());
    fprintf(fP,"transfer_output_files = %s/output, %s/log\n",  label.Data(), label.Data());
    fprintf(fP,"request_memory = 2000\n");
  } else {
    fprintf(fP,"log = %s/%s.log\n",log_dir.Data(),label.Data());
  }
  fprintf(fP,"notification = never\n");
  fprintf(fP,"rank=memory\n");
  fprintf(fP,"queue\n");

  fclose(fP);

  return 0;
}

//______________________________________________________________________________
vector<int> 
CWB::Toolbox::getCondorJobList(TString condor_dir, TString label) {
// 
// get the job list from the dag condor file 
//
// Input: condor_dir  - dag,sub condor directory 
//        label       - label used for dag file = 'condor_dir'/'label'.dag
//   

  char ifile[1024];
  sprintf(ifile,"%s/%s.dag",condor_dir.Data(),label.Data());
  ifstream in;
  in.open(ifile);
  if(!in.good()) {cout << "Error Opening File : " << ifile << endl;gSystem->Exit(1);}

  vector<int> jobs;

  char istring[1024];
  while(1) {
    in.getline(istring,1024);
    if (!in.good()) break;
    TObjArray* token = TString(istring).Tokenize(TString(' '));
    TObjString* stoken =(TObjString*)token->At(0);
    TString jobLabel = stoken->GetString();
    if(jobLabel.CompareTo("JOB")!=0) continue;
    stoken =(TObjString*)token->At(1);
    int jobID=TString(stoken->GetString()).ReplaceAll("A","").Atoi(); 
    jobs.push_back(jobID);
    if(token) delete token;
  }

  in.close();

  vector<int> jobList;
  Int_t *id = new Int_t[jobs.size()];
  Int_t *jd = new Int_t[jobs.size()];
  for(int i=0;i<(int)jobs.size();i++) jd[i]=jobs[i];
  TMath::Sort((int)jobs.size(),jd,id,false);
  for(int i=0;i<(int)jobs.size();i++) jobList.push_back(jd[id[i]]);
  jobs.clear();
  delete [] id;
  delete [] jd;

  return jobList;
}

//______________________________________________________________________________
vector<float> 
CWB::Toolbox::getJobBenchmark(TString ifName, int stageID, TString bench) {
//
// extract benchmark info from bench status history line
// return in vector vbench : vbench[0]=job number - vbench[1]=bench value
// 
  return getJobBenchmark(ifName, stageID, -1, -1, bench); 
}

//______________________________________________________________________________
vector<float> 
CWB::Toolbox::getJobBenchmark(TString ifName, int stageID, int resID, int factorID, TString bench) {
//
// extract benchmark info from bench status history line
// return in vector vbench : vbench[0]=job number - vbench[1]=bench value
// 
// These are the format of the bench status history line
// GPS:X1-JOB:X2-STG:X3-FCT:X4-JET:X5-SET:X6-MEM:X7-JFS:X8
// GPS:X1-JOB:X2-STG:X3-FCT:X4-RES:X5-THR:X6
// GPS:X1-JOB:X2-STG:X3-FCT:X4-RES:X5-CSIZE:X6-PSIZE:X7
//
// ifName   - root file name containig the history object
// stageID  - stage ID (defined in cwb.hh)
// resID    - resolution level ID
// factorID - factor ID (defined in user_parameters.C)
//
// bench is :  
//   GPS : gps time (sec)
//   JOB : job number 
//   STG : stage number (defined in cwb.hh)
//   FCT : factors index array 
//   JET : Job Elapsed Time (sec)
//   SET : Stage Elapsed Time (sec)
//   MEM : Memory usage (MB)
//   JFS : Job File Size - temporary file (bytes)
//   THR : Threshold @ Coherence Stage
// CSIZE : cluster size per lag @ coherence stage   
// PSIZE : pixels  size per lag @ coherence stage   
// 

  vector<float> vbench(2);
  vbench[0]=-1;		// job number
  vbench[1]=-1;		// bench value

  TFile *ifile = TFile::Open(ifName);
  if(ifile==NULL) {
    cout << "Failed to open " << ifName.Data() << endl;
    return vbench;
  }

  CWB::History* ihistory = (CWB::History*)ifile->Get("history");
  if(ihistory==NULL) {
    cout << "Error : history is not present!!!" << endl;
    return vbench;
  }

  // build stage tag in the bench line status 
  char stageName[64]; sprintf(stageName,"STG:%d-",stageID); 

  // build resolution tag in the bench line status 
  char resName[64]; sprintf(resName,"RES:%d-",resID); 

  // build factor tag in the bench line status 
  char factorName[64]; sprintf(factorName,"FCT:%d-",factorID); 

  int log_size = ihistory->GetLogSize(CCAST("FULL"));
  for(int i=0;i<log_size;i++) {
    TString log = ihistory->GetLog(CCAST("FULL"),i);
    //cout << "log " << log.Data() << endl;
    if(!log.Contains(stageName)) continue;
    if(resID>=0 && !log.Contains(resName)) continue;
    if(factorID>=0 && !log.Contains(factorName)) continue;
    if(log.Contains(bench+TString(":"))) {
       TObjArray* token = log.Tokenize('\n');
       for(int j=0;j<token->GetEntries();j++) {
         TString line = ((TObjString*)token->At(j))->GetString();
         // extract bench line which contains the "bench:" string
         if(line.Contains(bench+TString(":"))) {
           //cout << j << " " << (((TObjString*)token->At(j))->GetString()).Data() << endl;
           // extract tokens separated by "-"
           TObjArray* ltoken = line.Tokenize('-');
           for(int k=0;k<ltoken->GetEntries();k++) {
             TString stat = ((TObjString*)ltoken->At(k))->GetString();
             //cout << k << " " << stat.Data() << endl;
             TObjArray* stoken = stat.Tokenize(':');
             if(stoken->GetEntries()!=2) continue;
             TString stat_name  = ((TObjString*)stoken->At(0))->GetString();
             TString stat_value = ((TObjString*)stoken->At(1))->GetString();
             //cout << stat_name.Data() << " " << stat_value.Data() << endl;
             // extract value associated to bench "bench:value"
             if(stat_name==bench) {
               float value = stat_value.Atof();
               if(value>vbench[1]) vbench[1]=value;
               //cout << stat_name.Data() << " " << stat_value.Data() << endl;
             }
             if(stat_name=="JOB") {
               vbench[0] = stat_value.Atoi();    // Get JOB ID from bench status line
             }
             delete stoken;
           }
         delete ltoken;
         }
       }
       delete token;
    }
  }

  ifile->Close();

  // if jobID is not present in the bench line then it is extracted from file name
  if(vbench[0]==-1) vbench[0] = getJobId(ifName);    

  return vbench;
}

//______________________________________________________________________________
TString
CWB::Toolbox::DAG2LSF(char* dagFile, char* data_label, char* nodedir, char* data_dir,
              char* condor_dir, char* log_dir, char* output_dir, char* work_dir) {
//
// Input: dagFile	- Dag File 
//	  data_label    - data label
//	  nodedir       - temporary dir directory
//	  data_dir      - data directory
//	  condor_dir    - condor directory
//	  log_dir	- log dir
//	  output_dir	- output dir
//	  work_dir	- working dir
//
// Output: input dagFile is used to produce the script LSF file 
//         if number of jobs>0 return lsf file name otherwise empty string
//

  // get user name
  UserGroup_t* uinfo = gSystem->GetUserInfo();
  TString uname = uinfo->fUser;

  // get home wat path
  TString cwb_home_wat="";
  if(gSystem->Getenv("HOME_WAT")!=NULL) {
    cwb_home_wat=TString(gSystem->Getenv("HOME_WAT"));
  }                                                   
  if(cwb_home_wat=="") {                              
    cout << "CWB::Toolbox::DAG2LSF : Error : HOME_WAT not defined !!!" << endl;
    gSystem->Exit(1);                                                        
  }                                                                          

  // get LSF queue name
  TString lsf_queue="";
  if(gSystem->Getenv("LSF_QUEUE")!=NULL) {
    lsf_queue=TString(gSystem->Getenv("LSF_QUEUE"));
  }                                                 
  if(lsf_queue=="") {                               
    cout << "CWB::Toolbox::DAG2LSF : Error : LSF_QUEUE not defined !!!" << endl;
    gSystem->Exit(1);                                                         
  }                                                                           

  // --------------------------------------------------
  // create LSF job script file                        
  // --------------------------------------------------

  // open dag file
  ifstream in;                                     
  in.open(dagFile);                                  
  if(!in.good()) {
    cout << "CWB::Toolbox::DAG2LSF - Error Opening File : " << dagFile << endl;
    gSystem->Exit(1);
  }

  // create LSF file
  TString lsfFile = dagFile;
  lsfFile.ReplaceAll(".dag",".lsf");
  ofstream out;                                   
  out.open(lsfFile.Data(),ios::out);                       
  if(!out.good()) {
    cout << "CWB::Toolbox::DAG2LSF - Error Opening File : " << lsfFile << endl;
    gSystem->Exit(1);
  }                                                                                                                 

  out << "#!/bin/bash" << endl << endl; 

  int lsf_job_cnt=0;
  char istring[1024];
  while(1) {         
    int PID=0;       
    TString CWB_UFILE="";
    TString CWB_STAGE="";
    in.getline(istring,1024);
    if (!in.good()) break;   
    if(!TString(istring).BeginsWith("VARS")) continue;
    TObjArray* token = TString(istring).Tokenize(TString(' '));
    for(int j=2;j<token->GetEntries();j++) {                   
      TString item = ((TObjString*)token->At(j))->GetString(); 
      if(item.BeginsWith("PID=")) {                            
        item.ReplaceAll("PID=","");                            
        item.ReplaceAll("\"","");                              
        PID = item.Atoi();                                     
        continue;                                              
      }                                                        
      if(item.BeginsWith("CWB_UFILE=")) {                      
        item.ReplaceAll("CWB_UFILE=","");                      
        item.ReplaceAll("\"","");                              
        CWB_UFILE=item;;                                       
        continue;                                              
      }                                                        
      if(item.BeginsWith("CWB_STAGE=")) {                      
        item.ReplaceAll("CWB_STAGE=","");                      
        item.ReplaceAll("\"","");                              
        CWB_STAGE=item;;                                       
        continue;                                              
      }                                                        
    }                                                          
    //cout << "PID : " << PID << " CWB_UFILE : " << CWB_UFILE << " CWB_STAGE : " << CWB_STAGE << endl;
    char lsf_cmd[2048];                                                                               
    char lsf_label[1024]; 	// label used in the execution node
    if(CWB_STAGE=="CWB_STAGE_FULL") {
      sprintf(lsf_label,"%s",data_label);	
    } else {
      sprintf(lsf_label,"%s_%s",data_label,CWB_STAGE.Data());	
    }
    sprintf(lsf_cmd,"bsub -q %s -J A%d -g /%s/%s \
                     \\\n -f \"%s > %s/%d_%s.ufile\" \
                     \\\n -f \"%s/%s.tgz > %s/%d_%s.tgz\" \
                     \\\n -o %d_%s.out -f \"%s/%s/%d_%s.out < %d_%s.out\" \
                     \\\n -e %d_%s.err -f \"%s/%s/%d_%s.err < %d_%s.err\" \
                     \\\n -f \"%s/%s/%d_%s.tgz < %s/%d_%s.tgz\" \
                     \\\n -Ep \"rm %d_%s.out\" \
                     \\\n -Ep \"rm %d_%s.err\" \
                     \\\n -Ep \"rm %s/%d_%s.tgz\" \
                     \\\n \"/bin/bash %s/tools/cwb/scripts/cwb_net_lsf.sh %d %s %s %s %s %s\"",
                     lsf_queue.Data(), PID, uname.Data(), data_label,
                     CWB_UFILE.Data(), nodedir, PID, lsf_label,
                     condor_dir, lsf_label, nodedir, PID, lsf_label,
                     PID, lsf_label, work_dir, log_dir, PID, lsf_label, PID, lsf_label,
                     PID, lsf_label, work_dir, log_dir, PID, lsf_label, PID, lsf_label,
                     work_dir, output_dir, PID, lsf_label, nodedir, PID, lsf_label,
                     PID, lsf_label,
                     PID, lsf_label,
                     nodedir, PID, lsf_label,
                     cwb_home_wat.Data(), PID, CWB_STAGE.Data(), nodedir, data_label, data_dir, CWB_UFILE.Data());
    out << lsf_cmd << endl << endl;                                                                   
    //cout << lsf_cmd << endl;                                                                        
    if(token) delete token;                                                                           
    //if(PID%1000==0) cout << PID << endl;                                                               
    lsf_job_cnt++;
  }                                                                                                   
  out.close();                                                                                        
  in.close();                                                                                         

  return lsf_job_cnt>0 ? lsfFile : "";
}

//______________________________________________________________________________
vector<int> 
CWB::Toolbox::getMergeJobList(TString merge_dir, TString label, int version) {
//
// return vector job numbers extracted from the merged file list
// the job number ID is extracted from the file name : XXX_jobID.root
//
// Input: merge_dir     - merge directory
//        label         - label of the merge name list
//        version       - version of the merge name list
//              	- merge name list = 'merge_dir'/merge_'label'.M'version.lst
//

  char ifile[1024];
  sprintf(ifile,"%s/merge_%s.M%d.lst",merge_dir.Data(),label.Data(),version);
  vector<TString> jfname;
  return getMergeJobList(ifile,jfname);
}

//______________________________________________________________________________
vector<int> 
CWB::Toolbox::getMergeJobList(TString ifname) {
//
// return vector job numbers extracted from the ifname file list
// the job number ID is extracted from the file name : XXX_jobID.root
//
// Input: ifname	- input file name which contains the list of files
//

  vector<TString> jfname;
  return getMergeJobList(ifname,jfname);
}

//______________________________________________________________________________
vector<int> 
CWB::Toolbox::getMergeJobList(TString ifname, vector<TString>& jobFileList) {
//
// return vector job numbers extracted from the ifname file list
// the job number ID is extracted from the file name : XXX_jobID.root
//
// Input: ifname	- input file name which contains the list of files
//
// Output: jobFileList  - the job file names list  
//

  ifstream in;
  in.open(ifname.Data());
  if(!in.good()) {
     cout << "CWB::Toolbox::getMergeJobList : Error Opening File : " << ifname << endl;
     gSystem->Exit(1);
  }

  vector<int> job;
  vector<TString> jfname;

  char istring[1024];
  while(1) {
    in.getline(istring,1024);
    if (!in.good()) break;
    TObjArray* token = TString(istring).Tokenize(TString('_'));
    TObjString* stoken =(TObjString*)token->At(token->GetEntries()-1);
    int jobID=TString(stoken->GetString()).ReplaceAll("job","").ReplaceAll(".root","").Atoi();
    if(jobID==0) continue;	// file root is a merge file not a job file
    job.push_back(jobID);
    jfname.push_back(istring);
    if(token) delete token;
  }

  in.close();

  jobFileList.clear();
  vector<int> jobList;
  Int_t *id = new Int_t[job.size()];
  Int_t *jd = new Int_t[job.size()];
  for(int i=0;i<(int)job.size();i++) jd[i]=job[i];
  TMath::Sort((int)job.size(),jd,id,false);
  for(int i=0;i<(int)job.size();i++) {
    jobList.push_back(jd[id[i]]);
    jobFileList.push_back(jfname[id[i]]);
  }
  job.clear();
  jfname.clear();
  delete [] id;
  delete [] jd;

  return jobList;
}

/*--------------------------------------------------------------------------
  return the list of segments with length=segLen contained in the interval 
  ilist time_max-time_min
  the edges of the segments are a multiple of segLen 
--------------------------------------------------------------------------*/
//______________________________________________________________________________
vector<waveSegment> 
CWB::Toolbox::getSlagJobList(vector<waveSegment> ilist, int seglen) {
//
// Extract Time MIN and MAX from the ilist 
// recompute MIN MAX to be a integer multiple of seglen
// Compute the list of segments with length seglen within the MIN-MAX range 
//
//
// Input:  ilist  - list of segments (start,stop)
//         seglen - segment length 
//
// Output: return the list of segment time ranges 
// 

  if(ilist.size()==0) {cout << "CWB::Toolbox::getSlagJobList - Error ilist size=0" << endl;gSystem->Exit(1);}
  if(seglen<=0) {cout << "CWB::Toolbox::getSlagJobList - Error seglen<=0" << endl;gSystem->Exit(1);}

  waveSegment SEG;
  vector<waveSegment> jlist;

  int start=ilist[0].start;
  int stop=ilist[ilist.size()-1].stop;
  start=seglen*TMath::Nint(double(start/seglen));
  stop=seglen*TMath::Nint(double(stop/seglen));

  int njob=(stop-start)/seglen;
  for(int n=0;n<njob;n++) {
    SEG.start=start+n*seglen;
    SEG.stop=SEG.start+seglen;
    jlist.push_back(SEG);
  }

  return jlist;
}

//______________________________________________________________________________
waveSegment 
CWB::Toolbox::getMaxSeg(vector<waveSegment> list) {
//
//
// Input: list - segment list
//
// Return the segment with the maximum length
//

  waveSegment SEG={0,0,0};
  for(int i=0;i<(int)list.size();i++) {
    if((list[i].stop-list[i].start)>(SEG.stop-SEG.start)) SEG=list[i];
  }
  return SEG;
}

//______________________________________________________________________________
slag 
CWB::Toolbox::getSlag(vector<slag> slagList, int jobid) {
//
//
// Input: slagList - slag list
//        jobid    - job id 
//
// Return SLAG structure with slagList.jobId == jobid
//

  slag SLAG;SLAG.jobId=-1;
  for(int j=0;j<(int)slagList.size();j++) if(slagList[j].jobId==jobid) {SLAG=slagList[j];break;}
  return SLAG;
}

//______________________________________________________________________________
vector<slag> 
CWB::Toolbox::getSlagList(vector<slag> islagList, vector<TString> ifos,
                   double segLen, double segMin, double segEdge, 
                   int nDQF, dqfile* iDQF, CWB_CAT dqcat) {
//
//
// Input: islagList - vector list of slag structures
//        ifos      - vector list of ifo names
//        segLen    - Segment length [sec]
//        segMin    - Minimum Segment Length after dqcat [sec]
//        segEdge   - wavelet boundary offset [sec]
//        nDQF      - size of iDQF array 
//        iDQF      - DQ structure array
//        dqcat     - dq cat
//
// if dqcat=CWB_CAT1 -> return the list of slags with dq len > segMin+2*segEdge
// if dqcat=CWB_CAT2 -> return the list of slags with dq len > segMin
//

  // dqcat must be CWB_CAT1 or CWB_CAT2
  if((dqcat!=CWB_CAT1)&&(dqcat!=CWB_CAT2)) {
    cout << "CWB::Toolbox::getSlagList : dqcat must be CWB_CAT1 or CWB_CAT2 !!!" << endl;
    gSystem->Exit(1);
  }

  int nIFO=0;
  slag SLAG;
  int lsize=islagList.size();
  int nRejected=0;
  int livetime1=0;
  int livetime2=0;
  int segLength=0;
  waveSegment SEG,MSEG;
  vector<waveSegment> dqList;
  vector<waveSegment> dq1List;
  vector<waveSegment> jobList;
  vector<waveSegment> segList;
  vector<slag> oslagList;
  int jobID=0;int segID[20];int slagID=0;

  if(lsize>0) nIFO=islagList[0].segId.size(); 

  dqfile* DQF = new dqfile[nDQF];
  for(int i=0;i<nDQF;i++) DQF[i]=iDQF[i];

  // get zero lag merged dq cat 1 list
  dq1List=readSegList(nDQF, DQF, CWB_CAT1);
  // get number/list of the available super lag jobs
  jobList=getSlagJobList(dq1List, segLen);
  cout << "input list size " << lsize << endl;
  cout << "jobList    size " << jobList.size() << endl;
  dqList.clear();

  // read dq lists from files  and store in ilist
  vector<waveSegment>* ilist = new vector<waveSegment>[nDQF];
  for(int i=0;i<nDQF;i++) ilist[i]=readSegList(DQF[i]);

  // get range segment for this job intersecting with cat1 list (dqList)
  for(int j=0;j<lsize;j++) {
    if(j%1000==0) printf("%6d - Rejected slags = %d/%lu\n",j,nRejected,islagList.size());
    SLAG   = islagList[j];	// slag structure
    slagID = SLAG.slagId[0];
    jobID  = SLAG.jobId;
    for(int n=0; n<nIFO; n++) segID[n]=SLAG.segId[n];
    cout.precision(2);

    // shift dq files
    for(int i=0;i<nDQF;i++) DQF[i].shift=0.;
    setSlagShifts(SLAG, ifos, segLen, nDQF, DQF);

    // create shifted merged dq cat file list
    dqList.clear();
    bool first=true;
    for(int i=0;i<nDQF;i++) {
      if(DQF[i].cat<=dqcat) {
        for(int k=0;k<(int)ilist[i].size();k++) {	// apply time shifts
          ilist[i][k].start+=DQF[i].shift;
          ilist[i][k].stop+=DQF[i].shift;
        }
        if(first) {dqList=ilist[i];first=false;}
        else dqList = mergeSegLists(ilist[i],dqList);
        for(int k=0;k<(int)ilist[i].size();k++) {	// restore time shifts
          ilist[i][k].start-=DQF[i].shift;
          ilist[i][k].stop-=DQF[i].shift;
        }
      }
    }

    // create job file list
    if(dqList.size()==0) continue;

    segList.clear();
    SEG.start = jobList[segID[0]-1].start-segEdge;
    SEG.stop  = jobList[segID[0]-1].stop+segEdge;

    if(dqcat==CWB_CAT2) {

      // create shifted merged dq cat1 file list
      dq1List.clear();
      bool first=true;
      for(int i=0;i<nDQF;i++) {
        if(DQF[i].cat<=CWB_CAT1) {
          for(int k=0;k<(int)ilist[i].size();k++) {       // apply time shifts
            ilist[i][k].start+=DQF[i].shift;
            ilist[i][k].stop+=DQF[i].shift;
          }
          if(first) {dq1List=ilist[i];first=false;}
          else dq1List = mergeSegLists(ilist[i],dq1List);
          for(int k=0;k<(int)ilist[i].size();k++) {       // restore time shifts
            ilist[i][k].start-=DQF[i].shift;
            ilist[i][k].stop-=DQF[i].shift;
          }
        }
      }

      segList.push_back(SEG);
      segList = mergeSegLists(dq1List,segList); // merge detector segments with dq cat 1
      SEG = getMaxSeg(segList);                 // extract max length segment (only this segment is used for analysis)
      SEG.start+=segEdge;
      SEG.stop-=segEdge;
      segList.clear();
    }

    segList.push_back(SEG);
    segList = mergeSegLists(dqList,segList);    // merge detector segments with dq cat 
    double lenTHR;
    if(dqcat==CWB_CAT1) {
      MSEG = getMaxSeg(segList);                  // max length segment
      segLength=MSEG.stop-MSEG.start;
      lenTHR=segMin+2*segEdge;
    }
    if(dqcat==CWB_CAT2) {
      segLength=getTimeSegList(segList);  
      lenTHR=segMin;
    }
    livetime1+=segLength;
    if(segLength<lenTHR) {
      nRejected++;
    } else {
      livetime2+=segLength;
      oslagList.push_back(SLAG);
    }
  }
  printf("%6lu - Rejected slags = %d/%lu\n\n",islagList.size(),nRejected,islagList.size());
  printf("Slag livetime before dq cat%d is approximately %d sec \t= %.2f days\n",dqcat,livetime1,double(livetime1)/86400.);
  printf("Slag livetime after  dq cat%d is approximately %d sec \t= %.2f days\n",dqcat,livetime2,double(livetime2)/86400.);

  dqList.clear();
  dq1List.clear();
  for(int i=0;i<nDQF;i++) ilist[i].clear();
  delete [] ilist;
  delete [] DQF;

  return oslagList;
}

//______________________________________________________________________________
void
CWB::Toolbox::setSlagShifts(slag SLAG, vector<TString> ifos, double segLen, int nDQF, dqfile* DQF) {
//
// Apply slag shifts into the DQF structures (DQF[].shift)
//
// Input: SLAG   - slag structure
//        segLen - segment length
//        nDQF   - size of DQF array
//        DQF    - array of DQ structures
// 


  // shift dq files
  for(int i=0;i<nDQF;i++) {
    int ifoID=-1;
    for(int j=0;j<(int)ifos.size();j++) if(ifos[j]==DQF[i].ifo) ifoID=j;
    if(ifoID==-1) {
      cout << "CWB::Toolbox::setSlagShifts - " << DQF[i].ifo << " is not in the contained into ifos vector" << endl;
      gSystem->Exit(1);
    }
    if(ifoID>=(int)SLAG.segId.size()) {
      cout << "CWB::Toolbox::setSlagShifts - " << DQF[i].ifo << " index " << ifoID << "is not correct" << endl;
      gSystem->Exit(1);
    }
    DQF[i].shift+=-segLen*(SLAG.segId[ifoID]-SLAG.segId[0]);  // apply slag shift to dq
  }

  return;
}

//______________________________________________________________________________
vector<waveSegment> 
CWB::Toolbox::getSegList(slag SLAG, vector<waveSegment> jobList,
                  double segLen, double segMLS, double segEdge, 
                  vector<waveSegment> dqList) {
//
// create merged dq cat 1 file list
//
// Input: SLAG     - slag structure
//        jobList  - list of available super lag segments
//        segLen   - Segment length [sec]
//        segMLS   - Minimum Segment Length after DQ_CAT1 [sec]
//        segEdge  - wavelet boundary offset [sec]
//        dqList   - dq cat1 list
//
// Return the detector's slag segment ranges associated to SLAG infos & dq cat1 list
//

  waveSegment SEG,MSEG;
  vector<waveSegment> segList;

  int nIFO=SLAG.segId.size(); 

  // get range segment for this job intersecting with cat1 list (dqList)
  SEG.start = jobList[SLAG.segId[0]-1].start-segEdge;
  SEG.stop  = jobList[SLAG.segId[0]-1].stop+segEdge;
  segList.push_back(SEG);
  segList = mergeSegLists(dqList,segList);    // merge detector segments with dq cat 1
  MSEG    = getMaxSeg(segList);               // extract max length segment (only this segment is used for analysis)
  int segLength=MSEG.stop-MSEG.start;
  if(segLength<segMLS+2*segEdge) {
    cout << "CWB::Toolbox::getSegList - Segment length too small : " << segLength << endl;
    gSystem->Exit(1);
  }

  segList.clear();
  for(int n=0; n<nIFO; n++) {
    SEG.start = MSEG.start+segLen*(SLAG.segId[n]-SLAG.segId[0])+segEdge;
    SEG.stop  = MSEG.stop+segLen*(SLAG.segId[n]-SLAG.segId[0])-segEdge;
    segList.push_back(SEG);
  }
  jobList.clear();

  return segList;
}

//_______________________________________________________________________________________
double 
CWB::Toolbox::getLiveTime(vector<waveSegment>& jobList, vector<waveSegment>& dqList) {

  vector<double> dummy;
  return getLiveTime(jobList, dqList, dummy);
}

//_______________________________________________________________________________________
double 
CWB::Toolbox::getLiveTime(vector<waveSegment>& jobList, vector<waveSegment>& dqList, vector<double> shiftList) {
//
// return live time in sec
//
// NOTE : this procedure compute the livetime of the jobList after dq=cat2 when lags are applied
//        the returned livetime should be with a good approximation consistent 
//        with the one computed by the cWB algorithm
//
// Input:  jobList   - job list
//         dqList    - data quality list
//         shiftList - list of time shifts (sec), the size is the detector number
//                     shift[0]=shift_ifo_1, shift[1]=shift_ifo_2, ... shift[n-1]=shift_ifo_n
//                     if shift list size is 0 then no shifts are applied
//
// Procedure : for each time shift and for each job the data quality list inside the job interval
//             are circulary rotated by a time=shift, the new shifted list is saved into a new segment list   
//             finaly the new segment list are merged and the total livetime is returned
//

  if(shiftList.size()==0)  shiftList.push_back(0.);
  int nshift = shiftList.size();       		// extract the number of shifts

  waveSegment seg; 
  vector<waveSegment>* segList = new vector<waveSegment>[nshift];
  vector<waveSegment> mergeList=mergeSegLists(jobList,dqList);

  for(int i=0;i<nshift;i++) {				// loop over shifts
    double shift = shiftList[i];
    if(shift==0) {segList[i]=mergeList;continue;}	// skip if shift=0
    int jsize = jobList.size();
    for(int j=0;j<jsize;j++) {				// loop over job list
      double jstart = jobList[j].start;
      double jstop  = jobList[j].stop;
      double jlen   = jstop-jstart;
      int ksize = mergeList.size();
      for(int k=0;k<ksize;k++) {			// loop over dq list
        double kstart = mergeList[k].start;
        double kstop  = mergeList[k].stop;
        if((kstop<=jstart)||(kstart>=jstop)) continue;	// skip dq outside the job range
        // add shift
        kstart+=shift;
        kstop +=shift;
        // circular shift
        if(kstop<=jstop) {
          seg.start=kstart;seg.stop=kstop;segList[i].push_back(seg);
        } 
        if(kstart>=jstop) {
          seg.start=kstart-jlen;seg.stop=kstop-jlen;segList[i].push_back(seg);
        } 
        if((kstart<jstop)&&(kstop>jstop)) {
          seg.start=kstart;seg.stop=jstop;segList[i].push_back(seg);
          seg.start=jstart;seg.stop=kstop-jlen;segList[i].push_back(seg);
        } 
      } 
    }
  }

  // sort lists (needed by mergeSegLists)
  for(int i=0;i<nshift;i++) segList[i]=unionSegments(segList[i]);

  // merge segLists
  for(int i=1;i<nshift;i++) segList[0]=mergeSegLists(segList[0],segList[i]);

  double livetime=getTimeSegList(segList[0]);
  for(int i=0;i<nshift;i++) segList[i].clear();
  delete [] segList;

  return livetime;
}

//______________________________________________________________________________
int 
CWB::Toolbox::shiftBurstMDCLog(std::vector<std::string>& mdcList, vector<TString> ifos, double mdc_shift) {
//
// apply a time shift to the log entries contained in the mdcList
//
// Input: mdcList   - list of the log string (one entry for each injected MDC)
//        ifos      - array of input detector's names
//        mdc_shift - time shift
//  

  cout.precision(14);
  for(int i=0;i<(int)mdcList.size();i++) {
    char mdc_string[1024]="";
    //cout << "IN-" << mdcList[i] << endl;
    TObjArray* token = TString(mdcList[i]).Tokenize(' ');
    sprintf(mdc_string,"%s",(((TObjString*)token->At(0))->GetString()).Data());
    for(int j=1;j<9;j++) 
      sprintf(mdc_string,"%s %s",mdc_string,(((TObjString*)token->At(j))->GetString()).Data());
    double frTime = (((TObjString*)token->At(9))->GetString()).Atof();
    sprintf(mdc_string,"%s %10.6f",mdc_string,frTime+mdc_shift);
    //cout << "TIME " << frTime << endl;
    double mdcTime = (((TObjString*)token->At(10))->GetString()).Atof();
    sprintf(mdc_string,"%s %10.6f",mdc_string,mdcTime+mdc_shift);
    //cout << "TIME " << mdcTime << endl;
    for(int j=11;j<15;j++) 
      sprintf(mdc_string,"%s %s",mdc_string,(((TObjString*)token->At(j))->GetString()).Data());
    for(int j=15;j<token->GetEntries();j++) {
      TString stoken = ((TObjString*)token->At(j))->GetString();
      sprintf(mdc_string,"%s %s",mdc_string,stoken.Data());
      for(int n=0;n<(int)ifos.size();n++) {         
        if(ifos[n]==stoken) {
          double ifoTime = (((TObjString*)token->At(j+1))->GetString()).Atof();
          sprintf(mdc_string,"%s %10.6f",mdc_string,ifoTime+mdc_shift);
          //cout << "TIME " << ifos[n].Data() << " " << ifoTime << endl;
          j++; 
        }
      }
    }
    mdcList[i]=mdc_string;
    //cout << "OUT-" << mdcList[i] << endl;
    if(token) delete token;
  }

  return 0;
}

//______________________________________________________________________________
TString 
CWB::Toolbox::SetMDCLog(TString log, int pos, double val) {
//
// set entry at position pos contained in the log string
//
// Input: log      - log string 
//        pos      - entry position in the log string
//        val      - string value used to modify the log string  
//
// return the modified log string
//

   char sval[32];sprintf(sval,"%e",val); 
   return SetMDCLog(log, pos, sval);
}

//______________________________________________________________________________
TString 
CWB::Toolbox::SetMDCLog(TString log, int pos, TString val) {
//
// set entry at position pos contained in the log string
//
// Input: log      - log string 
//        pos      - entry position in the log string
//        val      - string value used to modify the log string  
//
// return the modified log string
//

  char log_string[1024]="";

  TObjArray* token = TString(log).Tokenize(' ');
  for(int n=0;n<token->GetEntries();n++) {
    if(n==pos) sprintf(log_string,"%s %s",log_string, val.Data());
    else       sprintf(log_string,"%s %s",log_string, (((TObjString*)token->At(n))->GetString()).Data());
  }
  if(token) delete token;

  return log_string;
}

//______________________________________________________________________________
int
CWB::Toolbox::GetMDCLogSize(TString log) {
//
// return the number of entries contained in the log string
//
// Input: log      - log string 
//

  TObjArray* token = TString(log).Tokenize(' ');
  int size = token->GetEntries();
  delete token;
  return size;
}

//______________________________________________________________________________
TString 
CWB::Toolbox::GetMDCLog(TString log, int pos) {
//
// return the entry at position pos contained in the log string
//
// Input: log      - log string 
//        pos      - entry position in the log string
//

  TObjArray* token = TString(log).Tokenize(' ');
  for(int n=0;n<token->GetEntries();n++) {
    if(n==pos) return (((TObjString*)token->At(n))->GetString()).Data();
  }
  delete token;
  return "";
}

//______________________________________________________________________________
double 
CWB::Toolbox::getMDCShift(mdcshift mshift, double time) {
//
// compute the MDC time shift to be applied to MDC data
// 
//        mshift       - mdc shift structure
//        time         - time 
//
// the GPS range P=[mshift.startMDC,mshift.stopMDC] defines the MDC period to be used 
//
//                 startMDC        stopMDC
//                 ^               ^ 
// ................|xxxxxx P xxxxxx|..............
//
// the period P is replicated starting from mshift.offset
//
//    mshift.offset                      time
//    ^                                  ^
// ...|xxxxxx P xxxxxx|xxxxxx P xxxxxx|xxxxxx P xxxxxx|...
//                                    ^
//                                    mshift.offset+mlength*int(time/mlength)
//
//                 |.. return shift ..|
//
 
  if(mshift.startMDC<=0.) return 0.;
  double mlength = mshift.stopMDC-mshift.startMDC;
  if(mlength<=0) return 0.;
  return mshift.offset+mlength*int(time/mlength)-mshift.startMDC;
}

//______________________________________________________________________________
vector<TString> 
CWB::Toolbox::mergeCWBTrees(int nthreads, TString dir_name, bool simulation, TString odir, 
                            TString label, bool brms, bool bvar, bool bpsm) {
//
// merge list of tree contained in dir_name using threads 
//
//        nthreads     - number of threads
//        simulation   - true/false -> simulation/background
//        odir         - output directory where to store the merged file
//        label        - label used for the output file name
//        brms         - true -> merge rms tree 
//        bvar         - true -> merge variability tree 
//        bpsm         - true -> merge probability skymap tree 
//

  vector<TString> fileList = getFileListFromDir(dir_name,".root","","wave_");
  mergeCWBTrees(nthreads, fileList, simulation, odir, label, brms, bvar, bpsm); 
  return fileList;
}

//______________________________________________________________________________
void 
CWB::Toolbox::mergeCWBTrees(int nthreads, vector<TString> fileList, bool simulation, 
                            TString odir, TString label, bool brms, bool bvar, bool bpsm) {
//
// merge list of tree contained in the root files listed in fileList using threads 
//
//        nthreads     - number of threads
//        ...          - see descriptions of parameters in the n the overload methods  
//

  if(nthreads>MAX_THREADS) {
    cout << "CWB::Toolbox::mergeCWBTrees : Error - nthreads must be <= : " << MAX_THREADS << endl;
    exit(1);
  }

  TThread *thread[MAX_THREADS];
  MergeParms mergeParms[MAX_THREADS];
  vector<TString> fList[MAX_THREADS];
  int lSize[MAX_THREADS];
  TString tLabel[MAX_THREADS];
  vector<TString> waveList;
  vector<TString> mdcList;
  vector<TString> liveList;
  vector<TString> rmsList;
  vector<TString> varList;

  // compute file lists size
  int listSize=1;
  if(fileList.size()<=nthreads) {
    nthreads=1;
    lSize[0]=fileList.size();
  } else {
    if(fileList.size()%nthreads!=0) {
      for(int i=0;i<nthreads;i++) lSize[i]=int(fileList.size()/nthreads);
      lSize[0]+=fileList.size()%nthreads;
    } else {
      for(int i=0;i<nthreads;i++) lSize[i]=int(fileList.size()/nthreads);
    }
  }

  // fill file lists merged files name
  int k=0;
  for(int i=0;i<nthreads;i++) {
    for(int j=0;j<lSize[i];j++) fList[i].push_back(fileList[k++]);
    tLabel[i]=TString::Format("%s.T%d",label.Data(),i);

    waveList.push_back(TString::Format("%s/wave_%s.root",odir.Data(),tLabel[i].Data()));
    if(simulation) {
      mdcList.push_back(TString::Format("%s/mdc_%s.root",odir.Data(),tLabel[i].Data()));
    } else {
      liveList.push_back(TString::Format("%s/live_%s.root",odir.Data(),tLabel[i].Data()));
      if(brms) rmsList.push_back(TString::Format("%s/rms_%s.root",odir.Data(),tLabel[i].Data()));
      if(bvar) varList.push_back(TString::Format("%s/var_%s.root",odir.Data(),tLabel[i].Data()));
    }
  }

  // fill merge structures
  for(int i=0;i<nthreads;i++) {
    mergeParms[i].threadID   = i;
    mergeParms[i].fileList   = fList[i];
    mergeParms[i].simulation = simulation;
    mergeParms[i].odir       = odir;
    mergeParms[i].label      = tLabel[i];
    mergeParms[i].brms       = brms;
    mergeParms[i].bvar       = bvar;
    mergeParms[i].bpsm       = bpsm;
  }

  // create & start threads
  for(int i=0;i<nthreads;i++) {
    char name[128]; sprintf(name,"Merge.T%d",i); 
    thread[i] = new TThread(name, MergeHandle, (void*) &mergeParms[i]);
    thread[i]->Run();
    printf("Starting Thread : %s\n",name);
  }

  // join threads
  for(int i=0;i<nthreads;i++) thread[i]->Join();

  // merge final output merged files
  TString fName;

  fName=TString::Format("wave_%s",label.Data());
  if(waveList.size()) {
    mergeTrees(waveList, "waveburst", odir, fName, true);
    for(int i=0;i<waveList.size();i++) gSystem->Exec("rm "+waveList[i]);
  }
  fName=TString::Format("mdc_%s",label.Data());
  if(mdcList.size()) {
    mergeTrees(mdcList, "mdc", odir, fName, true);
    for(int i=0;i<waveList.size();i++) gSystem->Exec("rm "+mdcList[i]);
  }
  fName=TString::Format("live_%s",label.Data());
  if(liveList.size()) {
    mergeTrees(liveList, "liveTime", odir, fName, true);
    for(int i=0;i<liveList.size();i++) gSystem->Exec("rm "+liveList[i]);
  }
  fName=TString::Format("rms_%s",label.Data());
  if(rmsList.size()) {
    mergeTrees(rmsList, "noise", odir, fName, true);
    for(int i=0;i<rmsList.size();i++) gSystem->Exec("rm "+rmsList[i]);
  }
  fName=TString::Format("var_%s",label.Data());
  if(varList.size()) {
    mergeTrees(varList, "variability", odir, fName, true);
    for(int i=0;i<varList.size();i++) gSystem->Exec("rm "+varList[i]);
  }

  // delete threads
  for(int i=0;i<nthreads;i++) delete thread[i];
}

//______________________________________________________________________________
vector<TString> 
CWB::Toolbox::mergeCWBTrees(TString dir_name, bool simulation, TString odir, 
                            TString label, bool brms, bool bvar, bool bpsm) {
//
// merge trees contained in dir_name
//
//        simulation   - true/false -> simulation/background
//        odir         - output directory where to store the merged file
//        label        - label used for the output file name
//        brms         - true -> merge rms tree 
//        bvar         - true -> merge variability tree 
//        bpsm         - true -> merge probability skymap tree 
//

  vector<TString> fileList = getFileListFromDir(dir_name,".root","","wave_");
  mergeCWBTrees(fileList, simulation, odir, label, brms, bvar, bpsm); 
  return fileList;
}

//______________________________________________________________________________
void 
CWB::Toolbox::mergeCWBTrees(vector<TString> fileList, bool simulation, TString odir, 
                            TString label, bool brms, bool bvar, bool bpsm) {
//
// merge list of tree contained in the root files listed in fileList
//
// Input: fileList     - list of tree root files 
//        simulation   - true/false -> simulation/background
//        odir         - output directory where to store the merged file
//        label        - label used for the output file name
//                       if simulation=true the following files are created
//                       - 'odir'/wave_'label'.root	// detected events
//                       - 'odir'/mdc_'label'.root	// injected events
//                       - 'odir'/merge_'label'.lst	// list of merged files
//                       if simulation=false the following files are created
//                       - 'odir'/wave_'label'.root	// detected events
//                       - 'odir'/live_'label'.root	// live times
//                       - 'odir'/merge_'label'.lst	// list of merged files
//                       - 'odir'/rms_'label'.root	// if brms=true
//                       - 'odir'/var_'label'.root	// if bvar=true
//        brms         - true -> merge rms tree 
//        bvar         - true -> merge variability tree 
//        bpsm         - true -> merge probability skymap tree 
//

  char s[1024];
  char f[1024];
  char cmd[4096];

  char c;
  if (simulation==1) c = 's';
  else c = 'b';
  if (simulation)
    cout << "CWB::Toolbox::mergeCWBTrees - Start merging " 
         << fileList.size() << " files (simulation) ..." << endl; 
  else
    cout << "CWB::Toolbox::mergeCWBTrees - Start merging " 
         << fileList.size() << " files (production) ..." << endl; 

  TChain wav("waveburst");
  if(!bpsm) wav.SetBranchStatus("Psm",false);   // include/exclude from merge the probability skymap
  wav.SetMaxTreeSize(MAX_TREE_SIZE);
  TChain mdc("mdc");
  TChain rms("noise");
  TChain var("variability");
  TChain liv("liveTime");
  int ZombieCnt = 0;

  CWB::History* history = NULL;

  int nfiles=0;
  bool bhistory=true;
  cout << "CWB::Toolbox::mergeCWBTrees - Add file to chain in progress ..." << endl; 
  for(int n=0;n<(int)fileList.size();n++) {
    sprintf(f,"%s",fileList[n].Data());
    //cout << f << endl;

/*  this stuff has been commented out to speed up the merging phase
    TFile f1(f,"READ");

    if(f1.IsZombie()) {
      ZombieCnt++;
      sprintf(cmd,"mv %s %s.zombie",f,f);
      cout<<cmd<<endl;
      gSystem->Exec(cmd);

      // Check if file txt exist
      char txtfile[1024];
      sprintf(txtfile,"%s",f);
      Long_t id,size,flags,mt;
      int estat = gSystem->GetPathInfo(txtfile,&id,&size,&flags,&mt);
      if (estat==0) {
        TString cmd2 = cmd;
        cmd2.ReplaceAll(".root",".txt");
        gSystem->Exec(cmd2.Data());
        cout<<"Zombie file!!! Renamed root and txt files as :"<<f<<".zombie"<<endl;
      }
      continue;
    }
*/

    if(bhistory) {
      TFile f1(f,"READ");
      history=(CWB::History*)f1.Get("history");
      if(history==NULL) history=new CWB::History(); 
      bhistory=false;
    }

    wav.Add(f);
    if (++nfiles%1000==0) cout << "CWB::Toolbox::mergeCWBTrees - " << nfiles 
                               << "/" << fileList.size() << " files" << endl;
    if(c=='s') mdc.Add(f);
    else { liv.Add(f); if(brms) rms.Add(f); if(bvar) var.Add(f);}
  }
 
  sprintf(s,"%s/wave_%s.root",odir.Data(),label.Data());
  cout << "CWB::Toolbox::mergeCWBTrees - Merging " << s << " in progress ..." << endl; 
  wav.Merge(s);

  // update history 
  if(history!=NULL) {
    history->SetHistoryModify(true);
    const CWB::HistoryStage* stage = history->GetStage(CCAST("MERGE"));
    if(stage==NULL) {
      history->AddStage(CCAST("MERGE"));
      if(!history->TypeAllowed(CCAST("WORKDIR"))) history->AddType(CCAST("WORKDIR"));
      char work_dir[512]="";
      sprintf(work_dir,"%s",gSystem->WorkingDirectory());
      history->AddHistory(CCAST("MERGE"), CCAST("WORKDIR"), work_dir);
    }
  }
  // write history to merge file
  if(history!=NULL) {
    TFile fwav(s,"UPDATE");
    history->Write("history");
    fwav.Close();
  }

  if(c=='s') {
    sprintf(s,"%s/mdc_%s.root",odir.Data(),label.Data());
    cout << "CWB::Toolbox::mergeCWBTrees - Merging " << s << " in progress ..." << endl; 
    mdc.Merge(s);
    // write history to merge file
    if(history!=NULL) {
      TFile fmdc(s,"UPDATE");
      history->Write();
      fmdc.Close();
    }
  }
  else {
    sprintf(s,"%s/live_%s.root",odir.Data(),label.Data());
    cout << "CWB::Toolbox::mergeCWBTrees - Merging " << s << " in progress ..." << endl; 
    liv.Merge(s);
    // write history to merge file
    if(history!=NULL) {
      TFile fliv(s,"UPDATE");
      history->Write();
      fliv.Close();
    }
    if(brms) {
      sprintf(s,"%s/rms_%s.root",odir.Data(),label.Data());
      cout << "CWB::Toolbox::mergeCWBTrees - Merging " << s << " in progress ..." << endl; 
      rms.Merge(s);
      if(history!=NULL) {
        TFile frms(s,"UPDATE");
        history->Write();
        frms.Close();
      }
    }
    if(bvar) {
      sprintf(s,"%s/var_%s.root",odir.Data(),label.Data());
      cout << "CWB::Toolbox::mergeCWBTrees - Merging " << s << " in progress ..." << endl; 
      var.Merge(s);
      if(history!=NULL) {
        TFile fvar(s,"UPDATE");
        history->Write();
        fvar.Close();
      }
    }
  }
  if(history!=NULL) delete history;

  cout << "CWB::Toolbox::mergeCWBTrees - Merged files : " << nfiles<<endl;
  //cout << "CWB::Toolbox::mergeCWBTrees - Zombie files : " << ZombieCnt<<endl;
  return;
}

//______________________________________________________________________________
void 
CWB::Toolbox::mergeTrees(vector<TString> fileList, TString treeName, 
                         TString odir, TString ofName, bool bhistory) {
//
// merge list of tree contained in the root files listed in fileList
//
// Input: fileList     - list of tree root files 
//        treeName     - name of tree
//        odir         - output directory where to store the merged file
//        ofName       - name used for the output file name
//                       if ofName not ends with ".root" then ofile_name = 'odir'/'ofName'.root
//                       else ofile_name = 'odir'/'ofName'  
//        bhistory     - if true then the history is copied into the merged file  
//

  char s[1024];
  char f[1024];
  char cmd[4096];

  cout << "CWB::Toolbox::mergeTrees - Start merging " 
       << fileList.size() << " files (simulation) ..." << endl; 

  TChain tree(treeName);
  tree.SetMaxTreeSize(MAX_TREE_SIZE);
  int ZombieCnt = 0;

  CWB::History* history = NULL;

  int nfiles=0;
  cout << "CWB::Toolbox::mergeTrees - Add file to chain in progress ..." << endl; 
  for(int n=0;n<(int)fileList.size();n++) {
    sprintf(f,"%s",fileList[n].Data());
    //cout << f << endl;
    TFile f1(f,"READ");

    if(f1.IsZombie()) {
      ZombieCnt++;
      sprintf(cmd,"mv %s %s.zombie",f,f);
      cout<<cmd<<endl;
      gSystem->Exec(cmd);

      // Check if file txt exist
      char txtfile[1024];
      sprintf(txtfile,"%s",f);
      Long_t id,size,flags,mt;
      int estat = gSystem->GetPathInfo(txtfile,&id,&size,&flags,&mt);
      if (estat==0) {
        TString cmd2 = cmd;
        cmd2.ReplaceAll(".root",".txt");
        gSystem->Exec(cmd2.Data());
        cout<<"Zombie file!!! Renamed root and txt files as :"<<f<<".zombie"<<endl;
      }
      continue;
    }

    if(bhistory) {
      TFile f1(f,"READ");
      history=(CWB::History*)f1.Get("history");
      if(history==NULL) history=new CWB::History(); 
      bhistory=false;
    }

    tree.Add(f);
    if (++nfiles%1000==0) cout << "CWB::Toolbox::mergeTrees - " << nfiles 
                               << "/" << fileList.size() << " files" << endl;
  }

  if(ofName.EndsWith(".root")) { 
    sprintf(s,"%s/%s",odir.Data(),ofName.Data());
  } else {
    sprintf(s,"%s/%s.root",odir.Data(),ofName.Data());
  }
  cout << "CWB::Toolbox::mergeTrees - Merging " << s << " in progress ..." << endl; 
  tree.Merge(s);

  // write history to merge file
  if(history!=NULL) {
    TFile froot(s,"UPDATE");
    history->Write("history");
    froot.Close();
  }
  if(history!=NULL) delete history;

  cout << "CWB::Toolbox::mergeTrees - Merged files : " << nfiles<<endl;
  cout << "CWB::Toolbox::mergeTrees - Zombie files : " << ZombieCnt<<endl;
  return;
}

//______________________________________________________________________________
vector<TString>
CWB::Toolbox::readFileList(TString ifName) {
//
// read file from ifName file
//
// Input: ifName     - input file name
// 
// return a vector of the file names contained in ifName
//

  ifstream in; 
  in.open(ifName.Data(),ios::in);
  if(!in.good()) {
    cout << "CWB::Toolbox::readFileList - Error Opening File : " << ifName.Data() << endl;
    gSystem->Exit(1);
  }

  vector<TString> fileList;

  char istring[8192];
  while (1) {
    in.getline(istring,8192);
    if (!in.good()) break;
    fileList.push_back(istring);
  }

  return fileList;
}

//______________________________________________________________________________
void
CWB::Toolbox::dumpFileList(vector<TString> fileList, TString ofName) {
//
// write file list to ofName file
//
// Input: fileList   - input file list
//        ofName     - output file name
// 

  char ofile_name[1024];
  if(ofName[0]!='/') {
    sprintf(ofile_name,"%s/%s",gSystem->WorkingDirectory(),ofName.Data());
  } else { 
    sprintf(ofile_name,"%s",ofName.Data());
  }
  cout << ofile_name << endl;
  ofstream out;
  out.open(ofile_name,ios::out);
  if (!out.good()) {cout << "CWB::Toolbox::dumpFileList - Error Opening File : " << ofName.Data() << endl;gSystem->Exit(1);}

  for (int i=0;i<(int)fileList.size();i++) {
    out << fileList[i] << endl;
  }

  out.close();

  return;
}

//______________________________________________________________________________
vector<waveSegment>
CWB::Toolbox::getSegList(int jobId, int nIFO, 
                  double segLen, double segMLS, double segEdge,
                  vector<waveSegment> dqList) {
//
// create merged dq cat 1 file list
//
// Input: jobId    - job number
//        nIfO     - detector number
//        segLen   - Segment length [sec]
//        segMLS   - Minimum Segment Length after DQ_CAT1 [sec]
//        segEdge  - wavelet boundary offset [sec]
//        dqList   - dq cat1 list
//
// Return the detector's slag segment ranges associated to SLAG infos & dq cat1 list
//


  if (jobId<1) {cout << "CWB::Toolbox::getSegList - Error : jobId must be > 0" << endl;gSystem->Exit(1);}

  vector<waveSegment> jobList = getJobList(dqList, segLen, segMLS, segEdge);        // get  standard job list

  if (jobId>(int)jobList.size()) {
    cout << "CWB::Toolbox::getSegList - Error : jobId is > max jobs " << jobList.size() << endl;
    gSystem->Exit(1);
  }

  vector<waveSegment> detSegs(nIFO);
  for(int i=0;i<nIFO;i++) detSegs[i] = jobList[jobId-1];

  jobList.clear();

  return detSegs;
}

//______________________________________________________________________________
char* 
CWB::Toolbox::readFile(TString ifName) {
//
// return the file content in a char buffer
//
// Input: ifName  - file name
//

  char* fileBuffer=NULL;

  Long_t id,size,flags,mt;
  int estat = gSystem->GetPathInfo(ifName.Data(),&id,&size,&flags,&mt);
  if (estat==0) {

    ifstream in; 
    in.open(ifName.Data(),ios::in);
    if(!in.good()) {
      cout << "CWB::Toolbox::readFile - Error Opening File : " << ifName.Data() << endl;
      gSystem->Exit(1);
    }

    fileBuffer = new char[2*size+1];
    bzero(fileBuffer,(2*size+1)*sizeof(char));

    char istring[8192];
    int fileLength = 0;
    while (1) {
      in.getline(istring,8192);
      if (!in.good()) break;
      //cout << istring << endl;
      int len = strlen(istring);
      istring[len]=0x0a;
      strncpy(fileBuffer+fileLength,istring,len+1);
      fileLength += len+1;
    }
  }

  return fileBuffer;
}

//______________________________________________________________________________
int
CWB::Toolbox::setCuts(TString ifName, TString idir, TString odir, TString trname, TString cuts, TString olabel) {
// 
// apply selection cuts to trname tree in ifname root file 
// and create an output root file with the selected tree entries 
// return selected entries
//
// ifname 	: input root file name
// idir		: input dir
// trname	: tree name
// cuts		: selection cuts (Ex: netcc[2]>0.5)
// odir		: output dir
// olabel	: output label used in the output file name (Ex: cc2_gt_0d5)
//		  ofname.root = ifname.Colabel.root
//

  if(cuts=="") {
    cout << "CWB::Toolbox::setCuts - cuts is empty : nothing to do" << endl;
    gSystem->Exit(1);
  } 

  TString ifname = idir+"/"+ifName;

  // open input root file
  cout<<"Opening File : " << ifname.Data() << endl;
  TFile ifile(ifname);
  if (ifile.IsZombie()) {
    cout << "CWB::Toolbox::setCuts - Error opening file " << ifname.Data() << endl;
    gSystem->Exit(1);
  } 
  TTree *itree = (TTree*)ifile.Get(trname.Data());
  if (itree==NULL) {
    cout << "CWB::Toolbox::setCuts - tree " << trname 
         << " is not present in file " << ifname.Data() << endl;
    gSystem->Exit(1);
  } 
  Int_t isize = (Int_t)itree->GetEntries();
  cout << "tree size : " << isize << endl;
  if (isize==0) {
    cout << "CWB::Toolbox::setCuts - tree " << trname 
         << " is empty in file " << ifname.Data() << endl;
    gSystem->Exit(1);
  } 

  // check cuts
  TTreeFormula formula("trCuts", cuts.Data(), itree);
  int err = formula.Compile(cuts.Data()); 
  if(err) {
    cout << "CWB::Toolbox::setCuts - wrong input cuts " << cuts << endl; 
    return -1;
  } 

  // get selected entries
  itree->SetEstimate(isize);
  itree->Draw("Entry$",cuts.Data(),"goff",isize);
  int nentries = itree->GetSelectedRows();
  cout << "CWB::Toolbox::setCuts - selected entries : " 
       << nentries << "/" << itree->GetEntries() << endl; 

  // create output root file name
  TString ofname = odir+"/"+ifName;
  ofname.ReplaceAll(".root",TString(".C_")+olabel+".root");
  cout << "CWB::Toolbox::setCuts - Create cuts file : " << ofname << endl;
  bool overwrite = checkFile(ofname,true);
  if(!overwrite) gSystem->Exit(0);

  // create output root file with selection cuts
  TFile ofile(ofname,"RECREATE");
  TTree* otree = itree->CopyTree(cuts.Data());
  otree->SetDirectory(&ofile);
  otree->Write();
  ofile.Close();

  // create setCuts history 
  CWB::History* history = (CWB::History*)ifile.Get("history");
  if(history==NULL) history=new CWB::History(); 
  if(history!=NULL) {
    TList* stageList = history->GetStageNames();    // get stage list
    TString pcuts="";
    for(int i=0;i<stageList->GetSize();i++) {       // get previous cuts 
      TObjString* stageObjString = (TObjString*)stageList->At(i);
      TString stageName = stageObjString->GetString();
      char* stage = const_cast<char*>(stageName.Data());
      if(stageName=="CUTS") pcuts=history->GetHistory(stage,const_cast<char*>("PARAMETERS"));
    }
    // update history
    history->SetHistoryModify(true);
    if(!history->StageAlreadyPresent(CCAST("CUTS"))) history->AddStage(CCAST("CUTS"));
    if(!history->TypeAllowed(CCAST("PARAMETERS"))) history->AddType(CCAST("PARAMETERS"));
    if(!history->TypeAllowed(CCAST("WORKDIR"))) history->AddType(CCAST("WORKDIR"));
    char work_dir[1024]="";
    sprintf(work_dir,"%s",gSystem->WorkingDirectory());
    history->AddHistory(CCAST("CUTS"), CCAST("WORKDIR"), work_dir);
    TString _cuts = (pcuts!="") ? pcuts+" "+cuts : cuts;
    history->AddHistory(CCAST("CUTS"), CCAST("PARAMETERS"), CCAST(_cuts.Data()));
    char logmsg[2048]; sprintf(logmsg,"Apply Selection Cuts : %s",cuts.Data());
    history->AddLog(CCAST("CUTS"), CCAST(logmsg));
  }

  // close input root file
  ifile.Close(); 	

  // write setCuts history to merge+cuts file
  if(history!=NULL) {
    TFile ofile(ofname,"UPDATE");
    history->Write("history");
    ofile.Close();
    delete history;
  }

  return nentries;
}

//______________________________________________________________________________
int
CWB::Toolbox::setIFAR(TString ifName, TString idir, TString odir, TString trname, 
                     TString sels, TString farFile, int irho, TString olabel, bool inclusive) {
// 
// add a new ifar (inverse false alarm rate) leaf to the selected entries in the merged wave root file
// for each selected entry a new leaf is created "ifar" which value is obtained from the farFile
// farFile is a text file which contains a corrispondence between rho and far (from background)
// return selected entries
//
// ifname 	: input root file name
// idir		: input dir
// odir		: output dir
// trname	: tree name
// sels		: selection sels (Ex: netcc[2]>0.5)
// farFile	: far file (text file with a list entries : rho far)
// irho         ; rho index 0/1 : rho[0]/rho[1] to be used for rho to far association
// olabel	: output label used in the output file name (Ex: cc2_gt_0d5)
//		  ofname.root = ifname.Folabel.root
// inclusive    : true(def)/false=inclusive/exclusive : used to compute the ifar 
//

  if(irho!=0 && irho!=1) {
    cout << "CWB::Toolbox::setIFAR - bad irho value : irho must be 0/1" << endl;
    gSystem->Exit(1);
  } 

  if(sels=="") {
    cout << "CWB::Toolbox::setIFAR - sels is empty : nothing to do" << endl;
    gSystem->Exit(1);
  } 

  // check if farFile exist
  bool exist = checkFile(farFile,false);
  if(!exist) gSystem->Exit(0);
  // read farFile
  vector<double>  xrho,xfar;
  int xsize = GetStepFunction(farFile, xrho, xfar);
  // remove elements with xfar=0
  vector<double>  trho,tfar;
  for(int i=0;i<xsize;i++) {
    if(xfar[i]!=0) {
      trho.push_back(xrho[i]);
      tfar.push_back(xfar[i]);
    }
  }
  xrho=trho; xfar=tfar; 
  // insert a new element at the beginning with xrho=-DBL_MAX and xfar[xfar[0]]
  std::vector<double>::iterator it;
  it =  xrho.begin();  xrho.insert(it,-DBL_MAX);
  it =  xfar.begin();  xfar.insert(it,xfar[0]);
  // insert a new element at the end with xrho=DBL_MAX/2 and xfar[xfar[xfar.size()-1]]
  xrho.push_back(DBL_MAX/2);
  xfar.push_back(xfar[xfar.size()-1]);
  xsize=xfar.size();
  // create step graph used to get far from rho
  TGraph gr;
  int cnt=0;
  for(int i=1;i<xsize;i++) {
    double dx = (xrho[i]-xrho[i-1])/100000.;
    gr.SetPoint(cnt++,xrho[i]-dx,xfar[i-1]);
    gr.SetPoint(cnt++,xrho[i]+dx,xfar[i]);
  }

  TString ifname = idir+"/"+ifName;

  // open input root file
  cout<<"Opening File : " << ifname.Data() << endl;
  TFile ifile(ifname);
  if (ifile.IsZombie()) {
    cout << "CWB::Toolbox::setIFAR - Error opening file " << ifname.Data() << endl;
    gSystem->Exit(1);
  } 
  TTree *itree = (TTree*)ifile.Get(trname.Data());
  if (itree==NULL) {
    cout << "CWB::Toolbox::setIFAR - tree " << trname 
         << " is not present in file " << ifname.Data() << endl;
    gSystem->Exit(1);
  } 
  Int_t isize = (Int_t)itree->GetEntries();
  itree->SetEstimate(isize);
  cout << "tree size : " << isize << endl;
  if (isize==0) {
    cout << "CWB::Toolbox::setIFAR - tree " << trname 
         << " is empty in file " << ifname.Data() << endl;
    gSystem->Exit(1);
  } 

  // check sels
  TTreeFormula fsel("sels", sels.Data(), itree);
  int err = fsel.Compile(sels.Data()); 
  if(err) {
    cout << "CWB::Toolbox::setIFAR - wrong input sels " << sels << endl; 
    return -1;
  } 

  // create output root file name
  TString ofname = odir+"/"+ifName;
  ofname.ReplaceAll(".root",TString(".S_")+olabel+".root");
  cout << "CWB::Toolbox::setIFAR - Create sels file : " << ofname << endl;
  bool overwrite = checkFile(ofname,true);
  if(!overwrite) gSystem->Exit(0);

  // create output root file with selection sels
  TFile ofile(ofname,"RECREATE");
  TTree *otree = (TTree*)itree->CloneTree(0);
  otree->SetMaxTreeSize(MAX_TREE_SIZE);

  // add ifar,bin leaf to otree 
  TBranch* branch;
  bool ifar_exists=false;
  bool bin_exists=false;
  TIter next(itree->GetListOfBranches());
  while ((branch=(TBranch*)next())) {
    if(TString("ifar").CompareTo(branch->GetName())==0) ifar_exists=true;
    if(TString("bin").CompareTo(branch->GetName())==0) bin_exists=true;
  }
  next.Reset();
  float ifar=0;
  if (ifar_exists) itree->SetBranchAddress("ifar",&ifar);
  else            otree->Branch("ifar",&ifar,"ifar/F");
  string* bin = new string();
  if (bin_exists) itree->SetBranchAddress("bin",&bin);
  else            otree->Branch("bin",bin);

  // get selected entries
  itree->Draw("Entry$",sels.Data(),"goff",isize);
  double* entry = itree->GetV1(); 
  int nsel = itree->GetSelectedRows();
  cout << "CWB::Toolbox::setIFAR - selected entries : "
       << nsel << "/" << itree->GetEntries() << endl;
  // write ifar for the selected entries
  float rho[2];
  itree->SetBranchAddress("rho",rho);
  for(int i=0;i<nsel;i++) {
    ifar=0; 
    itree->GetEntry(int(entry[i]));
    double far  = gr.Eval(rho[irho]);
    double xifar = far>0 ? 1./far : -1;
    if(inclusive) {			// inclusive
      if(xifar>ifar) {
        ifar=xifar;
        *bin=olabel.Data();
      }
    } else {				// exclusive
      ifar = xifar;
      *bin = olabel.Data();
    }
    otree->Fill();
  }
  // get not selected entries
  itree->Draw("Entry$",TString::Format("!(%s)",sels.Data()).Data(),"goff",isize);
  int _nsel = itree->GetSelectedRows();
  double* _entry = itree->GetV1(); 
  cout << "CWB::Toolbox::setIFAR - not selected entries : "
       << _nsel << "/" << itree->GetEntries() << endl;
  // write ifar for the not selected entries
  for(int i=0;i<_nsel;i++) {
    itree->GetEntry(int(_entry[i]));
    if(!ifar_exists) ifar=0;	// initialize ifar
    if(!bin_exists) *bin="";	// initialize bin
    otree->Fill();
  }
  delete bin;
  otree->Write();
  ofile.Close();

  // create setIFAR history 
  CWB::History* history = (CWB::History*)ifile.Get("history");
  if(history==NULL) history=new CWB::History(); 
  if(history!=NULL) {
    TList* stageList = history->GetStageNames();    // get stage list
    TString psels="";
    char* pfile=NULL;
    for(int i=0;i<stageList->GetSize();i++) {       // get previous sels 
      TObjString* stageObjString = (TObjString*)stageList->At(i);
      TString stageName = stageObjString->GetString();
      char* stage = const_cast<char*>(stageName.Data());
      if(stageName=="IFAR") psels=history->GetHistory(stage,CCAST("PARAMETERS"));
      if(stageName=="IFAR") pfile=history->GetHistory(stage,CCAST("FARFILE"));
    }
    // update history
    history->SetHistoryModify(true);
    if(!history->StageAlreadyPresent(CCAST("IFAR"))) history->AddStage(CCAST("IFAR"));
    if(!history->TypeAllowed(CCAST("PARAMETERS"))) history->AddType(CCAST("PARAMETERS"));
    if(!history->TypeAllowed(CCAST("WORKDIR"))) history->AddType(CCAST("WORKDIR"));
    if(!history->TypeAllowed(CCAST("FARFILE"))) history->AddType(CCAST("FARFILE"));
    char work_dir[1024]="";
    sprintf(work_dir,"%s",gSystem->WorkingDirectory());
    history->AddHistory(CCAST("IFAR"), CCAST("WORKDIR"), work_dir);
    TString _sels = (psels!="") ? psels+"\n"+sels : sels;
    history->AddHistory(CCAST("IFAR"), CCAST("PARAMETERS"), CCAST(_sels.Data()));
    char logmsg[2048]; sprintf(logmsg,"Created ifar on the bin set (%s) specified by the selection : %s",olabel.Data(),sels.Data());
    history->AddLog(CCAST("IFAR"), CCAST(logmsg));
    char* buffer = readFile(farFile);                // read farFile file
    if(buffer!=NULL) {
      char* farBuffer=NULL;
      if(pfile!=NULL) {
        farBuffer = new char[strlen(pfile)+strlen(buffer)+farFile.Sizeof()+10];
        sprintf(farBuffer,"%s\n%s\n%s\n",pfile,farFile.Data(),buffer);
      } else {
        farBuffer = new char[strlen(buffer)+farFile.Sizeof()+10];
        sprintf(farBuffer,"%s\n%s\n",farFile.Data(),buffer);
      } 
      history->AddHistory(CCAST("IFAR"), CCAST("FARFILE"), farBuffer);
      delete [] buffer;
      delete [] farBuffer;
    }
  }

  // close input root file
  ifile.Close(); 	

  // write setIFAR history to merge+sels file
  if(history!=NULL) {
    TFile ofile(ofname,"UPDATE");
    history->Write("history");
    ofile.Close();
    delete history;
  }

  return nsel;
}

//______________________________________________________________________________
int
CWB::Toolbox::setChunk(TString ifName, TString idir, TString odir, TString trname, int chunk) {
// 
// add a new chunk (data chunk id) leaf to the selected entries in the merged wave root file
// return error status
//
// ifname 	: input root file name
// idir		: input dir
// odir		: output dir
// trname	: tree name
// chunk	: chunk id
//

  if(chunk<=0) {
    cout << "CWB::Toolbox::setChunk - bad chunk chunk value : chunk must be >0" << endl;
    gSystem->Exit(1);
  } 

  TString ifname = idir+"/"+ifName;

  // open input root file
  cout<<"Opening File : " << ifname.Data() << endl;
  TFile ifile(ifname);
  if (ifile.IsZombie()) {
    cout << "CWB::Toolbox::setChunk - Error opening file " << ifname.Data() << endl;
    gSystem->Exit(1);
  } 
  TTree *itree = (TTree*)ifile.Get(trname.Data());
  if (itree==NULL) {
    cout << "CWB::Toolbox::setChunk - tree " << trname 
         << " is not present in file " << ifname.Data() << endl;
    gSystem->Exit(1);
  } 
  Int_t isize = (Int_t)itree->GetEntries();
  itree->SetEstimate(isize);
  cout << "tree size : " << isize << endl;
  if (isize==0) {
    cout << "CWB::Toolbox::setChunk - tree " << trname 
         << " is empty in file " << ifname.Data() << endl;
    gSystem->Exit(1);
  } 

  // create output root file name
  TString ofname = odir+"/"+ifName;
  char schunk[32];sprintf(schunk,"chunk%d",chunk);
  ofname.ReplaceAll(".root",TString(".K_")+schunk+".root");
  cout << "CWB::Toolbox::setChunk - Create sels file : " << ofname << endl;
  bool overwrite = checkFile(ofname,true);
  if(!overwrite) gSystem->Exit(0);

  // create output root file with selection sels
  TFile ofile(ofname,"RECREATE");
  TTree *otree = (TTree*)itree->CloneTree(0);
  otree->SetMaxTreeSize(MAX_TREE_SIZE);

  // add chunk leaf to otree 
  TBranch* branch;
  bool chunk_exists=false;
  TIter next(itree->GetListOfBranches());
  while ((branch=(TBranch*)next())) {
    if(TString("chunk").CompareTo(branch->GetName())==0) chunk_exists=true;
  }
  next.Reset();
  if (chunk_exists) itree->SetBranchAddress("chunk",&chunk);
  else              otree->Branch("chunk",&chunk,"chunk/I");

  // get selected entries
  itree->Draw("Entry$","","goff",isize);
  double* entry = itree->GetV1(); 
  int nsel = itree->GetSelectedRows();
  cout << "CWB::Toolbox::setChunk - selected entries : "
       << nsel << "/" << itree->GetEntries() << endl;
  // write chunk for the selected entries
  for(int i=0;i<nsel;i++) {
    itree->GetEntry(int(entry[i]));
    otree->Fill();
  }
  otree->Write();
  ofile.Close();

  // create setChunk history 
  CWB::History* history = (CWB::History*)ifile.Get("history");
  if(history==NULL) history=new CWB::History(); 
  if(history!=NULL) {
    TList* stageList = history->GetStageNames();    // get stage list
    TString psels="";
    char* pfile=NULL;
    for(int i=0;i<stageList->GetSize();i++) {       // get previous sels 
      TObjString* stageObjString = (TObjString*)stageList->At(i);
      TString stageName = stageObjString->GetString();
      char* stage = const_cast<char*>(stageName.Data());
      if(stageName=="CHUNK") psels=history->GetHistory(stage,CCAST("PARAMETERS"));
      if(stageName=="CHUNK") pfile=history->GetHistory(stage,CCAST("FARFILE"));
    }
    // update history
    history->SetHistoryModify(true);
    if(!history->StageAlreadyPresent(CCAST("CHUNK"))) history->AddStage(CCAST("CHUNK"));
    if(!history->TypeAllowed(CCAST("PARAMETERS"))) history->AddType(CCAST("PARAMETERS"));
    if(!history->TypeAllowed(CCAST("WORKDIR"))) history->AddType(CCAST("WORKDIR"));
    char work_dir[1024]="";
    sprintf(work_dir,"%s",gSystem->WorkingDirectory());
    history->AddHistory(CCAST("CHUNK"), CCAST("WORKDIR"), work_dir);
    char logmsg[2048]; sprintf(logmsg,"Created chunk leaf : %d",chunk);
    history->AddLog(CCAST("CHUNK"), CCAST(logmsg));
  }

  // close input root file
  ifile.Close(); 	

  // write setChunk history to merge+sels file
  if(history!=NULL) {
    TFile ofile(ofname,"UPDATE");
    history->Write("history");
    ofile.Close();
    delete history;
  }

  return 0;
}

//_______________________________________________________________________________________
void 
CWB::Toolbox::setUniqueEvents(TString ifwave, TString ofwave, int nIFO, int pp_irho) {                 
//                                                                                       
// select from the root ifwave file a unique reconstructed event for each injected       
// signal. The signal is selected from the reconstructed events which have the same      
// injected gps time, the selected one has the maximum rho                               
//                                                                                       
// ifwave       : name of simulation input wave file                                    
// ofwave       : name of simulation output wave file (unique reconstructed events)     
// nifo         : number of detectors                                                   
// nfactor      : number of simulation factors                                          
// factors      : array of simulation factors                                           
// pp_irho      : index used for rho = rho[pp_irho]                                     
//                                                                                       

  cout << "CWB::Toolbox::setUniqueEvents - Create unique events file : " << ofwave << endl;
  bool overwrite = checkFile(ofwave,true);
  if(!overwrite) gSystem->Exit(0);

  TFile* iwroot = TFile::Open(ifwave);
  if(iwroot==NULL||!iwroot->IsOpen()) {
    cout << "CWB::Toolbox::setUniqueEvents - Error : input file " << ifwave << " not found" <<  endl;
    exit(1);                                                                           
  }                                                                                    

  TTree* iwtree = (TTree*) iwroot->Get("waveburst");
  if(iwtree==NULL) {                                
    cout << "CWB::Toolbox::setUniqueEvents - Error : tree waveburst not found in file " << ifwave <<  endl;
    exit(1);                                                                                 
  }                                                                                          
  cout << endl;                                                                              
  cout << "CWB::Toolbox::setUniqueEvents : number of input events : " << iwtree->GetEntries() << endl;     
  cout << endl;                                                                              

  TFile* owroot = new TFile(ofwave,"RECREATE");
  if(owroot->IsZombie()) {                               
    cout << "CWB::Toolbox::setUniqueEvents - Error : output file " << ofwave << " not opened" <<  endl;
    exit(1);                                                                                  
  }                                                                                           

  TTree *owtree = (TTree*)iwtree->CloneTree(0);
  owtree->SetMaxTreeSize(MAX_TREE_SIZE);          

  // -------------------------------------------------------------------------
  // add new leave to owtree (fsize : number of fragments)
  // -------------------------------------------------------------------------
  TString leaf_name="fsize";
  TBranch* branch;
  bool replaceLeaf=false;
  TIter next(iwtree->GetListOfBranches());
  while ((branch=(TBranch*)next())) {
    if (leaf_name.CompareTo(branch->GetName())==0) {
      cout << "fragment leaf [" << leaf_name << "] already applied" << endl;
      char answer[1024];
      strcpy(answer,"");
      do {
        cout << "CWB::Toolbox::setUniqueEvents : Do you want to overwrite the previous leaf ? (y/n) ";
        cin >> answer;
      } while ((strcmp(answer,"y")!=0)&&(strcmp(answer,"n")!=0));
      if (strcmp(answer,"y")==0) {
         replaceLeaf=true;
      } else {
         gSystem->Exit(-1);
      }
    }
  }
  next.Reset();
  int fsize=1;
  if (replaceLeaf) owtree->SetBranchAddress("fsize",&fsize);
  else             owtree->Branch("fsize",&fsize,"fsize/I");
                                                 
  char sel[1024];                                
  sprintf(sel,"Entry$:rho[%i]:time[%i]:factor",pp_irho,nIFO);
  int iwsize=iwtree->GetEntries();       
  iwtree->Draw(sel,"","goff",iwsize);                               
  double* entry  = iwtree->GetV1();                               
  double* rho    = iwtree->GetV2();                                 
  double* time   = iwtree->GetV3();                                 
  double* factor = iwtree->GetV4();                                 

  // extract factors
  std::map <float, float> mfactors;    
  for(int i=0;i<iwsize;i++) {          // loop over the tree 
    iwtree->GetEntry(i);
    mfactors[factor[i]]=factor[i];                                              
  }          
  vector<float> factors;
  std::map<float, float>::iterator iter;
  for (iter=mfactors.begin(); iter!=mfactors.end(); iter++) {
    factors.push_back(iter->second);
  }
  int nfactor=factors.size(); 

  char cuts[1024];                                 
  for(int i=0;i<nfactor;i++) {          // loop over the factors 
    sprintf(cuts,"abs(factor-%g)/%g<0.0001",factors[i],factors[i]); 
    iwtree->Draw(sel,cuts,"goff",iwsize);                               
    int size=iwtree->GetSelectedRows();                           
    int unique_evt=0;                   // number of unique events for factor i
    if(size==1) {                                                              
      iwtree->GetEntry(entry[0]);                                              
      owtree->Fill();                                                          
      unique_evt++;
    } else if(size>1) {                                                        
      int* index = new int[size];                                              
      TMath::Sort(size,time,index,kFALSE);      // sort events in time                       
      int    j=index[0];                                                                     
      float  rhomax=rho[j];                                                                  
      int    imax=j;                                                                         
      double time1=time[j];                                                                  
      double time2=0.;                                                                       
      for(int k=1; k<size; k++) {                                                            
        j=index[k];                                                                          
        time2=time[j];                                                                       
        if(abs(time1-time2)>0.001 || k==size-1) {                                                         
          iwtree->GetEntry(entry[imax]);                                                     
          owtree->Fill();               // store the event with max rho[pp_irho]             
          imax=j;
          rhomax=rho[j];
          unique_evt++;
          fsize=1;
        } else {
          if(rho[j]>rhomax) {rhomax=rho[j];imax=j;}
          fsize++;
        }
        time1=time2;
      }
      delete index;
    }
    cout << " factor id = " << i << "\t reconstructed events : "
         << iwtree->GetSelectedRows() << "\t unique events : " 
         << unique_evt << "\tfactor = " << factors[i] << endl;
  }
  cout << endl;
  owtree->Write();
  owroot->Close();
  owroot->Close();

  // create setCuts history 
  CWB::History* history = (CWB::History*)iwroot->Get("history");
  if(history==NULL) history=new CWB::History(); 
  if(history!=NULL) {
    TList* stageList = history->GetStageNames();    // get stage list
    TString pcuts="";
    for(int i=0;i<stageList->GetSize();i++) {       // get previous cuts 
      TObjString* stageObjString = (TObjString*)stageList->At(i);
      TString stageName = stageObjString->GetString();
      char* stage = const_cast<char*>(stageName.Data());
      if(stageName=="CUTS") pcuts=history->GetHistory(stage,const_cast<char*>("PARAMETERS"));
    }
    // update history
    history->SetHistoryModify(true);
    if(!history->StageAlreadyPresent(CCAST("CUTS"))) history->AddStage(CCAST("CUTS"));
    if(!history->TypeAllowed(CCAST("PARAMETERS"))) history->AddType(CCAST("PARAMETERS"));
    if(!history->TypeAllowed(CCAST("WORKDIR"))) history->AddType(CCAST("WORKDIR"));
    char work_dir[1024]="";
    sprintf(work_dir,"%s",gSystem->WorkingDirectory());
    history->AddHistory(CCAST("CUTS"), CCAST("WORKDIR"), work_dir);
    TString cuts = (pcuts!="") ? pcuts+" (U)" : "(U)";
    history->AddHistory(CCAST("CUTS"), CCAST("PARAMETERS"), CCAST(cuts.Data()));
    char logmsg[2048]; sprintf(logmsg,"Apply Selection Cuts : %s","(U)");
    history->AddLog(CCAST("CUTS"), CCAST(logmsg));
  }

  // close input root file
  iwroot->Close(); 	

  // write setCuts history to merge+cuts file
  if(history!=NULL) {
    TFile ofile(ofwave,"UPDATE");
    history->Write("history");
    ofile.Close();
    delete history;
  }

  return;
}

//_______________________________________________________________________________________
int
CWB::Toolbox::CombineCBC(vector<TString> ifwave, vector<TString> ifmdc, TString ofwave, TString ofmdc,
                         int nIFO, float fthr, TString msearch, TString infos, int lag, int slag, float ifarthr) {
//
// Combine IMBHB & BBH searches
//
// ifwave       : name of simulation/background input wave files
//                ifwave[0] = IMBHB
//                ifwave[1] = BBH
// ifmdc        : simulation -> name of simulation input mdc files
//                ifmdc[0] = IMBHB
//                ifmdc[1] = BBH
//                background -> ifmdc.size()=0
// ofwave       : name of simulation output wave file (reconstructed events)
// ofmdc        : name of simulation output mdc file  (injected events)
// nIFO         : number of detectors
// fthr         : frequency threshold used to combine searches
//                IMBHB band -> freq<=fthr
//                BBH   band -> freq >fthr
// msearch      : is the master search. Used only by history -> save the history of the master search. Available values are IMBHB,BBH
// infos        : user bookkeeping informations to be saved to the history output root files
//
// lag,slag	: only for background (ifmdc.size()=0). Used to select the data where to apply the combined search
//                lag=slag=0 -> zero lag
// ifarthr      : Used to select events with ifar>ifarthr, ifarthr is in years 
//                ifarthr must be >=0. The condition with ifar>0 select all events that pass the post-production cuts
//
// NOTES
//
// This function implements the combine statistic for the IMBHB and BBH searches
// It can be used for simulation and background zero/non_zero lag
// Simulation:
//   The same injections (from XML files) should be made for both searches.
//     The algorithm select only the common injections (compare the MDC log strings used for simulation).
//   The wave files must contain unique reconstructed events (no more than 1 event for each injection).
//     The algorithm exit if entries are not unique.
//     Before to use this funcion apply the CWB::Toolbox::setUniqueEvents
//   The common injected events are written to the output MDC file ofmdc
//   The combined reconstructed events are written to the output WAVE file ofwave
// Background:
//   The combined reconstructed events are written to the output WAVE file ofwave
//
// Combined IFAR -> the IFAR of any event reconstructed by the BBH and IMBHB searches is defined in the table below
//
// -----------------------------------------------------------------------------------------------------------------
// Case 1 -> IMBHB band (IMBHB BBH), BBH band (         ) -> IFAR = IFAR(IMBHB)                 TRIALS=1
// Case 2 -> IMBHB band (         ), BBH band (IMBHB BBH) -> IFAR = IFAR(BBH  )                 TRIALS=1
// Case 3 -> IMBHB band (IMBHB    ), BBH band (      BBH) -> IFAR = max IFAR(IMBHB,BBH)/2       TRIALS=2
// Case 4 -> IMBHB band (      BBH), BBH band (IMBHB    ) -> IFAR = max IFAR(IMBHB,BBH)/2       TRIALS=2
// Case 5 -> IMBHB band (IMBHB    ), BBH band (         ) -> IFAR = max IFAR(IMBHB,BBH)         TRIALS=1
// Case 6 -> IMBHB band (      BBH), BBH band (         ) -> IFAR = max IFAR(IMBHB,BBH)/2       TRIALS=2
// Case 7 -> IMBHB band (         ), BBH band (IMBHB    ) -> IFAR = max IFAR(IMBHB,BBH)/2       TRIALS=2
// Case 8 -> IMBHB band (         ), BBH band (      BBH) -> IFAR = max IFAR(IMBHB,BBH)         TRIALS=1
// -----------------------------------------------------------------------------------------------------------------
//
// Note: 'Case' is stored into the output wave file in the 'combine' parameter
//       The selected search IMBHB/BBH is stored into the output wave file in the 'search' parameter 
//
 
  #define YEAR		(24.*3600.*365.)
  #define TRIALS	2		// number of trials used in the combine algorithm
  #define TIME_UNCERTAINTY    0.1	// (sec) used to select entries belonging to the same event
                                	// the injected time could be not the same for IMBHB and BBH searches
                                	// because it is the median of the whitened injected  event

  // check if output files exist
  cout << "CWB::Toolbox::CombineCBC - Create combine wave output file : " << ofwave << endl;
  bool overwrite = checkFile(ofwave,true);
  if(!overwrite) gSystem->Exit(1);

  // check if ifarthr>=0 
  if(ifarthr<0) {
    cout << "CWB::Toolbox::CombineCBC - Error : ifarthr must be >=0" << endl;
    exit(1);
  }

  // check if msearch is defined
  if(msearch!="IMBHB" && msearch!="BBH") {
    cout << "CWB::Toolbox::CombineCBC - Error : msearch not defined, valid values are: IMBHB,BBH" << endl;
    exit(1);
  }
  int msearch_id = (msearch=="IMBHB") ? 0 : 1;     // get master search id

  bool simulation = ifmdc.size()>0 ? true : false; // case simulation/zero_lag

  // open input wave files
  TFile* iwroot[2];
  TTree* iwtree[2];
  if(ifwave.size()!=2) {cout << "CWB::Toolbox::CombineCBC - Error in input wave files declaration" <<  endl; exit(1);}
  for(int n=0;n<2;n++) {
    iwroot[n] = TFile::Open(ifwave[n]);
    if(iwroot==NULL||!iwroot[n]->IsOpen()) {
      cout << "CWB::Toolbox::CombineCBC - Error : input wave file " << ifwave[n] << " not found" <<  endl;
      exit(1);
    }

    iwtree[n] = (TTree*) iwroot[n]->Get("waveburst");
    if(iwtree[n]==NULL) {
      cout << "CWB::Toolbox::CombineCBC - Error : tree waveburst not found in file " << ifwave[n] <<  endl;
      exit(1);
    }
    cout << endl;
    cout << "CWB::Toolbox::CombineCBC : number of input wave events : " << iwtree[n]->GetEntries()
         << " ( " << ((n==0)?"IMBHB":"BBH") << " )" << endl;
    cout << endl;
  }

  // open input mdc files
  TFile* imroot[2];
  TTree* imtree[2];
  if(simulation && ifmdc.size()!=2) {cout << "CWB::Toolbox::CombineCBC - Error in input mdc files declaration" <<  endl; exit(1);}
  if(simulation) {
    for(int n=0;n<2;n++) {
      imroot[n] = TFile::Open(ifmdc[n]);
      if(imroot==NULL||!imroot[n]->IsOpen()) {
        cout << "CWB::Toolbox::CombineCBC - Error : input mdc file " << ifmdc[n] << " not found" <<  endl;
        exit(1);
      }

      imtree[n] = (TTree*) imroot[n]->Get("mdc");
      if(imtree[n]==NULL) {
        cout << "CWB::Toolbox::CombineCBC - Error : tree mdc not found in file " << ifmdc[n] <<  endl;
        exit(1);
      }
      cout << endl;
      cout << "CWB::Toolbox::CombineCBC : number of input mdc events : " << imtree[n]->GetEntries()
           << " ( " << ((n==0)?"IMBHB":"BBH") << " )" << endl;
      cout << endl;
    }
  }

  // create temporary output wave files
  TString tfwave[2];
  TFile*  twroot[2];
  TTree*  twtree[2];
  if(simulation) {
    for(int n=0;n<2;n++) {
      tfwave[n] = ofwave;
      tfwave[n].ReplaceAll(".root",TString::Format("_tmp%d.root",n));
      twroot[n] = new TFile(tfwave[n],"RECREATE");
      if(twroot[n]->IsZombie()) {
        cout << "CWB::Toolbox::CombineCBC - Error : temporary output wave file " << tfwave[n] << " not opened" <<  endl;
        exit(1);
      }
      twtree[n] = (TTree*)iwtree[n]->CloneTree(0);
      twtree[n]->SetMaxTreeSize(MAX_TREE_SIZE);
    }
  }

  // create temporary output mdc files
  TString tfmdc[2];
  TFile*  tmroot[2];
  TTree*  tmtree[2];
  if(simulation) {
    for(int n=0;n<2;n++) {
      tfmdc[n] = ofmdc;
      tfmdc[n].ReplaceAll(".root",TString::Format("_tmp%d.root",n));
      tmroot[n] = new TFile(tfmdc[n],"RECREATE");
      if(tmroot[n]->IsZombie()) {
        cout << "CWB::Toolbox::CombineCBC - Error : temporary output mdc file " << tfmdc[n] << " not opened" <<  endl;
        exit(1);
      }
      tmtree[n] = (TTree*)imtree[n]->CloneTree(0);
      tmtree[n]->SetMaxTreeSize(MAX_TREE_SIZE);
    }
  }

  // -------------------------------------------------------------------------
  // add leaves 'search' and 'combine' to owtree
  // search -> IMBHB, BBH
  // combine   -> [1,8] is defined in the description of CombineCBC function
  // -------------------------------------------------------------------------
  // if 'search' or 'combine' leaves exist the combine is aborted -> the input root files are already combined
  for(int n=0;n<2;n++) {
    bool check_ifar=false;
    TBranch* branch;
    TIter next(iwtree[n]->GetListOfBranches());
    while ((branch=(TBranch*)next())) {
      if(TString(branch->GetName())=="search") {
        cout << "CWB::Toolbox::CombineCBC - Error : 'search' leaf already exist"; gSystem->Exit(1);
      }
      if(TString(branch->GetName())=="combine") {
        cout << "CWB::Toolbox::CombineCBC - Error : 'combine' leaf already exist"; gSystem->Exit(1);
      }
      if(TString(branch->GetName())=="ifar") check_ifar=true;
    }
    next.Reset();
    if(!check_ifar) {
      cout << "CWB::Toolbox::CombineCBC - Error : ifar is not defined in wave root file " << ifwave[n] << endl; 
      cout << endl;
      gSystem->Exit(1);
    }
  }
  string* SEARCH = new string();		// new tree leaf used to store the selected combine search
  int COMBINE=1;				// new tree leaf used to store the selected combine case [1,8]
  float IFAR;					// variable used to update the IFAR
  if(simulation) {
    for(int n=0;n<2;n++) {
      twtree[n]->Branch("search",SEARCH);
      twtree[n]->Branch("combine",&COMBINE,"combine/I");
      twtree[n]->SetBranchAddress("ifar",&IFAR);
    }
  }

  // define vectors used for zero lag case, contains the combine event parameters for zero lag background
  vector<double> vGPS;
  vector<string> vSEARCH;
  vector<float>  vIFAR;
  vector<int>    vCOMBINE;

  // extract factors
  vector<float> factors;
  std::map <float, float> mfactors;
  if(simulation) {
    for(int n=0;n<2;n++) {
      int size=imtree[n]->GetEntries();
      imtree[n]->Draw("factor","","goff",size);
      double* factor = imtree[n]->GetV1();
      for(int i=0;i<size;i++) {          // loop over the tree
        mfactors[factor[i]]=factor[i];
      }
    }
  }
  std::map<float, float>::iterator iter;
  for (iter=mfactors.begin(); iter!=mfactors.end(); iter++) {
    factors.push_back(iter->second);
  }
  int nfactor= simulation ? factors.size() : 1;	

  // loop over the factors
  int nREC=0;	// total number of recostructed events in combined search
  int nINJ=0;	// total number of injected     events in combined search;
  for(int i=0;i<nfactor;i++) {
    char wcuts[1024]="";
    char mcuts[1024]="";
    if(simulation) {		// simulation -> select factor where to apply the combined search
      sprintf(wcuts,"(abs(factor-%g)/%g<0.0001)",factors[i],factors[i]);
      sprintf(mcuts,"(abs(factor-%g)/%g<0.0001)",factors[i],factors[i]);
    } else {			// background -> select lag/slag where to apply the combined search
      sprintf(wcuts,"lag[%d]==0&&slag[%d]==0",nIFO,nIFO);
    }
    if(ifarthr>=0) {char cuts[1024];sprintf(cuts,"%s&&(ifar>%e)",wcuts,ifarthr*YEAR);strcpy(wcuts,cuts);}
    cout << endl << " Factor " << i+1 << "\twave cuts -> " << wcuts << endl;
    if(simulation) cout << " Factor " << i+1 << "\tmdc  cuts -> " << mcuts << endl;
    cout << endl;
    // extract Entry,time,ifar,frequency from wave tree
    int     iwsize[2];		// total number of entries in the wave tree
    double* iwentry[2]; 	// wave tree entry index
    double* iwtime[2];  	// injection GPS time of the first ifo
    double* iwifar[2];  	// ifar
    double* iwfreq[2];  	// reconstructed median frequency
    int     cwsize[2]={0,0};  	// number of wave entries with factor=factors[i]
    char wsel[1024];
    if(simulation) {
      sprintf(wsel,"Entry$:time[%i]:ifar:frequency[0]",nIFO);
    } else {			// zero lag
      sprintf(wsel,"Entry$:time[0]:ifar:frequency[0]");
    }
    for(int n=0;n<2;n++) {
      iwsize[n]=iwtree[n]->GetEntries();
      iwtree[n]->Draw(wsel,wcuts,"goff",iwsize[n]);
      iwentry[n] = iwtree[n]->GetV1();
      iwtime[n]  = iwtree[n]->GetV2();
      iwifar[n]  = iwtree[n]->GetV3();
      iwfreq[n]  = iwtree[n]->GetV4();
      cwsize[n]  = iwtree[n]->GetSelectedRows();
    }

    // extract Entry,time,factor from mdc tree
    int     imsize[2];		// total number of entries in the mdc tree
    double* imentry[2]; 	// mdc tree entry index
    double* imtime[2];  	// injection GPS time of the first ifo
    int     cmsize[2]={0,0};  	// number of mdc entries with factor=factors[i]
    char msel[1024];
    sprintf(msel,"Entry$:time[0]");
    if(simulation) cout << " Factor " << i+1 << "\tmdc  selections-> " << msel << endl << endl;
    if(simulation) {
      for(int n=0;n<2;n++) {
        imsize[n]=imtree[n]->GetEntries();
        imtree[n]->Draw(msel,mcuts,"goff",imsize[n]);
        imentry[n]  = imtree[n]->GetV1();
        imtime[n]   = imtree[n]->GetV2();
        cmsize[n]   = imtree[n]->GetSelectedRows();
      }
    }

    // merge arrays extracted from wave,mdc tree
    int size = cwsize[0]+cwsize[1]+cmsize[0]+cmsize[1];
    if(size==0) {cout << "CWB::Toolbox::CombineCBC - Error : no entries in wave/mdc tree" <<  endl; return 1;}

    int*   index = new int[size];       // used by sort function
    int*   entry = new int[size];	// entry id for root tree
    int*   ttype = new int[size];	// tree type.   0/1 -> wave/mdc
    int*   stype = new int[size];       // search type. 0/1 -> IMBHB/BBH
    double* time = new double[size];	// GPS time of the injection for the first ifo
    double* ifar = new double[size];	// inverse FAR in sec
    double* freq = new double[size];	// reconstructed median event frequency
    string* log  = new string[size];	// mdc log string, contains the input parameters for simulation
    string* ilog = new string();	// used by SetBranchAddress to exctract log from ttree entry

    int k=0;
    for(int n=0;n<2;n++) {		// loop over search type 0/1 -> IMBHB/BBH
      // fill wave
      for(int i=0;i<cwsize[n];i++) {
        entry[k] = iwentry[n][i];
        ttype[k] = 0;			// wave
        stype[k] = n;
        time[k]  = iwtime[n][i];
        ifar[k]  = iwifar[n][i];
        freq[k]  = iwfreq[n][i];
        log[k]   = "";
        k++;
      }
      if(!simulation) continue; 	// skip mdc if wave files are background
      // fill mdc
      imtree[n]->SetBranchAddress("log",&ilog);
      for(int i=0;i<cmsize[n];i++) {
        entry[k] = imentry[n][i];
        ttype[k] = 1;			// mdc
        stype[k] = n;
        time[k]  = imtime[n][i];
        ifar[k]  = 0;
        freq[k]  = 0;
        imtree[n]->GetEntry(entry[k]);
        log[k]   = *ilog;
        k++;
      }
    }

    if(!simulation) {
      cout << endl << "List of uncombined events found in " << msearch << " search" << endl << endl;
      for(int k=0;k<size;k++) {
        cout << k+1 << "\tsearch -> " << (stype[k]==0 ? "IMBHB" : "BBH") 
             << "\tGPS -> " << std::setprecision(14) << time[k]
             << "\tFREQ -> " << std::setprecision(4) << freq[k]
             << "\tIFAR -> " << std::setprecision(6)  << ifar[k]/YEAR << " (years)" << endl;
      }
      cout << endl;
    }

    // sort events wrt time
    TMath::Sort(size,time,index,kFALSE);

    // array to store the number of entries in wave/mdc
    int    esize[2]={0,0};			// wave/mdc -> 0/1

    // arrays used to store the event with the same GPS time. position 0/1 -> IMBHB/BBH
    float  eifar[2]={0,0};			// inverse far
    float  efreq[2]={0,0};			// reconstructed median frequency
    string  elog[2]={"",""};			// MDC log string used for simulation
    int    ewidx[2]={-1,-1};			// ttree entry wave index
    int    emidx[2]={-1,-1};			// ttree entry mdc  index

    // events found in wave/mdc/cases
    int    nWAVE=0;				// reconstructed events
    int    nMDC=0;				// injected events
    int    nX[8]={0,0,0,0,0,0,0,0};		// events found in combine cases

    // fill first entry
    int    j=index[0];
    double ptime=time[j];                       // previous GPS event time
    double ctime=0.;                            // current  GPS event time

    if(ttype[j]==0) {
      eifar[stype[j]]=ifar[j];
      efreq[stype[j]]=freq[j];
      ewidx[stype[j]]=j;
    }
    if(ttype[j]==1) {
      elog[stype[j]]=log[j];
      emidx[stype[j]]=j;
    }
    esize[ttype[j]]++;

    int iSEARCH=0;				// store the search type IMBHB/BBH -> 0/1
    int kstart = size==1 ? 0 : 1;		// manage the case with size=1;  
    int kstop  = size+1;			// manage the last event
    for(int k=kstart; k<kstop; k++) {  		// loop over the IMBHB,BBH wave/mdc events
      // event with the same GPS time belong to the same event
      // event is combined if the same injection has been performed in both searches: IMBHB and BBH -> esize[1]=2 and elog[0]=elog[1]
      // esize[0]=0 -> event not reconstructed in both searches
      // esize[0]=1 -> event reconstructed in only in one search -> iSEARCH=stype[j]
      // esize[0]=2 -> event reconstructed on both searches
      if(k<size) {
        j=index[k];
        ctime=time[j];
      } else {					// last event
        ctime=ptime+1;				// just to make ctime>ptime
      }
      if(fabs(ctime-ptime)>TIME_UNCERTAINTY) {	// New event, we process the previous

        if(esize[0]>2) {
           cout << "CWB::Toolbox::CombineCBC - Error : ";
           if(simulation) cout << "number of reconstructed events per injection is > 2" <<  endl; 
           else           cout << "number of reconstructed events within " << TIME_UNCERTAINTY << " (sec) is > 2" <<  endl; 
           cout << endl;
           exit(1);
        }
        if(esize[1]>ifmdc.size()) {
          cout << "CWB::Toolbox::CombineCBC - Error : number of injections per event is > " <<  ifmdc.size() << endl;
          exit(1);
        }

        // case simulation: if injection is present in BBH and IMBHB then write mdc entry in output wave/mdc files
        if(esize[1]==ifmdc.size()) {

          // MDC log IMBHB & BBH must be the same !!!
          if(elog[0].compare(elog[1])!=0) {
            cout << endl << "CWB::Toolbox::CombineCBC - Error : mdc log string are different" << endl << endl;
            cout << "IMBHB MDC log string: " << endl << elog[0].c_str() << endl << endl;
            cout << "BBH   MDC log string: " << endl << elog[1].c_str() << endl << endl;
            exit(1);
          }

          if(simulation) nMDC++;

          COMBINE=0;
          if(esize[0]==1) {	// event has been reconstructed only in one search
            iSEARCH=stype[ewidx[ewidx[0]!=-1?0:1]];					// store search type into output var 0/1 -> IMBHB/BBH
            if(iSEARCH==0 && efreq[0]<=fthr) {COMBINE=5;IFAR=eifar[0];}			// found by IMBHB search in IMBHB band
            if(iSEARCH==1 && efreq[1]<=fthr) {COMBINE=6;IFAR=eifar[1]/TRIALS;}		// found by BBH   search in IMBHB band
            if(iSEARCH==0 && efreq[0] >fthr) {COMBINE=7;IFAR=eifar[0]/TRIALS;}		// found by IMBHB search in BBH   band
            if(iSEARCH==1 && efreq[1] >fthr) {COMBINE=8;IFAR=eifar[1];}			// found by BBH   search in BBH   band
            nWAVE++;
          }
          if(esize[0]==2) {	// event has been reconstructed in both searches
            if(efreq[0]<=fthr && efreq[1]>fthr) COMBINE=3;
            if(efreq[1]<=fthr && efreq[0]>fthr) COMBINE=4;
            iSEARCH = eifar[0]>eifar[1] ? 0 : 1; 					// search type is the one with greater ifar
            IFAR = eifar[0]>eifar[1] ? eifar[0]/TRIALS : eifar[1]/TRIALS; 		// both freqs are in the wrong search band
            if(efreq[0]<=fthr && efreq[1]<=fthr) {COMBINE=1;IFAR=eifar[0];iSEARCH=0;}	// both freqs are <= fthr (IMBHB band)
            if(efreq[0] >fthr && efreq[1] >fthr) {COMBINE=2;IFAR=eifar[1];iSEARCH=1;}	// both freqs are  > fthr (BBH   band)
            nWAVE++;
          }
          nX[COMBINE-1]++;

          // save entry into mdc output file
          if(esize[1]==2) {
            imtree[iSEARCH]->GetEntry(entry[emidx[iSEARCH]]);
            tmtree[iSEARCH]->Fill();
            nINJ++;
          }
          // if event has been reconstructed in IMBHB or/and BBH searches then save entry into wave output file
          if(esize[0]>0) {
            *SEARCH = iSEARCH==0 ? "IMBHB" : "BBH";
            if(simulation) {
              iwtree[iSEARCH]->GetEntry(entry[ewidx[iSEARCH]]);
              twtree[iSEARCH]->Fill();
            } else {			// zero lag
              vGPS.push_back(ptime);
              vSEARCH.push_back(*SEARCH);
              vIFAR.push_back(IFAR);
              vCOMBINE.push_back(COMBINE);
	    }
            nREC++;
          }
        }
        // reset counters 
        esize[0]= 0;esize[1]= 0;
        ewidx[0]=-1;ewidx[1]=-1;
        emidx[0]=-1;emidx[1]=-1;
      }

      // update combine event parameters
      ptime=ctime;
      if(ttype[j]==0) {		// wave
        eifar[stype[j]]=ifar[j];
        efreq[stype[j]]=freq[j];
        ewidx[stype[j]]=j;
      }
      if(ttype[j]==1) {		// mdc
        elog[stype[j]]=log[j];
        emidx[stype[j]]=j;
      }
      esize[ttype[j]]++;
    }

    // print selected events statistic
    cout << endl;
    cout << "----------------------------------------------------------------------------------------------" << endl;
    cout << " factor id = " << i+1 << endl;
    cout << "             " << " -> injected      events IMBHB/BBH         : " << cmsize[0] << "/" << cmsize[1] << endl;
    cout << "             " << " -> reconstructed events IMBHB/BBH         : " << cwsize[0] << "/" << cwsize[1] << endl;
    cout << "             " << " -> combined reconstructed/injected events : " << nWAVE << "/" << nMDC << endl;
    cout << "----------------------------------------------------------------------------------------------" << endl;
    cout << " Case 1 -> TRIALS=1 IMBHB band (IMBHB BBH), BBH band (         )  = " << nX[0] << endl;
    cout << " Case 2 -> TRIALS=1 IMBHB band (         ), BBH band (IMBHB BBH)  = " << nX[1] << endl;
    cout << " Case 3 -> TRIALS=2 IMBHB band (IMBHB    ), BBH band (      BBH)  = " << nX[2] << endl;
    cout << " Case 4 -> TRIALS=2 IMBHB band (      BBH), BBH band (IMBHB    )  = " << nX[3] << endl;
    cout << " Case 5 -> TRIALS=1 IMBHB band (IMBHB    ), BBH band (         )  = " << nX[4] << endl;
    cout << " Case 6 -> TRIALS=2 IMBHB band (      BBH), BBH band (         )  = " << nX[5] << endl;
    cout << " Case 7 -> TRIALS=2 IMBHB band (         ), BBH band (IMBHB    )  = " << nX[6] << endl;
    cout << " Case 8 -> TRIALS=1 IMBHB band (         ), BBH band (      BBH)  = " << nX[7] << endl;
    cout << "----------------------------------------------------------------------------------------------" << endl;
    cout << endl;

    // release memory
    delete[] index;
    delete[] entry;
    delete[] ttype;
    delete[] stype;
    delete[] time;
    delete[] ifar;
    delete[] freq;
    delete[] log;
    delete   ilog;
  }
  delete SEARCH;
  cout << endl;
  cout << " Total number of REC/INJ = " << nREC << "/" << nINJ << endl;
  cout << endl;

  if(simulation) {

    // write temporary output mdc/wave file
    for(int n=0;n<2;n++) {
      twroot[n]->cd();
      twtree[n]->Write();
      twroot[n]->Close();
    }
    for(int n=0;n<2;n++) {
      tmroot[n]->cd();
      tmtree[n]->Write();
      tmroot[n]->Close();
    }

    // merge temporary wave/mdc root files into final root files
    TList *owlist = new TList;
    for(int n=0;n<2;n++) {
      twroot[n] = new TFile(tfwave[n],"READ");		// open output temporary wave root file
      twtree[n] = (TTree *)twroot[n]->Get("waveburst");
      owlist->Add(twtree[n]);
    }
    TFile* owroot = new TFile(ofwave,"RECREATE");		// open output wave root file
    if(owroot->IsZombie()) {
      cout << "CWB::Toolbox::CombineCBC - Error : output wave file " << ofwave << " not opened" <<  endl;
      exit(1);
    } 
    gROOT->cd(); 
    owroot->cd(); 
    TTree* owtree = TTree::MergeTrees(owlist);
    owtree->SetName("waveburst");
    owtree->Write();
    owroot->Close();
    for(int n=0;n<2;n++) twroot[n]->Close();

    // create output mdc file
    TList *omlist = new TList;
    for(int n=0;n<2;n++) {
      tmroot[n] = new TFile(tfmdc[n],"READ");		// open output temporary mdc root file
      tmtree[n] = (TTree *)tmroot[n]->Get("mdc");
      omlist->Add(tmtree[n]);
    }
    TFile* omroot = new TFile(ofmdc,"RECREATE");	// open output mdc root file	
    if(omroot->IsZombie()) {
      cout << "CWB::Toolbox::CombineCBC - Error : output mdc file " << ofmdc << " not opened" <<  endl;
      exit(1);
    }
    gROOT->cd(); 
    TTree* omtree = TTree::MergeTrees(omlist);
    omroot->cd(); 
    omtree->SetName("mdc");
    omtree->Write();
    omroot->Close();
    for(int n=0;n<2;n++) twroot[n]->Close();

    // remove temporary mdc/wave root files 
    for(int n=0;n<2;n++) {
      char cmd[1024]; sprintf(cmd,"rm %s",tfwave[n].Data());  
      //cout << cmd << endl;
      gSystem->Exec(cmd);
    }
    for(int n=0;n<2;n++) {
      char cmd[1024]; sprintf(cmd,"rm %s",tfmdc[n].Data());  
      //cout << cmd << endl;
      gSystem->Exec(cmd);
    }
  } else {	// write combined output wave file, only zero lag events are updated

    iwroot[msearch_id]->cd();
    int isize = iwtree[msearch_id]->GetEntries();

    TFile* owroot = new TFile(ofwave,"RECREATE");		// open output wave root file
    if(owroot->IsZombie()) {
      cout << "CWB::Toolbox::CombineCBC - Error : output wave file " << ofwave << " not opened" <<  endl;
      exit(1);
    } 
    TTree* owtree = (TTree*)iwtree[msearch_id]->CloneTree(0);
    owtree->SetMaxTreeSize(MAX_TREE_SIZE);

    float   oIFAR;                                 // variable used to update the IFAR
    int     oCOMBINE=0;                            // new tree leaf used to store the selected combine case
    string* oSEARCH = new string();                // new tree leaf used to store the selected combine search
    float   iIFAR;
    double* iTIME   = new double[2*nIFO]; 
    float*  iFREQ   = new float[nIFO]; 
    float*  iLAG    = new float[nIFO+1]; 
    float*  iSLAG   = new float[nIFO+1]; 

    owtree->SetBranchAddress("ifar",&oIFAR);
    owtree->Branch("combine",&oCOMBINE,"combine/I");
    owtree->Branch("search",oSEARCH);
    iwtree[msearch_id]->SetBranchAddress("ifar",&iIFAR);
    iwtree[msearch_id]->SetBranchAddress("time",iTIME);
    iwtree[msearch_id]->SetBranchAddress("frequency",iFREQ);
    iwtree[msearch_id]->SetBranchAddress("lag",iLAG);
    iwtree[msearch_id]->SetBranchAddress("slag",iSLAG);

    cout << "List of combined events found in " << msearch << " search" << endl;
    cout << endl;
    int k=1;
    for(int i=0;i<isize;i++) {
      iwtree[msearch_id]->GetEntry(i);
      if(iLAG[nIFO]==0 && iSLAG[nIFO]==0) {
        for(int j=0;j<vGPS.size();j++) {
          if(fabs(iTIME[0]-vGPS[j])<TIME_UNCERTAINTY) {
            *oSEARCH = vSEARCH[j];  
            oCOMBINE = vCOMBINE[j];  
            oIFAR    = vIFAR[j];  
            cout << k++ << "\tcombine -> " << oCOMBINE << "\tsearch -> " << oSEARCH->c_str()
                 << "\tGPS -> " << std::setprecision(14) << vGPS[j]
                 << "\tFREQ -> " << std::setprecision(4) << iFREQ[0]
                 << "\toIFAR -> " << std::setprecision(6)  << oIFAR/YEAR << " (years)" << endl;
          }
        }
      } else {
        *oSEARCH  = msearch;  
        oCOMBINE  = 0;  
        oIFAR     = iIFAR;  
      }
      owtree->Fill(); 
    }
    cout << endl;

    owroot->cd();
    owtree->Write();
    owroot->Close();

    delete   oSEARCH; 
    delete[] iTIME; 
    delete[] iLAG; 
    delete[] iSLAG; 
  }

  // create CombineCBC history
  // note: only the history of master search is saved
  CWB::History* history = (CWB::History*)iwroot[msearch_id]->Get("history");
  if(history==NULL) history=new CWB::History();
  if(history!=NULL) {
    TList* stageList = history->GetStageNames();    // get stage list
    TString pcuts="";
    for(int i=0;i<stageList->GetSize();i++) {       // get previous cuts
      TObjString* stageObjString = (TObjString*)stageList->At(i);
      TString stageName = stageObjString->GetString();
      char* stage = const_cast<char*>(stageName.Data());
      if(stageName=="CUTS") pcuts=history->GetHistory(stage,const_cast<char*>("PARAMETERS"));
    }
    // update history
    history->SetHistoryModify(true);
    if(!history->StageAlreadyPresent(CCAST("COMBINE"))) history->AddStage(CCAST("COMBINE"));
    if(!history->TypeAllowed(CCAST("PARAMETERS"))) history->AddType(CCAST("PARAMETERS"));
    if(!history->TypeAllowed(CCAST("WORKDIR"))) history->AddType(CCAST("WORKDIR"));
    char work_dir[1024]="";
    sprintf(work_dir,"%s",gSystem->WorkingDirectory());
    history->AddHistory(CCAST("COMBINE"), CCAST("WORKDIR"), work_dir);
    TString cuts = (pcuts!="") ? pcuts+" ("+infos+")" : "("+infos+")";
    history->AddHistory(CCAST("CUTS"), CCAST("PARAMETERS"), CCAST(cuts.Data()));
    char iparms[1024];
    sprintf(iparms,"IMBHB_ifwave=%s, BBH_ifwave=%s, fthr=%0.2f",ifwave[0].Data(),ifwave[1].Data(),fthr);
    history->AddHistory(CCAST("COMBINE"), CCAST("PARAMETERS"), CCAST(iparms));
    char logmsg[2048]; sprintf(logmsg,"Combine IMBHB and BBH searches : %s","(J)");
    history->AddLog(CCAST("COMBINE"), CCAST(logmsg));
  }

  // close input root file
  for(int n=0;n<2;n++) iwroot[n]->Close();

  // write CombineCBC history into wave/mdc files
  if(history!=NULL) {
    TFile owfile(ofwave,"UPDATE");
    history->Write("history");
    owfile.Close();
    if(simulation) {
      TFile omfile(ofmdc,"UPDATE");
      history->Write("history");
      omfile.Close();
    }
    delete history;
  }

  return 0;
}

//______________________________________________________________________________
int
CWB::Toolbox::setVeto(TString ifName, TString idir, TString odir, vector<TString> ifos, int nVDQF, dqfile* VDQF, 
                      int nDQF, dqfile* DQF, double segLen, double segMLS, double segEdge) {
//
// apply Data Quality and veto to waveburst in ifname root file
// and create an output root file adding a flag for each Data Quality
// return selected entries
//
// ifname       : input root file name
// idir         : input dir
// odir         : output dir
// ifos         : detector lists
// nVDQF        : number of Data Quality files
// VDQF         : Data Quality files structures
// nDQF         : 
// DQF          :
// segLen       : segment job lenght
// segMLS       : minimum segment lenght after DQ2
// segEdge      : scratch lenght  
//


  int nIFO = ifos.size();

  CWB::History* history = NULL;
  bool bhistory=true;

  // extract from VDQF the set of declared detector's names not included in the ifos list
  vector<TString> wifos;
  for(int j=0;j<nVDQF;j++) {
    if(VDQF[j].cat==CWB_CAT0) continue;
    if(VDQF[j].cat==CWB_CAT1) continue;
    bool bifo=true;
    for(int i=0;i<(int)wifos.size();i++) if(wifos[i].CompareTo(VDQF[j].ifo)==0) bifo=false;
    for(int i=0;i<nIFO;i++) if(ifos[i].CompareTo(VDQF[j].ifo)==0) bifo=false;
    if(bifo) wifos.push_back(VDQF[j].ifo);
  }
  for(int i=0;i<(int)wifos.size();i++) cout << wifos[i].Data() << " ";
  cout << endl;
//gSystem->Exit(0);

  // build XDQF structure used for internal's computations
  int nXDQF=0;
  dqfile* XDQF = new dqfile[nVDQF*nIFO]; 
  // netifo contains network with dimension > nIFO (for exclusion vetoes)
  int nNET=TMath::Factorial(5);
  vector<TString>* xifos = new vector<TString>[nNET];
  for(int j=0;j<nVDQF;j++) {
    if(VDQF[j].cat==CWB_CAT0) continue;
    if(VDQF[j].cat==CWB_CAT1) continue;
    bool bifo=false;
    for(int i=0;i<nIFO;i++) if(ifos[i].CompareTo(VDQF[j].ifo)==0) bifo=true;
    // change CAT if ifo not presents into ifos array
    if(!bifo) {  
      // insert CWB_EXC   (WARNING : miss cases with det>3 !!!)
      for(int i=0;i<nIFO;i++) {
        for(int k=0;k<(int)wifos.size();k++) {
          if(wifos[k].CompareTo(VDQF[j].ifo)==0) {
            XDQF[nXDQF]=VDQF[j];
            strcpy(XDQF[nXDQF].ifo,ifos[i].Data());
            XDQF[nXDQF].cat=CWB_EXC;
            xifos[nXDQF]=ifos;
            xifos[nXDQF].push_back(wifos[k]);
            nXDQF++;
          }
        }
      }
    } else {
      xifos[nXDQF]=ifos;
      XDQF[nXDQF]=VDQF[j];
      nXDQF++;
    }
  } 
  if(nXDQF==0) {
    cout << "CWB::Toolbox::setVeto - No veto found : setVeto terminated" << endl;
    gSystem->Exit(1);
  }
  for(int j=0;j<nXDQF;j++) {
    cout << "XDQF[" << j << "] " << CWB_CAT_NAME[XDQF[j].cat] << "_" << XDQF[j].ifo << " NETWORK : ";
    for(int i=0;i<(int)xifos[j].size();i++) cout << xifos[j][i].Data() << " ";
    cout << endl; 
  }

  // get  standard job list
  vector<waveSegment> cat1List = readSegList(nDQF, DQF, CWB_CAT1);
  vector<waveSegment> jobList  = getJobList(cat1List, segLen, segMLS, segEdge);  
  //vector<waveSegment> cat2List = readSegList(nDQF, DQF, CWB_CAT2);

  // build veto file label
  char ifo_label[10][1024];
  char veto_label[10][1024];
  int nset=0;
  for(int n=0;n<nXDQF;n++) {
    TString veto_name = CWB_CAT_NAME[XDQF[n].cat];
    bool bnew=true;
    for(int m=0;m<nset;m++) if(veto_name.CompareTo(veto_label[m])==0) bnew=false;
    if(bnew) strcpy(veto_label[nset++],veto_name.Data());
  }
  for(int m=0;m<nset;m++) strcpy(ifo_label[m],"");
  for(int n=0;n<nXDQF;n++) {
    TString veto_name = CWB_CAT_NAME[XDQF[n].cat];
    for(int m=0;m<nset;m++) {
      if(veto_name.CompareTo(veto_label[m])==0) {
        if(!TString(ifo_label[m]).Contains(XDQF[n].ifo)) {
          sprintf(ifo_label[m],"%s%s",ifo_label[m],XDQF[n].ifo);
        }
      }
    }
  }
  char veto_file_label[1024]="";
  for(int m=0;m<nset;m++) {
    if(strlen(veto_file_label)==0) {
      sprintf(veto_file_label,".V_%s%s",veto_label[m],ifo_label[m]);
    } else {
      sprintf(veto_file_label,"%s_%s%s",veto_file_label,veto_label[m],ifo_label[m]);
    }
  }
  strcpy(veto_file_label,TString(veto_file_label).ReplaceAll("veto_","").Data());
  strcpy(veto_file_label,TString(veto_file_label).ReplaceAll("1","").Data());
  sprintf(veto_file_label,"%s.root",veto_file_label);
  cout << "veto file label : " << veto_file_label << endl;

  TString efname = odir+"/"+ifName;
  efname.ReplaceAll(".root",veto_file_label);
  bool overwrite = checkFile(efname,true);
  if(!overwrite) gSystem->Exit(0);

  TString ifname = idir+"/"+ifName;
  TString ofname = odir+"/"+ifName;
  ofname.ReplaceAll(".root",".root.tmp.1");

//cout.precision(14);
//for(int k=0;k<jobList.size();k++)
//  cout << k << " " << jobList[k].start << " " <<  jobList[k].stop << endl;
//gSystem->Exit(0);

  vector<waveSegment> vlist;
  for(int n=0;n<nXDQF;n++) {
    // find ifo index
    int iIFO=-1;
    for(int i=0;i<nIFO;i++) if(ifos[i].CompareTo(XDQF[n].ifo)==0) iIFO=i;
    if (iIFO==-1) continue;  // sky ifo not presents into ifos array

//cout << n << " ifile " << ifname.Data() << endl;
//cout << n << " ofile " << ofname.Data() << endl;
//cout << endl;

    // get veto list
    if(XDQF[n].cat==CWB_CAT2) {
      // for cat2 is necessary to merge all CAT2 in XDQF
      // RDQF contains only ifo which are presents in ifos array
      int nRDQF=0;
      dqfile* RDQF = new dqfile[nVDQF]; 
      for(int j=0;j<nVDQF;j++) {
        for(int k=0;k<nIFO;k++) if(ifos[k].CompareTo(VDQF[j].ifo)==0) RDQF[nRDQF++]=VDQF[j];
      } 
//      for(int j=0;j<nRDQF;j++) cout << "RDQF[" << j << "] " << CWB_CAT_NAME[RDQF[j].cat] << "_" << RDQF[j].ifo.Data() << endl;
//gSystem->Exit(0);
      vlist=readSegList(nRDQF, RDQF, CWB_CAT2);
      delete [] RDQF;
    } else if(XDQF[n].cat==CWB_EXC) {
      // extract list (higher net conf) for exclusion veto 
      // RDQF contains only ifo which are presents in ifos array
      int nRDQF=0;
      dqfile* RDQF = new dqfile[nVDQF]; 
      for(int j=0;j<nVDQF;j++) {
        for(int k=0;k<(int)xifos[n].size();k++) if(xifos[n][k].CompareTo(VDQF[j].ifo)==0) RDQF[nRDQF++]=VDQF[j];
      } 
//      for(int j=0;j<nRDQF;j++) cout << "EXC -> RDQF[" << j << "] " << CWB_CAT_NAME[RDQF[j].cat] << "_" << RDQF[j].ifo.Data() << endl;
      vlist=readSegList(nRDQF, RDQF, CWB_CAT2);
//cout.precision(14);
//for(int i=0;i<vlist.size();i++) {
//  cout << i << " " << vlist[i].start << " " << vlist[i].stop <<  endl;
//}
//cout << CWB_CAT_NAME[XDQF[n].cat].Data() << " list size " << vlist.size() << endl;
    } else {
      vlist=readSegList(XDQF[n]);
    }
    double vlist_time = getTimeSegList(vlist);
/*
if(XDQF[n].cat==CWB_CAT3) {
cout.precision(14);
for(int i=0;i<vlist.size();i++) {
  cout << i << " " << vlist[i].start << " " << vlist[i].stop <<  endl;
}
cout << CWB_CAT_NAME[XDQF[n].cat].Data() << " list size " << vlist.size() << endl;
gSystem->Exit(1);
}
*/
    TString veto_name = CWB_CAT_NAME[XDQF[n].cat]+"_"+XDQF[n].ifo;
    cout << veto_name << " list duration " << int(vlist_time) << " list size " << vlist.size() << endl;

    // -------------------------------------------------------------------------
    // open root file
    // -------------------------------------------------------------------------
    cout<<"Opening BKG File : " << ifname.Data() << endl;
    TFile ifile(ifname);
    if (ifile.IsZombie()) {
      cout << "CWB::Toolbox::setVeto - Error opening file " << ifname.Data() << endl;
      gSystem->Exit(1);
    } 
    TTree *itree = (TTree*)ifile.Get("waveburst");
    if (itree==NULL) {
      cout << "CWB::Toolbox::setVeto - tree waveburst is not present in file " << ifname.Data() << endl;
      gSystem->Exit(1);
    } 
    Int_t isize = (Int_t)itree->GetEntries();
    cout << "tree size : " << isize << endl;
    if (isize==0) {
      cout << "CWB::Toolbox::setVeto - tree waveburst is empty in file " << ifname.Data() << endl;
      gSystem->Exit(1);
    } 
    itree->SetEstimate(isize);
    char selection[1024];sprintf(selection,"time[%d]",iIFO);
    itree->Draw(selection,"","goff",isize);
    double* time = itree->GetV1(); 
    // check if time is NaN
    for(int i=0;i<isize;i++) {
      if(TMath::IsNaN(time[i])) {
        cout.precision(14);
        cout << "CWB::Toolbox::setVeto - tree waveburst file " << ifname.Data() << " contains time NaN in ifo=" << ifos[iIFO] << " time=" << time[i] << endl;
        gSystem->Exit(1);
      }
    }
    // sort list
    Int_t *id = new Int_t[isize];
    bool *bveto = new bool[isize];
    TMath::Sort((int)isize,time,id,false);
    // add dummy veto to vlist to permit a full trigger loop
    waveSegment SEG;
    SEG.start=time[id[isize-1]];
    SEG.stop=time[id[isize-1]];
    vlist.push_back(SEG);
    int vsize = vlist.size();
/*
for(int h=0;h<isize;h++) {
cout << h << " " << iIFO << " " << time[id[h]] << endl;
}
gSystem->Exit(1);
*/

    if(bhistory) {
      history=(CWB::History*)ifile.Get("history");
      if(history==NULL) history=new CWB::History(); 
      bhistory=false;
    }

    // -------------------------------------------------------------------------
    // add new leave to itree
    // -------------------------------------------------------------------------
    TBranch* branch;
    bool replaceVeto=false;
    TIter next(itree->GetListOfBranches());
    while ((branch=(TBranch*)next())) {
      if (veto_name.CompareTo(branch->GetName())==0) {
        cout << "Veto [" << veto_name << "] already applied" << endl;
        char answer[1024];
        strcpy(answer,"");
        do {
          cout << "Do you want to overwrite the previous spurious ? (y/n) ";
          cin >> answer;
        } while ((strcmp(answer,"y")!=0)&&(strcmp(answer,"n")!=0));
        if (strcmp(answer,"y")==0) {
           replaceVeto=true;
        } else {
           gSystem->Exit(-1);
        }
      }
    }
    next.Reset();

    // -------------------------------------------------------------------------
    // create temporary root file
    // -------------------------------------------------------------------------
    UChar_t bVeto;
    TFile* ftmp = new TFile(ofname,"RECREATE");
    if (ftmp->IsZombie()) {
      cout << "CWB::Toolbox::setVeto - Error opening file " << ofname.Data() << endl;
      gSystem->Exit(1);
    } 
    TTree *trtmp = (TTree*)itree->CloneTree(0);
    trtmp->SetMaxTreeSize(MAX_TREE_SIZE);
    TString tVeto = veto_name;tVeto+="/b";
    if (replaceVeto) {
      trtmp->SetBranchAddress(veto_name.Data(),&bVeto);
    } else {
      trtmp->Branch(veto_name.Data(),&bVeto,tVeto.Data());
    }
    trtmp->Write();
  
    // -------------------------------------------------------------------------
    // insert flag into event tree
    // we add a dummy flag entry to the end of onoff array to avoid to discart good events when
    // no flag entries are present
    // -------------------------------------------------------------------------
    cout << "Start applying flag to time["<<iIFO<<"]..." << endl;

    ftmp->cd();

    int h=0;
    int pc = 0;
    int ipc = (int)((double)(isize+1)/10.); if(ipc==0) ipc=1;
    int count=0;
    double ttime=time[id[h]];
    for (int h=0;h<isize;h++) bveto[h]=0;
    for (int j=0;j<vsize;j++) {
      while ((ttime<=vlist[j].stop) && (h<isize)) {
        if ((ttime>vlist[j].start)&&(ttime<=vlist[j].stop)) {bveto[id[h]]=1;count++;} else bveto[id[h]]=0;
        if(++h<isize) ttime=time[id[h]];
        if (h%ipc==0) {cout << pc;if (pc<100) cout << " - ";pc+=10;cout.flush();}
      }
    }
    cout << pc << endl << endl;

    for (int h=0;h<isize;h++) {
      itree->GetEntry(h);
      bVeto=bveto[h];
      trtmp->Fill();
    }
    trtmp->Write();
    delete [] id;
    delete [] bveto;

    cout << "Writing new ofile "<< ofname.Data() << endl;
    Int_t osize = (Int_t)trtmp->GetEntries();
    cout << "osize : " << osize << endl;
    cout << "Flagged events:  " << count << " Percentage: "<< (double)count/osize << endl;

    if(history!=NULL) {
      // update history
      history->SetHistoryModify(true);
      if(!history->StageAlreadyPresent(CCAST("VETO"))) history->AddStage(CCAST("VETO"));

      // save pp configuration file
      if(!history->TypeAllowed(CCAST("PARAMETERS"))) history->AddType(CCAST("PARAMETERS"));
      TString ecwb_pparameters_name = TString(gSystem->Getenv("CWB_PPARAMETERS_FILE"));
      for(int i=0;i<gApplication->Argc()-1;i++) {  // skip last argument (net.C)
        if(TString(gApplication->Argv(i)).Contains(".C")) {
          char* parametersBuffer = readFile(gApplication->Argv(i));
          if(parametersBuffer!=NULL) {
            if(TString(gApplication->Argv(i))==ecwb_pparameters_name) {
              history->AddHistory(CCAST("VETO"), CCAST("PARAMETERS"), parametersBuffer);
            } else {
              //history->AddLog(job_stage, parametersBuffer);
            }
          }
          delete [] parametersBuffer;
        }
      }

      if(!history->TypeAllowed(CCAST("WORKDIR"))) history->AddType(CCAST("WORKDIR"));
      char work_dir[1024]="";
      sprintf(work_dir,"%s",gSystem->WorkingDirectory());
      history->AddHistory(CCAST("VETO"), CCAST("WORKDIR"), work_dir);
      //history->AddLog(CCAST("VETO"), CCAST(veto_name.Data()));
      char veto_log[1024];
      sprintf(veto_log,"%s Flagged events: %d/%d - Percentage : %f ", 
                       veto_name.Data(), count, osize, (double)count/osize ); 
      history->AddLog(CCAST("VETO"), CCAST(veto_log));
    }

    ftmp->Close();
    ifile.Close(); 
    cout << endl;

    ifname=ofname;
    if(n%2) ofname.ReplaceAll(".root.tmp.2",".root.tmp.1");
    else    ofname.ReplaceAll(".root.tmp.1",".root.tmp.2");
  }

  delete [] XDQF; 
  for(int i=0;i<nNET;i++) xifos[i].clear();
  delete [] xifos; 

  TString rfname = ofname;
  ofname = efname;

  TString lfname = ifName;			// out live file name
  lfname.ReplaceAll(".root",veto_file_label);
  lfname.ReplaceAll("wave_","live_");
  TString ilfname = idir+"/"+ifName;		// in live file name
  ilfname.ReplaceAll("wave_","live_");
  TString imfname = ilfname;			// in mdc file name
  imfname.ReplaceAll("live_","mdc_");
  TString mfname = lfname;			// out mdc file name
  mfname.ReplaceAll("live_","mdc_");
  TString ilstfname = idir+"/"+ifName;		// in list file name
  ilstfname.ReplaceAll("wave_","merge_");
  ilstfname.ReplaceAll(".root",".lst");
  TString lstfname = lfname;			// out list file name
  lstfname.ReplaceAll("live_","merge_");
  lstfname.ReplaceAll(".root",".lst");

  cout << " rfile    : " << rfname.Data() << endl;
  cout << " ifile    : " << ifname.Data() << endl;
  cout << " ofile    : " << ofname.Data() << endl;
  cout << " lfile    : " << lfname.Data() << endl;
  cout << " mfile    : " << mfname.Data() << endl;
  cout << " ilstfile : " << ilstfname.Data() << endl;
  cout << " lstfile  : " << lstfname.Data() << endl;

  char cmd[1024];
  sprintf(cmd,"mv %s %s",ifname.Data(),ofname.Data());  
  cout << cmd << endl;
  gSystem->Exec(cmd);
  int estat;
  Long_t id,size,flags,mt;
  estat = gSystem->GetPathInfo(ilfname,&id,&size,&flags,&mt);
  if (estat==0) {
    sprintf(cmd,"cd %s;ln -sf ../%s %s",odir.Data(),ilfname.Data(),lfname.Data());  
    cout << cmd << endl;
    gSystem->Exec(cmd);
  }
  estat = gSystem->GetPathInfo(imfname,&id,&size,&flags,&mt);
  if (estat==0) {
    sprintf(cmd,"cd %s;ln -sf ../%s %s",odir.Data(),imfname.Data(),mfname.Data());  
    cout << cmd << endl;
    gSystem->Exec(cmd);
  }
  estat = gSystem->GetPathInfo(ilstfname,&id,&size,&flags,&mt);
  if (estat==0) {
    sprintf(cmd,"cd %s;ln -sf ../%s %s",odir.Data(),ilstfname.Data(),lstfname.Data());  
    cout << cmd << endl;
    gSystem->Exec(cmd);
  }
  sprintf(cmd,"rm %s",rfname.Data());  
  cout << cmd << endl;
  gSystem->Exec(cmd);

  // write history to merge+veto file
  if(history!=NULL) {
    TFile fhist(ofname,"UPDATE");
    history->Write("history");
    fhist.Close();
    delete history;
  }

  vlist.clear();

  return 0;
}

//______________________________________________________________________________
bool
CWB::Toolbox::isFileExisting(TString fName) {
//
// check if file exists 
//
// Input: fName    - file name to be checked
//
// Return true if it exists
//

  // ----------------------------------------------------------------
  // Check if file exists
  // ----------------------------------------------------------------
  Long_t id,size=0,flags,mt;
  int estat = gSystem->GetPathInfo(fName.Data(),&id,&size,&flags,&mt);
  return (estat!=0) ? false : true;
}

//______________________________________________________________________________
bool
CWB::Toolbox::checkFile(TString fName, bool question, TString message) {
//
// check if file exists 
//
// Input: fName    - file name to be checked
//        question - true -> print question
//        message  - question message
//
// Return true if answer=y
//

  bool overwrite=true;

  // ----------------------------------------------------------------
  // Check if file exists
  // ----------------------------------------------------------------
  Long_t id,size=0,flags,mt;
  int estat = gSystem->GetPathInfo(fName.Data(),&id,&size,&flags,&mt);
  if ((estat!=0)&&(!question)) {
    cout << "CWB::Toolbox::checkFile - Error - File/Dir \"" << fName.Data() << "\" not exist" << endl;
    gSystem->Exit(1);
  }

  // -------------------------------------------------------------------------
  // Check if output file already exists and asks if you want to overwrite it
  // -------------------------------------------------------------------------
  if ((estat==0)&&(question)) {
    char answer[1024];
    strcpy(answer,"");
    do {
      cout << "File/Dir " << fName.Data() << " already exist" << endl;
      if(message.Sizeof()>1) cout << message.Data() << endl;
      cout << "Do you want to overwrite it ? (y/n) ";
      cin >> answer;
      cout << endl << endl;
    } while ((strcmp(answer,"y")!=0)&&(strcmp(answer,"n")!=0));
    if (strcmp(answer,"n")==0) overwrite=false;
  }

  return overwrite;
}

//______________________________________________________________________________
void
CWB::Toolbox::mkDir(TString dir, bool question, bool remove) {
//
// make dir 
//
// Input: dir      - directory name to be created
//        question - true -> ask confirm before creation
//        remove   - false/true(def) -> not_remove/remove existing files in the dir
//
// Return true if answer=y
//

  char cmd[1024];
  Long_t id,size=0,flags,mt;
  // Check if dir exist
  int estat = gSystem->GetPathInfo(dir.Data(),&id,&size,&flags,&mt);
  if((estat==0)&&(question==true)) {
    char answer[1024];
    strcpy(answer,"");
    do {
      cout << endl;
      cout << "dir \"" << dir.Data() << "\" already exist" << endl;
      cout << "Do you want to remove the files & recreate dir ? (y/n) ";
      cin >> answer;
      cout << endl;
    } while ((strcmp(answer,"y")!=0)&&(strcmp(answer,"n")!=0));
    if (strcmp(answer,"y")==0) {
      if(remove) sprintf(cmd,"rm %s/*",dir.Data());
      cout << cmd << endl;
      gSystem->Exec(cmd);
      sprintf(cmd,"mkdir -p %s",dir.Data());
      cout << cmd << endl;
      gSystem->Exec(cmd);
    }
  } else if((estat==0)&&(question==false)) {
    if(remove) {
      sprintf(cmd,"rm %s/*",dir.Data());
      cout << cmd << endl;
      gSystem->Exec(cmd);
    }
    sprintf(cmd,"mkdir -p %s",dir.Data());
    cout << cmd << endl;
    gSystem->Exec(cmd);
  } else {
    sprintf(cmd,"mkdir -p %s",dir.Data());
    cout << cmd << endl;
    gSystem->Exec(cmd);
  }

  return;
}

//______________________________________________________________________________
bool
CWB::Toolbox::rmDir(TString dir, bool question) {
//
// remove dir 
//
// Input: dir      - directory name to be removed
//        question - true -> ask confirm before remove
//
// Return true if answer=y
//

  bool banswer=false;
  char cmd[1024];
  Long_t id,size=0,flags,mt;
  int estat = gSystem->GetPathInfo(dir.Data(),&id,&size,&flags,&mt);
  // Check if dir exist
  if((estat==0)&&(question==true)) {
    char answer[1024];
    strcpy(answer,"");
    do {
      cout << endl;
      sprintf(cmd,"rm -rf %s",dir.Data());
      cout << cmd << endl;
      cout << "Do you want to remove the dir ? (y/n) ";
      cin >> answer;
      cout << endl;
    } while ((strcmp(answer,"y")!=0)&&(strcmp(answer,"n")!=0));
    if (strcmp(answer,"y")==0) {
      banswer=true;
      gSystem->Exec(cmd);
    }
  } else if((estat==0)&&(question==false)) {
    sprintf(cmd,"rm -rf %s",dir.Data());
    cout << cmd << endl;
    gSystem->Exec(cmd);
  }

  return banswer;
}

//______________________________________________________________________________
bool
CWB::Toolbox::question(TString question) {
//
// print question and wait for answer y/n 
//
//      question - question to be printed 
//
// Return true/false if answer=y/n
//

  char answer[1024];
  strcpy(answer,"");
  do {
    cout << endl;
    cout << question << " (y/n) ";
    cin >> answer;
    cout << endl;
  } while ((strcmp(answer,"y")!=0)&&(strcmp(answer,"n")!=0));
  return (strcmp(answer,"y")==0) ? true : false;
}

//______________________________________________________________________________
TString
CWB::Toolbox::addCWBFlags(TMacro macro, TString ofname) {
//
// Add Compiler cWB Flags to macro
//
//   macro  - input macro
//   ofname - if != "" than ofname is used to save macro file name with Compiler cWB Flags added
//                     otherwise the name is obtained substituting ".C" in macro.GetTitle() with "_CWBFlags.C"
//
//   return macro file name with Compiler cWB Flags added
//

  bool extension = ofname=="" ? true : false;
  if(ofname=="") ofname = macro.GetTitle();
  if(!ofname.EndsWith(".C"))
    {cout << "CWB::Toolbox::addCWBFlags - Error : macro : " << ofname << " must have extension .C" << endl;gSystem->Exit(1);}
  if(extension) ofname.ReplaceAll(".C","_CWBFlags.C");

  // write macro with cWB compiler flags
  ofstream out;
  out.open(ofname.Data(),ios::out);
  if(!out.good()) {cout << "CWB::Toolbox::addCWBFlags - Error : Opening File : " << ofname << endl;gSystem->Exit(1);}

  // add compiler cWB flags
  out << "// -------------------------------------------------------------------------" << endl;
  out << "// Compiler cWB Flags" << endl;
  out << "// -------------------------------------------------------------------------" << endl;
  if(gSystem->Getenv("_USE_HEALPIX")) {
    out << "#ifndef _USE_HEALPIX" << endl;
    out << "#define _USE_HEALPIX" << endl;
    out << "#endif" << endl;
  }
  if(gSystem->Getenv("_USE_EBBH")) {
    out << "#ifndef _USE_EBBH" << endl;
    out << "#define _USE_EBBH" << endl;
    out << "#endif" << endl;
  }
  if(gSystem->Getenv("_USE_LAL")) {
    out << "#ifndef _USE_LAL" << endl;
    out << "#define _USE_LAL" << endl;
    out << "#endif" << endl;
  }
  out << "// -------------------------------------------------------------------------" << endl;
  out << endl;

  // write macro code
  TList* fLines = macro.GetListOfLines();
  TObjString *obj;
  TIter next(fLines);
  while ((obj = (TObjString*) next())) {
    TString line = obj->GetName();
    out <<  line.Data() << endl;
  }
  out.close();

  return ofname;
}

//______________________________________________________________________________
double
CWB::Toolbox::getZeroLiveTime(int nIFO, TChain& liv, int dummy) {
//
// Calculate live time of zero lag
//
// Input: nIFO    : number of detectors
//        liv     : tree containing live time informatio
//        dummy   : not used
//
// Output: value of live time
//

  // check if slag is presents in liv tree
  TBranch* branch;                        
  bool slagFound=false;                   
  TIter next(liv.GetListOfBranches());    
  while ((branch=(TBranch*)next())) {     
    if (TString("slag").CompareTo(branch->GetName())==0) slagFound=true;
  }                                                                     
  next.Reset();                                                         
  if(!slagFound) {
    cout << "CWB::Toolbox::getZeroLiveTime : Error - live tree do not contains slag leaf" << endl;
    gSystem->Exit(1);
  }

  long int ntot = liv.GetEntries();
  liv.Draw("live",TString::Format("lag[%d]==0&&slag[%d]==0",nIFO,nIFO).Data(),"goff",ntot);
  long int nsel = liv.GetSelectedRows();
  double* live = liv.GetV1();

  double liveZero = 0.;
  for(long int i=0; i<nsel; i++) liveZero += live[i];

  return liveZero;
}

//______________________________________________________________________________
double
CWB::Toolbox::getLiveTime(int nIFO, TChain& liv, wavearray<double>& Trun,
                          wavearray<double>* Wlag, wavearray<double>* Wslag,
                          wavearray<double>& Tlag, wavearray<double>& Tdlag,
                          int lag_number, int slag_number, int dummy) {
//
// Calculate the live time of the total shift performed.
// If defining a lag_number/slag_number, calculate live time referred to this
//
// Input: nIFO        : detector number   
//        liv         : tree containing live time information
//
// Output:Trun        : wavearray containing live time for each run number 
//                      - the returned size is the maximum run number
//        Wlag        : wavearray containing lag progressive number
//        Wslag       : wavearray containing slag progressive number
//        TLag        : output list of time lag
//        TdLag       : output list of time index 
//        lag_number  : specific lag number  
//        slag_number : specific slag number 
//        dummy       : used to increase calculation speed
// 
// Examples : (lag_number==-1) && (slag_number==-1) -> all non zero lags       : excluded (lag==0&&slag==0)
//            (lag_number>=0)  && (slag_number>=0)  -> select (lag==lag_number) && (slag==slag_number)
//            (lag_number>=0)  && (slag_number==-1) -> select lag=lag_number   : excluded (lag==0&&slag==0)
//            (lag_number==-1) && (slag_number>=0)  -> select slag=slag_number : excluded (lag==0&&slag==0)
//                     
// Note: the variable dummy has been introduced to optimize 
//       speed of the method (the reason is unknown!!!)


  float  xlag[NIFO_MAX+1];
  float  xslag[NIFO_MAX+1];
  Int_t  xrun;
  double xlive;
  double liveTot = 0.;
  double Live=0.;
  bool   save;

  wavearray<double>  Qdlag;
  wavearray<double>* Qlag = new wavearray<double>[nIFO+1];
  wavearray<double>* Qslag = new wavearray<double>[nIFO+1];
  wavearray<double>  Qtlag;

  std::map <double, std::map <double, int> > QMap;       

  // check if slag is presents in liv tree
  TBranch* branch;
  bool slagFound=false;
  TIter next(liv.GetListOfBranches());
  while ((branch=(TBranch*)next())) {
    if (TString("slag").CompareTo(branch->GetName())==0) slagFound=true;
  }
  next.Reset();
  if(!slagFound) {
    cout << "CWB::Toolbox::getLiveTime : Error - live tree do not contains slag leaf" << endl;
    gSystem->Exit(1);
  }

  int run_max = liv.GetMaximum("run");

  Trun.resize(run_max);
  long int ntot = liv.GetEntries();

  // get max lag
  TLeaf* leaf = liv.GetLeaf("lag");
  if(!leaf) {
    cout << "CWB::Toolbox::getLiveTime : Error - lag leaf is not present" << endl;
    gSystem->Exit(1);
  }
  branch = leaf->GetBranch();
  int lag_max=0;
  for (long int i=0; i<ntot; ++i) {
    long int entryNumber = liv.GetEntryNumber(i);
    if (entryNumber < 0) break;
    branch->GetEntry(entryNumber);
    double val = leaf->GetValue(nIFO);
    if (val>lag_max) lag_max = (int)val;
  }

  liv.SetBranchAddress("slag",xslag);
  liv.SetBranchAddress("lag",xlag);
  liv.SetBranchAddress("run",&xrun);
  liv.SetBranchAddress("live",&xlive);
  liv.SetBranchStatus("*",false);
  liv.SetBranchStatus("live",true);
  liv.SetBranchStatus("lag",true);
  liv.SetBranchStatus("slag",true);
  liv.SetBranchStatus("run",true);

  int pc = 0;
  int ipc = double(ntot)/10.; if(ipc==0) ipc=1;
  for(long int i=0; i<ntot; i++) {
    liv.GetEntry(i);

    if(i%ipc==0) {if(ntot>100) {cout << pc<<"%";if (pc<100) cout << " - ";pc+=10;cout.flush();}}

    if((lag_number>=0)&&(slag_number>=0)) {
      Live = (xlag[nIFO]==lag_number&&xslag[nIFO]==slag_number) ? xlive : 0.; // lag/slag live time
    }
    if((lag_number>=0)&&(slag_number==-1)) {
      Live = ((xlag[nIFO]==lag_number)&&!((xlag[nIFO]==0)&&(xslag[nIFO]==0))) ? xlive : 0.; // lag live time
    }
    if((lag_number==-1)&&(slag_number>=0)) {
      Live = ((xslag[nIFO]==slag_number)&&!((xlag[nIFO]==0)&&(xslag[nIFO]==0))) ? xlive : 0.; // slag live time
    }
    if((lag_number==-1)&&(slag_number==-1)) {
      Live = !((xlag[nIFO]==0)&&(xslag[nIFO]==0)) ? xlive : 0.; // non-zero live time
    }

    liveTot += Live;
    Trun.data[xrun-1] += Live;

    save = true;

    int K=-1;
    if(QMap.find(xlag[nIFO])!=QMap.end()) {
      if(QMap[xlag[nIFO]].find(xslag[nIFO])!=QMap[xlag[nIFO]].end()) {
        K=QMap[xlag[nIFO]][xslag[nIFO]];
        Qtlag.data[K] += xlive;
        save=false;
      }
    }

    if((lag_number>=0)&&(slag_number>=0)) {
      save = (xlag[nIFO]==lag_number&&xslag[nIFO]==slag_number) ? save : false;
    }
    if((lag_number>=0)&&(slag_number==-1)) {
      save = ((xlag[nIFO]==lag_number)&&!((xlag[nIFO]==0)&&(xslag[nIFO]==0))) ? save : false;
    }
    if((lag_number==-1)&&(slag_number>=0)) {
      save = ((xslag[nIFO]==slag_number)&&!((xlag[nIFO]==0)&&(xslag[nIFO]==0))) ? save : false;
    }
    if((lag_number==-1)&&(slag_number==-1)) {
      save = !((xlag[nIFO]==0)&&(xslag[nIFO]==0)) ? save : false;
    }
    if(save) {
      for (int j=0; j<nIFO; j++) Qlag[j].append(xlag[j]);
      for (int j=0; j<nIFO; j++) Qslag[j].append(xslag[j]);
      Qlag[nIFO].append(xlag[nIFO]);
      Qslag[nIFO].append(xslag[nIFO]);
      Qtlag.append(xlive);
      QMap[xlag[nIFO]][xslag[nIFO]] = Qtlag.size()-1;
      double dlag=0;
      //for (int j=0; j<nIFO; j++) dlag+=(xslag[j]+xlag[j])*(xslag[j]+xlag[j]);
      dlag=xslag[nIFO]*lag_max+xlag[nIFO];
      Qdlag.append(sqrt(dlag));
    }
  }
  if(pc==100) cout << pc<<"%";;
  cout << endl << endl;

  // sort lag using Qdlag (time distance respect to zero lag)
  double* dlag = new double[Qdlag.size()];
  for(int i=0;i<(int)Qdlag.size();i++) dlag[i]=Qdlag.data[i];
  Int_t *id = new Int_t[Qdlag.size()];
  TMath::Sort((int)Qdlag.size(),dlag,id,false);
  for(int i=0;i<(int)Qdlag.size();i++) {
    Tlag.append(Qtlag[id[i]]);
    //Tdlag.append(Qdlag[id[i]]);
    Tdlag.append(i);
    for (int j=0; j<nIFO+1; j++) {
      Wlag[j].append(Qlag[j][id[i]]);
      Wslag[j].append(Qslag[j][id[i]]);
    }
  }
  delete [] id;
  delete [] dlag;

  for (int j=0; j<nIFO+1; j++) {
    Qlag[j].resize(0);
    Qslag[j].resize(0);
  }
  delete [] Qlag;
  delete [] Qslag;
  Qtlag.resize(0);
  Qdlag.resize(0);

  return liveTot;
}

//______________________________________________________________________________
vector<TString>
CWB::Toolbox::getFileListFromDir(TString dir_name, TString endString, TString beginString, TString containString, bool fast) {
//
// get the file list contained in a directory
//
// Input: dir_name	- input dir name
//        endString     - matched end string in file name
//        beginString   - matched initial string in file name
//        containString - matched contained string in file name
//        fast          - true -> skip check if exist (faster mode)
//
// Output: Return file name list contained in the directory "dir_name" 
//

  TString wdir = gSystem->WorkingDirectory();
  TSystemDirectory gdir("", dir_name);
  TList *dfiles = NULL;
  if(fast) {	// skip check if are dirs or files, faster when access to remote files		
    void *dir = gSystem->OpenDirectory(dir_name);
    if(dir) {
      const char *file = 0;
      dfiles  = new TList;
      dfiles->SetOwner();
      while ((file = gSystem->GetDirEntry(dir))) {
        dfiles->Add(new TSystemFile(file, dir_name));
      }
      gSystem->FreeDirectory(dir);
    }
  } else {	// call standard root method (check if are dirs or files)
    dfiles = gdir.GetListOfFiles();
  }
  if(!dfiles) {
    cout << "CWB::Toolbox::getFileListFromDir : Error - dir not found !!!" << endl;
    gSystem->Exit(1);
  }
  TIter dnext(dfiles);
  TSystemFile *dfile;
  TString fname;
  vector<TString> fileList;
  char path[1024];

  while ((dfile = (TSystemFile*)dnext())) {
    fname = dfile->GetName();
    sprintf(path,"%s/%s",dir_name.Data(),fname.Data());
    bool fsave=true;
    if ((endString!="")&&!fname.EndsWith(endString)) fsave=false;
    if ((beginString!="")&&!fname.BeginsWith(beginString)) fsave=false;
    if ((containString!="")&&!fname.Contains(containString)) fsave=false;
    if(fsave) fileList.push_back(path);
  }

  gSystem->ChangeDirectory(wdir);  // restore original dir

  delete dfiles;

  //cout << "CWB::Toolbox::getDirFileLists - nfiles : " << fileList.size() << endl;

  return fileList;
}

//______________________________________________________________________________
std::map<int,TString>
CWB::Toolbox::getJobFileMapFromDir(TString dir_name, TString endString, TString beginString, TString containString, bool fast) {
//
// get the job file map (jobId -> file) contained in a directory
//
// Input: dir_name	- input dir name
//        endString     - matched end string in file name
//        beginString   - matched initial string in file name
//        containString - matched contained string in file name
//        fast          - true -> skip check if exist (faster mode)
//
// Output: Return job file name map (jobId -> file) contained in the directory "dir_name" 
//

  vector<TString> fileList = getFileListFromDir(dir_name, endString, beginString, containString, fast); 

  std::map<int, TString> fileMap;
  for(int i=0;i<fileList.size();i++) {
    fileMap[getJobId(fileList[i])]=fileList[i];
  }
  return fileMap;
}

//______________________________________________________________________________
void
CWB::Toolbox::doPoissonPlot(int nIFO, wavearray<double>* Wlag, wavearray<double> Tlag, 
                            wavearray<double> Rlag, TString odir) {
//
// It produces Poisson Plot on background stage
//
// Input nIFO : detector number
//       Wlag : wavearray containing lag progressive number
//       Tlag : wavearray containing live time
//       Rlag : wavearray containing rate of each lag 
//       odir : output dir
//

// NOTE: Wlag,Tlag,Rlag are sorted using the time distance respect to zero lag -> Tlag[0] is the zero lag livetime
//       Uses lags with livetime which differ max 10% respect to lag zero livetime

  int    NC=0;
  double LT=0.,LAMBDA;
  double normalization=0;
  double enormalization=0;
  double mean=0;
  double emean=0;
  double chi2=0;
  int    ndf=0;
  int    myNDF =-1; //one degree of freedom less to account for the mean estimation
  double plevel=0;
  double ChiSquare=0.0;
  double dChiSquare=0.0;
  int    minCountsPerBin=5;

  double Tlag0 = Tlag.data[0];	// zero lag time
  int    Nlags = 0;		// number of lags which are used for poisson check

  // compute min/max false alarms per lag
  int min_fa_per_lag=1000000000;
  int max_fa_per_lag=0;
  for(int j=0; j<(int)Wlag[nIFO].size(); j++){
    if(Wlag[nIFO].data[j]>0.) { 
      if(fabs(Tlag.data[j]-Tlag0)/Tlag0>0.1) continue; 
      int fa_per_lag=Rlag.data[j]*Tlag.data[j];
      if(fa_per_lag<min_fa_per_lag) min_fa_per_lag=fa_per_lag;
      if(fa_per_lag>max_fa_per_lag) max_fa_per_lag=fa_per_lag;
      Nlags++;
    }
  }
  int nbin = max_fa_per_lag-min_fa_per_lag;

  //Histogram for the Poisson Test on the background estimation
  TH1I *h1 = new TH1I("h1","h1",nbin,min_fa_per_lag,max_fa_per_lag);

  for(int j=0; j<(int)Wlag[nIFO].size(); j++){
    if(Wlag[nIFO].data[j]>0.) { //Fill histogram for the Poisson Test on the background estimation
      if(fabs(Tlag.data[j]-Tlag0)/Tlag0>0.1) continue; 
      h1->Fill(Rlag.data[j]*Tlag.data[j]);
      NC+=Rlag.data[j]*Tlag.data[j];
      LT+=Tlag.data[j];
    }
  }
  LAMBDA= LT>0 ? NC/LT : 0;

  cout << "Total Number of bkg coinc.: " << NC << " total Live Time: " 
       << LT << " nb = " << LAMBDA << endl;

  if(NC==0) return;

  TCanvas *canvas = new TCanvas("BackgroundPoissonFit", "Background Poisson Fit",32,55,750,502);
  canvas->Range(-0.75,-1.23433,6.75,8.17182);
  canvas->SetBorderSize(1);
  canvas->SetFillColor(0);
  canvas->SetGridx();
  canvas->SetGridy();
  canvas->SetFrameFillColor(0);

  //Poisson Integer Function
  TF1 poisson("PoissonIFunction",PoissonIFunction,min_fa_per_lag,max_fa_per_lag,2);
  poisson.FixParameter(0,h1->GetEntries());
  poisson.SetParameter(1,h1->GetMean());

  TH1D* hfit = (TH1D*)h1->Clone("hfit");

  //loop over the populated bins (counts > minCountsPerBin) to calculate the p-level
  for (int k=1;k<=hfit->GetNbinsX();k++) {
    int n = k-1;
    if (poisson.Eval(n)>minCountsPerBin) {
      myNDF++; //myNDF starts from -1 !!
      dChiSquare=pow((hfit->GetBinContent(k)-poisson.Eval(n)),2)/poisson.Eval(n);
      ChiSquare+=dChiSquare;
      cout << "NCoinc = " << n << " Found = " << hfit->GetBinContent(k) << "   Expected = "
           << poisson.Eval(n) << "     Chi2_bin" << n << "= " << dChiSquare<<endl;
    }
  }

  if(myNDF<=0) {
    cout << "CWB::Toolbox::doPoissonPlot - Warning !!! - Poisson NDF <=0 " << endl;
    return;
  }

  // chi2
  double myPLevel=TMath::Prob(ChiSquare,myNDF);
  cout<<"myChiSquare  :"<<ChiSquare<<endl;
  cout << "NDF       : " << myNDF << endl;
  cout<<"myPlevel       :"<<myPLevel<<endl;
  int XMIN=kMaxInt,XMAX=kMinInt;
  int g=0;
  while((hfit->GetBinContent(g)==0)&&(g<hfit->GetNbinsX())){XMIN=g;g++;}
  g=0;
  while((hfit->GetBinContent(hfit->GetNbinsX()-g)==0)&&(g<hfit->GetNbinsX())){XMAX=hfit->GetNbinsX()-g;g++;}
  XMIN--;
  XMAX++;
  cout<<XMIN<<" "<<XMAX<<endl;
  hfit->Fit("PoissonIFunction","R");
  TF1* fpois = hfit->GetFunction("PoissonIFunction");
  fpois->SetFillColor(kGreen);
  fpois->SetFillStyle(3002);
  fpois->SetLineColor(kGreen);
  fpois->SetLineStyle(1);
  fpois->SetLineWidth(1);
  fpois->SetLineColor(kGreen);
  hfit->SetTitle("Poisson Fit (Black : Data - Green : Fit)");
  //hfit->SetTitle(name);
  hfit->SetLineWidth(4);
  hfit->GetXaxis()->SetTitle("#False Alarms/lag");
  hfit->GetYaxis()->SetTitle("#Counts");
  hfit->GetXaxis()->SetTitleSize(0.04);
  hfit->GetYaxis()->SetTitleSize(0.05);
  hfit->GetXaxis()->SetLabelSize(0.04);
  hfit->GetYaxis()->SetLabelSize(0.04);
  hfit->GetXaxis()->SetTitleOffset(0.97);
  hfit->GetYaxis()->SetTitleOffset(1.0);
  hfit->GetXaxis()->SetRange(XMIN,XMAX);
  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0);

  hfit->Draw();

  //Standard Chi2 & Probability
  chi2 = fpois->GetChisquare();
  ndf = fpois->GetNDF()-1;
  plevel = TMath::Prob(chi2,ndf);
  normalization = fpois->GetParameter(0);
  enormalization = fpois->GetParError(0);
  mean = fpois->GetParameter(1);
  emean = fpois->GetParError(1);

  cout << "Mean      : "<<mean<<endl;
  cout << "Norm     : "<<normalization<<endl;
  cout << "ChiSquare : " << chi2 << endl;
  cout << "NDF       : " << ndf << endl;
  cout << "Plevel    : " << plevel << endl;

  TLegend *legend = new TLegend(0.597,0.660,0.961,0.921,NULL,"brNDC");
  legend->SetLineColor(1);
  legend->SetTextSize(0.033);
  legend->SetBorderSize(1);
  legend->SetLineStyle(1);
  legend->SetLineWidth(1);
  legend->SetFillColor(10);
  legend->SetFillStyle(1001);

  legend->SetTextSize(0.03);
  char entry[1024];
  sprintf(entry,"# lags = %d ",Nlags);
  legend->AddEntry("",entry,"");
  //sprintf(entry,"Nc = %d @ rho>%3.1f",NC,T_out);
  sprintf(entry,"Nc = %d",NC);
  legend->AddEntry("*",entry,"");
  sprintf(entry,"Total Live time = %d s",int(LT));
  legend->AddEntry("*",entry,"");
  sprintf(entry,"Mean = %.3f +/-%.3f  FA/lag",mean,emean);
  legend->AddEntry("*",entry,"");
  sprintf(entry,"Chi2 (counts>%d) = %.3f",minCountsPerBin,ChiSquare);
  legend->AddEntry("*",entry,"");
  sprintf(entry,"NDF = %d",myNDF);
  legend->AddEntry("*",entry,"");
  sprintf(entry,"p-level = %.3f",myPLevel);
  legend->AddEntry("*",entry,"");
  legend->Draw();

  canvas->Update();
  char fname[1024];
  sprintf(fname,"%s/Lag_PoissonFit.png",odir.Data());
  canvas->Print(fname);

  delete h1;
  delete canvas;

  return;
}

//______________________________________________________________________________
char*
CWB::Toolbox::getEnvCWB() {
//
// return a buffer with all environment name variables and their values
//

  const int nCWB = 35;

  TString ecwb_name[nCWB] = {
                              "SITE_CLUSTER",
                              "_USE_ICC",
                              "_USE_CPP11",
                              "_USE_ROOT6",
                              "_USE_PEGASUS",
                              "_USE_HEALPIX",
                              "_USE_LAL",
                              "_USE_EBBH",
                              "HOME_WAT",
                              "HOME_WAT_FILTERS",
                              "HOME_WAT_INSTALL",
                              "HOME_FRLIB",
                              "HOME_HEALPIX",
                              "HOME_CFITSIO",
                              "HOME_CVODE",
                              "HOME_LAL",
                              "LAL_INC",
                              "LALINSPINJ_EXEC",
                              "CWB_GWAT",
                              "CWB_HISTORY",
                              "CWB_MACROS",
                              "CWB_ROOTLOGON_FILE",
                              "CWB_PARAMETERS_FILE",
                              "CWB_UPARAMETERS_FILE",
                              "CWB_PPARAMETERS_FILE",
                              "CWB_UPPARAMETERS_FILE",
                              "CWB_NETS_FILE",
                              "CWB_NETC_FILE",
                              "CWB_HTML_INDEX",
                              "CWB_HTML_HEADER",
                              "CWB_HTML_BODY_PROD",
                              "CWB_HTML_BODY_SIM",
                              "ROOT_VERSION",
                              "ROOTSYS",
                              "LD_LIBRARY_PATH"
                             };

  TString ecwb_value[nCWB];
  for(int i=0;i<nCWB;i++) {
    if(gSystem->Getenv(ecwb_name[i])==NULL) {
      ecwb_value[i]="";
      //cout << "Error : environment " << ecwb_name[i].Data() << " is not defined!!!" << endl;gSystem->Exit(1);
    } else {
      ecwb_value[i]=TString(gSystem->Getenv(ecwb_name[i]));
    }

  }

  int ecwb_csize=0;
  for(int i=0;i<nCWB;i++) ecwb_csize+=ecwb_name[i].Sizeof();
  for(int i=0;i<nCWB;i++) ecwb_csize+=ecwb_value[i].Sizeof();

  char* iBuffer = new char[2*ecwb_csize];
  bzero(iBuffer,(2*ecwb_csize)*sizeof(char));

  int len;
  int iLength = 0;
  for(int i=0;i<nCWB;i++) {
    len = ecwb_name[i].Sizeof();
    strncpy(iBuffer+iLength,ecwb_name[i].Data(),len);
    iLength += len-1;
    iBuffer[iLength]='=';
    iLength += 1;
    len = ecwb_value[i].Sizeof();
    strncpy(iBuffer+iLength,ecwb_value[i].Data(),len);
    iLength += len-1;
    iBuffer[iLength]=0x0a;
    iLength += 1;
  }

  return iBuffer;
}

//______________________________________________________________________________
bool
CWB::Toolbox::isLeafInTree(TTree* itree, TString leaf) {
//
// check if leaf is present in the input tree
//
// Input: itree       - pointer to the input tree
//        leaf        - name of the leaf
//
// Return true/false if leaf is present or not present  
//

  TBranch* branch;
  bool bleaf=false;
  TIter next(itree->GetListOfBranches());
  while ((branch=(TBranch*)next())) {
    if (TString(leaf.Data()).CompareTo(branch->GetName())==0) bleaf=true;
  }
  next.Reset();
  return bleaf;
}


//______________________________________________________________________________
void
CWB::Toolbox::makeSpectrum(wavearray<double>& psd, wavearray<double> x, double chuncklen, double scratchlen, bool oneside) {
//
// make PSD
// PSD is computed averaging N=x.size()/chuncklen energy FFT with length=chuncklen 
//
// Input: x           - input data sample
//        chuncklen   - FFT length
//        scratchlen  - input data scratch length
//        oneside     - true/false = one side/double side PSD
//
// Output : psd       - array with the psd values
//

  if(chuncklen<=0) {
    cout << "CWB::Toolbox::makeSpectrum : Error - chuncklen<=0" << endl;
    exit(1);
  }

  if(scratchlen<0) {
    cout << "CWB::Toolbox::makeSpectrum : Error - scratchlen<0" << endl;
    exit(1);
  }

  int scratchsize = scratchlen*x.rate();
  int blocksize = chuncklen*x.rate();
  scratchsize-=scratchsize%2;	// make it even
  blocksize-=blocksize%2;	// make it even
  double df=(double)x.rate()/(double)(blocksize);

  if(x.size()-2*scratchsize<blocksize) {
    cout << "CWB::Toolbox::makeSpectrum : Error - data size not enough to produce PSD" << endl;
    exit(1);
  }

  int loops = (x.size()-2*scratchsize)/blocksize;
  cout << "Rate: " << x.rate() << endl;

  double* window = new double[blocksize];
  blackmanharris(window, blocksize);

  psd.resize(blocksize/2);
  psd=0;

  wavearray<double> y(blocksize);
  y.rate(x.rate());

  for (int n=0;n<loops;n++) {

    int shift=n*blocksize;
    //cout << "shift: " << shift << endl;
    for (int i=0;i<blocksize;i++) y.data[i]=x.data[i+scratchsize+shift];
    for (int i=0;i<blocksize;i++) y.data[i]*=window[i];

    y.FFTW(1);
    for (int i=0;i<blocksize;i+=2) psd[i/2]+=pow(y.data[i],2)+pow(y.data[i+1],2);
  }

  for (int i=0;i<blocksize/2; i++) psd[i]=sqrt(psd[i]/(double)loops);
  if(oneside) {
    for (int i=0;i<blocksize/2; i++) psd[i]*=sqrt(2/df);  // one side psd
  } else {
    for (int i=0;i<blocksize/2; i++) psd[i]*=sqrt(1/df);  // double side psd
  }

  psd.start(0.);
  psd.stop(df*psd.size());

  delete [] window;
  return;
}

//______________________________________________________________________________
void
CWB::Toolbox::makeSpectrum(TString ofname, wavearray<double> x, double chuncklen, double scratchlen, bool oneside) {
//
// make PSD
// PSD is computed averaging N=x.size()/chuncklen energy FFT with length=chuncklen 
// the psd valued are saved to ofname file
//
// Input: ofname      - output psd file name (format : freq  psd)
//        x           - input data sample
//        chuncklen   - FFT length
//        scratchlen  - input data scratch length
//        oneside     - true/false = one side/double side PSD
//

  wavearray<double> psd;
  makeSpectrum(psd, x, chuncklen, scratchlen, oneside);
  double df = (psd.stop()-psd.start())/psd.size();

  ofstream out;
  out.open(ofname.Data(),ios::out);
  if (!out.good()) {cout << "CWB::Toolbox::makeSpectrum - Error : Opening File : " << ofname.Data() << endl;gSystem->Exit(1);}
  for (int i=0;i<psd.size();i++) {
    double freq = i*df+psd.start();
    out << freq << " " << psd[i] << endl;
  }
  out.close();

  return;
}

//______________________________________________________________________________
void
CWB::Toolbox::getSimNoise(wavearray<double>& u, TString fName, int seed, int run) {
//
// get colored gaussian noise
//
// Input:  fName       - input file with PSD values
//                       format : freq(Hz) singleside-PSD(1/sqrt(Hz))
//         seed        - base seed for random number generator
//                       if seed<0 input wavearray is used instead of gaussian random noise 
//                                 we pad freq<fabs(seed) with value A = amplitude at freq=fabs(seed)
//                                 'run' parameter is used as a multiplicative factor for A -> A=A*run
//         run         - auxiliary seed for random number generator
//                       SEED = seed+run
//         u.size()    - input wavearray size 
//
// Output: u           - wavearray filled with colored gaussian noise @ rate=16384 Hz
//

  #define OTIME_LENGHT   616		// minimum time lenght simulation
  #define TIME_SCRATCH	 184		// time scratch used by the FFT 
  #define LOW_CUT_FREQ   2.0            // output noise is 0 for freq<LOW_CUT_FREQ
  #define FREQ_RES       0.125		// input PSD is resampled with dF=FREQ_RES 
  #define SRATE          16384.		// output noise is produced @ rate=SRATE

  TRandom3 random;

  // input lenght must be <= OTIME_LENGHT
  double ilenght = u.size()/u.rate();
  int otime_lenght = (ilenght<OTIME_LENGHT) ? OTIME_LENGHT : int(ilenght);

  // read PSD
  ifstream in;
  in.open(fName.Data(),ios::in);
  if (!in.good()) {
    cout << "CWB::Toolbox::getSimNoise - Error Opening File : " << fName.Data() << endl;
    gSystem->Exit(1);
  }

  int size=0;
  char str[1024];
  while(true) {
    in.getline(str,1024);
    if (!in.good()) break;
    if(str[0] != '#') size++;
  }
  //cout << "size " << size << endl;
  in.clear(ios::goodbit);
  in.seekg(0, ios::beg);

  wavearray<double> ifr(size);
  wavearray<double> ish(size);

  int cnt=0;
  while (1) {
    in >> ifr.data[cnt] >> ish.data[cnt];
    if (!in.good()) break;
    if(ish.data[cnt]<=0)
      {cout << "CWB::Toolbox::getSimNoise - input sensitivity file : " << fName.Data() 
            << " contains zero at frequency : " << ifr.data[cnt] << " Hz " << endl;gSystem->Exit(1);}
    cnt++;
  }
  in.close();

  // convert frequency sample
  size=int((SRATE/2)/FREQ_RES);
  wavearray<double> ofr(size);
  wavearray<double> osh(size);
  for(int i=0;i<(int)ofr.size();i++) ofr[i]=i*FREQ_RES;
  convertSampleRate(ifr,ish,ofr,osh);
  ifr.resize(0);
  ish.resize(0);

  osh*=1./sqrt(osh.size()/(SRATE/2));  // normalization

  int time_factor = (otime_lenght+TIME_SCRATCH)/(2*osh.size()/double(SRATE));

  // change output time lenght
  wavearray<double> y;         	// temporary time series
  y.resize(2*osh.size()*time_factor);
  y=0.;
  for (int i=0;i<(int)osh.size();i++) {
    for (int j=0;j<time_factor;j++) {
      y.data[i*time_factor+j]=osh.data[i];
    }
  }
  y.rate(SRATE);
  y*=1./sqrt(time_factor);
  ofr.resize(0);
  osh.resize(0);

  wavearray<double> z;         	// temporary time series
  z.rate(y.rate());
  z.resize(y.size());
  double df=z.rate()/z.size();
  for (int i=0;i<(int)z.size()/2;i++) {
    double am = y.data[i];
    double frequency = df*i;
    if (frequency>=LOW_CUT_FREQ) {
      z.data[2*i]=am;  		// (A)
      z.data[2*i+1]=am;  	// (B)
    } else {
      z.data[2*i]=0;
      z.data[2*i+1]=0;
    }
  }
  z*=1./sqrt(2.);  		// because of (A & B)
  y.resize(0);

  int scratch=z.size()/z.rate()-otime_lenght;
  cout << "CWB::Toolbox::getSimNoise - scratch : " << scratch << " osize : " << z.size()/z.rate() << endl;
  if (scratch<0) {cout << "Error : bad data length -> " << z.size()/z.rate() << endl;gSystem->Exit(1);}

  wavearray<double> x;         // temporary time series
  x.rate(z.rate());
  x.resize(z.size());

  if(seed>=0) {
    // generate random gaussian noise -> FFT
    random.SetSeed(seed+run);
    for (int i=0;i<x.size()-scratch*x.rate();i++) x.data[i]=random.Gaus(0,1);
    random.SetSeed(seed+run+1);  // to be syncronized with the next frame
    for (int i=x.size()-scratch*x.rate();i<(int)x.size();i++) x.data[i]=random.Gaus(0,1);
  } else {
    // input wavearray is used instead of gaussian random noise
    if(u.size()>x.size()) {
      cout << "CWB::Toolbox::getSimNoise - Error : whitened data size " << u.size() 
           << " is greater than temporary datat size " << x.size() << endl;
      gSystem->Exit(1);
    }
    x=0;
    int iS = (x.size()-u.size())/2;
    for(int i=0;i<u.size();i++) x[i+iS]=u[i];
  }
  x.FFTW(1);
  x*=sqrt(x.size());  // Average x^2 = 1

  // gaussian white noise -> gaussian coloured noise
  u=z;
  u*=x;
  x.resize(0);
  z.resize(0);

  if(seed<0) {	
    // cwb whitened real noise is 0 for freq<16Hz
    // we pad freq<fabs(seed) with value = amplitude at freq=fabs(seed)*run
    // run parameter is used as a multiplicative factor 

    double df=u.rate()/u.size();
    int ipad = fabs(seed)/df;
    if(ipad>u.size()/2-1) ipad=u.size()/2-1;
    double fpad = sqrt(pow(u.data[2*ipad],2)+pow(u.data[2*ipad+1],2))*run;
    for(int i=0;i<(int)u.size()/2;i++) {
      double frequency = df*i;
      if(frequency<=fabs(seed)) {
        u.data[2*i]   = gRandom->Gaus(0,fpad);
        u.data[2*i+1] = gRandom->Gaus(0,fpad);
      }
    }
  }

  // change output frequency SRATE
  df=u.rate()/u.size();
  int osize = SRATE/df;
  u.resize(osize);
  u.rate(SRATE);

  u.FFTW(-1);

  scratch=(osize-otime_lenght*u.rate())/2;
  for (int i=0;i<otime_lenght*u.rate();i++) u.data[i]=u.data[scratch+i];

  u.resize(otime_lenght*u.rate());

  return;
}

//______________________________________________________________________________
wavearray<double>
CWB::Toolbox::GetDetectorPSD(TString fName, double fWidth, double dFreq) {
//
// read PSD from file
//
// input : fName  - input file name [freq(Hz)  psd]
//         fWidth - bandwidth (Hz)
//         dFreq  - frequency resolution (Hz)
//
// return PSD sampled at df=dFreq up to fWidth Hz
//

  ifstream in;
  in.open(fName.Data(),ios::in);
  if (!in.good()) 
    {cout << "CWB::Toolbox::ReadDetectorPSD - Error Opening File : " 
          << fName.Data() << endl;gSystem->Exit(1);}

  cout << "CWB::Toolbox::ReadDetectorPSD - Read File : " << fName.Data() << endl;                                   

  int size=0;
  char str[1024];
  while(true) {  
    in.getline(str,1024);
    if (!in.good()) break;
    if(str[0] != '#') size++;
  }                          
  //cout << "size " << size << endl;
  in.clear(ios::goodbit);           
  in.seekg(0, ios::beg);            

  wavearray<double> ifr(size);
  wavearray<double> ish(size);

  int cnt=0;
  while (1) {
    in >> ifr.data[cnt] >> ish.data[cnt];
    if (!in.good()) break;               
    if(ish.data[cnt]<=0)                 
      {cout << "CWB::Toolbox::ReadDetectorPSD - input sensitivity file : " << fName.Data()
            << " contains zero at frequency : " << ifr.data[cnt] << " Hz " << endl;gSystem->Exit(1);}
    cnt++;                                                                                  
  }                                                                                         
  in.close();                                                                               

  // convert frequency sample
  size=int(fWidth/dFreq);
  wavearray<double> ofr(size); 
  wavearray<double> osh(size); 

  for(int i=0;i<(int)ofr.size();i++) ofr[i]=i*dFreq;
  convertSampleRate(ifr,ish,ofr,osh);               

  osh.rate(size*dFreq);

  return osh;
}            

//______________________________________________________________________________
void
CWB::Toolbox::convertSampleRate(wavearray<double> iw, wavearray<double> ow) {
//
// convert sample rate data
// data are resampled using a linear approximation
//
// input  : iw     - input samples wavearray, sampled @ iw.rate() 
// output : ow     - output samples wavearray, sampled @ ow.rate() 
//                   ow.rate() is an input parameter
//

  if(iw.size()<=0) {cout << "CWB::Toolbox::convertSampleRate : Error -  input size <=0" << endl;gSystem->Exit(1);}
  if(iw.rate()<=0) {cout << "CWB::Toolbox::convertSampleRate : Error -  input rate <=0" << endl;gSystem->Exit(1);}
  if(ow.size()<=0) {cout << "CWB::Toolbox::convertSampleRate : Error - output size <=0" << endl;gSystem->Exit(1);}
  if(ow.rate()<=0) {cout << "CWB::Toolbox::convertSampleRate : Error - output rate <=0" << endl;gSystem->Exit(1);}

  // initialize input array
  wavearray<double> ix(iw.size());
  double idx=1./iw.rate();
  for (int i=0;i<(int)ix.size();i++) ix[i]=i*idx;

  // initialize output array
  wavearray<double> ox(ow.size());
  double odx=1./ow.rate();
  for (int i=0;i<(int)ox.size();i++) ox[i]=i*odx;

  // smooth amplitudes
  TGraph* grin = new TGraph(ix.size(), ix.data, iw.data);
  TGraphSmooth* gs = new TGraphSmooth("normal");
  TGraph* grout = gs->Approx(grin,"linear", ox.size(), ox.data);

  // save amplitudes
  for (int i=0;i<(int)ox.size();i++) {
    grout->GetPoint(i,ox[i],ow[i]);
//    cout << i << " " << ox[i] << " " << ow[i] << endl;
  }

  delete grin;
  delete gs;

  return;
}

//______________________________________________________________________________
void
CWB::Toolbox::convertSampleRate(wavearray<double> ix, wavearray<double> iy, wavearray<double> ox, wavearray<double>& oy) {
//
// convert sample rate data
//
// data are resampled using a linear approximation
// ox are used to resample input amplitudes iy sampled at times ix
//
// input  : ix     - input times wavearray 
//          iy     - input amplitudes wavearray 
//          ox     - output times wavearray 
// output : oy     - output amplitutes wavearray 
//

  if(ix.size()<=0) {cout << "CWB::Toolbox::convertSampleRate : Error -  input size <=0" << endl;gSystem->Exit(1);}
  if(ox.size()<=0) {cout << "CWB::Toolbox::convertSampleRate : Error - output size <=0" << endl;gSystem->Exit(1);}

  // smooth amplitudes
  TGraph* grin = new TGraph(ix.size(), ix.data, iy.data);
  TGraphSmooth* gs = new TGraphSmooth("normal");
  TGraph* grout = gs->Approx(grin,"linear", ox.size(), ox.data);

  // save amplitudes
  for (int i=0;i<(int)ox.size();i++) {
    grout->GetPoint(i,ox[i],oy[i]);
  }

  delete grin;
  delete gs;

  return;
}

//______________________________________________________________________________
int
CWB::Toolbox::setMultiplicity(TString ifName, TString idir, TString odir, int nIFO, double Tgap) {
//
// Create Multiplicity plot
//
// Input: ifName : input file name
//        idir   : input dir
//        odir   : output dir
//        nIFO   : detector number
//        Tgap   : maximum time gap to associate two events as same event
//

  gRandom->SetSeed(1);

  CWB::History* history = NULL;
  bool bhistory=true;

  TString efname = odir+"/"+ifName;
  efname.ReplaceAll(".root",".N.root");
  bool overwrite = checkFile(efname,true);
  if(!overwrite) gSystem->Exit(0);

  TString ifname = idir+"/"+ifName;
  TString ofname = odir+"/"+ifName;
  ofname.ReplaceAll(".root",".root.tmp.1");

  for(int n=0;n<nIFO;n++) {

    // -------------------------------------------------------------------------
    // open root file
    // -------------------------------------------------------------------------
    cout<<"Opening BKG File : " << ifname.Data() << endl;
    TFile ifile(ifname);
    if (ifile.IsZombie()) {
      cout << "CWB::Toolbox::setMultiplicity - Error opening file " << ifname.Data() << endl;
      gSystem->Exit(1);
    } 
    TTree *itree = (TTree*)ifile.Get("waveburst");
    if (itree==NULL) {
      cout << "CWB::Toolbox::setMultiplicity - tree waveburst is not present in file " << ifname.Data() << endl;
      gSystem->Exit(1);
    } 
    Int_t tsize = (Int_t)itree->GetEntries();
    cout << "tree size : " << tsize << endl;
    if (tsize==0) {
      cout << "CWB::Toolbox::setMultiplicity - tree waveburst is empty in file " << ifname.Data() << endl;
      gSystem->Exit(1);
    } 
    itree->SetEstimate(tsize);
    char selection[1024];sprintf(selection,"time[%d]:Entry$",n);
    char cut[1024]="";
//    if(n>0) sprintf(cut,"Mm[%d]==1",n-1);
    int isize = itree->Draw(selection,cut,"goff",tsize);
    double* time = itree->GetV1(); 
    double* entry = itree->GetV2(); 
    // sort list
    Int_t *id = new Int_t[isize];
    TMath::Sort((int)isize,time,id,false);

    if(bhistory) {
      history=(CWB::History*)ifile.Get("CWB::History");
      if(history==NULL) history=new CWB::History(); 
      bhistory=false;
    }

    // -------------------------------------------------------------------------
    // add new leave to itree
    // -------------------------------------------------------------------------

    TBranch* branch;
    if(n==0) {
      TIter next(itree->GetListOfBranches());
      while((branch=(TBranch*)next())) {
        if(TString("Msize").CompareTo(branch->GetName())==0) {
          cout << "Multiplicity already applied" << endl;
          char answer[1024];
          strcpy(answer,"");
          do{
            cout << "Do you want to overwrite the previous values ? (y/n) ";
            cin >> answer;
          } while((strcmp(answer,"y")!=0)&&(strcmp(answer,"n")!=0));
          if(strcmp(answer,"n")==0) {
            gSystem->Exit(-1);
          }
        }
      }
      next.Reset();
    }

    // -------------------------------------------------------------------------
    // create temporary root file
    // -------------------------------------------------------------------------

    int* Msize = new int[nIFO];
    int* Mid = new int[nIFO];
    UChar_t* Mm = new UChar_t[nIFO];

    char cMsize[32];sprintf(cMsize,"Msize[%1d]/I",nIFO);
    char cMid[32];sprintf(cMid,"Mid[%1d]/I",nIFO);
    char cMm[32];sprintf(cMm,"Mm[%1d]/b",nIFO);

    TFile* ftmp = new TFile(ofname,"RECREATE");
    if (ftmp->IsZombie()) {
      cout << "CWB::Toolbox::setMultiplicity - Error opening file " << ofname.Data() << endl;
      gSystem->Exit(1);
    } 
    TTree *trtmp = (TTree*)itree->CloneTree(0);
    trtmp->SetMaxTreeSize(MAX_TREE_SIZE);
    if (n>0) {
      trtmp->SetBranchAddress("Msize",Msize);
      trtmp->SetBranchAddress("Mid",Mid);
      trtmp->SetBranchAddress("Mm",Mm);
    } else {
      trtmp->Branch("Msize",Msize,cMsize);
      trtmp->Branch("Mid",Mid,cMid);
      trtmp->Branch("Mm",Mm,cMm);
    }
    trtmp->Write();
  
    // -------------------------------------------------------------------------
    // insert flag into event tree
    // we add a dummy flag entry to the end of onoff array to avoid to discart good events when
    // no flag entries are present
    // -------------------------------------------------------------------------
    cout << "Start applying flag to time["<<n<<"]..." << endl;

    ftmp->cd();

    int pc = 0;
    int ipc = (int)((double)(isize+1)/10.);
    int *oMsize = new int[tsize];
    int *oMid = new int[tsize];
    bool *oMm = new bool[tsize];
    for(int j=0;j<tsize;j++) {
      oMsize[j]=1;
      oMid[j]=1;
      oMm[j]=false;
    }
    vector<int> mlist;
    int xMsize=1;
    int xMid=1;
    int ientry=int(entry[id[0]]);
    oMsize[ientry]=xMsize; oMid[ientry]=xMid; oMm[ientry]=false;
    mlist.push_back(id[0]);
    for (int j=1;j<isize;j++) {
      ientry=int(entry[id[j]]);
      oMsize[ientry]=1; 
      oMm[ientry]=false; 
      if((time[id[j]]-time[id[j-1]])<Tgap) {
        xMsize++;
        mlist.push_back(id[j]);
      } else {
        int xMm=mlist[int(gRandom->Uniform(0,mlist.size()))];  // select random multiplicity event
        oMm[int(entry[xMm])]=true;  // set master multiplicity event
        for(int m=0;m<(int)mlist.size();m++) oMsize[int(entry[mlist[m]])]=xMsize;
        mlist.clear();
        xMsize=1;xMid++;
        mlist.push_back(id[j]);
      } 
      oMid[ientry]=xMid; 

      if(ipc!=0) if(j%ipc==0) {cout << pc;if (pc<100) cout << " - ";pc+=10;cout.flush();}
    }
    cout << pc << endl << endl;

    int* iMsize = new int[nIFO];
    int* iMid = new int[nIFO];
    UChar_t* iMm = new UChar_t[nIFO];
    if(n>0) {
      itree->SetBranchAddress("Msize",iMsize);
      itree->SetBranchAddress("Mid",iMid);
      itree->SetBranchAddress("Mm",iMm);
    }
    int nevt_master=0;
    for (int j=0;j<tsize;j++) {
      itree->GetEntry(j);
      for(int k=0;k<nIFO;k++) {
        Msize[k]=iMsize[k];
        Mid[k]=iMid[k];
        Mm[k]=iMm[k];
      }
      Msize[n]=oMsize[j];
      Mid[n]=oMid[j];
      Mm[n]=oMm[j];
      if(Mm[n]) nevt_master++;
      trtmp->Fill();
/*
if(n==(nIFO-1)) {
  if(Mm[n]) trtmp->Fill();
} else {
  trtmp->Fill();
} 
*/
    }
    trtmp->Write();
    delete [] id;
    delete [] oMsize;
    delete [] oMid;
    delete [] oMm;
    delete [] Msize;
    delete [] Mid;
    delete [] Mm;
    delete [] iMsize;
    delete [] iMid;
    delete [] iMm;

    cout << "Writing new ofile "<< ofname.Data() << endl;
    Int_t osize = (Int_t)trtmp->GetEntries();
    cout << "osize : " << osize << endl;
    cout << "Master events:  " << nevt_master << " Percentage: "<< 100.*double(nevt_master)/double(tsize) << endl;

    if((history!=NULL)&&(n==nIFO-1)) {
      // update history
      history->SetHistoryModify(true);
      if(!history->StageAlreadyPresent(CCAST("MULTIPLICITY"))) history->AddStage(CCAST("MULTIPLICITY"));

      if(!history->TypeAllowed(CCAST("WORKDIR"))) history->AddType(CCAST("WORKDIR"));
      char work_dir[512]="";
      sprintf(work_dir,"%s",gSystem->WorkingDirectory());
      history->AddHistory(CCAST("MULTIPLICITY"), CCAST("WORKDIR"), work_dir);
      char hTgap[1024];sprintf(hTgap,"Tgap = %f",Tgap);
      history->AddLog(CCAST("MULTIPLICITY"), CCAST(hTgap));
    }

    ftmp->Close();
    ifile.Close(); 
    cout << endl;

    ifname=ofname;
    if(n%2) ofname.ReplaceAll(".root.tmp.2",".root.tmp.1");
    else    ofname.ReplaceAll(".root.tmp.1",".root.tmp.2");
  }

  TString rfname = ofname;
  ofname = efname;

  // define output live merge file
  TString lfname = ifName;
  lfname.ReplaceAll(".root",".N.root");
  lfname.ReplaceAll("wave_","live_");
  TString ilfname = idir+"/"+ifName;
  ilfname.ReplaceAll("wave_","live_");

  // define output list merge file
  TString mfname = ifName;
  mfname.ReplaceAll(".root",".N.root");
  mfname.ReplaceAll("wave_","merge_");
  mfname.ReplaceAll(".root",".lst");
  TString imfname = idir+"/"+ifName;
  imfname.ReplaceAll("wave_","merge_");
  imfname.ReplaceAll(".root",".lst");

  cout << " rfile " << rfname.Data() << endl;
  cout << " ifile " << ifname.Data() << endl;
  cout << " ofile " << ofname.Data() << endl;
  cout << " lfile " << lfname.Data() << endl;
  cout << " mfile " << mfname.Data() << endl;

  char cmd[1024];
  sprintf(cmd,"mv %s %s",ifname.Data(),ofname.Data());  
  cout << cmd << endl;
  gSystem->Exec(cmd);
  // create symbolic link to live root file
  sprintf(cmd,"cd %s;ln -s ../%s %s",odir.Data(),ilfname.Data(),lfname.Data());  
  cout << cmd << endl;
  gSystem->Exec(cmd);
  // create symbolic link to list merge file
  sprintf(cmd,"cd %s;ln -s ../%s %s",odir.Data(),imfname.Data(),mfname.Data());  
  cout << cmd << endl;
  gSystem->Exec(cmd);
  sprintf(cmd,"rm %s",rfname.Data());  
  cout << cmd << endl;
  gSystem->Exec(cmd);

  // write history to merge+veto file
  if(history!=NULL) {
    TFile fhist(ofname,"UPDATE");
    history->Write();
    fhist.Close();
    delete history;
  }

  return 0;
}

//______________________________________________________________________________
wavecomplex
CWB::Toolbox::getRate(double rho, double Tgap, int nIFO, TChain& wav, wavearray<int>& Wsel,
                      wavearray<double>* Wlag, wavearray<double>* Wslag, wavearray<double> Tlag) {
//
// Calculate FAR rate at a certain rho treshold
//
// Input rho   : rho threshold
//       Tgap  : 
//       nIFO  : detector number
//       wav   : input tree
//       Wsel  :
//       Wlag  : wavearray containing progressive lag number
//       Wslag : wavearray containing progressive slag number
//       Tlag  : wavearray containing live time
//

//  cout << "rhoThr       : " <<  rho  << endl;
  gRandom->SetSeed(1);
  netevent W(&wav,nIFO);
  // slag is not defined in the wat-5.4.0
  float* iWslag = new float[nIFO+1];
  W.fChain->SetBranchAddress("slag",iWslag);

  int size = Wlag[nIFO].size();

  double liveTot=0;
  for(int i=0;i<(int)Tlag.size();i++) liveTot+=Tlag[i];
//  cout << "liveTot : " << liveTot << endl;

  // compute the number of indipendent slags,lags
  int n;
  bool save;
  wavearray<int> slag;
  wavearray<int> lag;
  for(int i=0; i<size; i++) {
    int islag=int(Wslag[nIFO].data[i]);
    save = true;
    n = slag.size();
    for(int j=0; j<n; j++) {
      if(slag[j]==islag) {save=false; j=n;}
    }
    if(save) slag.append(islag);

    int ilag=int(Wlag[nIFO].data[i]);
    save = true;
    n = lag.size();
    for(int j=0; j<n; j++) {
      if(lag[j]==ilag) {save=false; j=n;}
    }
    if(save) lag.append(ilag);
  }

//  for(int i=0;i<(int)lag.size();i++) cout << i << " LAG " << lag[i] << endl;
//  for(int i=0;i<(int)slag.size();i++) cout << i << " SLAG " << slag[i] << endl;

  // construct the array to associate slag,lag to index
  int* slag2id = new int[slag.max()+1];
  int* lag2id = new int[lag.max()+1];
  for(int i=0;i<(int)slag.size();i++) slag2id[slag[i]]=i;
  for(int i=0;i<(int)lag.size();i++) lag2id[lag[i]]=i;

  // construct the matrix to associate lag,slag to live times
  int** lag2live = (int**)malloc((slag.size())*sizeof(int*));
  for (int i=0;i<(int)slag.size();i++) lag2live[i] = new int[lag.size()];

  for(int i=0; i<size; i++) {
    int islag=int(Wslag[nIFO].data[i]);
    int ilag=int(Wlag[nIFO].data[i]);
    lag2live[slag2id[islag]][lag2id[ilag]]=Tlag[i];
  }

  // list of rejected lags
  bool** lagRej = (bool**)malloc((slag.size())*sizeof(bool*));
  for (int i=0;i<(int)slag.size();i++) lagRej[i] = new bool[lag.size()];
  for(int i=0;i<(int)slag.size();i++) for(int j=0;j<(int)lag.size();j++) lagRej[i][j]=false;

  int ntrg = wav.GetEntries();

  int *id = new int[ntrg];
  int *wslag = new int[ntrg];
  int *wlag = new int[ntrg];
  int *entry = new int[ntrg];
  double** time = (double**)malloc(nIFO*sizeof(double*));
  for (int i=0;i<nIFO;i++) time[i] = new double[ntrg];
  bool *trgMaster = new bool[ntrg];
  for (int i=0;i<ntrg;i++) trgMaster[i] = false;
  bool *wrej = new bool[ntrg];
  for (int i=0;i<ntrg;i++) wrej[i] = false;

  int isel=0;
  for(int i=0; i<ntrg; i++) {
    if(!Wsel[i]) continue;
    W.GetEntry(i);
    if(W.rho[1]>rho) {
      bool skip=false; 
      for(int j=0;j<nIFO;j++) if(!TMath::IsNaN(W.time[j])) time[j][isel]=W.time[j]; else skip=true;
      if(skip) continue; 
      //wslag[isel]=W.slag[nIFO];
      wslag[isel]=iWslag[nIFO];
      wlag[isel]=W.lag[nIFO];
      entry[isel]=i;
      isel++;
    }
  }
//  cout << "selTrg       : " <<  isel  << endl;

  int rejTrg=0;
  double rejLive=0;
  vector<double> livlist;
  vector<int> idlist;
  for(int i=0;i<nIFO;i++) {
    if(isel==0) continue;
    TMath::Sort((int)isel,time[i],id,false);
    int islag = slag2id[wslag[id[0]]];    
    int ilag = lag2id[wlag[id[0]]];    
    double live = lag2live[islag][ilag];
    livlist.push_back(live);
    idlist.push_back(id[0]);
    for (int j=1;j<isel;j++) {
      if(wrej[id[j]]) continue;  // skip trigger if already rejected
      if((time[i][id[j]]-time[i][id[j-1]])>Tgap) {
/*
        // find minimum live time in the multiplicity list
        double liveMin=kMaxInt; int M=0;
        for(int m=0;m<(int)livlist.size();m++) {
          if(liveMin>livlist[m]) {liveMin=livlist[m];M=m;}
        }
*/
        int M=int(gRandom->Uniform(0,(int)idlist.size()));  // select random multiplicity event
        trgMaster[idlist[M]]=true;
        for(int m=0;m<(int)idlist.size();m++) {
          // select only the trigger in the multiplicity list with minimum live time  
          // check if live time has been already rejected
          int islag = slag2id[wslag[idlist[m]]];    
          int ilag = lag2id[wlag[idlist[m]]];    
          if(!trgMaster[idlist[m]]) {
             if(!lagRej[islag][ilag]) rejLive+=livlist[m];
             lagRej[islag][ilag]=true;
             wrej[idlist[m]]=true;   // Tag the rejected trigger
             Wsel[entry[idlist[m]]]=2;   // Tag the rejected trigger with 2
             rejTrg++;
          }
        }
        livlist.clear();
        idlist.clear();
      }
      int islag = slag2id[wslag[id[j]]];    
      int ilag = lag2id[wlag[id[j]]];    
      live = lag2live[islag][ilag];
      livlist.push_back(live);
      idlist.push_back(id[j]);
    }
/*
    // find minimum live time in the multiplicity list
    double liveMin=kMaxInt; int M=0;
    for(int m=0;m<(int)livlist.size();m++) {
      if(liveMin>livlist[m]) {liveMin=livlist[m];M=m;}
    }
*/
    int M=int(gRandom->Uniform(0,(int)idlist.size()));  // select random multiplicity event
    trgMaster[idlist[M]]=true;
    for(int m=0;m<(int)idlist.size();m++) {
      // select only the trigger in the multiplicity list with minimum live time  
      // check if live time has been already rejected
      int islag = slag2id[wslag[idlist[m]]];    
      int ilag = lag2id[wlag[idlist[m]]];    
      if(!trgMaster[idlist[m]]) {
        if(!lagRej[islag][ilag]) rejLive+=livlist[m];
        lagRej[islag][ilag]=true;
        wrej[idlist[m]]=true;   // Tag the rejected trigger
        Wsel[entry[idlist[m]]]=2;   // Tag the rejected trigger with 2
        rejTrg++;
      }
    }
    livlist.clear();
    idlist.clear();
  }
//  cout << "rejTrg       : " <<  rejTrg  << endl;

  // reject triggers which belong to the rejected lags
  for (int j=0;j<isel;j++) {
    if(wrej[j]) continue;   // only the triggers not already rejected
    int islag = slag2id[wslag[j]];    
    int ilag = lag2id[wlag[j]];    
    if(lagRej[islag][ilag]) {
      Wsel[entry[j]]=2;   // Tag the rejected trigger with 2
      rejTrg++;
    }
  }  
//  cout << "rejTrg final : " <<  rejTrg  << endl;
//  cout << "rejLive      : " <<  rejLive << endl;

/*
int xrejTrg=0;
for (int j=0;j<isel;j++) {
  int islag = slag2id[wslag[j]];    
  int ilag = lag2id[wlag[j]];    
  if(lagRej[islag][ilag]) {
    xrejTrg++;
  }
}
cout << "xrejTrg final: " <<  xrejTrg  << endl;
int nlagRej=0;
for(int i=0;i<(int)slag.size();i++) for(int j=0;j<(int)lag.size();j++) if(lagRej[i][j]) nlagRej++;
cout << "nlagRej final: " <<  nlagRej  << endl;
double xlivRej=0;
for(int i=0;i<(int)slag.size();i++) for(int j=0;j<(int)lag.size();j++) if(lagRej[i][j]) xlivRej+=lag2live[i][j];
cout << "xlivRej final: " <<  xlivRej  << endl;
*/
/*
  for(int j=0; j<size; j++){
    printf("slag=%i\t lag=%i\t ",int(Wslag[nIFO].data[j]),int(Wlag[nIFO].data[j]));
    //for (int jj=0; jj<nIFO; jj++) printf("%d:slag=%5.2f\tlag=%5.2f\t",jj,Wslag[jj].data[j],Wlag[jj].data[j]);
    printf("live=%9.2f\n",Tlag.data[j]);
  }
*/ 

  for (int i=0;i<nIFO;i++) delete [] time[i];
  delete [] id;
  delete [] wslag;
  delete [] wlag;
  for (int i=0;i<(int)slag.size();i++) delete [] lag2live[i];
  delete [] slag2id;
  delete [] lag2id;
  delete [] iWslag;
  delete [] trgMaster;
  delete [] wrej;
  delete [] entry;

  double rate = (isel-rejTrg)/(liveTot-rejLive);
  double erate = sqrt(isel-rejTrg)/(liveTot-rejLive);
  return wavecomplex(rate,erate);
}


//_______________________________________________________________________________________
int 
CWB::Toolbox::GetStepFunction(TString fName, vector<double>&  x, vector<double>&  y, 
          vector<double>& ex, vector<double>& ey) {      

  GetStepFunction("", fName, 0, x, y, ex, ey);

  // remove the first and the last elements added by the GetStepFunction method
  int size=x.size();

  // erase the first element (adde)
  x.erase(x.begin()+size-1);
  y.erase(y.begin()+size-1);
  ex.erase(ex.begin()+size-1);
  ey.erase(ey.begin()+size-1);

  // erase the last element
  x.erase(x.begin());
  y.erase(y.begin());
  ex.erase(ex.begin());
  ey.erase(ey.begin());

  return x.size();
}

//_______________________________________________________________________________________
double 
CWB::Toolbox::GetStepFunction(TString option, TString fName, double V,
         vector<double>&  x, vector<double>&  y, vector<double>& ex, vector<double>& ey) {      
//                                                                                       
// read x,y coodinated from fName
// y[x] is approximated with a step function 
// 
// V       : input value
//
// fName   : ascii file with the list of x,y coordinateds
//           x1	y1
//           x2	y2
//           .....
//           xN	yN
//
// option  : "xmin" - return minimum x value   
//           "xmax" - return maximim x value
//           "ymin" - return minimum y value
//           "ymax" - return maximim y value
//           "y"    - return the y value corresponding to x=V
//
// ERRORS  : the same options can be used for errors : use 'e' in front of the option
//           for example to get error of 'x' use 'ex'
//           the input file must have 4 columns format : x y ex ey
//

  bool ferror=false;	// if true return errors

  // get option
  option.ToUpper();
  if(option!="") { 
    if((option!="XMIN")&&(option!="XMAX")&&
       (option!="YMIN")&&(option!="YMAX")&&
       (option!="Y")&&(option!="EX")&&(option!="EY")) {                    
      cout<<"GetLinearInterpolation : Error : option --value bad value -> "
          <<"must be y/xmin/xmax/ymin/ymax"<<endl;
      exit(1);                                                                              
    }                                                                                       
    if(option.BeginsWith("E")) {ferror=true;option.Remove(0,1);}
  }                                                                                         

  double xmin=DBL_MAX;
  double ymin=DBL_MAX;  
  double xmax=-DBL_MAX;
  double ymax=-DBL_MAX;  

  double exmin=DBL_MAX;
  double eymin=DBL_MAX;  
  double exmax=-DBL_MAX;
  double eymax=-DBL_MAX;  

  ifstream in;
  in.open(fName.Data());
  if (!in.good()) {cout << "GetGraphValue : Error Opening File : " << fName << endl;exit(1);}
  double ix,iy;    
  double iex,iey;    
  char line[1024];                                                                            
  while(1) {
    in.getline(line,1024);
    if (in.eof()) break;                                                                   
    std::stringstream linestream(line);
    if(ferror) {	// get error value
      if(linestream >> ix >> iy >> iex >> iey) {
        ex.push_back(iex); 
        ey.push_back(iey); 
        if(iex<exmin) exmin=iex;
        if(iey<eymin) eymin=iey;
        if(iex>exmax) exmax=iex;
        if(iey>eymax) eymax=iey;
      } else {
        linestream.str(line);
        linestream.clear();
        if(!(linestream >> ix >> iy)) {
          cout << "GetGraphValue : Wrong line format : must be 'x y' " << endl;
          exit(1);
        }
        ex.push_back(0); ey.push_back(0); 
        exmin=0;eymin=0;exmax=0;eymax=0;
      }
    } else {
      if(!(linestream >> ix >> iy)) {
        cout << "GetGraphValue : Wrong line format : must be 'x y' " << endl;
        exit(1);
      }
      ex.push_back(0); ey.push_back(0); 
      exmin=0;eymin=0;exmax=0;eymax=0;
    }
    x.push_back(ix); 
    y.push_back(iy); 
    if(ix<xmin) xmin=ix;
    if(iy<ymin) ymin=iy;
    if(ix>xmax) xmax=ix;
    if(iy>ymax) ymax=iy;
  }
  in.close();

  if(option=="XMIN") return ferror ? exmin : xmin;
  if(option=="XMAX") return ferror ? exmax : xmax;
  if(option=="YMIN") return ferror ? eymin : ymin;
  if(option=="YMAX") return ferror ? eymax : ymax;

  if(x.size()==0) cout << "GetGraphValue : Error - input file empty" << endl;

  // Sort respect to x values
  int size = x.size();
  wavearray<int> I(size);
  wavearray<double> X(size);
  for(int i=0;i<size;i++) X[i]=x[i];
  TMath::Sort(size,X.data,I.data,false);

  // insert a new element at the beginning with x=-DBL_MAX and y[I[0]]
  // necessary for TGraph
  std::vector<double>::iterator it;
  it =  x.begin();  x.insert(it,-DBL_MAX);
  it =  y.begin();  y.insert(it,y[I[0]]);
  it = ex.begin(); ex.insert(it,ex[I[0]]);
  it = ey.begin(); ey.insert(it,ey[I[0]]);
  size = x.size();

  // Re-Sort respect to x values
  I.resize(size);
  X.resize(size);
  for(int i=0;i<size;i++) X[i]=x[i];
  TMath::Sort(size,X.data,I.data,false);

  x.push_back(x[I[size-1]]); y.push_back(-DBL_MAX);
  ex.push_back(ex[I[size-1]]); ey.push_back(ey[I[size-1]]);  
  I.resize(size+1); I[size]=I[size-1]; size++;

  if(!ferror) {				// return y
    if(V<x[I[0]]) return y[I[0]];
    if(V>=x[I[size-1]]) return y[I[size-1]];
    for(int i=1;i<size;i++) if((V>=x[I[i-1]]) && (V<x[I[i]])) return y[I[i-1]];
  } else {	
    if(option=="X") {			// return ex
      if(V<x[I[0]]) return ex[I[0]];
      if(V>=x[I[size-1]]) return ex[I[size-1]];
      for(int i=1;i<size;i++) if((V>=x[I[i-1]]) && (V<x[I[i]])) return ex[I[i-1]];
    } else if(option=="Y") {		// return ey
      if(V<x[I[0]]) return ey[I[0]];
      if(V>=x[I[size-1]]) return ey[I[size-1]];
      for(int i=1;i<size;i++) if((V>=x[I[i-1]]) && (V<x[I[i]])) return ey[I[i-1]];
    }
  }
}


//______________________________________________________________________________
void 
CWB::Toolbox::resampleToPowerOfTwo(wavearray<double>& w) {
//
// convert to a power of 2 rate > original rate
//

  // compute the new rate 
  int R=1;while (R < 2*(int)w.rate()) R*=2; 

  if(R==2*w.rate()) return;

  Meyer<double> BB(256);
  WSeries<double> ww(BB);
  wavearray<double> yy(2*w.size());
  yy.rate(2*w.rate());
  yy.start(w.start());

  ww.Forward(yy,BB,1);
  ww=0.;
  ww.putLayer(w,0);
  ww.Inverse();
  yy.resample(ww,R,32);
  ww.Forward(yy,BB,1);
  ww.getLayer(w,0);

  fprintf(stdout,"--------------------------------\n");
  fprintf(stdout,"After resampling: rate=%f, size=%d, start=%f\n", w.rate(),(int)w.size(),w.start());
  fprintf(stdout,"--------------------------------\n");

  return;
}

//______________________________________________________________________________
void 
CWB::Toolbox::PrintProcInfo(TString str) {
//
// write date and memory usage at a certain point of the analysis
//

   gSystem->Exec("date");
   TString s;
   FILE *f = fopen(Form("/proc/%d/statm", gSystem->GetPid()), "r");
   s.Gets(f);
   Long_t total, rss;
   sscanf(s.Data(), "%ld %ld", &total, &rss);
   cout << str.Data() << " virtual : " <<  total * 4 / 1024 << " (mb)  rss  : " <<  rss * 4 / 1024 << " (mb)" << endl;
   fclose(f);
   return;
}

//______________________________________________________________________________
TString 
CWB::Toolbox::getTemporaryFileName(TString label, TString extension) {
//
//  get name of temporary file on the node
//
// Input: label associate to the analysis
//        file extension
//

  // create temporary file
  gRandom->SetSeed(0);
  int rnID = int(gRandom->Rndm(13)*1.e9);
  UserGroup_t* uinfo = gSystem->GetUserInfo();
  TString uname = uinfo->fUser;
  // create temporary dir
  gSystem->Exec(TString("mkdir -p /dev/shm/")+uname);

  char fName[1024]="";

  sprintf(fName,"/dev/shm/%s/%s_%d.%s",uname.Data(),label.Data(),rnID,extension.Data());

  return fName;
}

//______________________________________________________________________________
vector<TString>
CWB::Toolbox::sortStrings(vector<TString> ilist) {
//
// sort in alphabetic order the TString list 
//
// Input: TString list
// return the sorted TString list
//

  std::vector<std::string> vec(ilist.size()); 
  for(int i=0;i<ilist.size();i++) vec[i]=ilist[i].Data();
  std::sort( vec.begin(), vec.end() ) ;
  vector<TString> slist(ilist.size());
  for(int i=0; i<vec.size(); i++) slist[i]=vec[i]; 

  return slist;
}

//______________________________________________________________________________
int 
CWB::Toolbox::getJobId(TString file, TString fext) {
//
// return jobId from the file name
// file must have the following format *_job#id.'fext'
//

  if(!file.EndsWith("."+fext)) {
    cout << "CWB::Toolbox::getJobId : input file " 
         << file.Data() << " must terminate with ." << fext << endl; 
    gSystem->Exit(1);     
  }

  file.ReplaceAll("."+fext,"");
  TObjArray* token = file.Tokenize('_');
  TObjString* sjobId = (TObjString*)token->At(token->GetEntries()-1);
  TString check = sjobId->GetString(); check.ReplaceAll("job","");
  if(!check.IsDigit()) {
    cout << "CWB::Toolbox::getJobId : input file " << file.Data() 
         << " must terminate with _job#id." << fext << endl; 
    cout << "#id is not a digit" << endl; 
    gSystem->Exit(1);     
  }
  int jobId = TString(sjobId->GetString()).ReplaceAll("job","").Atoi();
  if(token) delete token;

  return jobId;
}

//______________________________________________________________________________
TString 
CWB::Toolbox::getParameter(TString options, TString param) {
// get params from option string                                                     
// options format : --par1 val1 --par2 val2 ...

  TObjArray* token;

  // check if input param format is correct
  token = param.Tokenize(TString(" "));
  if((token->GetEntries())>0) {
    TObjString* otoken = (TObjString*)token->At(0);
    TString stoken = otoken->GetString();
    if(((token->GetEntries())>1) || (!stoken.BeginsWith("--"))) {
      cout << "CWB::Toolbox::getParameter : Bad param format \"" << param.Data() << "\"" << endl; 
      cout << "Correct format is : --par1" << endl;
      gSystem->Exit(1);     
    }
  }
  if(token) delete token;

  token = options.Tokenize(TString(" "));
  // check if input options format is correct
  if((token->GetEntries())%2==1) {
    cout << "CWB::Toolbox::getParameter : Bad options format \"" << options.Data() << "\"" << endl; 
    cout << "Correct format is : --par1 val1 --par2 val2 ..." << endl;
    gSystem->Exit(1);     
  }
  for(int i=0;i<token->GetEntries();i+=2) {
    TObjString* otoken = (TObjString*)token->At(i);
    TString stoken = otoken->GetString();
    if(!stoken.BeginsWith("--")) {
      cout << "CWB::Toolbox::getParameter : Bad options format \"" << stoken.Data() << "\"" << endl; 
      cout << "Correct format is : --par1 val1 --par2 val2 ..." << endl;
      gSystem->Exit(1);     
    }
  }
  for(int i=0;i<token->GetEntries();i++) {
    TObjString* otoken = (TObjString*)token->At(i);
    TString stoken = otoken->GetString();
    // exctract value
    if(stoken==param) {
      if(i<token->GetEntries()-1) {
        otoken = (TObjString*)token->At(i+1);
        return otoken->GetString();
      }
    }
  }
  if(token) delete token;

  return "";
}

//______________________________________________________________________________
TString 
CWB::Toolbox::getFileName(FILE* fp) { 
//
// get file name from its pointer fp
//

  int MAXSIZE = 0xFFF;
  char proclnk[0xFFF];
  char filename[0xFFF];
  int fno;
  ssize_t r;

  if (fp != NULL)
  {
    fno = fileno(fp);
    sprintf(proclnk, "/proc/self/fd/%d", fno);
    r = readlink(proclnk, filename, MAXSIZE);
    if (r < 0) return "";
    filename[r] = '\0';
    return filename;
  }
  return "";
}

//______________________________________________________________________________
TString 
CWB::Toolbox::getFileName(char* symlink) { 
//
// get file name from its symbolic link
//

  int MAXSIZE = 0xFFF;
  char proclnk[0xFFF];
  char filename[0xFFF];
  ssize_t r;

  if (symlink != NULL)
  {
    sprintf(proclnk, "%s", symlink);
    r = readlink(symlink, filename, MAXSIZE);
    if (r < 0) return "";
    filename[r] = '\0';
    return filename;
  }
  return "";
}

//______________________________________________________________________________
void  
CWB::Toolbox::getUniqueFileList(TString ifile, TString ofile) {
//
// extract from ifile (contains a list of files) a unique file list and write it to ofile
//

  vector<std::string> ifileList;
  vector<std::string> ipathList;
  vector<std::string> ofileList;
  vector<std::string> opathList;

  ifstream in;
  in.open(ifile.Data());
  if(!in.good()) {
    cout << "CWB::Toolbox::getUniqueFileList : Error Opening Input File : " << ifile.Data() << endl;
    gSystem->Exit(1);
  }

  // read file list
  char istring[1024];
  while(1) {
    in.getline(istring,1024);
    if (!in.good()) break;
    TObjArray* token = TString(istring).Tokenize(TString('/'));
    // extract last entry -> file name
    TObjString* stoken =(TObjString*)token->At(token->GetEntries()-1);
    TString fName = stoken->GetString();
    //cout << fName.Data() << endl;
    ipathList.push_back(istring);
    ifileList.push_back(fName.Data());
  }
  in.close();

  // extract unique file list
  for(int i=0;i<(int)ifileList.size();i++) {
    bool check=false;
    for(int j=0;j<(int)ofileList.size();j++) {
      if(TString(ofileList[j])==TString(ifileList[i])) {check=true;break;}
    }
    if(!check) {
      ofileList.push_back(ifileList[i]);
      opathList.push_back(ipathList[i]);
    }
    //cout << i << " " << ipathList[i] << endl;
    //cout << i << " " << ifileList[i] << endl;
  }

  // write unique file list
  ofstream out;
  out.open(ofile.Data(),ios::out);
  if(!out.good()) {
    cout << "CWB::Toolbox::getUniqueFileList : Error Opening Output File : " << ofile.Data() << endl;
    gSystem->Exit(1);
  }

  for(int i=0;i<(int)ofileList.size();i++) {
    out << opathList[i] << endl;
    //cout << i << " " << opathList[i] << endl;
    //cout << i << " " << ofileList[i] << endl;
  }

  out.close();

  return;
}

//______________________________________________________________________________
wavearray<double> 
CWB::Toolbox::getHilbertTransform(wavearray<double> x) {
//
// return the Hilbert transform
//

  wavearray<double> y=x;

  y.FFTW(1);
  for(int i=1;i<(int)y.size()/2;i++) {
    double temp=y[2*i+1];             
    y[2*i+1]=y[2*i];                  
    y[2*i]=-1*temp;                   
  }                                   
  y.FFTW(-1);                         

  return y;
}          

//______________________________________________________________________________
wavearray<double> 
CWB::Toolbox::getHilbertEnvelope(wavearray<double> x) {
//
// return the envelope from analytic signal analysis
//

  wavearray<double> y=x;

  // compute hilbert transform 
  wavearray<double> h = getHilbertTransform(y);

  // compute envelope
  for(int i=0;i<(int)y.size();i++) {
    y[i]=sqrt(y[i]*y[i]+h[i]*h[i]); 
  }                                 

  return y;
}          

//______________________________________________________________________________
wavearray<double> 
CWB::Toolbox::getHilbertIFrequency(wavearray<double> x) {
//
// return the instantaneous frequency from analytic signal analysis
//

  wavearray<double> y=x;

  double twopi = TMath::TwoPi();
  double dt = 1./x.rate();

  // compute hilbert transform
  wavearray<double> h = getHilbertTransform(y);

  // compute phase (use p as temporary array for phase)
  wavearray<float> p(h.size());
  for(int i=0;i<(int)y.size();i++) {
    p[i]=TMath::ATan2(y[i],h[i]);
  }

  unWrapPhase(p);
  for(int i=1;i<(int)p.size()/2-1;i++) p[2*i]=(p[2*i-1]+p[2*i+1])/2;   // cut alias frequency 

  // compute frequency
  for(int i=1;i<(int)x.size();i++) {
    y[i] = (p[i]-p[i-1])/(twopi*dt);
    if(y[i]>x.rate()/2) y[i]=0;		// cut unphysical frequency
    if(y[i]<0) y[i]=0;		        // cut unphysical frequency
  }
  if(y.size()>1) y[0]=y[1]; else y[0]=0;
  if(y.size()>1) y[y.size()-1]=y[y.size()-2]; else y[y.size()-1]=0;

  return y;
}

//______________________________________________________________________________
wavearray<double>  
CWB::Toolbox::getWignerVilleTransform(wavearray<double> x) {
//
// return TF map wavearray NxN : Matrix[i,j]=N*j+i : i=time_index, j=freq_index  
//

  int  N = x.size();
  int dN = 2*N;

  wavearray<double> y = x;
  y.FFTW(1);
  y.resize(2*y.size());
  for(int i=(int)x.size();i<(int)y.size();i++) y[i]=0;
  y.FFTW(-1);

  wavearray<double> z(3*dN);
  z=0;for(int i=dN; i<2*dN; ++i) z[i] = y[i-dN];

  wavearray<double> wv(N*N);
  wv.start(x.start());
  wv.rate(x.rate());

  y.resize(dN);
  for(int n=1; n<=N; ++n ) {

    for(int i=0;  i<N; ++i)    y[i] = z[dN+2*n+i]*z[dN+2*n-i];
    for(int i=-N; i<0; ++i) y[dN+i] = z[dN+2*n+i]*z[dN+2*n-i];

    y.FFTW(1);

    for(int m=0; m<N; m++) wv[m*N+(n-1)] = y[2*m];
  }

  return wv;
} 

//______________________________________________________________________________
void 
CWB::Toolbox::unWrapPhase(wavearray<float>& p) {
//
// Phase unwrapping ensures that all appropriate multiples of 2*pi have been included in phase evolution
// p : phase array
// return in p the unwrapped phase
//

  int N = p.size();

  // ported from matlab (Dec 2002)
  wavearray<float> dp(N);     
  wavearray<float> dps(N);    
  wavearray<float> dp_corr(N);
  wavearray<float> cumsum(N);
  float cutoff = M_PI;               /* default value in matlab */
  int j;

  // incremental phase variation 
  // MATLAB: dp = diff(p, 1, 1);
  for (j = 0; j < N-1; j++) dp[j] = p[j+1] - p[j];
    
  // equivalent phase variation in [-pi, pi]
  // MATLAB: dps = mod(dp+dp,2*pi) - pi;
  for (j = 0; j < N-1; j++)
      dps[j] = (dp[j]+M_PI) - floor((dp[j]+M_PI) / (2*M_PI))*(2*M_PI) - M_PI;

  // preserve variation sign for +pi vs. -pi
  // MATLAB: dps(dps==pi & dp>0,:) = pi;
  for (j = 0; j < N-1; j++)
    if ((dps[j] == -M_PI) && (dp[j] > 0))
       dps[j] = M_PI;

  // incremental phase correction
  // MATLAB: dp_corr = dps - dp;
  for (j = 0; j < N-1; j++)
    dp_corr[j] = dps[j] - dp[j];
      
  // Ignore correction when incremental variation is smaller than cutoff
  // MATLAB: dp_corr(abs(dp)<cutoff,:) = 0;
  for (j = 0; j < N-1; j++)
    if (fabs(dp[j]) < cutoff)
      dp_corr[j] = 0;

  // Find cumulative sum of deltas
  // MATLAB: cumsum = cumsum(dp_corr, 1);
  cumsum[0] = dp_corr[0];
  for (j = 1; j < N-1; j++)
    cumsum[j] = cumsum[j-1] + dp_corr[j];

  // Integrate corrections and add to P to produce smoothed phase values
  // MATLAB: p(2:m,:) = p(2:m,:) + cumsum(dp_corr,1);
  for (j = 1; j < N; j++) p[j] += cumsum[j-1];

}

//______________________________________________________________________________
void
CWB::Toolbox::getSineFittingParams(double a, double b, double c, double rate, 
                                   double& amplitude, double& omega, double& phase) {
//
// return the amplitude,phase,omega of the 3 numbers Sine fitting 
// a,b,c are 3 consecutive discrete values sampled at rate='rate'
// return != 0 only if (b>a && b>c) or (b<a && b<c) 
//

  double dt = 1./rate;

  // The time of maximum is the time of b - cp.phase/cp.omega

  if(((a<b && b>=c && b>0)||(a>b && b<=c && b<0)) && (fabs((a+c)/(2*b))<1)){

    omega     = acos((a+c)/(2*b))/dt;
    phase     = atan2(a*a+c*a-2*b*b, 2*a*b*sin(omega*dt));
    amplitude = (a!=0)?(a/cos(phase)):(b/sin(omega*dt));

    return;  
  } else {
    amplitude = 0;
    return;
  }
}

//______________________________________________________________________________
TString                                                                         
CWB::Toolbox::WriteFrameFile(wavearray<double> x, TString chName, TString frName, TString frLabel, TString frDir) {
//                                                                                                                     
// Write x array to frame file                                                                                             
//                                                                                                                     
// Input: x         - wavearray data                                                                             
//        chName    - channel name                                                    
//        frName    - frame name                                                    
//        frLabel   - label used for output file name                                                                  
//                    file name path = frDir/frLabel-gps-length.gwf                      
//        frDir     - output directory                                                                                 
//                                                                                                                     
// return file name 
//

  // check input wavearray data
  if(x.size()==0 || x.rate()==0) {
    cout << "CWB::Toolbox::WriteFrameFile - Error : wavearray size=0 or rate=0 " << endl;
    exit(1);                                                                                              
  }                                                                                                       

  double gps = x.start();
  double length = (x.size()/x.rate());

  // check gps integer
  if(fmod(gps,1)!=0) {
    cout << "CWB::Toolbox::WriteFrameFile - Error : start gps time is not integer" << endl;
    exit(1);                                                                                              
  }                                                                                                       
  // check length integer
  if(fmod(length,1)!=0) {
    cout << "CWB::Toolbox::WriteFrameFile - Error : length is not integer" << endl;
    exit(1);                                                                                              
  }                                                                                                       

  if(frDir=="") frDir=".";

  char frFile[1024];
  sprintf(frFile,"%s/%s-%lu-%lu.gwf",frDir.Data(),frLabel.Data(),(int)gps,(int)length);
  cout << frFile << endl;                                                                        
  FrFile *ofp = FrFileONew(frFile,1);   // gzip compression                                      

  /*----------------------- Create a new frame ---------*/

  FrameH* simFrame = FrameNew(const_cast<char*>(frLabel.Data()));
  simFrame->frame = 0;                                     
  simFrame->run = -1;                                      
  simFrame->dt = length;                                   
  simFrame->GTimeS = gps;                                  
  simFrame->GTimeN = 0;                                    

  cout << "Size (sec) " << x.size()/x.rate() << endl;
  FrProcData* proc = FrProcDataNew(simFrame,const_cast<char*>(chName.Data()),x.rate(),x.size(),-64);
  if(proc == NULL) {cout << "CWB::Toolbox::WriteFrameFile - Cannot create FrProcData" << endl; exit(-1);}
  proc->timeOffset = 0;                                                                              
  proc->tRange = simFrame->dt;                                                                       
  proc->type = 1;   // Time Serie                                                                    

  for (int i=0;i<(int)proc->data->nData;i++) proc->data->dataD[i] = x[i];

  int err=FrameWrite(simFrame,ofp);
  if (err) {cout << "CWB::Toolbox::WriteFrameFile - Error writing frame" << endl;exit(1);}
  FrameFree(simFrame);                                                                

  if (ofp!=NULL) FrFileOEnd(ofp);

  return frFile; 
}                  

//______________________________________________________________________________
wavearray<double> 
CWB::Toolbox::GetPhase(wavearray<double> hi, wavearray<double> hr, wavearray<double>& fi, wavearray<double>& fr, 
                                                                   wavearray<double>& s,  wavearray<double>& tt) {
//
// measurement of phase sync  between two waveforms
// hi - whitened CBC PE template
// hr - cWB whitened reconstructed signal
// The cWB and CBC waveforms are synchronized and the following data is calculated
//     as a function of number of cycles from the nominal merger time
// fi - frequency evolution of the CBC template
// fr - frequency evolution of the CWB waveform
// tt - waveform evolution time [sec]
//  s - accumulated cWB SNR
//  p - phase difference between the cWB and CBC waveforms
//  q - phase difference between the CBC L and H detectors
//

   double R,a;
   wavearray<double> x,y,p,q,r;
   wavearray<double> wi,wr;
   wavearray<double> HI,HR,WI,WR;
   double thi,thr,T;
   double t,dt,am;
   int j=0;
   int l=0;

   thi=hi.start();

   T=hi.start();
   if(thr>T) T=thr;
   if(thi>T) T=thi;
  
   hr.start(hr.start()-T); hi.start(hi.start()-T);
   R = int(hr.rate()+0.01);

   am = 0;
   for(int i=0; i<hi.size(); i++) {
      t = hi.start()+i/R;
      if(t<0) continue;
      if(fabs(hi.data[i])>am) {am=fabs(hi.data[i]); T=t;}
   }
   int Lt = int(R*(int(T)+1)+0.01);

   HR=hr; HI=hi;

   WI=hi; WR=hr; 
   wi=WI; wr=WR;

   for(int m=HI.size()-1; m>0; m--) {
      if(hi.data[m]>0) hi.data[m]=1; else hi.data[m]=-1;
      if(hr.data[m]>0) hr.data[m]=1; else hr.data[m]=-1;
      if(wi.data[m]>0) wi.data[m]=1; else wi.data[m]=-1;
      if(wr.data[m]>0) wr.data[m]=1; else wr.data[m]=-1;
   }
   
   int nl,nr,mm;

   tt=HI; tt=0;
   fi=HI; fi=0;
   fr=HI; fr=0;
   p=HI;  p=0;
   q=HI;  q=0;
   s=HI;  s=0;
   double Er,Ei;
   j=0; dt=1; am=0; l=1; nl=0;
   for(int m=HI.size()-1; m>0; m--) {
      t = wi.start()+m/wi.rate();
      if(t<=0 || t>T) continue;
      if(j==0 && am==0.) am=wi.data[m];
      if(wi.data[m]*am <= 0) {
         nr = -wi.data[m]*wr.data[m];
         am=wi.data[m];
         tt.data[j]=t;
         p.data[j] /= l; 
         q.data[j] /= l; l=0;
         fi.data[j]=0.5/dt; dt=0.;
         mm=m; while(wr.data[mm++]*wr.data[m]>0) fr.data[j]+=1;
         mm=m; while(wr.data[mm--]*wr.data[m]>0) fr.data[j]+=1;
         fr.data[j] = 0.5*HI.rate()/fr.data[j];
         if(j>0) s.data[j]=s.data[j]+s.data[j-1];
         j++;
      }
      a = wr.data[m]*wi.data[m];
      if(a<0) p.data[j]+=a*nr;
      q.data[j]+=hi.data[m]; l++;
      s.data[j]+=HR.data[m]*HR.data[m]; 
      dt += 1./HI.rate();
   }

   fi.resize(j); fi.rate(2);
   fr.resize(j); fr.rate(2);
   tt.resize(j); tt.rate(2);

   p.resize(j); p.rate(2); p.data[0]=1.; p.data[1]=-2.;
   q.resize(j); q.rate(2);
   s.resize(j); s.rate(2);

   for(int j=0;j<s.size();j++) s[j]=sqrt(s[j]);

   for(int m=p.size()-2; m>0; m--) {
      if(p.data[m]-p.data[m+1] > 1) mm=2; else mm=0;
      p.data[m]-=mm;
   }

   for(int j=0;j<p.size();j++) {
     if(p[j]<-1) p[j]+=2;
     if(p[j]>1)  p[j]-=2;
   }

   return p;
}

//_______________________________________________________________________________________
Int_t
CWB::Toolbox::getTableau10BlindColor( Int_t index ) {
//
// @brief Returns a color out of the Tableau color blind set
// @details Uses a numerical index to determine which color to return
//
// @param index  0-9 for the different colors
// @return Int_t Reference to the ROOT color number

    switch ( index ) {
        case 0:  return TColor::GetColor(  0, 107, 164);
        case 1:  return TColor::GetColor(255, 128,  14);
        case 2:  return TColor::GetColor(171, 171, 171);
        case 3:  return TColor::GetColor( 89,  89,  89);
        case 4:  return TColor::GetColor( 95, 158, 209);
        case 5:  return TColor::GetColor(200,  82,   0);
        case 6:  return TColor::GetColor(137, 137, 137);
        case 7:  return TColor::GetColor(163, 200, 236);
        case 8:  return TColor::GetColor(255, 188, 121);
        case 9:  return TColor::GetColor(207, 207, 207);
        default: return kBlack;
    }
}

//_______________________________________________________________________________________
Int_t
CWB::Toolbox::getTableau10BlindColor( TString name ) {
//
// @brief Returns a color out of the Tableau color blind set
// @details Uses a string representation to determine which color to return
// to increase readability of the code. Names represent similar colors from
// LaTeX' xcolor package with svgnames or x11names (if there is a number at
// end of the name).
//
// @param name   Name of the color
// @return Int_t Reference to the ROOT color number
 
    if      ( name == "DeepSkyBlue4" ) return getTableau10BlindColor(0);
    else if ( name == "DarkOrange1" )  return getTableau10BlindColor(1);
    else if ( name == "DarkGray" )     return getTableau10BlindColor(2);
    else if ( name == "DimGray" )      return getTableau10BlindColor(3);
    else if ( name == "SkyBlue3" )     return getTableau10BlindColor(4);
    else if ( name == "Chocolate3" )   return getTableau10BlindColor(5);
    else if ( name == "Gray" )         return getTableau10BlindColor(6);
    else if ( name == "SlateGray1" )   return getTableau10BlindColor(7);
    else if ( name == "SandyBrown" )   return getTableau10BlindColor(8);
    else if ( name == "LightGray" )    return getTableau10BlindColor(9);
    else return kBlack;
}

//_______________________________________________________________________________________
void
CWB::Toolbox::getTableau10BlindColorPalette( const int nSteps, Int_t *palette ) {
//
// @brief Generate a palette for color blinds
// @details Use the tableau 10 blind palette to generate a palette for
// ROOT's 2D histograms.
//
// @param const int   Number of steps to generate
// @param Int_t *     Pointer to the array storing the palette colors

    // define the colors
    const int nColors = 4;
    TColor * col[nColors] = {
        gROOT->GetColor(getTableau10BlindColor("SlateGray1")),
        gROOT->GetColor(getTableau10BlindColor("DeepSkyBlue4")),
        gROOT->GetColor(getTableau10BlindColor("SandyBrown")),
        gROOT->GetColor(getTableau10BlindColor("Crimson"))
    };
    Double_t stop[nColors] = {0., .25, .65, 1.0};

    // get color values
    Double_t r[nColors], g[nColors], b[nColors];
    for ( int c = 0; c < nColors; c++ ) {
        r[c] = col[c]->GetRed();
        g[c] = col[c]->GetGreen();
        b[c] = col[c]->GetBlue();
    }

    // generate palette
    Int_t FI = TColor::CreateGradientColorTable(nColors, stop, r, g, b, nSteps);
    for ( int i = 0; i < nSteps; i++ ) palette[i] = FI+i;
}

