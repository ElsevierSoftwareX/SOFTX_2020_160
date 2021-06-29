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


#include "frame.hh"
#include "Meyer.hh"

////////////////////////////////////////////////////////////////////////////////
/* BEGIN_HTML
The frame class is a partial c++ wrapper to the c frame library<br> 
framelib : (http://lappweb.in2p3.fr/virgo/FrameL/)<br>
it is designed to implements only the functions used by cWB<br> 
<br>
the following are some examples which illustrate how to use the class

<h3><a name="example1">Example 1</a></h3>
This example shows how to read a frame from frame file list<br>
using the frfile structure
<br><br>
<pre>
    root[] #define FRLIST_NAME  "input/H1_LDAS_C02_L2.frl"	// frame files list name
    root[] #define CHANNEL_NAME "H1:LDAS-STRAIN"		// channel to be read
    root[] double start = 942449664;				// start time
    root[] double stop = start+600;				// stop time
    root[] wavearray<double> x;						// array 
    root[] CWB::frame fr(FRLIST_NAME);				// open frame file list
    root[] int nfiles=fr.getNfiles();				// get number of files in frame list
    root[] cout << "nfiles : " << nfiles << endl;		
    root[] frfile FRF = fr.getFrList(start, stop, 0);		// fill the frfile structure
    root[] fr.readFrames(FRF,CHANNEL_NAME,x);			// read data in x array 
</pre>
<h3><a name="example2">Example 2</a></h3>
This example shows how to read a frame from frame file list using the read operator ">>"
<br><br>
<pre>
    root[] #define FRLIST_NAME  "input/H1_LDAS_C02_L2.frl"	// frame files list name
    root[] #define CHANNEL_NAME "H1:LDAS-STRAIN"		// channel to be read
    root[] wavearray<double> x;						// array 
    root[] x.start(942449664);					// start time
    root[] x.stop(x.start()+600);				// stop time
    root[] frame ifr(FRLIST_NAME,CHANNEL_NAME,"READ");		// open frame file list
    root[] ifr >> x;						// read data in x array 
    root[] cout << x.start() << " " << x.rate() << " " << x.size() << endl;	// print data info
</pre>
<h3><a name="example3">Example 3</a></h3>
This example shows how to write data to frame file using the write operator "<<"
<br><br>
<pre>
    root[] #define CHANNEL_NAME "H1:LDAS-STRAIN"		// channel to be read
    root[] #define FRNAME "TEST"				// output frame name
    root[] wavearray<double> x;						// array 
    root[] x.start(942449664);					// set start time
    root[] x.stop(x.start()+600);				// set stop time
    root[] x.rate(16384);					// set rate (Hz)
    root[] for(int i=0;i&lt;x.size();i++) x[i]=i;			// fill array
    root[] frame ofr("TEST.gwf","CHNAME","WRITE");		// open output frame file
    root[] ofr.setFrName(FRNAME);				// set frame name
    root[] ofr << x;						// read data in x array 
    root[] ofr.close();						// close output frame file
</pre>
<h3><a name="example4">Example 4</a></h3>
This example shows how to resample read data
<br><br>
<pre>
    root[] #define FRLIST_NAME  "input/H1_LDAS_C02_L2.frl"	// frame files list name
    root[] #define CHANNEL_NAME "H1:LDAS-STRAIN"		// channel to be read
    root[] wavearray<double> x;						// array 
    root[] x.start(942449664);					// start time
    root[] x.stop(x.start()+600);				// stop time
    root[] frame ifr(FRLIST_NAME,CHANNEL_NAME,"READ");		// open frame file list
    root[] ifr.setSRIndex(10);					// read data are resampled @ 1024 Hz
    root[] ifr >> x;						// read data in x array 
    root[] cout << x.start() << " " << x.rate() << " " << x.size() << endl;	// print data info
</pre>


END_HTML */
////////////////////////////////////////////////////////////////////////////////

ClassImp(CWB::frame)

//______________________________________________________________________________
CWB::frame::frame() : frtree_List(NULL), chName(""), frName("Frame"), 
                      rfName(""), frFile(NULL), fOption(""), srIndex(-1), 
                      verbose(false), frRetryTime(60), xstart(0.), xstop(0.) {
// default constructor
// initialize only the frame class parameters 
//
}

//______________________________________________________________________________
CWB::frame::frame(TString ioFile, TString chName, Option_t* option, bool onDisk, TString label, unsigned int mode) : 
  frtree_List(NULL), chName(""), frName("Frame"), rfName(""), 
  frFile(NULL), fOption(""), srIndex(-1), verbose(false), frRetryTime(60), xstart(0.), xstop(0.) {
//
// initialize frame class parameters & open frame file list
//
// see the 'CWB::frame::open' method description
//
  open(ioFile, chName, option, onDisk, label, mode);
}

//______________________________________________________________________________
CWB::frame::~frame() {
//
// destructor
//

  // remove temporary root file
  if(rfName!="") gSystem->Exec(TString("rm "+rfName).Data());
  close();
}

//______________________________________________________________________________
wavearray<double>& operator >> (CWB::frame& fr, wavearray<double>& x) {
//
// operator : read frame into wavearray x
//

  fr.readFrames(x);
  return x;
}

//______________________________________________________________________________
CWB::frame& operator << (CWB::frame& fr, wavearray<double>& x) {
//
// operator : write wavearray x to frame
//

  if(fr.getOption()=="READ") {
    cout << "CWB::frame::operator(<<) : Not allowed in READ mode" << endl;
    exit(1);
  }

  if(fr.getChName()=="") {
    cout << "CWB::frame::operator(<<) : channel not defined";
    exit(1);
  }

  if(fr.getFrName()=="") {
    cout << "CWB::frame::operator(<<) : frame name not defined" << endl;
    exit(1);
  }

  fr.writeFrame(x, fr.getFrName(), fr.getChName()); 
  return fr;
}

//______________________________________________________________________________
CWB::frame& operator >> (wavearray<double>& x, CWB::frame& fr) {
//
// operator : write wavearray x to frame
//

  fr << x;
  return fr;
}

//______________________________________________________________________________
void
CWB::frame::writeFrame(wavearray<double> x, TString frName, TString chName) {
//
// write in frame the data contained in the wavearray x
// frame is saved in the frFile (object must be opened in write mode)
//
// x             : in  - data array which contains the data 
// chName        : in  - name of channel to be extracted
// frName        : in  - frame name  
//

  if (frFile==NULL) {
    cout << "CWB::frame::writeFrame : output frame file close" << endl;
    exit(1);
  }
  if(x.rate()<=0) {
    cout << "CWB::frame::writeFrame : input rate must be > 0" << endl;
    exit(1);
  }
  if(x.size()==0) {
    cout << "CWB::frame::writeFrame : input size must be > 0" << endl;
    exit(1);
  }

  /*----------------------- Create a new frame ---------*/
  FrameH* frFrame = FrameNew(const_cast<char*>(frName.Data()));
  frFrame->frame = 0;
  frFrame->run = -1;
  frFrame->dt = x.size()/x.rate();
  frFrame->GTimeS = x.start();
  frFrame->GTimeN = 0;

  cout << "Size (sec) " << x.size()/x.rate() << endl;
  FrProcData* proc = FrProcDataNew(frFrame,const_cast<char*>(chName.Data()),x.rate(),x.size(),-64);
  if(proc == NULL) {
    cout << "CWB::frame::writeFrame : Cannot create FrProcData" << endl;
    exit(1);
  }
  proc->timeOffset = 0;
  proc->tRange = frFrame->dt;
  proc->type = 1;   // Time Serie

  for (int i=0;i<(int)proc->data->nData;i++) proc->data->dataD[i] = x[i];

  int err=FrameWrite(frFrame,frFile);
  if (err) {
    cout << "CWB::frame::writeFrame : Error writing frame" << endl;
    exit(1);
  }
  FrameFree(frFrame);

  return;
}

//______________________________________________________________________________
void
CWB::frame::open(TString ioFile, TString chName, Option_t* option, bool onDisk, TString label, unsigned int mode) {
//
// open frame file in read or write mode
//
// ioFile        : in  - if option="READ"  -> ioFile = input  frame file list name 
//                 out - if option="WRITE" -> ioFile = output frame file name 
// chName        : in  - name of the data channel (def="")
// option        : in  - READ/WRITE mode  (def="")
// onDisk        : in  - true  -> tmp tree on root file 
//                       false (def) -> tmp tree is stored in ram
//                       true must be used only if the input file frame list is huge
// label         : in  - this is a string used to select files listed in the input frame file list
//                       def = ".gwf" 
// mode          : in  - if mode=0 (default) the frame file list is a list of frame file paths 
//
//                       has the following format :
//                       DIR_NAME/XXXXX-GPS-LEN.gwf where :
//                       DIR_NAME : is the directory which contains the frame file.
//                                  The strings 'file://localhost' or 'framefile=' at the 
//                                  beginning of the path will be automatically removed
//                       GPS : is the start frame time in sec
//                       LEN : is the frame length time in sec
//
//                       if mode=1 the frame file list is a list of line with the following format :
//
//                       FILE_PATH GPS LEN 
//

  if(rfName!="") gSystem->Exec(TString("rm "+rfName).Data());
  close();

  this->chName=chName; 
  this->frName="Frame"; 

  // if onDisk store tree on root file otherwise it is stored in ram
  if(onDisk) {
    // create temporary file
    gRandom->SetSeed(0);
    int rnID = int(gRandom->Rndm(13)*1.e9);
    UserGroup_t* uinfo = gSystem->GetUserInfo();
    TString uname = uinfo->fUser;
    // create temporary dir
    gSystem->Exec(TString("mkdir -p /dev/shm/")+uname);
    char fName[1024]="";
    sprintf(fName,"/dev/shm/%s/cwbframe_%d.root",uname.Data(),rnID);
    rfName=fName;
  }

  fOption = option;
  fOption.ToUpper();
  if((fOption!="READ")&&(fOption!="WRITE")) fOption="READ";

  if(fOption=="READ")  nfiles = frl2FrTree(ioFile, rfName, label, mode);
  if(fOption=="WRITE") {
    if(!fNameCheck(ioFile)) exit(1);
    frFile = FrFileONew(const_cast<char*>(ioFile.Data()),1);   // gzip compression
    if(frFile==NULL) {
      cout << "CWB::frame::frame : Error opening file " << ioFile.Data() << endl;
      exit(1);
    }
  }
  return;
}

//______________________________________________________________________________
void
CWB::frame::close() {
//
// close frame file and reset parameters
//

  if(frtree_List!=NULL) frtree_List->Reset();
  frtree_List=NULL;
  if (frFile!=NULL) FrFileOEnd(frFile);
  frFile=NULL;
  fOption="";
  chName="";
  frName="";
}

//______________________________________________________________________________
int 
CWB::frame::frl2FrTree(TString iFile, TString rfName, TString label, unsigned int mode) {
//
// convert input frame file list to a tree 
// note : this method can be used only in READ mode
//
// iFile         : in  - input frame file list name 
// rfName        : in  - name of auxiliary file name used to store the temporary tree
//                       if rfName=="" the tree is saved & sorted to the private frtree_List
// label         : in  - this is a string used to select files listed in the input frame file list
//                       def = ".gwf" 
// mode          : in  - 0/1 (def=0) 
//

  if(fOption!="READ") {
    cout << "CWB::frame::frl2FrTree : allowed only in READ mode" << endl;
    exit(1);
  }

  TFile* ofile = NULL;
  if(rfName.Sizeof()>1) ofile=new TFile(rfName.Data(),"RECREATE");

  TTree* frtree = new TTree("frl","frl");
  double gps;
  double length;
  double start;
  double stop;
  char atlas_path[1024];
  frtree->Branch("gps",&gps,"gps/D");
  frtree->Branch("length",&length,"length/D");
  frtree->Branch("start",&start,"start/D");
  frtree->Branch("stop",&stop,"stop/D");
  frtree->Branch("path",atlas_path,"path/C");

  if(rfName.Sizeof()>1) frtree->SetDirectory(ofile);
  else                  frtree->SetDirectory(0);

  char ifile_name[1024];
  sprintf(ifile_name,"%s",iFile.Data());
  cout << "CWB::frame::frl2FrTree - " << ifile_name << endl;
  FILE *in;
  char s[512];
  char f[512];
  TString idir_name;
  cout.precision(14);
  int cnt=0;
  // read frame file list
  if( (in=fopen(ifile_name,"r"))!=NULL ) {
      while(fgets(s,512,in) != NULL) {
          sprintf(f,"%s",s);
          //cout << f << endl;
          TString line(f);
          // remove header created with gw_data_find
          line.ReplaceAll("framefile=","");
          line.ReplaceAll("file://localhost","");
          line.ReplaceAll("gsiftp://ldr.aei.uni-hannover.de:15000","");
          //cout << line.Data() << endl;
          if (line.Contains(label)) {
            // select the first token of the string
            TObjArray* token0 = TString(line).Tokenize(TString(" "));
            TObjString* sline = (TObjString*)token0->At(0);
            TString path = sline->GetString().Data();

            //cout << f << endl;
            sprintf(atlas_path,"%s",path.Data());

            if(mode==0) {
              TObjArray* token1 = TString(path).Tokenize(TString("/"));
              TObjString* path_tok1 = (TObjString*)token1->At(token1->GetEntries()-1);
              //cout << path_tok1->GetString().Data() << endl;

              TObjArray* token2 = TString(path_tok1->GetString()).Tokenize(TString("-"));
              TObjString* gps_tok = (TObjString*)token2->At(token2->GetEntries()-2);
              TString sgps = gps_tok->GetString().Data();
              //cout << sgps.Data() << endl;
              gps = sgps.Atof();

              TObjString* length_tok = (TObjString*)token2->At(token2->GetEntries()-1);
              TString slength = length_tok->GetString().Data();
              slength.ReplaceAll(label,"");
              //cout << slength.Data() << endl;
              length = slength.Atof();

              delete token1;
              delete token2;
            } else {
              TObjArray* token1 = TString(line).Tokenize(TString(" "));

              if(token1->GetEntries()<3) {
                cout << "CWB::frame::frl2FrTree : bad format " << ifile_name << endl;
                exit(1);
              }

              TObjString* gps_tok = (TObjString*)token1->At(1);
              TString sgps = gps_tok->GetString().Data();
              //cout << sgps.Data() << endl;
              gps = sgps.Atof();

              TObjString* length_tok = (TObjString*)token1->At(2);
              TString slength = length_tok->GetString().Data();
              //cout << slength.Data() << endl;
              length = slength.Atof();

              delete token1;
            }
            delete token0;

            start = gps;

            stop = start+length;

            // skip file if it is out of time range
            if(xstart && stop<xstart) continue;  
            if(xstop  && start>xstop) continue;

            frtree->Fill();

            //cout << "Frame info : " << gps << " " << atlas_path << " " << length << " " << start << " " << stop << endl;
            //if (++cnt%100==0) cout << cnt << endl;
            ++cnt;
          }
      }
  } else {
    cout << "CWB::frame::frl2FrTree : Error opening file " << ifile_name << endl;
    exit(1);
  }

  if(ofile!=NULL) {
    frtree->Write();
    ofile->Close();
  } else { // copy tree to internal frtree
    if(frtree_List!=NULL) delete frtree_List;
    frtree_List = (TTree*)frtree->CloneTree();
    sortFrTree();
    delete frtree;
  }
  return cnt;
}

//______________________________________________________________________________
int 
CWB::frame::sortFrTree(TString iFile, TString rfName) {
//
// sort the tree contained in iFile according the GPS time and save it to the auxiliary file rfName
// note : this method can be used only in READ mode
//
// iFile         : in  - name of input root file containing tree to be sorted 
// rfName        : in  - name of auxiliary file name used to store the temporary tree
//                       if rfName=="" the tree is saved & sorted to the private frtree_List
//

  if(fOption!="READ") {
    cout << "CWB::frame::sortFrTree : allowed only in READ mode" << endl;
    exit(1);
  }
  TFile f(iFile.Data());
  TTree *tree = (TTree*)f.Get(LST_TREE_NAME);
  Int_t nentries = (Int_t)tree->GetEntries();
  tree->Draw("gps","","goff");
  Int_t *index = new Int_t[nentries];
  TMath::Sort(nentries,tree->GetV1(),index,false);

  TFile f2(rfName.Data(),"recreate");
  TTree *tsorted = (TTree*)tree->CloneTree(0);
  for (Int_t i=0;i<nentries;i++) {
    //if (i%1000==0) cout << i << endl;
    tree->GetEntry(index[i]);
    tsorted->Fill();
  }
  tsorted->Write();
  f2.Close();
  delete [] index;

  return nentries;
}

//______________________________________________________________________________
int 
CWB::frame::sortFrTree() {
//
// sort the the auxiliary tree according the GPS time
// note : this method can be used only in READ mode
//

  if(fOption!="READ") {
    cout << "CWB::frame::sortFrTree : allowed only in READ mode" << endl;
    exit(1);
  }
  Int_t nentries = (Int_t)frtree_List->GetEntries();
  frtree_List->Draw("gps","","goff");
  double* gps = frtree_List->GetV1(); 
  Int_t *index = new Int_t[nentries];
  TMath::Sort(nentries,gps,index,false);

  double pgps = -1;
  TTree *sorted_frtree_List = (TTree*)frtree_List->CloneTree(0);
  for (Int_t i=0;i<nentries;i++) {
    //if (i%1000==0) cout << i << endl;
    int j = index[i];		// j = sorted index
    if(gps[j]==pgps) continue;  // skip duplicated entries
    pgps=gps[j];
    frtree_List->GetEntry(j);
    sorted_frtree_List->Fill();
  }
  delete [] index;

  // copy sorted tree to the original not sorted tree
  delete frtree_List;
  frtree_List = (TTree*)sorted_frtree_List->CloneTree();
  delete sorted_frtree_List;

  return nentries;
}

//______________________________________________________________________________
frfile 
CWB::frame::getFrList(int istart, int istop, int segEdge) {
//
// return the list of frfile structures of the frame files 
// cointained in the range [istart,istop]
//
// istart        : in  - start time sec
// istop         : in  - stop  time sec
// segEdge       : in  - wavelet boundary offset [sec]
//

  frfile frf;
  if(rfName!="") {
    frf = getFrList(rfName, istart, istop, segEdge);
  } else {
    frf = getFrList(istart, istop, segEdge, NULL);
  }
  return frf;
}

//______________________________________________________________________________
frfile 
CWB::frame::getFrList(TString rfName, int istart, int istop, int segEdge) {
//
// return the frame list of frames contained in the range [start-segEdge,stop-segEdge]
// note : this method can be used only in READ mode
//
// rfName        : in  - name of auxiliary file name used to store the temporary tree
//                       use the tree contained in rfName root file
// istart        : in  - start time sec
// istop         : in  - stop  time sec
// segEdge       : in  - wavelet boundary offset [sec]
//

  if(fOption!="READ") {
    cout << "CWB::frame::getFrList : allowed only in READ mode" << endl;
    exit(1);
  }

  frfile frf;

  TFile *ifile = TFile::Open(rfName.Data());
  if (ifile==NULL) {cout << "No " << rfName.Data() << " file !!!" << endl;exit(1);}
  TTree* itree = (TTree *) gROOT->FindObject(LST_TREE_NAME);
  if (itree==NULL) {cout << "No ffl tree name match !!!" << endl;exit(1);}
  frf=getFrList(istart, istop, segEdge, itree);
  delete itree;
  delete ifile;

  return frf;
}

//______________________________________________________________________________
vector<frfile> 
CWB::frame::getFrList(int istart, int istop) {
//
// return the list of frfile structures of the frame files 
// cointained in the range [istart,istop]
// note : this method can be used only in READ mode
//
// istart        : in  - start time sec
// istop         : in  - stop  time sec
//

  if(fOption!="READ") {
    cout << "CWB::frame::getFrList : allowed only in READ mode" << endl;
    exit(1);
  }

  vector<frfile> frlist;
  frfile frf;

  TTree* itree=frtree_List;
  if(itree==NULL) {cout << "CWB::frame::getFrList - Error : No tree defined !!!" << endl;exit(1);}

  int size = itree->GetEntries();
  //cout << size << endl;

  double gps;
  double length;
  char path[1024];
  itree->SetBranchAddress("gps",&gps);
  itree->SetBranchAddress("length",&length);
  itree->SetBranchAddress("path",path);

  double start = istart!=0 ? istart : 0;
  double stop  = istop!=0  ? istop  : 2000000000;
  char cut[1024];
  sprintf(cut,"(gps>=%f && gps<=%f)", start,stop);
  //cout << cut << endl;
  itree->SetEstimate(size);
  itree->Draw("Entry$",cut,"goff");
  int isel = itree->GetSelectedRows();
  double* entry = itree->GetV1();

  for (int i=0;i<isel;i++) {
    itree->GetEntry(i);
    frf.start=gps;
    frf.stop=gps+length;
    frf.length=length;
    frf.file.push_back(path);
    frlist.push_back(frf);
    frf.file.clear();
  }

  return frlist;
}

//______________________________________________________________________________
frfile 
CWB::frame::getFrList(int istart, int istop, int segEdge, TTree* itree) {
//
// return the frame list of frames contained in the range [start-segEdge,stop-segEdge]
// note : this method can be used only in READ mode
//
// istart        : in  - start time sec
// istop         : in  - stop  time sec
// segEdge       : in  - wavelet boundary offset [sec]
// itree         : in  - if itree=NULL use the private frtree_List
//

  if(fOption!="READ") {
    cout << "CWB::frame::getFrList : allowed only in READ mode" << endl;
    exit(1);
  }

  frfile frf;

  if(itree==NULL) itree=frtree_List;
  if(itree==NULL) {cout << "CWB::frame::getFrList - Error : No tree defined !!!" << endl;exit(1);}

  int size = itree->GetEntries();
  //cout << size << endl;

  double gps;
  double length;
  char path[1024];
  itree->SetBranchAddress("gps",&gps);
  itree->SetBranchAddress("length",&length);
  itree->SetBranchAddress("path",path);

  double start=istart-segEdge;
  double stop= istop+segEdge;
  char cut[1024];
  sprintf(cut,"((%f>=gps && %f<=gps+length) || (%f>=gps && %f<=gps+length) || (%f<=gps && %f>=gps+length))",
          start,start,stop,stop,start,stop);
  //cout << cut << endl;
  itree->SetEstimate(size);
  itree->Draw("Entry$",cut,"goff");
  int isel = itree->GetSelectedRows();
  double* entry = itree->GetV1();


  // CHECK
  if (isel==0) {
    cout << endl; 
    cout << "CWB::frame::getFrList - Error : No files matched in the range " 
         << start << ":" << stop << " !!!" << endl << endl;
    sprintf(cut,"gps>=%f && gps<=%f",start-500,stop+500);
    cout << "Dump frame segments available in the input list in the range : " << cut << endl << endl;
    itree->SetScanField(0);
    itree->Scan("Entry$:gps:length",cut,"goff");
    exit(1);
  }
  if (isel==1) {
    itree->GetEntry(entry[0]);
    if (start < gps) {
      cout << endl;
      cout << "CWB::frame::getFrList - Error : " << " No files matched !!!" << endl;
      cout << "start buffer (" << start << ") < begin gps frame (" << gps << ")" << endl;
      cout << "check input frame list" << endl << endl;
      exit(1);
    }
    if (stop > gps+length) {
      cout << endl;
      cout << "CWB::frame::getFrList - Error : " << " No files matched !!!" << endl;
      cout << "stop buffer (" << stop << ") > end gps frame (" << gps+length << ")" << endl;
      cout << "check input frame list" << endl << endl;
      exit(1);
    }
  }
  if (isel>1) {
    itree->GetEntry(entry[0]);
    if (start < gps) {
      cout << endl;
      cout << "CWB::frame::getFrList - Error : " << " No files matched !!!" << endl;
      cout << "start buffer (" << start << ") < begin gps frame (" << gps << ")" << endl;
      cout << "check input frame list" << endl << endl;
      exit(1);
    }
    itree->GetEntry(entry[isel-1]);
    if (stop > gps+length) {
      cout << endl;
      cout << "CWB::frame::getFrList - Error : " << " No files matched !!!" << endl;
      cout << "stop buffer (" << stop << ") > end gps frame (" << gps+length << ")" << endl;
      cout << "check input frame list" << endl << endl;
      exit(1);
    }
    // check if frames are contiguous
    itree->GetEntry(entry[0]);
    double stop0=gps+length;
    for(int i=1;i<isel;i++) {
      itree->GetEntry(entry[i]);
      if(gps!=stop0) {
        cout << "CWB::frame::getFrList - Error : " << " Requested frames are not contiguous !!!" << endl;
        cout << "start : " << gps << " is != from the end of the previous segment " << stop0 << endl;
        exit(1);
      }
      stop0=gps+length;
    }
  }

  itree->GetEntry(entry[isel-1]);
/*
  cout << isel << endl;
  cout << int(length) << endl;
  cout << int(start) << endl;
  cout << int(stop) << endl;
*/
  frf.start=start;
  frf.stop=stop;
  frf.length=length;
  frf.file.clear();
  for (int i=0;i<isel;i++) {
    itree->GetEntry(entry[i]);
    frf.file.push_back(path);
  }

  return frf;
}

//______________________________________________________________________________
waveSegment
CWB::frame::getFrRange(TTree* itree) {
//
// return the begin and end range of the frame list
// note : this method can be used only in READ mode
//
// itree         : in  - if itree=NULL use the private frtree_List
//

  if(fOption!="READ") {
    cout << "CWB::frame::getFrRange : allowed only in READ mode" << endl;
    exit(1);
  }

  waveSegment SEG={0,0.,0.};

  if(itree==NULL) itree=frtree_List;
  if(itree==NULL) {cout << "CWB::frame::getFrRange - No tree defined !!!" << endl;exit(1);}

  int size = itree->GetEntries();

  double gps;
  double length;
  itree->SetBranchAddress("gps",&gps);
  itree->SetBranchAddress("length",&length);

  double min=4000000000.;
  double max=0.;
  for (int i=0;i<size;i++) {
    itree->GetEntry(i);
    if(min>gps) min=gps;
    if(max<gps+length) max=gps+length;
  }

  SEG.start=min;
  SEG.stop=max;

  return SEG;
}

//______________________________________________________________________________
void 
CWB::frame::readFrames(wavearray<double>& w) {
//
// read in w the channel data contained in the frames
// the range of data to be read is defined in the wavearray object
// user must initialize the object with : 
// - w.start(GPS)
// - w.stop(GPS+LEN)
//
// note : this method can be used only in READ mode
//
// w             : out - data array which contains the data read
//

  if(getOption()!="READ") {
    cout << "CWB::frame::readFrame  : allowed only in READ mode" << endl;
    exit(1);
  }
  if(getChName()=="") {
    cout << "CWB::frame::readFrame  : channel not defined" << endl;
    exit(1);
  }
  frfile frf = getFrList(w.start(), w.stop(), 0);
  readFrames(frf,const_cast<char*>(getChName().Data()),w);
  return;
}

//______________________________________________________________________________
void 
CWB::frame::readFrames(char* filename, char *channel, wavearray<double> &w) {
//
// read in w the channel data contained in the frames defined in the filename
// this method is obsolete, was used by the wat-5.3.9 pipeline infrastructure 
// the input file must contains the following informations :
//
// - nfiles              // number of files to be read 
// - frame length        // dummy, not used anymore
// - start time          // (sec)
// - stop time           // (sec)
// - rate                // dummy, not used anymore
// - frame file path     // this field is repeated nfiles times
//
// note : this method can be used only in READ mode
//
// filename      : in  - name of frame file
// channel       : in  - name of channel to be extracted
// w             : out - data array which contains the data read
//

  if(fOption!="READ") {
    cout << "CWB::frame::readFrames : allowed only in READ mode" << endl;
    exit(1);
  }

  int buffersize=1024;
  char *buffer=new char[buffersize];

  double rf_start, rf_stop;
  double frlen;
  int nframes;
  double rf_rate;

  vector<frfile> frlist;
  frfile frf;

  FILE *framelist=fopen(filename,"r");

  fgets(buffer,buffersize,framelist);
  nframes=(int)atoi(buffer);

  if(nframes<=0) {
      fprintf(stderr,"nframes=%d\n",nframes);
      exit(1);
  }

  fgets(buffer,buffersize,framelist);
  frlen=(double)atof(buffer);//dummy, not used anymore
  fgets(buffer,buffersize,framelist);
  rf_start=(double)atof(buffer);
  fgets(buffer,buffersize,framelist);
  rf_stop=(double)atof(buffer);
  fgets(buffer,buffersize,framelist);
  rf_rate=(double)atof(buffer);//dummy, not used anymore

  if(verbose) {
    fprintf(stdout,"-------------\n");
    fprintf(stdout,"filename=%s channel=%s nframes=%d frlen=%f start=%f stop=%f rate=%f \n",
            filename, channel, nframes, frlen, rf_start, rf_stop, rf_rate);
  }

  frf.start=rf_start;
  frf.stop=rf_stop;
  frf.length=frlen;
  for(int i=0;i<nframes;i++) {
    fgets(buffer,buffersize,framelist);
    if(verbose) fprintf(stdout,"framefile=%s\n",buffer);
    frf.file.push_back(TString(buffer));
  }
  delete [] buffer;
  fclose(framelist);

  readFrames(frf, channel, w);

  return;
}

//______________________________________________________________________________
void 
CWB::frame::readFrames(frfile frf, char *channel, wavearray<double> &w) {
//
// read in w the channel data contained in the frames defined in the frf structure
// the range of data to be read is defined in the frf structure
//
// note : this method can be used only in READ mode
//
// frf           : in  - frame file structure 
// channel       : in  - name of channel to be extracted
// w             : out - data array which contains the data read
//

  if(fOption!="READ") {
    cout << "CWB::frame::readFrames : allowed only in READ mode" << endl;
    exit(1);
  }

  FrVect *v;
  FrFile *ifile;

  double rf_start, rf_stop, fstart, fend;
  double frlen;
  int nframes;
  double rf_rate;
  unsigned int wcounter=0;
  double rf_time=0;
  double rf_step=-1.0;

  nframes=frf.file.size();
  rf_start=frf.start;
  rf_stop=frf.stop;
  frlen=frf.length;

  // BIG FRAMES : when frlen>8000 (Virgo VSR1-V2) FrFileIGetV crash (GV)
  bool bigFrame = false;
  if (frlen>1024) {
    if(verbose) cout << "CWB::frame::readFrames - bigFrame !!! : " << frlen << endl;
    bigFrame = true;
  }

  double fend0=0;
  for(int i=0;i<nframes;i++) {

    if(verbose) fprintf(stdout,"CWB::frame::readFrames - framefile=%s channel=%s\n",frf.file[i].Data(),channel);

    // Read Loop 
    int ntry=0;
    ifile=NULL;
    while(ifile==NULL&&ntry<3) {
      ifile=FrFileINew(const_cast<char*>(frf.file[i].Data()));
      if(ifile!=NULL || frRetryTime==0) break;
      ntry++;
      cout << "CWB::frame::readFrames - Failed to read file - "  
           << "retry after " << ntry*frRetryTime << " sec ..." 
           << "\t[" << ntry << "/3]" << endl;
      if(ifile==NULL) gSystem->Sleep(1000*ntry*frRetryTime);  	// wait 
    }

    if(ifile==NULL) {
      fprintf(stderr,"CWB::frame::readFrames - Cannot read file %s\n",frf.file[i].Data());
      fflush(stderr);
      fflush(stdout);
      exit(1);
      // break;
    }

    if(ifile->error) {
      fprintf(stderr,"CWB::frame::readFrames - Cannot read file %s, error=%d\n",frf.file[i].Data(),ifile->error);
      fflush(stderr);
      fflush(stdout);
      exit(1);
      // break;
    }
    fstart=(double)FrFileITStart(ifile);
    fend=(double)FrFileITEnd(ifile);

    // check if frames are contiguous
    if(fend0!=0 && fstart!=fend0) {
      cout << "CWB::frame::readFrames - Error : " << " Requested frames are not contiguous !!!" << endl;
      cout << "start : " << fstart << " is != from the end of the previous segment " << fend0 << endl;
      exit(1);
    }
    fend0=fend;

    // BIG FRAMES CHECK (GV)
    if (bigFrame) {
      if (fstart<rf_start) fstart=rf_start;
      if (fend>rf_stop) fend=rf_stop;
    }
    frlen=fend-fstart;
    if(verbose) cout << "frlen : " << frlen << endl;
    if (frlen<=0) continue;

    rf_time=fstart;

    v=FrFileIGetV(ifile,channel, fstart, frlen);

    // BIG FRAMES CHECK (GV)
    if (bigFrame&&(v!=NULL)) {
      double time_offset = fstart-v->GTime;
      if (time_offset!=0) {
        fprintf(stdout,"CWB::frame::readFrames - BigFrame Check Mismatch : v->GTime=%.8f != fstart=%.8f \n",v->GTime,fstart);
        FrVectFree(v);
        v=FrFileIGetV(ifile,channel, fstart+time_offset, frlen);
      }
    }

    if(v==NULL) {
      fprintf(stderr,"CWB::frame::readFrames - FrFileIGetV failed: channel %s not present\n",channel);
      exit(1);
    }
    if((v->type!=FR_VECT_4R && v->type!=FR_VECT_8R)&&
       (v->type!=FR_VECT_2S && v->type!=FR_VECT_4S && v->type!=FR_VECT_8S)&&
       (v->type!=FR_VECT_2U && v->type!=FR_VECT_4U && v->type!=FR_VECT_8U)) {
      fprintf(stderr,"CWB::frame::readFrames - Wrong vector type %d\n",v->type);fflush(stderr);
      exit(1);
    }
    rf_step=v->dx[0];
    rf_rate=1./rf_step;
    w.resize(int((rf_stop-rf_start)*rf_rate));
    w.rate(rf_rate);
    w.start(rf_start);
    w.stop(rf_stop);
    long samples=long(frlen*rf_rate);

    if(verbose) {
      fprintf(stdout,"CWB::frame::readFrames - samples=%ld v->nDim=%d v->nData=%ld v->type=%d 8R=%d v->GTime=%.8f v->localtime=%d rate=%.1f step=%.8f\n",
            samples,(int)v->nDim,v->nData,v->type,FR_VECT_8R,v->GTime, v->localTime, rf_rate, rf_step);
      for(long j=0;j<v->nDim;j++) {
        fprintf(stdout,"v->nx[%ld]=%ld v->dx[%ld]=%.8f v->startX[%ld]=%.8f\n",
                j,v->nx[j],j,v->dx[j],j,v->startX[j]);
      }
    }

    // BIG FRAMES CHECK (GV)
    if (bigFrame&&(v->GTime!=fstart)) {fprintf(stdout,"CWB::frame::readFrames - Error : v->GTime=%.8f != fstart=%.8f \n",v->GTime,fstart);exit(1);}

    for(long j=0;j<samples;j++) {
      rf_time=fstart+j*rf_step;
      if(rf_time>=rf_start && rf_time<rf_stop) {
        if(wcounter<w.size()) {
          if(v->type==FR_VECT_2S) w.data[wcounter]=double(v->dataS[j]);
          if(v->type==FR_VECT_4S) w.data[wcounter]=double(v->dataI[j]);
          if(v->type==FR_VECT_8S) w.data[wcounter]=double(v->dataL[j]);
          if(v->type==FR_VECT_4R) w.data[wcounter]=double(v->dataF[j]);
          if(v->type==FR_VECT_8R) w.data[wcounter]=v->dataD[j];
          if(v->type==FR_VECT_2U) w.data[wcounter]=double(v->dataUS[j]);
          if(v->type==FR_VECT_4U) w.data[wcounter]=double(v->dataUI[j]);
          if(v->type==FR_VECT_8U) w.data[wcounter]=double(v->dataUL[j]);
          wcounter++;
        } else {
          fprintf(stderr,"w overflow\n");
          exit(1);
        }
      }
    }

    FrVectFree(v);
    FrFileIEnd(ifile);
  }
  if(verbose) fprintf(stdout,"CWB::frame::readFrames - wcounter=%d w.size()=%d\n",(int)wcounter,(int)w.size());

  if(srIndex>=0 && w.rate()!=(1<<srIndex)) {
    w.Resample(1<<srIndex);	// resampling
    if(verbose) {
      fprintf(stdout,"--------------------------------\n");
      fprintf(stdout,"After resampling: rate=%f, size=%d, start=%f\n", w.rate(),(int)w.size(),w.start());
      fprintf(stdout,"--------------------------------\n");
    }
  }

  return;
}

//______________________________________________________________________________
int  
CWB::frame::dumpFrList(frfile frf, TString ofName, double sRate) {
//
// dump to file ofName the list of frames defined in the frf structure
// note : this method can be used only in READ mode
// - list file used for the old cWB analysis
// - output file format
//
//   size    : number of frame files 
//   length  : stop-start
//   start   : start time
//   stop    : stop time
//   sRate   : rate
//   file1   :  
//   .....   :  
//   fileN   :  
//
// frf           : in  - frfile structure
// ofName        : in  - name of the output file
// sRate         : in  - sample rate of frame data
//

  if(fOption!="READ") {
    cout << "CWB::frame::dumpFrList : allowed only in READ mode" << endl;
    exit(1);
  }

  char ofile_name[1024];
  sprintf(ofile_name,"%s",ofName.Data());
  cout << ofile_name << endl;
  ofstream out;
  out.open(ofile_name,ios::out);
  if (!out.good()) {cout << "Error Opening File : " << ofName.Data() << endl;exit(1);}
  out.precision(14);

  out << frf.file.size() << endl;
  out << frf.length << endl;
  out << frf.start << endl;
  out << frf.stop << endl;
  out << int(sRate) << endl;
  for (int i=0;i<(int)frf.file.size();i++) {
    out << frf.file[i];
  }

  out.close();

  return 0;
}

//______________________________________________________________________________
bool  
CWB::frame::fNameCheck(TString fName) {
//
// check if file has the correct format
// xxxx-GPS-LENGTH.yyy
//

  // remove file extension
  TObjArray* token0 = fName.Tokenize(TString("."));
  TObjString* body_tok = (TObjString*)token0->At(0);
  TString body = body_tok->GetString().Data();

  TObjArray* token = body.Tokenize(TString("-"));

  if(token->GetEntries()<2) {
    cout << "CWB::frame::fNameCheck : File Name Format Error - " << fName.Data() << endl;  
    cout << "GPS & LENGTH must be integers : xxxx-GPS-LENGTH.yyy" << endl;
    return false;
  }
  TObjString* gps_tok = (TObjString*)token->At(token->GetEntries()-2);
  TString sgps = gps_tok->GetString().Data();

  TObjString* length_tok = (TObjString*)token->At(token->GetEntries()-1);
  TString slength = length_tok->GetString().Data();

  if((slength.IsDigit()) && (sgps.IsDigit())) return true;

  cout << "CWB::frame::fNameCheck : File Name Format Error - " << fName.Data() << endl;  
  cout << "GPS & LENGTH must be integers : xxxx-GPS-LENGTH.yyy" << endl;
  return false;
}
