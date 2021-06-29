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


/**********************************************************
 * Package:      frame Class Library
 * File name:    frame.hh
 * Author:       Gabriele Vedovato (vedovato@lnl.infn.it)
 **********************************************************/


#ifndef CWBFRAME_HH
#define CWBFRAME_HH

#include "TROOT.h"
#include "TMath.h"
#include "TSystem.h"
#include "TString.h"
#include "TTree.h"
#include "TFile.h"
#include "TRandom.h"
#include "TRandom3.h"
#include "TNamed.h"

#include <string>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <vector>
#include <string.h>

#include "wavearray.hh"
#include "network.hh"

#include "FrameL.h"

#define LST_TREE_NAME "frl"

#ifndef TOOLBOX_HH
struct frfile {
  int start;
  int stop;
  int length;
  vector<TString> file;
};
#endif

using namespace std;

namespace CWB {

class frame : public TNamed {

public:
  
  frame();  
  frame(TString ioFile, TString chName="", Option_t* option = "", bool onDisk=false, TString label=".gwf", unsigned int mode=0);
  ~frame();  

  void open(TString ioFile, TString chName="", Option_t* option = "", bool onDisk=false, TString label=".gwf", unsigned int mode=0);

  // read in w the channel data contained in the frames defined in the filename
  // filename      : in  - name of frame file
  // channel       : in  - name of channel to be extracted
  // w             : out - data array which contains the data read
  void readFrames(char* filename, char *channel, wavearray<double> &w);

  // read in w the channel data contained in the frames defined in the frf structure
  // frf           : in  - frame file structure
  // channel       : in  - name of channel to be extracted
  // w             : out - data array which contains the data read
  void readFrames(frfile frf, char *channel, wavearray<double> &w);

  // read in w the channel data contained in the frames
  // w             : out - data array which contains the data read
  void readFrames(wavearray<double> &w);

  // write in frame the channel data contained in the wavearray x
  void writeFrame(wavearray<double> x, TString frName, TString chName);

  // dump to file ofName the list of frames defined in frf structures (list file used for the old cWB analysis)
  int dumpFrList(frfile frf, TString ofName, double sRate=16384.);

  // return the frame list of frames contained in the range [start-segEdge,stop-segEdge]
  frfile getFrList(int istart, int istop, int segEdge);

  // return the frame list of frames contained in the range [start,stop]
  vector<frfile> getFrList(int istart=0, int istop=0);

  // return the begin and end range of the frame list  
  waveSegment getFrRange() {return getFrRange(NULL);}

  // return the number of frame files read from input list
  int     getNfiles() {return nfiles;}

  // set resample index : used to resample data after read from files
  //                      resample rate = pow(2,srIndex) 
  void    setSRIndex(int srIndex) {this->srIndex=srIndex;}

  // get resample index
  int     getSRIndex() {return srIndex;}

  // set data channel name
  void    setChName(TString chName) {this->chName=chName;}

  // get data channel name
  TString getChName() {return chName;}

  // set frame name
  void    setFrName(TString frName) {this->frName=frName;}

  // get frame name
  TString getFrName() {return frName;}

  // get frame file option (READ/WRITE)
  TString getOption() {return fOption;}

  void    close();

  // if true then print more infos
  void    setVerbose(bool verbose=true) {this->verbose=verbose;}

  // frame reading retry time (sec) : 0 -> disable
  // retry time = frRetryTime*(num of trials) : max trials = 3
  void    setRetryTime(int frRetryTime=60) {this->frRetryTime=frRetryTime;}

  bool fNameCheck(TString fName);

  // set the range interval for input files
  // if xstart!=0 then only files with gps>=xstart are stored in the tree
  // if xstop!=0  then only files with gps<=xstop are stored in the tree
  void setTimeRange(int xstart=0, int xstop=0) {
       this->xstart = xstart>0 ? this->xstart=xstart : 0;
       this->xstop  = xstop>0  ? this->xstop=xstop   : 0;
       }

private:

  // return the frame list of frames contained in the range [start-segEdge,stop-segEdge]
  // if itree=NULL use the private frtree_List
  frfile getFrList(int istart, int istop, int segEdge, TTree* itree);

  // return the frame list of frames contained in the range [start-segEdge,stop-segEdge]
  // use the tree contained in rfName root file 
  frfile getFrList(TString rfName, int istart, int istop, int segEdge=0);

  // return the begin and end range of the frame list  
  // if itree=NULL use the private frtree_List
  waveSegment getFrRange(TTree* itree);

  // convert frl (frame file list) file to tree and save to rfName
  // if rfName=="" the tree is saved & sorted to the private frtree_List
  int  frl2FrTree(TString iFile, TString rfName="", TString label=".gwf", unsigned int mode=0);

  // sort iFile Tree and save the sorted tree to rfName
  int  sortFrTree(TString iFile, TString rfName);

  int sortFrTree();

  TTree* frtree_List;   //!auxiliary tree used to store frame file infos
  TString chName;	// data channel name
  TString frName;	// frame name
  TString rfName; 	// auxiliary root file name containing the tree
  int nfiles;		// number of frame files read from input list
  FrFile *frFile;	//!frame file pointer
  TString fOption;	// define the READ/WRITE mode
  int srIndex;          // resample rate = pow(2,srIndex) - 0 : disabled resampling
  bool verbose;		// if true then print more infos
  int frRetryTime;	// frame reading retry time (sec) : 0 -> disable

  double xstart;	// start time range used to fill frtree_List 
  double xstop;	        // stop  time range used to fill frtree_List

  ClassDef(frame,3)
};  

} // end namespace

// get operator 
wavearray<double>& operator >> (CWB::frame& fr, wavearray<double>& x);
// put operator
CWB::frame& operator << (CWB::frame& fr, wavearray<double>& x);
// put operator 
CWB::frame& operator >> (wavearray<double>& x, CWB::frame& fr);

#endif
