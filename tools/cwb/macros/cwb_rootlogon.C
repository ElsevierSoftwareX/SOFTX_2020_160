// cwb rootlogon macro : define the cwb environment - must defined in .rootrc

{
  #include "Riostream.h"
  // Warning : This is used only for rootlogon.C
  //           Must be uncommented only for ROOT6 
  //           For the others macros this definition is provides by CLing dictionary 
  
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,0,0)
  #define _USE_ROOT6	
#else
  if(gSystem->Getenv("_USE_ROOT6")!=NULL) {
    cout << endl;
    cout << "Error : environment _USE_ROOT6 is defined but ROOT version is < 6.00.00 !!!" << endl;
    cout << endl;
    exit(1);
  }
#endif

  bool libSTAT=true;	// lib status 

  // ------------------------------------------------------------------------------------
  // check operative system
  // ------------------------------------------------------------------------------------
  TString OS = "";    
  if(TString(gSystem->GetBuildArch()).Contains("linux")) OS="Linux";
  if(TString(gSystem->GetBuildArch()).Contains("macos")) OS="Darwin";
  if(OS=="") {
    cout << "Error : Operative System not supported !!!	" << endl;exit(1);}
  else
    cout << endl << "OS : " << OS.Data() << endl << endl;   

  // ------------------------------------------------------------------------------------
  // check if environmental variables are defined
  // ------------------------------------------------------------------------------------
  if(gSystem->Getenv("_USE_HEALPIX")!=NULL) {
    if(gSystem->Getenv("HOME_CFITSIO")==NULL) {
      cout << "Error : environment HOME_CFITSIO is not defined !!!" << endl;
      exit(1);
    }
    if(gSystem->Getenv("HOME_HEALPIX")==NULL) {
      cout << "Error : environment HOME_HEALPIX is not defined !!!" << endl;
      exit(1);
    }
  }
  if(gSystem->Getenv("_USE_LAL")!=NULL) {
    if(gSystem->Getenv("HOME_LAL")==NULL) {
      cout << "Error : environment HOME_LAL is not defined !!!	" << endl;
      exit(1);
    }
  }
  if(gSystem->Getenv("_USE_EBBH")!=NULL) {
    if(gSystem->Getenv("HOME_CVODE")==NULL) {
      cout << "Error : environment HOME_CVODE is not defined !!!" << endl;
      exit(1);
    }
  }
  if(gSystem->Getenv("HOME_WAT_INSTALL")==NULL) {
    cout << "Error : environment HOME_WAT_INSTALL is not defined !!!" << endl;exit(1);}
  if(gSystem->Getenv("HOME_WAT")==NULL) {
    cout << "Error : environment HOME_WAT is not defined !!!	" << endl;exit(1);}
  if(gSystem->Getenv("HOME_FRLIB")==NULL) {
    cout << "Error : environment HOME_FRLIB is not defined !!!	" << endl;exit(1);}
  if(gSystem->Getenv("CWB_STFT")==NULL) {
    cout << "Error : environment CWB_STFT is not defined !!!	" << endl;exit(1);}
  if(gSystem->Getenv("CWB_GWAT")==NULL) {
    cout << "Error : environment CWB_GWAT is not defined !!!	" << endl;exit(1);}
  if(gSystem->Getenv("CWB_TOOLBOX")==NULL) {
    cout << "Error : environment CWB_TOOLBOX is not defined !!!	" << endl;exit(1);}
  if(gSystem->Getenv("CWB_HISTORY")==NULL) {
    cout << "Error : environment CWB_HISTORY is not defined !!!	" << endl;exit(1);}
  if(gSystem->Getenv("CWB_BICO")==NULL) {
    cout << "Error : environment CWB_BICO is not defined !!!	" << endl;exit(1);}
  if(gSystem->Getenv("CWB_FRAME")==NULL) {
    cout << "Error : environment CWB_FRAME is not defined !!!	" << endl;exit(1);}
  if(gSystem->Getenv("CWB_FILTER")==NULL) {
    cout << "Error : environment CWB_FILTER is not defined !!!	" << endl;exit(1);}
  if(gSystem->Getenv("HOME_CWB")==NULL) {
    cout << "Error : environment HOME_CWB is not defined !!!	" << endl;exit(1);}
  if(gSystem->Getenv("CWB_MACROS")==NULL) {
    cout << "Error : environment CWB_MACROS is not defined !!!" << endl;exit(1);}

  cout << endl;
  printf("ROOT/WAT/CWB initialization starting...\n");
  cout << endl;

  TString wat_dir     = gSystem->Getenv("HOME_WAT");  
  TString wat_install = gSystem->Getenv("HOME_WAT_INSTALL");  
  TString lib_path;

  // include paths
  // printf("Set Include Paths...\n");
  char sys_cmd[256];
  sprintf(sys_cmd,".include %s/inc",wat_install.Data());
  gROOT->ProcessLine(sys_cmd);
  sprintf(sys_cmd,".include %s/src",gSystem->Getenv("HOME_FRLIB"));
  gROOT->ProcessLine(sys_cmd);

  // load libraries
  printf("Load Aux Libraries...\n");

  // auxiliaries root libraries 
  if(gSystem->Load("libPhysics.so")<0) libSTAT=false;
  if(gSystem->Load("libFFTW.so")<0) libSTAT=false;
  if(gSystem->Load("libHtml.so")<0) libSTAT=false;
  if(gSystem->Load("libTreeViewer.so")<0) libSTAT=false;
  if(OS=="Darwin") {
     if(gSystem->Load("libFITSIO.so")<0) libSTAT=false;
  } else {
     if(gSystem->Load("libpng.so")<0) libSTAT=false;
  }
     
  // lal suite library
  if(gSystem->Getenv("_USE_LAL")) {
    TString lal_dir=TString(gSystem->Getenv("HOME_LAL"));
    if(lal_dir!="") { 
      printf("Loading LAL Suite     : %s ...\n",lal_dir.Data());
      if(gSystem->Load(lal_dir+"/lib/liblal.so")<0) libSTAT=false;
      if(gSystem->Load(lal_dir+"/lib/liblalsupport.so")<0) libSTAT=false;
      if(gSystem->Load(lal_dir+"/lib/liblalframe.so")<0) libSTAT=false;
      if(gSystem->Load(lal_dir+"/lib/liblalmetaio.so")<0) libSTAT=false;
      if(gSystem->Load(lal_dir+"/lib/liblalsimulation.so")<0) libSTAT=false;
      if(gSystem->Load(lal_dir+"/lib/liblalinspiral.so")<0) libSTAT=false;
      if(gSystem->Load(lal_dir+"/lib/liblalburst.so")<0) libSTAT=false;
    } else {
      printf("Loading LAL Suite ...\n");
      if(gSystem->Load("liblal.so")<0) libSTAT=false;
      if(gSystem->Load("liblalsupport.so")<0) libSTAT=false;
      if(gSystem->Load("liblalframe.so")<0) libSTAT=false;
      if(gSystem->Load("liblalmetaio.so")<0) libSTAT=false;
      if(gSystem->Load("liblalsimulation.so")<0) libSTAT=false;
      if(gSystem->Load("liblalinspiral.so")<0) libSTAT=false;
      if(gSystem->Load("liblalburst.so")<0) libSTAT=false;
    }
    if(gSystem->Load("libmetaio.so")<0) libSTAT=false;
  }

  // cvode library
  if(gSystem->Getenv("_USE_EBBH")) {
    TString cvode_dir=TString(gSystem->Getenv("HOME_CVODE"));
    printf("Loading cvode         : %s ...\n",cvode_dir.Data());
    if(OS=="Linux") {
      if(gSystem->Load(cvode_dir+"/lib/libsundials_cvode.so")<0) libSTAT=false;
      if(gSystem->Load(cvode_dir+"/lib/libsundials_nvecserial.so")<0) libSTAT=false;
    }
    if(OS=="Darwin") {
      if(gSystem->Load(cvode_dir+"/lib/libsundials_cvode.dylib")<0) libSTAT=false;
      if(gSystem->Load(cvode_dir+"/lib/libsundials_nvecserial.dylib")<0) libSTAT=false;
    }
  }

  TString libgomp="";
  // healpix & cfitsio libraries
  if(gSystem->Getenv("_USE_HEALPIX")) {
    // cfitsio library need to be standard
     if(OS=="Linux") {
	TString cfitsio_dir=TString(gSystem->Getenv("HOME_CFITSIO"));
	printf("Loading cfitsio       : %s ...\n",cfitsio_dir.Data());
	if(gSystem->Load(cfitsio_dir+"/libcfitsio.so")<0) libSTAT=false;
     }

    // healpix library
    TString healpix_dir=TString(gSystem->Getenv("HOME_HEALPIX"));
    printf("Loading HEALPix       : %s ...\n",healpix_dir.Data());
    if(OS=="Darwin") {
      if(gSystem->Load(healpix_dir+"/src/cxx/osx/lib/libcxxsupport.dylib")<0)  libSTAT=false;
      if(gSystem->Load(healpix_dir+"/src/cxx/osx/lib/libfftpack.dylib")<0)     libSTAT=false;
      if(gSystem->Load(healpix_dir+"/src/cxx/osx/lib/libc_utils.dylib")<0)     libSTAT=false;
      if(gSystem->Load(healpix_dir+"/src/cxx/osx/lib/libhealpix_cxx.dylib")<0) libSTAT=false;
    }
    if(OS=="Linux") {
      if(healpix_dir.Contains("2.20a")) {
        printf("\nHEALPix %s error : obsolete version (must be >= 3.00) !!!\n\n",healpix_dir.Data());
        gSystem->Exit(1);
      }

      if(gSystem->Load(healpix_dir+"/src/cxx/shared/lib/libcxxsupport.so")<0)  libSTAT=false;
      if(gSystem->Load(healpix_dir+"/src/cxx/shared/lib/libfftpack.so")<0)     libSTAT=false;
      if(gSystem->Load(healpix_dir+"/src/cxx/shared/lib/libc_utils.so")<0)     libSTAT=false;

      Long_t id,size,flags,mt; int estat; TString libname; bool blib=false;
      libname = healpix_dir+"/src/cxx/shared/lib/libpsht.so";
      estat = gSystem->GetPathInfo(libname,&id,&size,&flags,&mt);
      if(estat==0) {blib=true;if(gSystem->Load(libname)<0) libSTAT=false;}     // for HEALPix <  3.30
      libname = healpix_dir+"/src/cxx/shared/lib/libsharp.so";
      estat = gSystem->GetPathInfo(libname,&id,&size,&flags,&mt);
      if(estat==0) {blib=true;if(gSystem->Load(libname)<0) libSTAT=false;}     // for HEALPix >= 3.30
      if(blib==false) {
        printf("\nLoading HEALPix %s error !!!\n\n",healpix_dir.Data());
        cout << "check : " << healpix_dir+"/src/cxx/shared/lib/libpsht.so" << endl; 
        cout << "check : " << healpix_dir+"/src/cxx/shared/lib/libsharp.so" << endl << endl; 
        gSystem->Exit(1);
      }

      if(gSystem->Load(healpix_dir+"/src/cxx/shared/lib/libhealpix_cxx.so")<0) libSTAT=false;

      // gomp library (used in healpix : libpsht.so - SphericalHarmonic)
      if(gSystem->AccessPathName("/usr/lib/gcc/x86_64-redhat-linux/8/libgomp.so")==0)
	libgomp="/usr/lib/gcc/x86_64-redhat-linux/8/libgomp.so";
      if(gSystem->AccessPathName("/usr/lib/gcc/x86_64-redhat-linux/4.8.2/libgomp.so")==0) 
	libgomp="/usr/lib/gcc/x86_64-redhat-linux/4.8.2/libgomp.so";
      if(gSystem->AccessPathName("/usr/lib/gcc/x86_64-linux-gnu/8/libgomp.so")==0)
 	libgomp="/usr/lib/gcc/x86_64-linux-gnu/8/libgomp.so";
      if(gSystem->AccessPathName("/usr/lib/gcc/x86_64-linux-gnu/4.6/libgomp.so")==0)   
 	libgomp="/usr/lib/gcc/x86_64-linux-gnu/4.6/libgomp.so";
      if(gSystem->AccessPathName("/usr/lib/gcc/x86_64-linux-gnu/4.7/libgomp.so")==0)   
 	libgomp="/usr/lib/gcc/x86_64-linux-gnu/4.7/libgomp.so";
      if(gSystem->AccessPathName("/usr/lib/gcc/x86_64-linux-gnu/4.8/libgomp.so")==0)   
 	libgomp="/usr/lib/gcc/x86_64-linux-gnu/4.8/libgomp.so";
      if(gSystem->AccessPathName("/usr/lib/gcc/x86_64-linux-gnu/4.9/libgomp.so")==0)   
 	libgomp="/usr/lib/gcc/x86_64-linux-gnu/4.9/libgomp.so";
      if(libgomp!="") {
        if(gSystem->Load(libgomp)<0) libSTAT=false;
      } else {
        printf("Loading gomp          : library not found !!! (used by SphericalHarmonic in HEALPix)\n");
      }
    }
  }

  printf("Load cWB Libraries...\n");

  // wavelet library & macros
  lib_path = wat_install+"/lib/wavelet.so";
  printf("Loading WAT           : %s ...\n",lib_path.Data());
  if(gSystem->Load(lib_path)<0) libSTAT=false;

  printf("Loading Macros        : %s ...\n",wat_dir.Data());
  // these macros do not work with ROOT6, must be fixed !!!
  // if(OS=="Darwin") gROOT->LoadMacro(wat_dir+"/wat/macro/Plot.C");  // note: this macro produce an issue with treeviewer 
  // gROOT->LoadMacro(wat_dir+"/wat/macro/Spectrum.C");		// note: this macro could crash root (not safe)
  // gROOT->LoadMacro(wat_dir+"/wat/macro/WTSpectrum.C");	// note: this macro could crash root (not safe)
  gROOT->LoadMacro(wat_dir+"/wat/macro/Histogram.C");
  gROOT->LoadMacro(wat_dir+"/wat/macro/AddPulse.C");
  //gROOT->LoadMacro(wat_dir+"/wat/macro/AddSignals.C");
  gROOT->LoadMacro(wat_dir+"/wat/macro/readAscii.C");
  gROOT->LoadMacro(wat_dir+"/wat/macro/readtxt.C");
 
  // frame library 
  TString fr_dir=TString(gSystem->Getenv("HOME_FRLIB"));
  printf("Loading Frame         : %s ...\n",fr_dir.Data());
  if(OS=="Linux") {
     if(gSystem->Load(fr_dir+"/"+OS+"/libFrame.so")<0) libSTAT=false;
     if(gSystem->Load(fr_dir+"/"+OS+"/libFrameROOT.so")<0) libSTAT=false;
  }
  if(OS=="Darwin") {
     if(gSystem->Load(fr_dir+"/"+OS+"/libFrame.dylib")<0) libSTAT=false;
     if(gSystem->Load(fr_dir+"/"+OS+"/libFrameROOT.dylib")<0) libSTAT=false;
  }

  // eBBH library
  if(gSystem->Getenv("_USE_EBBH")) {
    lib_path = wat_install+"/lib/eBBH.so";
    printf("Loading eBBH          : %s ...\n",lib_path.Data());
    if(gSystem->Load(lib_path)<0) libSTAT=false;
  }

  // stft library
  lib_path = wat_install+"/lib/STFT.so";
  printf("Loading STFT          : %s ...\n",lib_path.Data());
  if(gSystem->Load(lib_path)<0) libSTAT=false;

  // gwat library
  lib_path = wat_install+"/lib/gwat.so";
  printf("Loading gwat          : %s ...\n",lib_path.Data());
  if(gSystem->Load(lib_path)<0) libSTAT=false;

  // toolbox library
  lib_path = wat_install+"/lib/Toolbox.so";
  printf("Loading Toolbox       : %s ...\n",lib_path.Data());
  if(gSystem->Load(lib_path)<0) libSTAT=false;

  // history library
  lib_path = wat_install+"/lib/History.so";
  printf("Loading History       : %s ...\n",lib_path.Data());
  if(gSystem->Load(lib_path)<0) libSTAT=false;

  // bicoherence library
  lib_path = wat_install+"/lib/Bicoherence.so";
  printf("Loading Bicoherence   : %s ...\n",lib_path.Data());
  if(gSystem->Load(lib_path)<0) libSTAT=false;

  // filter library
  lib_path = wat_install+"/lib/Filter.so";
  printf("Loading Filter        : %s ...\n",lib_path.Data());
  if(gSystem->Load(lib_path)<0) libSTAT=false;

  // cwb frame library
  lib_path = wat_install+"/lib/frame.so";
  printf("Loading CWB FRAME     : %s ...\n",lib_path.Data());
  if(gSystem->Load(lib_path)<0) libSTAT=false;

  // cwb library
  lib_path = wat_install+"/lib/cwb.so";
  printf("Loading cwb           : %s ...\n",lib_path.Data());
  if(gSystem->Load(lib_path)<0) libSTAT=false;

  // wavegraph library
  lib_path = wat_install+"/lib/wavegraph.so";
  printf("Loading wavegraph     : %s ...\n",lib_path.Data());
  if(gSystem->Load(lib_path)<0) libSTAT=false;

  // declare ACLiC includes environment 
  gSystem->AddIncludePath("-I$HOME_WAT_INSTALL/inc");
  gSystem->AddIncludePath("-I$ROOTSYS/include");
  gSystem->AddIncludePath("-I"+wat_install+"/inc");
  gSystem->AddIncludePath("-I"+fr_dir+"/src");
  if(gSystem->Getenv("_USE_HEALPIX")) {
    gSystem->AddIncludePath("-I$HOME_HEALPIX/src/cxx/Healpix_cxx");
    gSystem->AddIncludePath("-I$HOME_HEALPIX/src/cxx/cxxsupport");
  }
  if(gSystem->Getenv("_USE_LAL")) {
    TString lal_dir=TString(gSystem->Getenv("HOME_LAL"));
    if(lal_dir!="") { 
      gSystem->AddIncludePath("-I"+lal_dir+"/include");
    }
  }
  if(gSystem->Getenv("_USE_EBBH")) {
    TString cvode_dir=TString(gSystem->Getenv("HOME_CVODE"));
    if(cvode_dir!="") { 
      gSystem->AddIncludePath("-I"+cvode_dir+"/include");
    }
  }

  // declare ACLiCFlag options 
  TString fopts = gSystem->GetFlagsOpt();
  fopts.Append(" -D_USE_ROOT -fPIC -Wno-deprecated -mavx -Wall -Wno-unknown-pragmas");
  fopts.Append(" -fexceptions -O2 -D__STDC_CONSTANT_MACROS");
  if(gSystem->Getenv("_USE_HEALPIX")) fopts.Append(" -D_USE_HEALPIX");
  if(gSystem->Getenv("_USE_EBBH"))    fopts.Append(" -D_USE_EBBH");
  if(gSystem->Getenv("_USE_LAL"))     fopts.Append(" -D_USE_LAL");
#ifdef _USE_ROOT6
  fopts.Append(" -D_USE_ROOT6");
#endif
  if(OS=="Darwin") fopts.Append(" -fno-common -dynamiclib -undefined dynamic_lookup");  
  else             fopts.Append(" -fopenmp");
  gSystem->SetFlagsOpt(fopts.Data());

  // set the offset for TimeDisplay, the seconds declared in xaxis
  // are refered to "1980-01-06 00:00:00 UTC Sun" -> GPS = 0
  gStyle->SetTimeOffset(315964790); 

  gStyle->SetPalette(1,0);
  gStyle->SetNumberContours(256);

  gROOT->ForceStyle(0);

  if(gSystem->Getenv("HOME_WAT")) {
    printf("\ncWB library path      : %s\n",gSystem->Getenv("HOME_WAT"));
  }

  if(gSystem->Getenv("CWB_CONFIG")) {
    printf("\ncWB config path       : %s\n",gSystem->Getenv("CWB_CONFIG"));
  }

  // Print CWB logo
  PrintLogoCWB(GetLALVersion());

  // set prompt
  // Warning : prompt do not works in ROOT 6.00.02 (to be be fixed in the next version)
  ((TRint*)gROOT->GetApplication())->SetPrompt("cwb [%d] ");

  // check library
  if(!libSTAT) {
    cout << "Error Loading Libraries ..." << endl << endl;
    gSystem->Exit(1);
  }
}  

