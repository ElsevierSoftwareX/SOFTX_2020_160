#ifdef __CINT__ 

#pragma link off all globals; 
#pragma link off all classes; 
#pragma link off all functions; 

#pragma link C++ global FACTORS_MAX;
#pragma link C++ global DQF_MAX;

#pragma link C++ enum MDC_TYPE;
#pragma link C++ enum MDC_COORDINATES;
#pragma link C++ enum MDC_DISTRIBUTION;
#pragma link C++ enum MDC_DRAW;

#pragma link C++ struct waveform+;
#pragma link C++ struct mdcpar+;
#pragma link C++ struct source+;
#pragma link C++ struct mdcid+;

#pragma link C++ class CWB::mdc+;
#pragma link C++ class vector<waveform>+;
#pragma link C++ class vector<mdcpar>+;
#pragma link C++ class vector<source>+;

#pragma link C++ enum CWB_STAGE;
#pragma link C++ enum CWB_PLUGIN;
#pragma link C++ enum CWB_JOBF_OPTIONS;
#pragma link C++ enum CWB_OUTF_OPTIONS;

#pragma link C++ class cwb+;
#pragma link C++ class cwb1G+;
#pragma link C++ class cwb2G+;
#pragma link C++ class CWB::config-;
#pragma link C++ class CWB::ced;
#pragma link C++ class livetime;
#pragma link C++ class variability;
#pragma link C++ class wavenoise;


#endif // __CINT__
