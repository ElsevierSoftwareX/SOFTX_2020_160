#ifdef __CINT__ 

#pragma link off all globals; 
#pragma link off all classes; 
#pragma link off all functions; 

#pragma link C++ class CWB::Filter;

#ifdef _USE_WAT
#pragma link C++ function operator>>(CWB::Filter&,wavearray<double>&);
#pragma link C++ function operator>>(wavearray<double>&,CWB::Filter&);
#endif

#endif // __CINT__
