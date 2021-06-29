#ifdef __CINT__ 

#pragma link off all globals; 
#pragma link off all classes; 
#pragma link off all functions; 

//#pragma link C++ namespace CWB;

#pragma link C++ function operator>>(CWB::frame&,wavearray<double>&);
#pragma link C++ function operator>>(wavearray<double>&,CWB::frame&);
#pragma link C++ function operator<<(CWB::frame&,wavearray<double>&);

#pragma link C++ class CWB::frame+;
//#pragma link C++ struct frfile+;
#pragma link C++ class vector<frfile>+;

#endif // __CINT__
