#ifdef __CINT__ 

#pragma link off all globals; 
#pragma link off all classes; 
#pragma link off all functions; 

#pragma link C++ enum GWAT_DRAW;

//#pragma link C++ class gwavearray<Long64_t>+;
#pragma link C++ class gwavearray<int>+;
#pragma link C++ class gwavearray<unsigned int>+;
#pragma link C++ class gwavearray<long long>+;
#pragma link C++ class gwavearray<short>+;
#pragma link C++ class gwavearray<long>+;
#pragma link C++ class gwavearray<float>+;
#pragma link C++ class gwavearray<double>+;

//#pragma link C++ class gWSeries<int>+;
//#pragma link C++ class gWSeries<unsigned int>+;
//#pragma link C++ class gWSeries<long long>+;
//#pragma link C++ class gWSeries<short>+;
//#pragma link C++ class gWSeries<long>+;
//#pragma link C++ class gWSeries<float>+;
#pragma link C++ class gWSeries<double>+;

#pragma link C++ class gskymap+;
#pragma link C++ class gnetwork+;

#endif // __CINT__
