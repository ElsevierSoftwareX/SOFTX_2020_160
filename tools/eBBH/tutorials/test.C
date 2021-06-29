{  gSystem->Load("/Users/valentin/sundials/lib/libsundials_cvode.dylib");
   gSystem->Load("/Users/valentin/sundials/lib/libsundials_nvecserial.dylib");
   
   gSystem->Load("libMatrix.so");
   gSystem->Load("/Users/valentin/LIGO/trunk/wat/lib/wavelet.so");
   
   gSystem->Load("lib/eBBH.so");
   
   wavearray<double> hp, hx;
   getEBBH(16.395803, 19.401641, 25.317980, 0.24891893, hp, hx);
}
