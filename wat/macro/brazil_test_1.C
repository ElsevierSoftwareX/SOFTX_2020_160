{
watplot C;                                        // create default canvas
wavearray<double> ts(16384);   // create TS 16384 samples long
ts.rate(1024);                                 // set TS rate at 1024 Hz
AddGaus(ts,1);                               // add Gauss noise with unity var.
C.plot(ts,"APL");                            // plot TS
Biorthogonal<double> B(32);     // create bi-orthogonal wavelet 32 coeff.
WSeries<double> tf;                  // create TF object
tf.Forward(ts,B,1);                      // perform 1 step wavelet transform
wavearray<double> y;                 // create temporary dummy TS
tf.getLayer(y,0);                          // extract down-sampled TS into y
C.plot(ts,"same",2);                     // plot down-sampled TS
tf.getLayer(y,1); y=0;                 // extract HF data and zero 
tf.putLayer(y,1);                         // zero HF wavelet data
tf.Inverse();                                 // low-pass TS
C.plot(ts, "APL",1,2.,ts.stop()-2.,true,0.,0.,true,1.);       // plot x PSD
C.plot(tf, "same",2,2.,ts.stop()-2.,true,0.,0.,true,1.);  // plot down-sampled PSD
WDM<double> wdm(32, 64, 4, 8);  // define a WDM transform (32 bands) 
tf.Forward(ts,wdm);                            // perform WDM transform
C.plot(tf,0,2,ts.stop()-2,"colz");      // plot WDM series
CWB::mdc::AddSGBurst(ts,5, 200, 0.01,0);
tf.Forward(ts,wdm);                            // perform WDM transform
C.plot(tf,0,2,ts.stop()-2,"colz");      // plot WDM series
 double th=2;                            // threshold at 2 sigma
 int count0 = 0;
 double* p00 = tf.data;
 double* p90 = tf.data+tf.size()/2;
 for(int i=0; i<tf.size()/2; i++) {
    double a = p00[i]*p00[i]+p90[i]*p90[i]<th*2 ? 0. : 1.;
    p00[i] *= a; p90[i] *= a; count0 += int(a);
 }
C.plot(tf,0,2,ts.stop()-2,"colz");      // plot WDM series
}
