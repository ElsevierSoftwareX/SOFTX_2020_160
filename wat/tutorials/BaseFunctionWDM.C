{  
   // this macro extract the WDM Base Wave Function and compute time. frequency RMS 
   // rms are used to compute the time,frequency pixel erros for chirp mass in netcluster::mchirp

   #define nLAYERS 	32	// layers used in TF WDM transform
   #define RATE		512
   #define DURATION	2	// sec

   #define TIME_PIXEL_INDEX	100
   #define FREQ_PIXEL_INDEX	20

   // define wavearray time series
   int Rate = RATE;           // sampling rate (Hz)
   int Duration = DURATION;   // duration (seconds)
   int N = Rate*Duration;     // number of samples
   wavearray<double> ts(N);   // time series container
   ts.rate(Rate);	      // set rate
   ts=0;		      // set array to zero 

   // produce the TF map:
   WDM<double> wdm(nLAYERS, nLAYERS, 6, 10);   	// define a WDM transform (32 bands)
   WSeries<double> tf;           		// TF map container

   cout << endl;
   cout << "ts size = " << ts.size() << " ts rate = " << ts.rate() << endl;
   tf.Forward(ts, wdm);          		// apply the WDM to the time series 
   int levels = tf.getLevel();
   //cout << "levels = " << levels << endl;  
   cout << "tf size = " << tf.size() << endl;

   double dF = tf.resolution();             	// frequency bin resolution (hz)
   double dT = 1./(2*dF);                   	// time bin resolution (sec)

   cout << "rate(hz) : " << RATE << "\t layers : " << nLAYERS
        << "\t dF(hz) : " << dF << "\t dT(ms) : " << dT*1000. << endl;

   // compute tfmap index 
   int itime = TIME_PIXEL_INDEX;
   int ifreq = FREQ_PIXEL_INDEX;
   int index = (levels+1)*itime+ifreq;

   // compute pixel time and frequency resolution
   double time = itime*dT;
   double freq = (ifreq>0) ? ifreq*dF : dF/4;
   cout << endl;
   cout << "PIXEL TIME = " << time << " sec " << endl;
   cout << "PIXEL FREQ = " << freq << " Hz  " << endl;
   cout << endl;

   // get WDM Base Function
   wavearray<double> x;
   int j00 = wdm->getBaseWave(index,x,false);
   x.resize(N);
   x.rate(RATE);

   // plot WDM base function
   gwavearray<double> gx(&x);
   gx.Draw(GWAT_TIME);

   cout << endl;

   // ------------------------------------------------------
   // compute the mean and rms time
   // ------------------------------------------------------

   double dt = 1./x.rate();
   //cout << "sample time - dt (ms) = " << dt*1000 << endl;
   double tee=0.;
   double tavr=0.;
   for(int i=0;i<x.size();i++) {
     double t = i*dt;
     tavr+=x[i]*x[i]*t;
     tee+=x[i]*x[i];
   }
   tavr/=tee;
   cout << "mean time      : " << tavr << endl;

   double trms=0.;
   for(int i=0;i<x.size();i++) {
     double t = i*dt;
     trms+=x[i]*x[i]*pow(t-tavr,2);
   }
   trms/=tee;
   trms=sqrt(trms);
   cout << "rms  time      : " << trms*1000 << " (ms) " << endl;

   cout << endl;

   // ------------------------------------------------------
   // compute the mean and rms frequency
   // ------------------------------------------------------

   x.FFT(1);

   double df = (double)RATE/x.size();
   //cout << "frequency resolution : df (Hz) = " << df << endl;
   double fee=0.;
   double favr=0.;
   for(int i=0;i<x.size()/2;i++) {
     double f = i*df;
     double e = x[2*i]*x[2*i]+x[2*i+1]*x[2*i+1];
     favr+=e*f;
     fee+=e;
   }
   favr/=fee;
   cout << "mean frequency : " << favr << endl;

   double frms=0.;
   for(int i=0;i<x.size()/2;i++) {
     double f = i*df;
     double e = x[2*i]*x[2*i]+x[2*i+1]*x[2*i+1];
     frms+=e*pow(f-favr,2);
   }
   frms/=fee;
   frms=sqrt(frms);
   cout << "rms  frequency : " << frms << " (Hz) " << endl;
  
}
