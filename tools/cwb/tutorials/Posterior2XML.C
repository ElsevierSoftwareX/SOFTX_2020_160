// --------------------------------------------------------------------------------------------------------------------------
//
// sampleFile      : input posterior sample file produced by PE
// xmlFile         : output xml file
// options
// --gps_start_time  : gps start time 
//                     if(gps_start_time<=0) the posterior time is used
// -- gps_stop_time   : gps stop time 
// -- time_step       : injection time step
// -- seed            : if (<0) the parameters are read sequentially from posterior sample 
//                    : if (>0) the parameters are read randomly (uniformely with seed) from posterior sample 
// --waveform         : if declared it overwrite the approximant declared in the posteriors
// --source           : if declared it is added to the SimInspiralTable->source (max length=20) 
// --ninjections      : number of xml entries (only for gps_stop_time<=0)
// --clb_file         : if defined then it is used to dump the calibration parameters (same format as for posteriors sample)
//
// --------------------------------------------------------------------------------------------------------------------------

{
  // Example: create xml & clb files from posterior samples file for chunk 19

  TString options = "";
  options += "--gps_start_time 1185217218 ";	// start gps time chunk 19
  options += "--gps_stop_time 1185937218 ";	// stop  gps time chunk 19
  options += "--time_step 150 ";		// time step between injections
  options += "--seed 1 ";			// seed>=0 -> injections are selected randomly from the posterior samples 
  options += "--clb_file posterior_samples_K19_TS150s.clb ";	// name of output calibration file (optional)
  options += "--source GWXXYYZZ ";		// set source name 

//  options += "--waveform IMRPhenomPv2 ";	// force approximant
//  options += "--ninjections 1000 ";		// number of output injection (only if --gps_start_time<=0)

  TString posteriorFile = "posterior_samples.dat";		// input posterior samples file
  TString xmlFile       = "posterior_samples_K19_TS150s.xml";	// output xml file

  CWB::mdc::Posterior2XML(posteriorFile, xmlFile, options);

  exit(0);
}
