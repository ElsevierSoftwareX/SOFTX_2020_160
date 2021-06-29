#!/usr/bin/python

import commands, os, sys
from glue.segments import *
from glue.segmentsUtils import *
import getpass

#general setup
user=getpass.getuser()
title="TITLE"
accounting_group="ligo.prod.o2.burst.allsky.cwbonline"

label="LABEL"
online_dir="/home/%s/online/%s"%(user,label)

ifos=["L1","H1"]

channelname={}                                                                        
DQ_channel={}
DQ_channel_rate={}
bitmask={}
for ifo in ifos:
    channelname[ifo]="%s:GDS-CALIB_STRAIN"%(ifo)
    DQ_channel[ifo]=["%s:GDS-CALIB_STATE_VECTOR"%(ifo),"%s:DMT-DQ_VECTOR"%(ifo)]
    DQ_channel_rate[ifo]=[16,16]
    bitmask[ifo]=[3,1] 

inj_name=["BURST","CBC","STOCH","DETCHAR"]
inj_bitmask=[128,64,32,256]

#Job setup
run_offset=0
job_offset=8
job_timeout=90*60
#min_seg_duration=30
seg_duration=180
moving_step=60
look_back=5*60
debug=0
sleep=3
max_jobs=20
#science_segment_offset=[30]

#run_start=1126256840
#run_end=1126260736

#General cWB parameters
levelR = 3
l_high   = 10

#Background information

bkg_delay=5*60
bkg_nlags=500
bkg_job_duration=600
bkg_job_minimum=600
bkg_split=1
bkg_njobs=20
bkg_framesize=4

#Directories

run_dir=online_dir+"/RUN_cWB"
jobs_dir="JOBS"
seg_dir="SEGMENTS"
summaries_dir="SUMMARIES"
zerolag_par="%s/config/user_parameters.C"%(run_dir)
bkg_par="%s/config/user_parameters_bkg.C"%(run_dir)
pp_par="%s/config/user_pparameters.C"%(run_dir)
#pe_par="%s/config/user_parameters_pe.C"%(run_dir)
#pe_plugin="%s/tools/online/cWB_Plugin_PE_on.C"%(os.environ['HOME_WAT'])
prod_plugins=["%s/tools/cwb/plugins/CWB_Plugin_QLWveto.C"%(os.environ['HOME_WAT'])]

bkg_run_dir="%s/TIME_SHIFTS"%(online_dir)
postprod_dir="POSTPRODUCTION"
considered_segments_file="considered.txt"
processed_segments_file="processed.txt"
running_segments_file="running.txt"
missing_segments_file="missing.txt"
run_segments_file="run.txt"
job_segments_file="jobs.txt"

if (os.environ['SITE_CLUSTER']=="CIT"):
   frames_dir=["/dev/shm/llhoft/L1_O2","/dev/shm/llhoft/H1_O2"]
   bkg_dir=["/ifocache/llcache/llhoft/L1_O2/L-L1_O2_llhoft-?????/L-L1_O2_llhoft-","/ifocache/llcache/llhoft/H1_O2/H-H1_O2_llhoft-?????/H-H1_O2_llhoft-"]
   log_path="/usr1/%s"%(user)
   web_dir="/home/%s/public_html/online/%s"%(user,label)
   web_link="https://ldas-jobs.ligo.caltech.edu/~%s/online/%s"%(user,label)
   accounting_group_user="marco.drago"
   condor_requirements="""Requirements = (TARGET.Online_Burst_cWB =?= True)
+Online_Burst_cWB=True"""
if (os.environ['SITE_CLUSTER']=="ATLAS"):
   frames_dir=["/dev/shm/llhoft/O2L1","/dev/shm/llhoft/O2H1"]
   bkg_dir=["/atlas/data/llcache/O2L1/L1-llhoft-?????/L1-llhoft-","/atlas/data/llcache/O2H1/H1-llhoft-?????/H1-llhoft-"]
   log_path="/local/user/%s"%(user)
   web_dir="/home/%s/WWW/LSC/online/%s"%(user,label)
   web_link="https://atlas3.atlas.aei.uni-hannover.de/~%s/LSC/online/%s"%(user,label)
   log_dir="/atlas/user/atlas5/%s/%s"%(user,label)
   accounting_group_user="marco.drago"

#Data Quality
apply_veto=False

#E-mails
emails=["marco.drago@aei.mpg.de","klimenko@phys.ufl.edu","gabriele.vedovato@lnl.infn.it","francesco.salemi@aei.mpg.de","shubhanshu.tiwari@gssi.infn.it","clazzaro@pd.infn.it","maria.tringali@unitn.it","filipe.dasilva@ufl.edu"]

#Library information
version="1"
version_wat="2G"
search="r"
optim="false"
code_version="wat6.2.6"
hoft_version="C00"

#gracedb
gracedb_group="Burst"
gracedb_analysis="CWB"
gracedb_search="AllSky"

#pp threshold
id_rho=0
#th_rho_off=8.
th_far_off=3.17e-08
th_rho_lum=6.
id_cc=0
th_cc=0.6

#TCuts
#Cuts_file=""
#Cuts_list=["",""]
#Cuts_name=["",""]

#string of cwb parameters
cwb_par=""" 
  strcpy(analysis,"%s");

  nIFO = %i;
  cfg_search = '%s';
  optim=%s;
  
  
  // frequency
  fLow  = 16.;       // low frequency of the search
  fHigh = 2048.;      // high frequency of the search

  levelR = %i;
  l_low    = 4;       // low frequency resolution level
  l_high   = %i;      // high frequency resolution level

  strcpy(wdmXTalk,"wdmXTalk/OverlapCatalog16-1024.bin");

  healpix=7;
  nSky=196608;

  bpp   = 0.001;
  subnet = 0.5;
  subcut = 0.0;
  netRHO = 5.0;
  netCC = 0.5;
  Acore = 1.7;
  Tgap  = 0.2;
  Fgap  = 128.0;
  delta = 0.5;
  cfg_gamma = -1.0;
  LOUD = 300;

  pattern = 5;

  //simulation
  nfactor = 1;
  simulation = 0;

  """%(version_wat,len(ifos),search,optim,levelR,l_high)
