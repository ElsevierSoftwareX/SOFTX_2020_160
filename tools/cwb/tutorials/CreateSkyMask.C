//
// Create SkyMask for CWB analysis
// Author : Gabriele Vedovato

//#define OFILE_NAME "L1H1V1_EarthSkyMask_2DetMask_FxThr_0d01.txt" 
//#define Fx_THR 0.01

#define OFILE_NAME "L1H1V1_EarthSkyMask_2DetMask_FxThr_0d02.txt" 
#define Fx_THR 0.02

size_t setSkyMask(double, char*, bool*, int);

gskymap* gSM;

void CreateSkyMask() {

  skymap sm(0.4,0,180,0,360);
  int L = sm.size();

  detector L1(const_cast<char*>("L1"));
  detector H1(const_cast<char*>("H1"));
  detector V1(const_cast<char*>("V1"));

  gnetwork sp_l1h1; sp_l1h1.add(&L1);sp_l1h1.add(&H1); 
  gnetwork sp_l1v1; sp_l1v1.add(&L1);sp_l1v1.add(&V1);
  gnetwork sp_v1h1; sp_v1h1.add(&V1);sp_v1h1.add(&H1);

  for (int l=0;l<L;l++) {
    double phi = sm.getPhi(l);
    double theta = sm.getTheta(l);

    double Fp_l1h1 = sp_l1h1.GetAntennaPattern(phi,theta,0,true);
    double Fx_l1h1 = sp_l1h1.GetAntennaPattern(phi,theta,0,false);

    double Fp_l1v1 = sp_l1v1.GetAntennaPattern(phi,theta,0,true);
    double Fx_l1v1 = sp_l1v1.GetAntennaPattern(phi,theta,0,false);

    double Fp_v1h1 = sp_v1h1.GetAntennaPattern(phi,theta,0,true);
    double Fx_v1h1 = sp_v1h1.GetAntennaPattern(phi,theta,0,false);

    //if(Fx_l1h1<Fx_THR || Fx_l1v1<Fx_THR || Fx_v1h1<Fx_THR) sm.set(l,0); else sm.set(l,1);
    if(Fx_l1h1<Fx_THR || Fx_l1v1<Fx_THR || Fx_v1h1<Fx_THR) sm.set(l,1); else sm.set(l,0);
  }

  int n=0;
  for (int l=0;l<L;l++) if(sm.get(l)==1) n++;
  cout << 100*(double)n/(double)L << endl;
  double REJECTED_SKY_PIXEL_PERCENTAGE = (double)n/(double)L;

  //sm+=-sm.max();
  //sm*=-1;

#ifdef OFILE_NAME
  ofstream out;
  out.open(OFILE_NAME, ios::out);
  if (!out.good()) {cout << "Error Opening File : " << OFILE_NAME << endl;exit(1);}
  for (int l=0;l<L;l++) out << l << " " << sm.get(l) << endl;
  out.close();

  bool* mask = new bool[L];
  setSkyMask((double)REJECTED_SKY_PIXEL_PERCENTAGE, (char*)OFILE_NAME, (bool*)mask, (int)L);
  for (int l=0;l<L;l++) sm.set(l,(double)mask[l]);
#endif

  gSM = new gskymap(sm);
  gSM->Draw();
  return;
}


// read skyMask
size_t setSkyMask(double f, char* file, bool* mask, int L) {
  int i;
  size_t n = 0;
  size_t l;
  char   str[1024];
  FILE* in;
  char* pc;
  double a;

  wavearray<double> skyHole(L);     // static sky mask describing "holes"
  wavearray<double> probability(L); // sky probability

  if(!L) return 0;
  if(!file) return 0;
  if(!strlen(file)) return 0;

  if( (in=fopen(file,"r"))==NULL ) return 0;

  while(fgets(str,1024,in) != NULL){

     if(str[0] == '#') continue;
     if((pc = strtok(str," \t")) == NULL) continue;
     if(pc) i = atoi(pc);                                       // sky index
     if((pc = strtok((char*)NULL," \t")) == NULL) continue;
     if(pc && i>=0 && i<int(L)) {
       skyHole.data[i] = atof(pc);                        // probability
//       nSkyStat.set(i, atof(pc));
       n++;
     }
  }
  a = skyHole.mean()*skyHole.size();
  skyHole *= a>0. ? 1./a : 0.;
  if(f==0.) { skyHole = 1.; return n; }

  double* p  = skyHole.data;
  double** pp = (double **)malloc(L*sizeof(double*));
  for(l=0; l<L; l++) pp[l] = p + l;

  probability.waveSort(pp,0,L-1);

  a = double(L);
  for(l=0; l<L; l++) {
    a -= 1.;
    *pp[l] = a/L<f ? 0. : 1.;
    //if(*pp[l] == 0.) nSkyStat.set(pp[l]-p,*pp[l]);
    if(*pp[l] == 0.) mask[pp[l]-p]=false; else mask[pp[l]-p]=true;
  }
  free(pp);
  return n;
}


