#define XIFO 4

#pragma GCC system_header

#include <vector>
#include <algorithm>
#include <stdexcept>
#include <iterator>

void
wavegraph_demo() {

  // load wavegraph library
  if(!gROOT->GetClass("wavegraph")) {
    TString wat_path = gSystem->Getenv("HOME_WAT");
    TString wavegraph_path = wat_path+"/tools/wavegraph/lib/wavegraph.so";
    printf("Loading wavegraph          : %s ...\n",wavegraph_path.Data());
    if(gSystem->Load(wavegraph_path)) {
      cout << "error loading wavegraph library" << endl;
      gSystem->Exit(1);
    }
  }

  // load data
  TFile inputfile("wavegraph_demo_data.root");
  if(!inputfile.IsOpen()) {
    cout << "File .root not exist" << endl; 
    gSystem->Exit(1);
  }
  inputfile.ls();

  // WSeries<double>* test=(WSeries<double>*)inputfile.Get("V1:7");
  // cout << "time bins="<<num_of_time_bins(test)<<" freq bins="<<num_of_freq_bins(test)<<"\n";
  // double mytest=get_map00(test,3231,13);                                                                                                                         
  // cout << mytest << "\n";
  // exit(0);

  // store data in STL vectors
  const int nscales=6;
  Int_t my_scales[nscales]={3,4,5,6,7,8};

  std::vector<int> scales;
  for (int n=0; n<nscales; n++)
    scales.push_back(my_scales[n]);

  std::vector< WSeries<double>* > data;
  ostringstream convert_str;
  for (int n=0; n != nscales; ++n) {
    convert_str << "V1:" << my_scales[n];
    std::string label=convert_str.str();
    cout << "" << label << "\n";
    data.push_back((WSeries<double>*)inputfile.Get(label.c_str()));
    cout <<  "pointer to data[" << n << "] is " << data[n] << endl;
    label.clear();
    convert_str.clear();
    convert_str.str("");
  }  

  // load graph
  wavegraph graph;
  graph.create("wavegraph_demo.txt");

  if ( !graph.is_topologically_sorted() ) {
    std::cerr << "the graph is not topologically sorted\n"; 
    exit(0);
  }

  // apply clustering
  std::vector<cluster> clusters = graph.clustering(400.0,data,8.0);
  std::cout << "found " << clusters.size() << " clusters\n";
 
  // make plots
  watplot WTS(const_cast<char*>("WTS"));
  WSeries<double>* WS = data.at(4); 

  int layers = WS->maxLayer()+1;  // numbers of frequency bins (first & last bins have df/2)
  int slices = WS->sizeZero();    // number of time bins
  float df = WS->resolution();    // frequency bin resolution (hz)
  float dt = 1./(2*df);           // time bin resolution (sec)

  // draw TF maps

  double start = WS->start();
  double stop  = WS->start()+slices*dt;
  double flow  = 0;
  double fhigh = (layers-1)*df;
  cout.precision(14);
  cout << "start " << start << " stop " << stop << " flow " << flow << " fhigh " << fhigh << endl;
  WTS.plot(*WS, 2, start, stop, const_cast<char*>("COLZ"));

  // set frequency range
  WTS.hist2D->GetYaxis()->SetRangeUser(flow, fhigh);

  // draw clusters
  int colors[8]={kWhite, kYellow, kOrange, kMagenta, kGreen, kBlue, kRed, kCyan};
  std::vector<TPolyLine> clusters_lines (clusters.size());
  std::vector<TPolyMarker> clusters_points (clusters.size());

  for (int n=0; n<clusters.size(); ++n) {
    
    int N = clusters[n].size();
    wavearray<double>* t= new wavearray<double>(N);
    wavearray<double>* f= new wavearray<double>(N);
    wavearray<double>* scale= new wavearray<double>(N);
    
    for (int k = N; k-- > 0;) {
      (*t)[k]=(clusters[n])[k].time;
      (*f)[k]=(clusters[n])[k].freq;
      (*scale)[n]=(clusters[n])[k].log2scale;
      //cout << "t=" << t[n] << " sec, f=" << f[n] << " Hz, scale=" << scale[n] << endl;
    }
    clusters_lines[n].SetPolyLine(N,t->data,f->data);    
    
    clusters_lines[n].SetLineColor(colors[n%8]);
    clusters_lines[n].SetLineStyle(2);
    clusters_lines[n].SetLineWidth(4);  
    clusters_lines[n].Draw("same");
    
    clusters_points[n].SetPolyMarker(N,t->data,f->data);
    clusters_points[n].SetMarkerSize(8);
    clusters_points[n].SetMarkerStyle(kDot);
    clusters_points[n].SetMarkerColor(colors[n%8]);
    clusters_points[n].Draw("same");
    
    //delete t;
    //delete f;
    //delete scale;

  }

  // write plots and data to file
  char fname[1024];
  sprintf(fname,"wavegraph_demo_plots.root");
  WTS.canvas->cd();
  WTS.canvas->Print(fname);

}
