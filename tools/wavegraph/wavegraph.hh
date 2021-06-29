/*
# Copyright (C) 2019 Eric Chassande-Mottin, Philippe Bacon, Gayathri V, Archana Pai
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/


/**********************************************************
 * Package:      wavegraph Class Library
 * File name:    wavegraph.hh
 * Authors:      Eric Chassande-Mottin, Eric O. Le Bigot
 **********************************************************/

#ifndef WAVEGRAPH_HH
#define WAVEGRAPH_HH

#include <math.h>
#include <vector>
#include <set>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <algorithm>
#include <stdexcept>
#include <limits>

#include "wseries.hh"
#include "wavepath.hh"

// alias to functions from cWB WAT

// scale exponent = log_2(scale)
inline int get_scale(WSeries<double> *WS) {return WS->maxLayer();}

// number of freq bins (first & last bins have df/2)
inline int num_of_freq_bins(WSeries<double> *WS) {return WS->maxLayer()+1;}

// number of time bins
inline int num_of_time_bins(WSeries<double> *WS) {return WS->sizeZero();}

// get time (sec) and frequency (Hz) resolution
inline float freq_res(WSeries<double> *WS) {return (float)WS->rate()/WS->getLevel()/2;}
inline float time_res(WSeries<double> *WS) {return WS->getLevel()/(float)WS->rate();}

// get map00/90 value from index
inline float get_map00(WSeries<double> *WS, int index) 
{
#ifdef NDEBUG
  if ((index < 0) | (index > WS->maxIndex())) {
    std::cout << "debug clustering: requested index " << index << " is <0 or >" << WS->maxIndex() << "\n";
    throw std::string("Error: index is out of range");
  }
#endif
  return float(WS->pWavelet->pWWS[index]);
}
inline float get_map90(WSeries<double> *WS, int index) {return float(WS->pWavelet->pWWS[index+WS->maxIndex()+1]);}

// get map00/90 value at given time and frequency
inline float get_map00(WSeries<double> *WS, int time, int freq) {return get_map00(WS, time*(WS->maxLayer()+1)+freq);}
inline float get_map90(WSeries<double> *WS, int time, int freq) {return get_map90(WS, time*(WS->maxLayer()+1)+freq);}

class wavegraph : public TNamed {
  
public:
 
  wavegraph() {};
  ~wavegraph() {};
 
  void create(const std::string& srcfile);

  void clear(){stride=0; span=0; scalemin=0; scalemax=0; nscales=0; samp_freq=0; graph.clear();}

  void print();

  void add_node(const int nodeidx, const int timeidx, const int freqidx, const int scaleidx, 
		const double value_avg, const double value_stdev,
                const bool endnode, const nodeids&  ancestors);

  void add_node(const wavenode node);

  int size() {return graph.size();}

  wavenode get_node(int id) {return graph.at(id);}

  nodeids get_ancestors(int id){return graph.at(id).ancestors;} 

  int nodeid(int id) {return graph[id].nodeid;}

  int time(int id) {return graph[id].timeix;}

  int freq(int id) {return graph[id].freqix;}

  int scale(int id) {return graph[id].scaleix;}

  int is_endnode(int id) {return graph[id].is_endnode;}

  int get_span(){return span;}

  int get_stride(){return stride;}

  bool is_topologically_sorted();

  void heaviest_path(std::vector<double>& total_weight, std::vector<int>& total_length, 
                     std::vector<int>& best_ancestor, const std::vector<double> weights);

  bool is_compatible(const std::vector< WSeries<double>* >& data);

  std::vector<cluster> clustering(const double threshold, const std::vector< WSeries<double>* >& data, const double strip_edges, const int path_halfwidth, const double penal_factor, const std::vector<double>& energy_thresholds);
 
private:
  
  nodes graph;    // graph object
  int span;       // width of the graph (measured by a number of samples at the largest scale)
  int stride;     // number of samples at largest scale between a block and the next one
  int scalemin;   // smallest scale used to compute the graph
  int scalemax;   // largest scale used to compute the graph
  int nscales;    // number of scales used to compute the graph
  int samp_freq;  // sampling frequency in Hertz used to compute the graph

  ClassDef(wavegraph,0)

};
 
#endif 
