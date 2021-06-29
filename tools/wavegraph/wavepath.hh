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
 * File name:    wavepath.hh
 * Authors:      Eric Chassande-Mottin, Eric O. Le Bigot
 **********************************************************/

#ifndef WAVEPATH_HH
#define WAVEPATH_HH

#include <math.h>
#include <vector>
#include <set>
#include <iostream>
#include <string>
#include <algorithm>

#include "TNamed.h"

typedef std::set<int> nodeids;

typedef struct {
  int nodeid;  // node id (should match node rank in the graph -- stored for debugging purpose)
  int timeix;  // time index
  int freqix;  // freq index
  int scaleix; // scale index
  double value_avg;   // node average value 
  double value_stdev; // value standard deviation
  bool is_endnode;    // flag is true if node is an end node
  nodeids ancestors;  // ids of the node ancestors
} wavenode;

typedef struct {
  int scaleix; // scale index
  int timeix;  // time index
  int freqix;  // freq index
  int log2scale;  // scale
  int nodeid;     // node ID
  double time;    // time (sec)
  double freq;    // freq (Hz)
} pixel;

typedef std::vector<wavenode> nodes;
typedef std::vector<pixel> cluster;

class wavepath : public TNamed {
public:
 
  wavepath() {};
  ~wavepath() {};
 
  int length() {return path.size();}

  int nodeid(int id) {return path[id].nodeid;}

  int time(int id) {return path[id].timeix;}

  int freq(int id) {return path[id].freqix;}

  int scale(int id) {return path[id].scaleix;}

  double get_weight() {return weight;}
  double get_weight() const {return weight;}

  int is_endnode(int id) {return path[id].is_endnode;}

  void clear(){offset=0; weight=0.0; path.clear();}

  void init(const wavenode node, const int offset, const int refscaleidx, const double path_weight);

  void add_node(const wavenode node);

  void print();

  cluster convert_to_cluster(const int scalemin, const double fs, const int path_width);

private:
  
  nodes path;      // list of nodes in the path
  int offset;      // time origin of the path in the full data segment
                   // offset corresponds to an number of pixels in the reference scale plane
  int refscaleidx; // index of the reference scale
  double weight;   // total weight carried by the pixels in the path
 
  ClassDef(wavepath,0);

};

bool compare_paths(const wavepath& path1, const wavepath& path2);
std::vector<cluster> select_clusters_distinct_in_time(std::vector<cluster>& clusters,const double precision);

#endif 
