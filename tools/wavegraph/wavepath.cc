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


#include "wavepath.hh"

ClassImp(wavepath)

bool compare_paths(const wavepath& path1, const wavepath& path2) {
  return path1.get_weight() > path2.get_weight();
}

////////////////////////////////////////////////////////////////////////////////
/* BEGIN_HTML 

<p>Wavegraph is a graph-based clustering algorithm for coherent WaveBurst (cWB).
Wavegraph generates structured clusters of pixels in the time-frequency-scale
energy maps computed by cWB. The structure of the clusters is dictated by an
auxiliary user-provided directed graph. Admissible clusters are one-dimensional
paths in this graph. Wavegraph finds paths of connected pixels that maximizes a
figure-of-merit related to the cumulated sum of the time-frequency-scale energy
along the path plus a regularization term.</p>

<p>Wavegraph clustering analyses the data segment by successive blocks. The
block duration corresponds to that of the input graph. The algorithm collects
paths from the analysis of each blocks and then it sorts them and selects
distinct path that do not overlap significantly in time.</p>

<p>A graph connects nodes together. A node is associated to a
(time-frequency-scale) pixel. The node object contains only geometrical
information. For instance it includes the list of nodes to which it is connected
(its ancestors). The pixel object contains physical information (center of time,
frequency and scale bin in physical units) but it does not include any
connectivity information.</p>

<p>The wavepath class defines methods for paths in the graph. Wavepath are
ordered lists of nodes with some auxiliary information.</p>

<p>The cluster object is also an ordered list of pixels.</p>

<p>We use a fiducial scale as reference for time measurement. XXX TIME
CONVERSION FORMULA HERE! XXX</p>

END_HTML */
////////////////////////////////////////////////////////////////////////////////


void 
wavepath::init(const wavenode node, const int path_offset, const int path_refscaleidx, const double path_weight) {
  // Initialize a path and append a first node
  // 
  // Input: node         = node to append
  //        weight       = initial path weight
  //        offset       = time origin of the block (number of samples in the reference scale plane
  //                       from the beginning of the data segment)
  //        refscale_idx = index of the reference scale
  //         

  weight=path_weight;
  offset=path_offset;
  refscaleidx=path_refscaleidx;
  add_node(node);
  
}

void 
wavepath::add_node(const wavenode node) {
  // Append one node to a path
  //
  // Input: node = node to append

  path.push_back(node);

}

void
wavepath::print() {
  // Print path to standard output

  for (unsigned int rank = length(); rank-- > 0; ) {
    std::cout << "node " << nodeid(rank);
    std::cout << " (t=" << time(rank)+offset*pow(2,refscaleidx-scale(rank)) << ",f=" << freq(rank) << ",s=" << scale(rank) << ") ->\n";
  }
  std::cout << "weight=" << weight << "\n";

}

cluster
wavepath::convert_to_cluster(const int scalemin,const double deltaf_refscale, const int path_width) {
  // Convert a path object to a cluster object
  //
  // Input: scalemin        = value of the minimum scale
  //        deltaf_refscale = frequency resolution at reference scale in Hertz
  //        path_width      = width of the path in time
  //
  // The time offset of the path/cluster is computed in the reference scale plane
  //
  // Output: cluster object
  
  cluster out;
  for (unsigned int rank = length(); rank-- > 0; ) {
    int refscale_to_thisscale = pow(2,refscaleidx-scale(rank));
    double deltaf_thisscale = deltaf_refscale * refscale_to_thisscale;
    
    // init and add pixel in the center
    pixel path_center = {
      scale(rank),                                // scale index
      time(rank)+offset*refscale_to_thisscale,    // time index
      freq(rank),                                 // freq index
      0,                  // log2(scale)
      nodeid(rank),       // node ID
      0.0,  // time
      0.0}; // freq
    path_center.log2scale=scalemin+path_center.scaleix;       // log2(scale)
    path_center.time=path_center.timeix/(2*deltaf_thisscale); // time (sec)
    path_center.freq=path_center.freqix*deltaf_thisscale;     // freq (Hz)

    out.push_back(path_center);

    // add pixels on the side if width > 1
    // XXX The node ID of the side pixels is that of the path center and do not match with the graph. XXX
    // XXX This will lead to issues in the chi^2 test. XXX
    for (int n=1; n<path_width; ++n) {
      pixel path_side=path_center;
      path_side.timeix++;
      path_side.time=path_side.timeix/(2*deltaf_thisscale); 
      out.push_back(path_side);

      path_side=path_center;
      path_side.timeix--;
      path_side.time=path_side.timeix/(2*deltaf_thisscale); 
      out.push_back(path_side);
    }

  }
  return out;
}

class is_coincident_with {
public:
  is_coincident_with(double ref_time, double precision): _ref_time(ref_time),_precision(precision){}

  bool operator()(cluster const & x) const {
    return abs(_ref_time-(x.back()).time) <= _precision;
  }
private:
  double _ref_time;
  double _precision;
};

std::vector<cluster>
select_clusters_distinct_in_time(std::vector<cluster>& sorted_clusters,const double precision) {
  // Select the loudest univocal clusters that do not overlap significantly with any others
  //
  // Input: sorted_clusters  = list of clusters sorted by descending weights
  //        precision = duration of the coincidence window according to which two clusters
  //                    are said to be overlapping
  //
  // Output: list of univocal clusters
  
  std::vector<cluster> selected;

#ifdef NDEBUG
  std::cout << "debug flag: select univocal clusters\n";
  std::cout << "total number of clusters is " << sorted_clusters.size() << "\n";
#endif

  while (! sorted_clusters.empty()) {
    cluster first=sorted_clusters.front();
    selected.push_back(first);
    double ref_time=(first.back()).time;
    sorted_clusters.erase(std::remove_if(sorted_clusters.begin(),sorted_clusters.end(),is_coincident_with(ref_time,precision)),sorted_clusters.end());
#ifdef NDEBUG
    std::cout << "remaining number of clusters is " << sorted_clusters.size() << "\n";
#endif
  }
  
#ifdef NDEBUG
  std::cout << "retained " << selected.size() << " univocal clusters\n";
#endif
  
  return selected;

}
