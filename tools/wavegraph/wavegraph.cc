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


#include "wavegraph.hh"

// Predicate function used to trim lines of text from double white spaces
bool bothspaces(char lhs, char rhs) { return (lhs == rhs) && (lhs == ' '); }

// Symbols for admissible keys in the header of the graph file
enum graph_header_code {
  xstride,
  xspan,
  xscalemin,
  xscalemax,
  xnscales,
  xsamp_freq
};


graph_header_code header_hash (std::string const& str) {
// Hash function for the header of the graph file
  if (str == "stride") return xstride;
  if (str == "span") return xspan;
  if (str == "scalemin") return xscalemin;
  if (str == "scalemax") return xscalemax;
  if (str == "nscales") return xnscales;
  if (str == "samp_freq") return xsamp_freq;
}

ClassImp(wavegraph)

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

<p>The wavegraph class includes methods for building and searching the graph. 
The wavegraph object is an ordered lists of nodes with some auxiliary information.</p>

END_HTML */
////////////////////////////////////////////////////////////////////////////////

void 
wavegraph::add_node(const int nodeidx, const int timeidx, const int freqidx, const int scaleidx, 
		    const double value_avg, const double value_stdev,
                    const bool endnode, const nodeids&  node_ancestors) {
  // Add a node to graph object
  //
  // Input: nodeidx        = node index (should match the node rank in the graph -- this is for debugging purpose)
  //        timeidx        = time index of the node
  //        freqidx        = frequency index of the node
  //        scaleidx       = scale index of the node
  //        value_avg      = average value of the node
  //        value_stdev    = value standard deviation of the node
  //        endnode        = flag saying if the node is an endnode or not
  //        node_ancestors = list of indices of the node ancestors

  wavenode node = {
    nodeidx,     // node id in graph
    timeidx,     // time index
    freqidx,     // freq index
    scaleidx,    // scale index
    value_avg,   // node average value
    value_stdev, // value standard deviation
    endnode,     // flag is true if node is end node
    node_ancestors  // ancestors of the node
  };

  add_node(node);

}

void 
wavegraph::add_node(const wavenode node) {
  // Add a node to graph object
  //
  // Input: node  = node to append to graph

  graph.push_back(node);

}

void 
wavegraph::create(const std::string& srcfile){
// Create the graph from the given source file
//
// Input: srcfile = file name
//
// The file format is as follows
//
// # Comments
// %% stride=XX span=XX scalemin=XX scalemax=XX nscales=XX samp_freq=XX
// 0 II TT FF SS E AA BB CC ...  << first node
// 1 II TT FF SS E AA BB CC ...  << second node
// ... etc ...
//
// II, TT, FF and SS: ID and time, freq and scale indices (integers)
// E: true if the node is an end node (boolean)
// AA, BB, CC, ...: node ancestors (there may be none)
//

#ifdef NDEBUG
  std::cout << "debug flag1 create: open file\n";
#endif

  std::ifstream source;
  source.open(srcfile.c_str(), std::ifstream::in);

  clear();

#ifdef NDEBUG
  std::cout << "debug flag2 create: read file\n";
#endif

  std::string line;    
  int count_nodes=0;

  // read line
  while ( std::getline(source,line) ){

    // ignore if empty line
    if (line.empty())
      continue;

    // ignore if commented
    if (line[0] == '#')
       continue;

    // parse header
    if ((line[0] == '%')&(line[1] == '%')) {

      // trim leading characters and multiple spaces
      line.erase(0,3);
      std::string::iterator line_end = std::unique(line.begin(), line.end(), bothspaces);
      line.erase(line_end, line.end()); 

      std::string delimiter_fields = " ";
      std::string delimiter_values = "=";
      size_t pos_field = 0;
      while (true) {

	// extract field
	pos_field = line.find(delimiter_fields);
	std::string field = line.substr(0, pos_field);
	if (field.empty())
	  continue;

	//std::cout << "line: " << line << "\n";
	//std::cout << "field: " << field << "\n";

	// read token and values from field
	size_t pos_value = field.find(delimiter_values);
	std::string token = field.substr(0, pos_value);
	field.erase(0, pos_value + delimiter_values.length());
	std::string value = field.substr(0, pos_value);

	// set variables
	switch (header_hash(token)) {
	case xstride:
	  stride=atoi(value.c_str());
	  break;
	case xspan:
	  span=atoi(value.c_str());
	  break;
	case xscalemin:
	  scalemin=atoi(value.c_str());
	  break;
	case xscalemax:
	  scalemax=atoi(value.c_str());
	  break;
	case xnscales:
	  nscales=atoi(value.c_str());
	  break;
	case xsamp_freq:
	  samp_freq=atoi(value.c_str());
	  break;
	default:
	  std::cout << "Invalid token: " << token << "\n";
	  throw std::invalid_argument("Error: Graph file has an invalid header");
	}

	// break if all done
	if (pos_field == std::string::npos)
	  break;

	// erase processed field from header line
	line.erase(0, pos_field + delimiter_fields.length());
      }

#ifdef NDEBUG
    std::cout << "header: stride= " << stride << " span=" << span;
    std::cout << " scalemin=" << scalemin <<  " scalemax=" << scalemax;
    std::cout << " nscales=" << nscales << " samp_freq=" << samp_freq << "\n";
#endif

      // go to next line
      continue;
    }

#ifdef NDEBUG
    std::cout << line << "\n";
#endif

    // parse data
    nodeids ancestors;
    std::stringstream streamline(line);

    // read nodeid
    int nodeid;
    streamline >> nodeid;

    if (nodeid != count_nodes)
       throw std::invalid_argument("Error: Invalid file");
      
    // read time index
    int time_idx;
    streamline >> time_idx;

    // read freq index
    int freq_idx;
    streamline >> freq_idx;

    // read scale index
    int scale_idx;
    streamline >> scale_idx;

    // read average graph value -- not used by wavegraph
    double value_avg;
    streamline >> value_avg;

    // read standard deviation of graph value -- not used by wavegraph
    double value_stdev;
    streamline >> value_stdev;

    // read end_node flag
    bool endnode_flag;
    streamline >> endnode_flag;

    // read node ancestors sequentially
    while (true) {
      int node;
      bool read_ok = static_cast<bool> (streamline >> node);
      if (!read_ok) 
         break;
      ancestors.insert(node);
    }

    // add node to graph
    add_node(nodeid, time_idx, freq_idx, scale_idx, value_avg, value_stdev, endnode_flag, ancestors);
    ancestors.clear();
    
    count_nodes++;
  }

#ifdef NDEBUG
  print();
#endif

  source.close();
}

void 
wavegraph::print() {
  // Print graph to standard output

  std::cout << "stride=" << stride << " span=" <<  span;
  std::cout << " scalemin=" << scalemin << " scalemax=" << scalemax;
  std::cout << " nscales=" << nscales << "samp_freq" << samp_freq << "\n";

  for (int node_id=0; node_id != graph.size(); ++node_id) {
    std::cout << node_id << " (" << nodeid(node_id) << "): ";
    
    if (get_ancestors(node_id).empty()) {
      std::cout << " [no ancestor] \n";
      continue;
    }

    nodeids::iterator anc;
    for (anc = get_ancestors(node_id).begin(); anc != get_ancestors(node_id).end(); ++anc) 
      std::cout << *anc << " "; 

    if (is_endnode(node_id))
      std::cout << "*";

    std::cout << "\n";
  }
}
  
bool 
wavegraph::is_topologically_sorted() {
// Test whether the node order defined by the node IDs is a topological order.
// The first nodes have no ancestor. The graph must be acyclic.
//
// Output: result of the test
    
  std::set<int> marked_nodes;
    
  bool result=true;
  
  // loop on nodes
  for (int nodeid=0; nodeid != graph.size(); ++nodeid) {

    // std::set<node_type>::iterator node;
    // std::cout << "looking to " << nodeid << ":\n";
    // for (node = marked_nodes.begin(); node != marked_nodes.end(); ++node) {
    // 	std::cout << *node << " ";
    // }
    // std::cout << "\n";
      
    // if node has no ancestor, mark the node and continue
    if (get_ancestors(nodeid).empty()) {
      marked_nodes.insert(nodeid);
      continue;
    }

    // loop on the node ancestors
    nodeids::iterator anc;
    for (anc = get_ancestors(nodeid).begin(); anc != get_ancestors(nodeid).end(); ++anc) {
      // if one of the ancestors is NOT in marked_nodes, the graph is NOT topologically sorted
      if ( marked_nodes.find(*anc) == marked_nodes.end() )
	return false;
    }
    // mark current node
    marked_nodes.insert(nodeid);

  }

  // the graph is topologically sorted
  return true;
};

void 
wavegraph::heaviest_path(std::vector<double>& total_weight, std::vector<int>& total_length, 
                         std::vector<int>& best_ancestor, const std::vector<double> weights) {
// Calculation of the "heaviest path" from the entry nodes (nodes with no
// ancestor) to each node of the graph, that has the largest total weight ( =
// sum of the weights of the node in the path). Uses dynamic programming.
//  
//  Input:  weights       = mapping from the node to its weight
//  Output: total_weight  = mapping from each node to the total weight of the heaviest path
//  Output: total_length  = mapping from each node to the total length of the heaviest path
//          best_ancestor = mapping from each node to the node ancestor in the heaviest path
//                          (or -1 if none)
//

  for (int nodeid=0; nodeid != graph.size(); ++nodeid) {

#ifdef NDEBUG
    std::cout << "debug flag1 heaviest_path: looking at node " << nodeid << ":";
#endif

    nodeids ancestors=get_ancestors(nodeid);
    
    // case with no ancestors
    if (ancestors.empty()) {
#ifdef NDEBUG
      std::cout << "[no ancestors]\n";
#endif
      total_weight[nodeid] = weights[nodeid];
      total_length[nodeid] = 1;
      best_ancestor[nodeid] = -1;
      continue;
    }
     
    // case with one ancestor
    nodeids::const_iterator first = ancestors.begin();      
    if (ancestors.size() == 1) {
#ifdef NDEBUG
      std::cout << "[one ancestor] select " << *first <<"\n";
#endif
      total_weight.at(nodeid) = total_weight.at(*first)+weights.at(nodeid);
      total_length.at(nodeid) = total_length.at(*first)+1;
      best_ancestor.at(nodeid) = *first;
      continue;
    }

#ifdef NDEBUG
    std::cout << "[many ancestors] find max ";
#endif

    // case with many ancestors
    double best_weight = total_weight.at(*first);
    int best_length = total_length.at(*first);
    int best_anc = *first;
    for (nodeids::const_iterator anc = ancestors.begin(); anc != ancestors.end(); ++anc) {
      
      if (total_weight.at(*anc) > best_weight) {
	best_weight = total_weight.at(*anc);
	best_length = total_length.at(*anc);
	best_anc = *anc;
      }
    }

#ifdef NDEBUG
    std::cout << "select " << best_anc << "\n";
#endif

    total_weight.at(nodeid) = best_weight+weights.at(nodeid);
    total_length.at(nodeid) = best_length+1;
    best_ancestor.at(nodeid) = best_anc;

  }
}

bool wavegraph::is_compatible(const std::vector< WSeries<double>* >& data)
{

  // assert that the graph and cWB data cube have the same scale axis
  if (log2(get_scale(data.front())) != scalemin) {

#ifdef NDEBUG
    std::cout << "data scalemin= " << log2(get_scale(data.front())) << ", graph scalemin= " << scalemin << "\n";
#endif
 
    std::cerr << "Error: min scale of graph does not match that of data\n";
    return false;
  }

  if (log2(get_scale(data.back())) != scalemax) {

#ifdef NDEBUG
    std::cout << "data scalemax= " << log2(get_scale(data.back())) << ", graph scalemax= " << scalemax << "\n";
#endif

    std::cerr << "Error: max scale of graph does not match that of data\n";
    return false;
  }

  if (data.size() != nscales) {
#ifdef NDEBUG
    std::cout << "data nscales= " << data.size() << ", graph nscales= " << nscales << "\n";
#endif

    std::cerr << "Error: scale axis of graph does not match that of data\n";
    return false;
  }

  if (nscales != scalemax-scalemin+1) {
    std::cerr << "Error: scale axis sampling is not contiguous\n";
    return false;
  }

  if (data.front()->rate() != samp_freq) {
    std::cerr << "Error: sampling frequency of the graph does not match that of data\n";
    return false;
  }
 
  return true;

}

double median(std::vector<double> & vec)
{
    size_t n = vec.size() / 2;
    nth_element(vec.begin(), vec.begin()+n, vec.end());
    return vec[n];
}

std::vector<cluster> 
wavegraph::clustering(const double threshold, const std::vector< WSeries<double>* >& data, const double strip_edges, const int path_halfwidth, const double penal_factor, const std::vector<double>& energy_thresholds){
// Compute clusters using wavegraph clustering. Main steps in pseudo-code:
//
// while data are available (avoid edges as required by strip_edges) 
//   fill the weight table of the graph with current data block
//   compute the heaviest path for this block
//   if its total weight > threshold
//     trace back (the list of nodes of) this path and store
//   endif
//   go to next block
// endwhile
// sort paths by their total weight
// convert paths into clusters
// retain univocal clusters that do not overlap in time
//
// Input:  threshold   =     selection threshold on path total weight
//         data        =     list of Wilson Daubechies Meyer transforms computed by cWB
//                           Each WSeries in the list corresponds to a given scale (or "resolution").
//                           The list should go contiguously from the smallest to the largest scale.
//         strip_edges =     remove strip_edges seconds at the segment edges to prevent contamination
//                           by edge artifacts
//         path_halfwidth  = half width of the path (path weight is computed summing weights of +/- path_halfwidth pixels)
//         penal_factor    = leading factor of length penalization term
//         energy_thresholds = vector of thresholds optionally applied to pixel energy
//                             one threshold per scale (no selection if threshold is 0)

  if (! is_compatible(data))
     throw std::invalid_argument("Error: graph and data are incompatible");

  // compute initial and end offsets to strip the edges from analysis
  int strip_scalemax = 0;
  if (strip_edges > 0)
    strip_scalemax = strip_edges*(2*data.back()->resolution());
  
#ifdef NDEBUG
  std::cout << "debug flag0 clustering: main loop\n";
  std::cout << "strip_max = " << strip_scalemax << "\n";
#endif
  
  // compute median noise level for each frequency and for each scale
  std::vector < std::vector<double> > noise_levels;
  for (int k=0; k<data.size(); k++) {
    std::vector<double> noise_level(num_of_freq_bins(data[k]));
    for (int m=0; m<num_of_freq_bins(data[k]); m++) {
      std::vector<double> non_zero_data;
      for (int n=0; n<num_of_time_bins(data[k]); n++) {
	double value = get_map00(data[k],n,m);
	if (value > 0)
	  non_zero_data.push_back(value);
      }
      if (non_zero_data.empty())
	noise_level.push_back(0.0);
      else
	noise_level.push_back(median(non_zero_data));
    }
    noise_levels.push_back(noise_level);
  }

#ifdef NDEBUG
  std::cout << "noise level estimates:\n";
  for (int k=0; k<noise_levels.size(); k++) {
    std::cout << "scale=" << k << ": ";
    for (int m=0; m<noise_levels[k].size(); m++) {
      std::cout << noise_levels[k].at(m) << " ";
    }
    std::cout << "\n";
  }
#endif 

  //
  // main loop: this loop runs by blocks of pixels in the max scale plane. 
  // The max scale is the last element in the data vector [hence data.back()]
  //
  std::vector<wavepath> selected_paths;

  int end_scalemax = num_of_time_bins(data.back())-get_span()-strip_scalemax;
  for (int block_offset_scalemax = strip_scalemax; block_offset_scalemax <= end_scalemax; block_offset_scalemax += get_stride()) {

#ifdef NDEBUG
    std::cout << "debug flag1 clustering: feed weight table\n";
#endif
    
    // fill the weight table of the graph with current data block
    std::vector<double> weights(size(),0.0);
    for (int nodeid=0; nodeid != size(); ++nodeid) {
      int node_offset = block_offset_scalemax * pow(2,scalemax-(scalemin+scale(nodeid)));
      int time_index = node_offset + time(nodeid);
#ifdef NDEBUG      
      //std::cout << "debug flag2 clustering: accessing pixel t=" << node_offset+time(nodeid);
      //std::cout << " f=" << freq(nodeid) << " scale=" << scale(nodeid) << " with block_offset=" << block_offset_scalemax << "\n";
      if ((node_offset+time(nodeid) >= num_of_time_bins(data.at(scale(nodeid)))) | (freq(nodeid) >= num_of_freq_bins(data.at(scale(nodeid))))) {
	std::cout << "debug flag2 clustering: accessing pixel t=" << node_offset+time(nodeid);
	std::cout << " f=" << freq(nodeid) << " scale=" << scale(nodeid) << " with block_offset=" << block_offset_scalemax << "/" << end_scalemax << "\n";
	std::cout << "at scale index= " << scale(nodeid);
	std::cout << ": num of time bins= " << num_of_time_bins(data.at(scale(nodeid)));
	std::cout << ", num of freq bins= " << num_of_freq_bins(data.at(scale(nodeid))) << "\n";
	throw std::string("Error: requested position exceed allocated limits");
      }	
#endif
      double pixel_data = get_map00(data.at(scale(nodeid)),time_index,freq(nodeid));
      weights[nodeid]= pixel_data > energy_thresholds.at(scale(nodeid)) ? pixel_data : 0.0 ; // apply energy threshold
      // if non zero path_width, sum over 
      for (int n=1; n<=path_halfwidth; ++n)  {
	pixel_data = get_map00(data.at(scale(nodeid)),time_index-n,freq(nodeid));
	weights[nodeid]+= pixel_data > energy_thresholds.at(scale(nodeid)) ? pixel_data : 0.0 ; // apply energy threshold
	pixel_data = get_map00(data.at(scale(nodeid)),time_index+n,freq(nodeid));
	weights[nodeid]+= pixel_data > energy_thresholds.at(scale(nodeid)) ? pixel_data : 0.0 ; // apply energy threshold
      }
      weights[nodeid]-=penal_factor*(2*path_halfwidth+1)*noise_levels.at(scale(nodeid)).at(freq(nodeid));
    }
    
#ifdef NDEBUG
    std::cout << "debug flag2 clustering: dynamic programming\n";
#endif

    // compute the heaviest path for this block
    std::vector<double> total_weights(size());
    std::vector<int> total_lengths(size());;
    std::vector<int> best_ancestors(size());;
    heaviest_path(total_weights,total_lengths,best_ancestors,weights);
  
#ifdef NDEBUG
    std::cout << "debug flag3 clustering: select and reconstruct best paths\n";
#endif

    // select path in block with largest total_weight
    int best_node_id=-1;
    double best_weight=0;
    for (int node_id=0; node_id != size(); ++node_id) {

      //   if its total weight > path_length * threshold
      if (is_endnode(node_id) & (total_weights[node_id] >= total_lengths[node_id]*(2*path_halfwidth+1)*threshold)) {
	if (total_weights[node_id] > best_weight) {
	  best_weight = total_weights[node_id];
	  best_node_id = node_id;
	}
      }
    }

    // trace back (the list of nodes of) this path and store
    if (best_node_id>0) {
      wavepath path;
      path.clear();
      
      int node=best_node_id;
      // paths use the maximum scale as the reference scale
      path.init(get_node(node),block_offset_scalemax,nscales-1,total_weights[node]);
      
#ifdef NDEBUG
      std::cout << "path: " << node;
#endif
      while (true) {
	node = best_ancestors.at(node);
#ifdef NDEBUG
	std::cout << "->" << node;
#endif
	if (node == -1)
	  break;
	path.add_node(get_node(node));
      }
#ifdef NDEBUG
      std::cout << "\n";
      path.print();
#endif 
 
      // store the best path (largest total weight)
      selected_paths.push_back(path);  
    }

#ifdef NDEBUG
    std::cout << selected_paths.size() << " paths selected\n";
#endif 

    // go to next block
    //block_offset_scalemax+=get_stride();
  }
  
#ifdef NDEBUG
  std::cout << "debug flag4 clustering: sort selected paths\n";
#endif
  
  // sort paths by their total weight
  std::sort(selected_paths.begin(), selected_paths.end(), compare_paths);
  
#ifdef NDEBUG
  std::cout << selected_paths.size() << " selected paths\n";
  // wavepath best_path=selected_paths.front();
  // std::cout << "best path is:\n";
  // best_path.print();
#endif

  // convert paths to clusters
  std::vector<cluster> clusters;
  std::vector<wavepath>::const_iterator this_path;
  for (this_path = selected_paths.begin(); this_path != selected_paths.end(); ++this_path)
    clusters.push_back(((wavepath) *this_path).convert_to_cluster(log2(get_scale(data.front())),data.back()->resolution(),path_halfwidth));

#ifdef NDEBUG
  std::cout << clusters.size() << " selected clusters\n";
#endif

  // retain univocal clusters that are distinct in time
  std::vector<cluster> univocal_clusters;
  univocal_clusters=select_clusters_distinct_in_time(clusters,get_span()/(2*data.back()->resolution()));

#ifdef NDEBUG
  std::cout << univocal_clusters.size() << " univocal clusters\n";
  if ( !univocal_clusters.empty() ) {
    cluster best_cluster=univocal_clusters.front();
    std::cout << "best cluster is:\n";
    cluster::const_iterator this_pixel;
    for (this_pixel=best_cluster.begin(); this_pixel != best_cluster.end(); ++this_pixel) {
      std::cout << "timeix=" << ((pixel) *this_pixel).timeix <<  " freqix=" << ((pixel) *this_pixel).freqix  << " scaleix=" << ((pixel) *this_pixel).scaleix << " : ";
      std::cout << " time=" << ((pixel) *this_pixel).time <<  " freq=" << ((pixel) *this_pixel).freq  << " log_2 scale=" << ((pixel) *this_pixel).log2scale << "\n";
    }
  }
#endif
  
  return univocal_clusters;

}
