#include "util.h"
#include <vector>
#include <map>
#include <set>
using namespace std;

namespace jtf{
  namespace util{

    static size_t mapping(const vector<size_t> * ptr, size_t index){
      return (ptr?(*ptr)[index]:index);
    }

    void extract_chain_from_edges(
        const vector<pair<size_t,size_t> > & possible_edges,
        vector<deque<pair<size_t,size_t> > > & chain,
        const vector<size_t> * point_mapping)
    {
      chain.clear();
      map<size_t,vector<size_t> > points_with_edges;
      for(size_t t = 0; t < possible_edges.size(); ++t){
          const pair<size_t,size_t> & edge = possible_edges[t];
          points_with_edges[mapping(point_mapping,edge.first)].push_back(t);
          points_with_edges[mapping(point_mapping,edge.second)].push_back(t);
        }

      set<pair<size_t,size_t> > edge_set;
      for(size_t ei = 0; ei < possible_edges.size(); ++ei){
          const pair<size_t,size_t> & one_edge = possible_edges[ei];
          if(one_edge.first < one_edge.second)
            edge_set.insert(one_edge);
          else
            edge_set.insert(make_pair(one_edge.second, one_edge.first));
        }
      while(!edge_set.empty()){
          deque<pair<size_t,size_t> > one_chain;
          one_chain.push_back(*edge_set.begin());
          edge_set.erase(edge_set.begin());
          while(!edge_set.empty() &&
                points_with_edges[mapping(point_mapping,one_chain.back().second)].size() == 2) // regular point needs to link another edge
            {
              const size_t edge_set_num_prev = edge_set.size();
              for(set<pair<size_t,size_t> >::iterator sit = edge_set.begin();
                  sit != edge_set.end(); ++sit){
                  if(mapping(point_mapping,sit->first) ==
                     mapping(point_mapping,one_chain.back().second)) {
                      one_chain.push_back(*sit);
                      edge_set.erase(sit);
                      break;
                    }
                  if(mapping(point_mapping,sit->second) ==
                     mapping(point_mapping,one_chain.back().second)){
                      one_chain.push_back(make_pair(sit->second, sit->first));
                      edge_set.erase(sit);
                      break;
                    }
                }
              if(edge_set.size() == edge_set_num_prev) break;
            }
          while(!edge_set.empty() &&
                points_with_edges[mapping(point_mapping,one_chain.front().first)].size() == 2)
            {
              const size_t edge_set_num_prev = edge_set.size();
              for(set<pair<size_t,size_t> >::iterator sit = edge_set.begin();
                  sit != edge_set.end(); ++sit){
                  if(mapping(point_mapping,sit->first) ==
                     mapping(point_mapping,one_chain.front().first)){
                      one_chain.push_front(make_pair(sit->second, sit->first));
                      edge_set.erase(sit);
                      break;
                    }
                  if(mapping(point_mapping,sit->second) ==
                     mapping(point_mapping,one_chain.front().first)){
                      one_chain.push_front(*sit);
                      edge_set.erase(sit);
                      break;
                    }
                }
              if(edge_set.size() == edge_set_num_prev) break;
            }
          chain.push_back(one_chain);
        }
    }
  }
}
