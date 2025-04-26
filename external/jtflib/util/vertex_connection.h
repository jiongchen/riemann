#ifndef COMMON_VERTEX_CONNECTION_H
#define COMMON_VERTEX_CONNECTION_H

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/kruskal_min_spanning_tree.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/static_assert.hpp>
#include <iostream>

enum DIRECTION{UNDIRECT,DIRECT};

template<DIRECTION TRAITS_TYPE>
class direction_traits;

template<>
class direction_traits<DIRECT>{
public:
  typedef boost::adjacency_list < boost::vecS, boost::vecS, boost::directedS,
  boost::no_property, boost::property < boost::edge_weight_t, double > > Graph;
};

template<>
class direction_traits<UNDIRECT>{
public:
  typedef boost::adjacency_list < boost::vecS, boost::vecS, boost::undirectedS,
  boost::no_property, boost::property < boost::edge_weight_t, double > > Graph;
};

template<DIRECTION DIRECTION_TYPE>
class vertex_connection{
  BOOST_STATIC_ASSERT((DIRECTION_TYPE == UNDIRECT || DIRECTION_TYPE == DIRECT));
public:
  typedef typename direction_traits<DIRECTION_TYPE>::Graph Graph;
  typedef typename boost::graph_traits<Graph>::vertex_descriptor Vertex;
  typedef typename boost::graph_traits<Graph>::vertex_descriptor vertex_descriptor;

  static vertex_connection<DIRECTION_TYPE> * create(
      const std::map<std::pair<size_t,size_t>,double> & edge_weight_)
  {
    std::auto_ptr<vertex_connection<DIRECTION_TYPE> > vc(
          new vertex_connection<DIRECTION_TYPE> );

    if(vc->init(edge_weight_)) // fail
      return 0;
    return vc.release();
  }

  int get_shortest_path(const size_t begin_vertex,
                        const size_t end_vertex,
                        std::vector<size_t> &path) const
  {
    typedef std::map<size_t,size_t>::const_iterator mcit;
    const Graph & g = *graph_;

    std::vector<double> distance(num_vertices(g));
    mcit begin_vertex_name_iter = smp_.find(begin_vertex);
    mcit end_vertex_name_iter = smp_.find(end_vertex);
    if(begin_vertex_name_iter == smp_.end() || end_vertex_name_iter == smp_.end()) {
//      std::cerr << "# [error] begin/end vertex does not belong to this graph."
//                << std::endl;
      return __LINE__;
    }

    vertex_descriptor s = vertex(begin_vertex_name_iter->second, g);
    std::vector<vertex_descriptor> p(num_vertices(g));
    dijkstra_shortest_paths(g, s, boost::predecessor_map(&p[0]).distance_map(&distance[0]));

    vertex_descriptor vi = vertex(end_vertex_name_iter->second,g);
    path.clear();

    while(vi != s){
      if(p[vi] == vi){
        return __LINE__;
      }
      path.push_back(vertex_name_[vi]);
      vi = p[vi];
    }
    path.push_back(vertex_name_[vi]);

    std::reverse(path.begin(),path.end());
    return 0;
  }

  // find the shortes path by using dijkstra, path records the whole path from begin_vertex to end_vertex
  // if can not find a path, return non-zeros

private:
  int init(const std::map<std::pair<size_t,size_t>,double> & edge_weight_)
  {
    typedef std::map<std::pair<size_t,size_t>,double>::const_iterator mcit;
    std::set<size_t> vertex_;
    for(mcit it = edge_weight_.begin(); it != edge_weight_.end(); ++it){
      vertex_.insert(it->first.first);
      vertex_.insert(it->first.second);
    }

    vertex_name_.resize(vertex_.size());
    copy(vertex_.begin(),vertex_.end(),vertex_name_.begin());

    edges_.reserve(edge_weight_.size());//(edge_weight.size());
    weights_.reserve(edge_weight_.size());

    for(size_t t = 0; t < vertex_name_.size(); ++t) smp_[vertex_name_[t]] = t;

    for(mcit it = edge_weight_.begin(); it != edge_weight_.end(); ++it){
      edges_.push_back(std::make_pair(smp_[it->first.first],smp_[it->first.second]));
      weights_.push_back(it->second);
    }

    graph_.reset(new Graph(edges_.begin(),edges_.end(),weights_.begin(),
                           vertex_name_.size()));
    return 0;
  }


private:
  std::vector<std::pair<size_t,size_t> > edges_;
  std::vector<double> weights_;
  std::auto_ptr<Graph> graph_;
  std::map<size_t,size_t> smp_;
  std::vector<size_t> vertex_name_;
};

#endif // VERTEX_CONNECTION_H
