#ifndef DUAL_GRAPH_H
#define DUAL_GRAPH_H

#include <boost/graph/adjacency_list.hpp>
#include <jtflib/mesh/mesh.h>

namespace riemann {

using jtf::mesh::edge2cell_adjacent;
using boost::adjacency_list;
using boost::vecS;
using boost::undirectedS;
using boost::edge_index_t;
using boost::edge_weight_t;
using boost::property;
using boost::property_map;
using Graph=adjacency_list<vecS, vecS, undirectedS, property<edge_index_t, size_t>, property<edge_weight_t, int>>;
using mati_t=zjucad::matrix::matrix<size_t>;

struct graph_t {
  size_t vert_num;
  std::vector<size_t> u;
  std::vector<size_t> v;
  std::vector<size_t> first;
  std::vector<size_t> next;
};

int build_tri_mesh_dual_graph(const mati_t &tris, std::shared_ptr<edge2cell_adjacent> &ec, std::shared_ptr<Graph> &g, const char *dotfile=nullptr);
int get_minimum_spanning_tree(const std::shared_ptr<const Graph> &g, graph_t &mst, const char *dotfile=nullptr);

}
#endif
