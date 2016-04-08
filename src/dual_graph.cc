#include "dual_graph.h"

#include <boost/graph/graph_traits.hpp>
#include <boost/graph/graphviz.hpp>
#include <boost/graph/kruskal_min_spanning_tree.hpp>

using namespace std;
using namespace zjucad::matrix;

namespace riemann {

int build_tri_mesh_dual_graph(const mati_t &tris, shared_ptr<edge2cell_adjacent> &ec,
                              shared_ptr<Graph> &g, const char *dotfile) {
  ec.reset(edge2cell_adjacent::create(tris));
  if ( ec.get() == nullptr ) {
    cerr << "[Error] can not create edge2cell\n";
    return EXIT_FAILURE;
  }
  g.reset(new Graph(tris.size(2)));
  if ( g.get() == nullptr ) {
    cerr << "[Error] new graph fail\n";
    return EXIT_FAILURE;
  }
  property_map<Graph, edge_weight_t>::type weight = get(boost::edge_weight, *g);
  for (size_t i = 0; i < ec->edges_.size(); ++i) {
    pair<size_t, size_t> faces = ec->query(ec->edges_[i].first, ec->edges_[i].second);
    if ( ec->is_boundary_edge(faces) )
      continue;
    boost::graph_traits<Graph>::edge_descriptor e; bool inserted;
    tie(e, inserted) = boost::add_edge(faces.first, faces.second, *g);
    weight[e] = 1;
  }
  if ( dotfile != nullptr ) {
    ofstream ofs(dotfile);
    write_graphviz(ofs, *g);
    ofs.close();
  }
  return EXIT_SUCCESS;
}

int get_minimum_spanning_tree(const shared_ptr<const Graph> &g, graph_t &mst, const char *dotfile) {
  if ( g.get() == nullptr ) {
    cerr << "[Error] graph object dosen't exist\n";
    return EXIT_FAILURE;
  }
  typedef boost::graph_traits<Graph>::edge_descriptor Edge;
  vector<Edge> spanning_tree;
  kruskal_minimum_spanning_tree(*g, std::back_inserter(spanning_tree));

  mst.u.resize(2*spanning_tree.size());
  mst.v.resize(2*spanning_tree.size());
  mst.first.resize(boost::num_vertices(*g));
  mst.next.resize(2*spanning_tree.size());
  std::fill(mst.first.begin(), mst.first.end(), -1);
  std::fill(mst.next.begin(), mst.next.end(), -1);

  size_t k = 0;
  for (vector<Edge>::iterator ei = spanning_tree.begin(); ei != spanning_tree.end(); ++ei) {
    mst.u[k] = source(*ei, *g);
    mst.v[k] = target(*ei, *g);
    mst.next[k] = mst.first[mst.u[k]];
    mst.first[mst.u[k]] = k;
    ++k;
    mst.u[k] = target(*ei, *g);
    mst.v[k] = source(*ei, *g);
    mst.next[k] = mst.first[mst.u[k]];
    mst.first[mst.u[k]] = k;
    ++k;
  }
  if ( dotfile != nullptr ) {
    Graph tree(boost::num_vertices(*g));
    for (vector<Edge>::iterator ei = spanning_tree.begin(); ei != spanning_tree.end(); ++ei) {
      boost::add_edge(source(*ei, *g), target(*ei, *g), tree);
    }
    ofstream ofs(dotfile);
    write_graphviz(ofs, tree);
    ofs.close();
  }
  return EXIT_SUCCESS;
}

}
