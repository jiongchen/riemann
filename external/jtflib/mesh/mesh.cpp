#include <memory>
#include <iostream>
#include <set>
#include <deque>
#include <stack>
#include <numeric>

#include <jtflib/util/container_operation.h>
#include <jtflib/util/util.h>

#include <zjucad/matrix/io.h>
#include "mesh.h"
#include "util.h"
#include "config.h"

using namespace std;
using namespace zjucad::matrix;

namespace jtf{
  namespace mesh{
	   size_t edge2cell_adjacent::get_edge_idx(size_t vi, size_t vj) const
    {
      std::pair<size_t,size_t> edge(vi,vj);
      if(vi > vj)
        swap(edge.first,edge.second);
      std::vector<std::pair<size_t,size_t> >::const_iterator i =
          lower_bound(edges_.begin(), edges_.end(),edge);
      if(i->first == edge.first && i->second == edge.second)
        return static_cast<size_t>(i - edges_.begin());
      else
        return -1;
    }

	  size_t  edge2cell_adjacent::get_edge_idx(const size_t *v) const
    {
      return get_edge_idx(v[0],v[1]);
    }

	   std::pair<size_t, size_t>  edge2cell_adjacent::query(size_t vi, size_t vj)  const
    {
      const size_t vidx = get_edge_idx(vi, vj);
      if( vidx >= edges_.size())
        return make_pair(-1,-1);
      return edge2cell_[vidx];
    }

	   std::pair<size_t, size_t> edge2cell_adjacent::query(const size_t *v) const
    {
      return query(v[0],v[1]);
    }

	   edge2cell_adjacent* edge2cell_adjacent::create(
        const matrixst & tri,
        const bool &show_debug_info)
    {
      unique_ptr<edge2cell_adjacent> ea(new edge2cell_adjacent);
      if(ea->init(tri, show_debug_info))
        return 0;
      return ea.release();
    }

	   edge2cell_adjacent*  edge2cell_adjacent::create_on_dirty_mesh(const matrixst & tri)
    {
      unique_ptr<edge2cell_adjacent> ea(new edge2cell_adjacent);
      if(ea->init_on_dirty_mesh(tri))
        return 0;
      return ea.release();
    }

	    int edge2cell_adjacent::init(const matrixst & mesh, const bool & show_debug_info)
    {
      typedef std::map< std::pair<size_t, size_t>, std::pair<size_t,size_t> > map_type;
      map_type adj_map;
      assert(mesh.size(1) == 3 || mesh.size(1) == 4);
      vector<size_t> edge(2);
      for(size_t cidx = 0; cidx < mesh.size(2); ++cidx) {
          for(size_t vidx = 0; vidx < mesh.size(1); ++vidx) {
              edge[0] = mesh(vidx,cidx);
              edge[1] = mesh((vidx + 1)%mesh.size(1),cidx);
              sort(edge.begin(),edge.end());
              bool keep_order = true;
              if(edge.front() != mesh(vidx, cidx)) keep_order = false;

              map_type::iterator it = adj_map.find(make_pair(edge[0],edge[1]));
              if(it != adj_map.end()) {
                  if(keep_order == false){
                      if(it->second.second != -1){
                          if(it->second.first == -1)
                            cerr << "# [error] face should be ordered adjacent to edge." << endl;
                          else
                            cerr << "# [error] this edge associate with >= 3 faces: " << endl
                                 << it->second.first << " " << it->second.second << " " << cidx << endl;
                          return __LINE__;
                        }
                      it->second.second = cidx;
                    }else{
                      if(it->second.first != -1){
                          if(it->second.second == -1)
                            cerr << "# [error] face should be ordered adjacent to edge." << endl;
                          else
                            cerr << "# [error] this edge associate with >= 3 faces: " << endl
                                 << it->second.first << " " << it->second.second << " " << cidx << endl;
                          return __LINE__;
                        }
                      it->second.first = cidx;
                    }
                } else {
                  pair<size_t,size_t> empty_cell(-1,-1);
                  if(keep_order == false) empty_cell.second = cidx;
                  else empty_cell.first = cidx;
                  adj_map.insert(make_pair(make_pair(edge[0],edge[1]),empty_cell));
                }
            }
        }

      const size_t edgenum = adj_map.size();
      edges_.resize(edgenum);
      edge2cell_.resize(edgenum);
      map_type::const_iterator mit = adj_map.begin();
      for(size_t eidx = 0; eidx < edgenum; ++eidx,++mit) {
          edges_[eidx] = mit->first;
          edge2cell_[eidx] = mit->second;
        }

      if(show_debug_info)
        cerr << "# create edge2cell_adjacent success with "
             << edgenum << " entries" << endl;
      return 0;
    }


    int edge2cell_adjacent::init_on_dirty_mesh(const matrixst & mesh)
    {
      typedef std::map< std::pair<size_t, size_t>, std::pair<size_t,size_t> > map_type;
      map_type adj_map;
      assert(mesh.size(1) == 3 || mesh.size(1) == 4);
      vector<size_t> edge(2);
      for(size_t cidx = 0; cidx < mesh.size(2); ++cidx) {
          for(size_t vidx = 0; vidx < mesh.size(1); ++vidx) {
              edge[0] = mesh(vidx,cidx);
              edge[1] = mesh((vidx + 1)%mesh.size(1),cidx);
              sort(edge.begin(),edge.end());

              map_type::iterator it = adj_map.find(make_pair(edge[0],edge[1]));
              if(it != adj_map.end()) {
                  if(it->second.second != -1 && it->second.first != -1){
                      cerr << "# [error] this edge associate with >= 3 faces: " << endl
                           << it->second.first << " " << it->second.second << " " << cidx << endl;
                      return __LINE__;
                    }
                  if(it->second.first == -1) it->second.first = cidx;
                  if(it->second.second == -1) it->second.second = cidx;
                } else {
                  pair<size_t,size_t> empty_cell(cidx,-1);
                  adj_map.insert(make_pair(make_pair(edge[0],edge[1]),empty_cell));
                }
            }
        }

      const size_t edgenum = adj_map.size();
      edges_.resize(edgenum);
      edge2cell_.resize(edgenum);
      map_type::const_iterator mit = adj_map.begin();
      for(size_t eidx = 0; eidx < edgenum; ++eidx,++mit) {
          edges_[eidx] = mit->first;
          edge2cell_[eidx] = mit->second;
        }

      cerr << "# create edge2cell_adjacent success with "
           << edgenum << " entries" << endl;
      return 0;
    }


    int edge2cell_adjacent_general::init(const matrixst & mesh)
    {
      typedef std::map< std::pair<size_t, size_t>, std::vector<size_t> > map_type;
      map_type adj_map;
      assert(mesh.size(1) == 3 || mesh.size(1) == 4);
      const size_t cell_n = mesh.size(2);
      vector<size_t> edge(2);
      vector<size_t> cidx_(1);
      for(size_t cidx = 0; cidx < cell_n; ++cidx) {
          for(size_t vidx = 0; vidx < mesh.size(1); ++vidx) {
              edge[0] = mesh(vidx,cidx);
              edge[1] = mesh((vidx + 1)%mesh.size(1),cidx);
              sort(edge.begin(),edge.end());
              map_type::iterator it = adj_map.find(make_pair(edge[0],edge[1]));
              if(it != adj_map.end()) {
                  it->second.push_back(cidx);
                } else {
                  cidx_[0] = cidx;
                  adj_map.insert(make_pair(make_pair(edge[0],edge[1]),cidx_));
                }

            }
        }

      const size_t edgenum = adj_map.size();
      edges_.resize(edgenum);
      edge2cell_.resize(edgenum);
      map_type::const_iterator mit = adj_map.begin();
      for(size_t eidx = 0; eidx < edgenum; ++eidx,++mit) {
          edges_[eidx] = mit->first;
          edge2cell_[eidx] = mit->second;
        }
      cerr << "# create edge2face_adjacent success with "
           << edgenum << " entries" << endl;
      return 0;
    }

    size_t edge2cell_adjacent_general::get_edge_idx(size_t vi, size_t vj) const
    {

      std::pair<size_t,size_t> edge;
      edge.first = vi;
      edge.second = vj;
      if(vi > vj)
        swap(edge.first,edge.second);
      std::vector<std::pair<size_t,size_t> >::const_iterator i =
          lower_bound(edges_.begin(), edges_.end(),edge);
      if(i->first == edge.first && i->second == edge.second)
        return static_cast<size_t>(i - edges_.begin());
      else
        return -1;
    }

    size_t edge2cell_adjacent_general::get_edge_idx(const size_t *v) const
    {
      return get_edge_idx(v[0],v[1]);
    }

    std::vector<size_t> edge2cell_adjacent_general::query(size_t vi, size_t vj)  const
    {
      const size_t vidx = get_edge_idx(vi, vj);
      if( vidx >= edges_.size()){
          vector<size_t> zero_vec;
          return zero_vec;
        }
      return edge2cell_[vidx];
    }

    std::vector<size_t> edge2cell_adjacent_general::query(const size_t *v) const
    {
      return query(v[0],v[1]);
    }

    edge2cell_adjacent_general* edge2cell_adjacent_general::create(const matrixst & mesh)
    {
      unique_ptr<edge2cell_adjacent_general> ea(new edge2cell_adjacent_general);
      if(ea->init(mesh))
        return 0;
      return ea.release();

    }


    face2tet_adjacent *face2tet_adjacent::create(
        const matrixst &tet,
        const std::string & strategy)
    {
      unique_ptr<face2tet_adjacent> fa(new face2tet_adjacent);
      if(strategy == "geometry"){
          if(fa->init_geometry(tet))
            return 0;
        }else if(strategy == "topology"){
          if(fa->init_topology(tet))
            return 0;
        }else{
          cerr << "# [error] wrong face2tet_adjacent strategy: " << strategy << endl;
          return 0;
        }
      return fa.release();
    }

    int face2tet_adjacent::init_geometry(const matrixst &tet)
    {
      strategy_ = "geometry";
      typedef std::map<std::vector<size_t>, std::pair<size_t, size_t> > map_type;
      map_type adj_map;

      assert(tet.size(1) == 4);
      const size_t tet_n = tet.size(2);
      vector<size_t> tri(3), sort_tri(3);
      for(size_t tet_i = 0; tet_i < tet_n; ++tet_i) {
          for(size_t tri_i = 0; tri_i < 4; ++tri_i) {
              for(size_t ni = 0; ni < 3; ++ni)
                tri[ni] = tet((tri_i+ni)%4, tet_i);

              sort_tri = tri;
              sort(sort_tri.begin(), sort_tri.end());

              rotate(tri.begin(), min_element(tri.begin(), tri.end()), tri.end());

              // 0th, 2th is positive face of this tet, 1th and 3th is
              // negative face of this tet
              const bool is_in_this_tet = ((tri == sort_tri)^(tri_i%2));

              map_type::iterator i = adj_map.find(sort_tri);
              if(i == adj_map.end()) {
                  if(is_in_this_tet)
                    adj_map[sort_tri] = make_pair(tet_i, -1);
                  else
                    adj_map[sort_tri] = make_pair(-1, tet_i);
                }
              else {
                  bool success = true;
                  if(is_in_this_tet) {
                      if(i->second.first == -1)
                        i->second.first = tet_i;
                      else
                        success = false;
                    }
                  else {
                      if(i->second.second == -1)
                        i->second.second = tet_i;
                      else
                        success = false;
                    }
                  if(!success) {
                      cerr << "# [error] improper tet face, maybe tets self-intesection happens: " << endl;
                      copy(tri.begin(), tri.end(), ostream_iterator<size_t>(cerr, " "));
                      cerr << endl;
                      cerr << "# [error] in tet " << tet_i << trans(tet(colon(), tet_i)) << endl;
                      cerr << "# [error] is_in_this_tet "<< (is_in_this_tet?"true ":"false ") << endl;
                      cerr << "# [error] tet_pair " << i->second.first << "," << i->second.second << endl;
                      const size_t other_tet_idx = (i->second.first==-1?i->second.second:i->second.first);
                      cerr << "# [error] in tet " << other_tet_idx << trans(tet(colon(), other_tet_idx)) << endl;
                      return __LINE__;
                    }
                }
            }
        }
      const size_t face_num = adj_map.size();
      cerr << "# create face2tet_adjacent success with " << face_num << " entries." << endl;
      faces_.resize(face_num);
      face2tet_.resize(face_num);
      map_type::const_iterator itr = adj_map.begin();
      for(size_t i = 0; i < face_num; ++i, ++itr) {
          faces_[i] = itr->first;
          face2idx_[itr->first] = i;
          face2tet_[i] = itr->second;
        }

      return 0;
    }

    size_t face2tet_adjacent::get_face_idx(size_t a, size_t b, size_t c) const
    {
      vector<size_t> tri(3);
      tri[0] = a; tri[1] = b; tri[2] = c;
      sort(tri.begin(), tri.end());
      std::map<std::vector<size_t>,size_t>::const_iterator cit =
          face2idx_.find(tri);
      if(cit == face2idx_.end()) return -1;
      return cit->second;
    }

    size_t face2tet_adjacent::get_face_idx(const size_t *abc) const
    {
      return get_face_idx(abc[0], abc[1], abc[2]);
    }


    std::pair<size_t, size_t> face2tet_adjacent::query(size_t a, size_t b, size_t c) const
    {
      size_t idx = get_face_idx(a, b, c);
      if(idx >= faces_.size())
        return make_pair(-1, -1);
      return face2tet_[idx];
    }

    std::pair<size_t, size_t> face2tet_adjacent::query(const size_t *abc) const
    {
      return query(abc[0], abc[1], abc[2]);
    }


    int face2tet_adjacent::init_topology(const matrixst &tet)
    {
      strategy_ = "topology";
      typedef std::map<std::vector<size_t>, std::pair<size_t, size_t> > map_type;
      map_type adj_map;

      assert(tet.size(1) == 4);
      const size_t tet_n = tet.size(2);
      vector<size_t> tri(3), sort_tri(3);
      for(size_t tet_i = 0; tet_i < tet_n; ++tet_i) {
          for(size_t tri_i = 0; tri_i < 4; ++tri_i) {
              for(size_t ni = 0; ni < 3; ++ni)
                tri[ni] = tet((tri_i+ni)%4, tet_i);

              sort_tri = tri;
              sort(sort_tri.begin(), sort_tri.end());
              map_type::iterator it = adj_map.find(sort_tri);
              if(it == adj_map.end())
                adj_map.insert(make_pair(sort_tri,make_pair(tet_i,-1)));
              else{
                  //assert(it->second.second == -1);
                  if(it->second.second != -1){
                      cerr << "# [error] strange face " ;
                      for(size_t pi = 0; pi < sort_tri.size(); ++pi)
                        cerr << sort_tri[pi] << " ";
                      cerr << endl;
                      cerr << "# [error] this face associated with at least three tets "
                           << it->second.first << " " << it->second.second
                           << " " << tet_i << endl;
                      return __LINE__;
                    }
                  it->second.second = tet_i;
                }
            }
        }
      const size_t face_num = adj_map.size();
      cerr << "# create face2tet_adjacent success with " << face_num << " entries." << endl;
      faces_.resize(face_num);
      face2tet_.resize(face_num);
      map_type::const_iterator itr = adj_map.begin();
      for(size_t i = 0; i < face_num; ++i, ++itr) {
          faces_[i] = itr->first;
          face2tet_[i] = itr->second;
        }

      {
        matrixst outside_face, outside_face_idx;
        get_outside_face(*this, outside_face,true);
        get_outside_face_idx(*this, outside_face_idx);
        for(size_t fi = 0; fi < outside_face_idx.size(); ++fi){
            copy(outside_face(colon(),fi).begin(),
                 outside_face(colon(),fi).end(),
                 faces_[outside_face_idx[fi]].begin());
          }
        for(size_t fi = 0; fi < faces_.size(); ++fi){
            vector<size_t> one_face = faces_[fi];
            sort(one_face.begin(), one_face.end());
            face2idx_[one_face] = fi;
          }
      }

      return 0;
    }

    bool is_surface_mesh_manifold(
        const matrixst & surface,
        const dyn_one_ring_face_at_point * dorfap_ptr)
    {
      unique_ptr<edge2cell_adjacent> ea(edge2cell_adjacent::create(surface,false));
      if(!ea.get()) return false;

      if(dorfap_ptr){
          typedef dyn_one_ring_face_at_point dorfap_type;
          for(dorfap_type::pcit it = dorfap_ptr->p2f_.begin();
              it != dorfap_ptr->p2f_.end(); ++it){
              size_t adjecent_face_num = 0;
              set<size_t> adjacent_points;

              const dorfap_type::linked_face_ptr_types & lfpt = it->second;
              adjecent_face_num = lfpt.size();

              for(dorfap_type::linked_face_const_iterator lit = lfpt.begin();
                  lit != lfpt.end(); ++lit){
                  const dorfap_type::face_list_type::iterator fit = *lit;
                  const dorfap_type::face_type & face = *fit;
                  for(size_t i = 0; i < face.size(); ++i){
                      if(face[i] != it->first)
                        adjacent_points << face[i];
                    }
                }
              if(adjacent_points.size() != adjecent_face_num &&
                 adjacent_points.size() != adjecent_face_num + 1){
                  cerr << "# [error] non-manifold point" << endl;
                  return false;
                }
            }
        }
      return true;
    }

    int get_boundary_edge(const edge2cell_adjacent& ea, matrixst &edge)
    {
      vector<size_t> edge_;
      for(size_t ei = 0; ei < ea.edges_.size(); ++ei)
        {
          const pair<size_t,size_t> &etri = ea.edge2cell_[ei];
          if(ea.is_boundary_edge((etri)))
            {
              edge_.push_back(ea.edges_[ei].first);
              edge_.push_back(ea.edges_[ei].second);
            }
        }
      edge.resize(2,edge_.size()/2);
      copy(edge_.begin(),edge_.end(),edge.begin());
      return 0;
    }

    int get_boundary_edge_idx(const edge2cell_adjacent& ea, matrixst &edge_idx)
    {
      vector<size_t> edge_idx_;
      for(size_t ei = 0; ei < ea.edges_.size(); ++ei)
        {
          const pair<size_t,size_t> &etri = ea.edge2cell_[ei];
          if(ea.is_boundary_edge((etri)))
            edge_idx_.push_back(ei);
        }
      edge_idx.resize(edge_idx_.size());
      copy(edge_idx_.begin(),edge_idx_.end(),edge_idx.begin());
      return 0;
    }

    void get_outside_face(const face2tet_adjacent &m, matrixst &face,
                          bool check_order, const matrixd * node)
    {
      vector<size_t> face0;
      for(size_t i = 0; i < m.faces_.size(); ++i) {
          const pair<size_t, size_t> &tet_i = m.face2tet_[i];
          if(!m.is_outside_face(tet_i)) continue;

          const vector<size_t> &tri = m.faces_[i];
          if(tet_i.first >= 0)
            face0.insert(face0.end(), tri.begin(), tri.end());
          else {
              face0.push_back(tri[0]);
              face0.push_back(tri[2]);
              face0.push_back(tri[1]);
            }
        }
      face.resize(3, face0.size()/3);
      copy(face0.begin(), face0.end(), face.begin());

      if(check_order){
          reorder_face(face,node);
        }
    }

    void get_outside_face_idx(const face2tet_adjacent &m,
                              matrixst &face_idx)
    {
      vector<size_t> face0;
      for(size_t i = 0; i < m.faces_.size(); ++i) {
          if(m.is_outside_face(m.face2tet_[i]))
            face0.push_back(i);
        }
      face_idx.resize(face0.size());
      copy(face0.begin(), face0.end(), face_idx.begin());
    }

    inline size_t adjust_top_z_face(matrixst & face, const matrixd & node)
    {
      vector<double> face_bary_center_z (face.size(2),0);
      for(size_t fi = 0; fi < face.size(2); ++fi){
          for(size_t pi = 0; pi < face.size(1); ++pi)
            face_bary_center_z[fi] += node(2, face(pi, fi));
          face_bary_center_z[fi] /= face.size(1);
        }
      auto it = std::max_element(face_bary_center_z.begin(), face_bary_center_z.end());

      const size_t top_z_face_idx = static_cast<size_t>(it-face_bary_center_z.begin());
      matrix<double> face_normal(3,1);
      cal_face_normal(node(colon(),face(colon(), top_z_face_idx)), face_normal);
      if(face_normal[2] < 0) // normal is with the same direction of -z
        std::reverse(face(colon(), top_z_face_idx).begin()+1, face(colon(), top_z_face_idx).end());
      return top_z_face_idx;
    }

    int reorder_face(matrixst & face, const matrixd * node)
    {
      unique_ptr<jtf::mesh::edge2cell_adjacent> ea(
            jtf::mesh::edge2cell_adjacent::create_on_dirty_mesh(face));
      if(!ea.get()){
          cerr << "# [error] can not build edge2cell_adjacent." << endl;
          return __LINE__;
        }

      vector<bool> is_face_visited(face.size(2), false);

      stack<size_t> face_stack;
      if(node == nullptr)
        face_stack.push(0);
      else{
          size_t seed = adjust_top_z_face(face, *node);
          face_stack.push(seed);
        }

      while(!face_stack.empty()){
          const size_t face_idx = face_stack.top();
          face_stack.pop();
          is_face_visited[face_idx] = true;
          for(size_t pi = 0; pi < face.size(1); ++pi){
              const size_t & p0 = face(pi,face_idx);
              const size_t & p1 = face((pi+1)%face.size(1), face_idx);
              const size_t edge_idx = ea->get_edge_idx(p0,p1);
              if(edge_idx == -1){
                  cerr << "# [error] strange can not find edge " << p0 << " " << p1 << endl;
                  return __LINE__;
                }
              const pair<size_t,size_t> & tri_pair = ea->edge2cell_[edge_idx];
              if(ea->is_boundary_edge(tri_pair)) continue;
              const size_t & other_tri_idx = tri_pair.first + tri_pair.second - face_idx;
              assert(other_tri_idx < face.size(2));
              if(is_face_visited[other_tri_idx]) continue;
              matrixst::const_iterator cit =
                  find(face(colon(),other_tri_idx).begin(),
                       face(colon(),other_tri_idx).end(), p0);
              assert(cit != face(colon(), other_tri_idx).end());
              const size_t idx = static_cast<size_t>(cit - face(colon(),other_tri_idx).begin());
              if(face((idx+1)%face.size(1),other_tri_idx) == p1)
                reverse(face(colon(), other_tri_idx).begin(), face(colon(),other_tri_idx).end());
              face_stack.push(other_tri_idx);
            }
        }

      return 0;
    }

    void N_ring_face_at_point::add_all_faces(
        const matrixst &face,
        const edge2cell_adjacent &ea,
        const size_t N)
    {
      if(N == 0) {
          cerr << "# [error] strange, no need to find 0 ring face for point." << endl;
          return;
        }

      one_ring_face_at_point orfap;
      orfap.add_all_faces(face, ea);

      for(one_ring_face_at_point::p2f_type::const_iterator
          cit = orfap.p2f_.begin(); cit != orfap.p2f_.end(); ++cit){

          p2f_[cit->first].push_back(cit->second);

          for(size_t i = 1; i < N; ++i){
              const vector<size_t> & boundary_faces_vec = p2f_[cit->first].back();
              vector<size_t> new_ring_face;

              // the former layer of faces
              for(size_t fi = 0; fi < boundary_faces_vec.size(); ++fi){

                  for(size_t pi = 0; pi < face.size(1); ++pi){
                      if(boundary_faces_vec[fi] == -1) continue; // ignore boundary out face
                      const size_t & point_idx = face(pi, boundary_faces_vec[fi]);
                      const vector<size_t> & one_ring_face_vec = orfap.p2f_[point_idx];
                      for(size_t ffi = 0; ffi < one_ring_face_vec.size(); ++ffi){
                          if(find(boundary_faces_vec.begin(), boundary_faces_vec.end(),
                                  one_ring_face_vec[ffi]) == boundary_faces_vec.end())
                            new_ring_face.push_back(one_ring_face_vec[ffi]);
                        }
                    }
                }
              p2f_[cit->first].push_back(new_ring_face);
            }
        }
    }

    void one_ring_face_at_point::add_face(
        const matrixst &one_face,
        const edge2cell_adjacent & ea)
    {
      assert(one_face.size() == 3);
      for(size_t t = 0 ; t < one_face.size(); ++t){
          const size_t edge_idx =
              ea.get_edge_idx(one_face[t],
                              one_face[(t+1)%one_face.size()]);
          if(edge_idx == -1) {
              cerr << "# [error] this edge is invalid. < "
                   << one_face[t] << "," << one_face[(t+1)%one_face.size()]
                   << ">." << endl;
            }
          const pair<size_t,size_t> & face_pair = ea.edge2cell_[edge_idx];
          vector<size_t> & f0 = p2f_[one_face[t]];
          f0.push_back(face_pair.first);
          f0.push_back(face_pair.second);

		  vector<size_t> & f1 = p2f_[one_face[t]];
          f1.push_back(face_pair.first);
          f1.push_back(face_pair.second);
        }
    }

    int one_ring_face_at_point::sort_int_loop(
        const matrixst &face,
        const matrixd & node)
    {
      for(p2f_type::iterator pcit = p2f_.begin(); pcit != p2f_.end(); ++pcit){
          vector<size_t> & loop_face = pcit->second;
          std::set<pair<size_t,size_t> > segments;
          for(size_t t = 0 ;t < loop_face.size()/2; ++t){
              segments.insert(make_pair(loop_face[2*t],
                              loop_face[2*t+1]));
            }
          deque<size_t>  temp_loop;
          temp_loop.push_back(segments.begin()->first);
          temp_loop.push_back(segments.begin()->second);
          segments.erase(segments.begin());

          while(!segments.empty()){
              const size_t segments_prev_size = segments.size();

              for(set<pair<size_t,size_t> >::iterator sit = segments.begin();
                  sit != segments.end(); ++sit){
                  if(temp_loop.back() != -1){
                      if(sit->first == temp_loop.back() ||
                         sit->second == temp_loop.back()){
                          temp_loop.push_back(sit->second + sit->first
                                              - temp_loop.back());
                          segments.erase(sit);
                          break;
                        }
                    }
                  if(temp_loop.front() != -1){
                      if(sit->first == temp_loop.front() ||
                         sit->second == temp_loop.front()){
                          temp_loop.push_front(sit->first + sit->second
                                               - temp_loop.front());
                          segments.erase(sit);
                          break;
                        }
                    }
                }

              if(segments.size() == segments_prev_size){
                  cerr << "# [error] strange: possible edges can not be glued any more." << endl;
                  std::cerr << "# [error] wrong face loop: ";
                  return __LINE__;
                }
            }

          loop_face.resize(temp_loop.size());
          copy(temp_loop.begin(),temp_loop.end(),loop_face.begin());
        }
      fix_boundary_loop();
      return 0;
    }

    void one_ring_face_at_point::fix_boundary_loop()
    {
      for(p2f_type::iterator pcit = p2f_.begin(); pcit != p2f_.end(); ++pcit){
          vector<size_t> & loop_face = pcit->second;
          if(loop_face.front() != loop_face.back()){
              if(loop_face.front() == -1) loop_face.push_back(-1);
              else if(loop_face.back() == -1) loop_face.insert(loop_face.begin(),-1);
              else {
                  std::cerr << "# [error] strang one ring face." << std::endl;
                }
            }
        }
    }

    void one_ring_face_at_point::sort_int_loop_with_normal_info(
        const matrixst &face,
        const matrixd & node,
        const edge2cell_adjacent &ea,
        const matrixd & normal)
    {
      sort_int_loop(face,node);
      size_t face_pair[2];
      size_t edge[2];
      matrixd average_face_normal = zjucad::matrix::zeros<double>(3,1);
      for(p2f_type::iterator ptit = p2f_.begin(); ptit != p2f_.end(); ++ptit){
          vector<size_t> & loop_face = ptit->second;
          // if(count(loop_face.begin(),loop_face.end(),-1) == 1) continue; // sharp edge
          if(loop_face.front() == -1 && loop_face.back() == -1 && loop_face.size() == 3) continue; //sharp edge
          assert(loop_face.front() == loop_face.back());
          average_face_normal = zjucad::matrix::zeros<double>(3,1);
          for(size_t t = 0; t < loop_face.size(); ++t){
              if(loop_face[t] != -1){
                  average_face_normal += normal(zjucad::matrix::colon(),loop_face[t]);
                }
            }
          average_face_normal /= zjucad::matrix::norm(average_face_normal);
          assert(loop_face.size() > 2);
          if(loop_face.front() == -1){
              face_pair[0] = loop_face[1];
              face_pair[1] = loop_face[2];
            }else{
              face_pair[0] = loop_face.front();
              face_pair[1] = loop_face[1];
            }
          find_common_edge(face(zjucad::matrix::colon(),face_pair[0]),
              face(zjucad::matrix::colon(),face_pair[1]),
              &edge[0]);
          size_t other_point = face(0,face_pair[1]) +
              face(1,face_pair[1]) +
              face(2,face_pair[1]) -
              edge[0] - edge[1];
          matrixd edge_out =
              node(zjucad::matrix::colon(),edge[0] + edge[1] - ptit->first)
              - node(zjucad::matrix::colon(),ptit->first);
          matrixd edge_around =
              node(zjucad::matrix::colon(),other_point) -
              node(zjucad::matrix::colon(),edge[0] + edge[1] - ptit->first);
          if(zjucad::matrix::dot(zjucad::matrix::cross(edge_out,edge_around),
                                 average_face_normal) < 0){
              reverse(loop_face.begin(),loop_face.end());
            }
        }
    }

    void dyn_one_ring_face_at_point::add_all_faces(
        const matrixst &faces)
    {
      boost::unordered_map<size_t, vector<size_t> > point_to_face_idx;
      for(size_t t = 0; t < faces.size(2); ++t){
          for(size_t i = 0; i < faces.size(1); ++i){
              point_to_face_idx[faces(i,t)].push_back(t);
            }
        }

      faces_.clear();
      vector<size_t> one_face(faces.size(1));
      vector<face_list_type::iterator> face_to_iter(faces.size(2));

      for(size_t t = 0; t < faces.size(2); ++t){
          copy(faces(zjucad::matrix::colon(),t).begin(),
               faces(zjucad::matrix::colon(),t).end(), one_face.begin());
          sort(one_face.begin(), one_face.end());
          faces_ << one_face;
          face_list_type::iterator it_last = faces_.end();
          face_to_iter[t] = --it_last;
        }

      size_t show_idx = 0;
      const size_t per_100 = point_to_face_idx.size() / 100;
      for(boost::unordered_map<size_t, vector<size_t> >::const_iterator msvcit
          = point_to_face_idx.begin(); msvcit != point_to_face_idx.end();
          ++msvcit, ++show_idx){
          if(point_to_face_idx.size() > 10000){
              if(show_idx % per_100 == 0){
                  cerr << "# [info] construct dyn_one_ring_face_at_pint "
                       << show_idx / per_100 <<  "%" << endl;
                }
            }
          linked_face_ptr_types & lfpt = p2f_[msvcit->first];
          for(size_t t = 0; t < msvcit->second.size(); ++t){
              lfpt << face_to_iter[msvcit->second[t]];
            }
        }
    }

    void dyn_one_ring_face_at_point::add_face(
        const matrixst &one_face)
    {
      if(one_face.size() != 3){
          cerr << "# [error] sorry, I can only handle triangle face." << endl;
          return;
        }
      vector<size_t> one_face_vec(one_face.size());
      copy(one_face.begin(), one_face.end(), one_face_vec.begin());

      sort(one_face_vec.begin(), one_face_vec.end());
      faces_.push_back(one_face_vec);
      face_list_type::iterator this_face_ptr = faces_.begin();
      std::advance(this_face_ptr, faces_.size() - 1);//faces_.end()-1;
      for(size_t i = 0; i < one_face_vec.size(); ++i){
          p2f_[one_face_vec[i]].push_back(this_face_ptr);
        }
    }

    int dyn_one_ring_face_at_point::remove_face(
        const size_t *face, const size_t point_num)
    {
      assert(face);
      if(point_num != 3){
          cerr << "# [error] sorry, currently I can only handle triangle face." << endl;
          return __LINE__;
        }
      zjucad::matrix::itr_matrix<const size_t*> one_face(point_num, 1, &face[0]);
      vector<size_t> one_face_vec(point_num);
      copy(one_face.begin(), one_face.end(), one_face_vec.begin());
      sort(one_face_vec.begin(), one_face_vec.end());

      face_list_type::iterator fit =
          find(faces_.begin(), faces_.end(), one_face_vec);
      if(fit == faces_.end()) {
          cerr << "# [error] can not find face: ";
          copy(one_face.begin(), one_face.end(),
               std::ostream_iterator<size_t>(cerr, ","));
          cerr << endl;
          return __LINE__;
        }

      // found this face
      for(size_t i = 0; i < one_face_vec.size(); ++i){
          point2face_type::iterator pit = p2f_.find(one_face_vec[i]);
          if(pit == p2f_.end()){
              cerr << "# [error] can not find point " << one_face_vec[i]
                      << " in point2face map" << endl;
              return __LINE__;
            }
          linked_face_ptr_types & linked_faces = pit->second;
          linked_face_ptr_types::iterator lit =
              find(linked_faces.begin(), linked_faces.end(), fit);

          if(lit == linked_faces.end()){
              cerr << "# [error] strange. can not find face: ";
              copy(one_face_vec.begin(), one_face_vec.end(),
                   ostream_iterator<size_t>(cerr, ","));
              cerr << endl;
              cerr << "# ------- in " << one_face_vec[i]
                      << " one ring faces list." << endl;
              return __LINE__;
            }
          linked_faces.erase(lit); // delete the face
        }
      faces_.erase(fit); // remove current face
      return 0;
    }

    bool dyn_one_ring_face_at_point::is_face_inside(
        const std::vector<size_t> & sorted_face) const
    {
      face_list_type::const_iterator fcit =
          find(faces_.begin(), faces_.end(), sorted_face);
      if(fcit == faces_.end()) return false;
      return true;
    }

    int dyn_one_ring_face_at_point::glue_two_faces_at_shared_edge(
        const std::pair<size_t,size_t> & edge,
        const size_t & point_from, const size_t & point_to)
    {
      // check the validity of both faces
      vector<size_t> face0, face1;
      face0 << edge.first << edge.second << point_from;
      face1 << edge.first << edge.second << point_to;

      sort(face0.begin(), face0.end());
      sort(face1.begin(), face1.end());

      if(!is_face_inside(face0)){
          cerr << "# [error] can not find face in map: ";
          copy(face0.begin(), face0.end(), ostream_iterator<size_t>(cerr, ","));
          cerr << endl;
          return __LINE__;
        }
      if(!is_face_inside(face1)){
          cerr << "# [error] can not find face in map: ";
          copy(face1.begin(), face1.end(), ostream_iterator<size_t>(cerr, ","));
          cerr << endl;
          return __LINE__;
        }

      if(remove_face(&face0[0], face0.size())) {
          cerr << "# [error] can not remove face ";
          copy(face0.begin(), face0.end(), ostream_iterator<size_t>(cerr, ","));
          cerr << endl;
          return __LINE__;
        }

      if(remove_face(&face1[0], face1.size())) {
          cerr << "# [error] can not remove face ";
          copy(face1.begin(), face1.end(), ostream_iterator<size_t>(cerr, ","));
          cerr << endl;
          return __LINE__;
        }

      // merge point from point_from to point_to
      point2face_type::iterator pit_from = p2f_.find(point_from);
      point2face_type::iterator pit_to = p2f_.find(point_to);
      if(pit_from == p2f_.end() || pit_to == p2f_.end()){
          cerr << "# [error] can not find point " << point_from << " or "
               << point_to << " in point2face map." << endl;
          return __LINE__;
        }

      // modify all faces which has point as "point_from" to "point_to"

      linked_face_ptr_types & linked_faces_to = pit_to->second;
      linked_face_ptr_types & linked_faces_from = pit_from->second;

      for(linked_face_ptr_types::iterator lfptit = linked_faces_from.begin();
          lfptit != linked_faces_from.end();){
          face_type one_face_to_sort = *(*lfptit);
          face_type & one_face_orig = *(*lfptit);
          face_type::iterator point_from_ptr =
              find(one_face_to_sort.begin(), one_face_to_sort.end(), point_from);
          if(point_from_ptr == one_face_to_sort.end()){
              cerr << "# [error] strange, can not find " << point_from
                   << " in face ";
              copy(one_face_to_sort.begin(), one_face_to_sort.end(),
                   ostream_iterator<size_t>(cerr, ","));
              cerr << endl;
              return __LINE__;
            }
          *point_from_ptr = point_to;
          sort(one_face_to_sort.begin(), one_face_to_sort.end()); // to keep it ordered
          // check whether this face occure on point_to's arrounding,
          // if so remove them both,
          linked_face_iterator lfit =
              find_face_in_linked_ptr(linked_faces_to, one_face_to_sort);
          if(lfit == linked_faces_to.end()) {
              one_face_orig = one_face_to_sort;
              ++lfptit;
              continue;
            }else{
              // find a face which is the same as the sorted "one_face", it means
              // these two face should be removed, since they will be glued inside
              ++lfptit;
              remove_face(&one_face_orig[0], one_face_orig.size());
              const face_type & to_face_orig = *(*lfit);
              remove_face(&to_face_orig[0], to_face_orig.size());
            }
        }
      unite(linked_faces_to, linked_faces_from);
      p2f_.erase(pit_from);

      return 0;
    }

    dyn_one_ring_face_at_point::linked_face_iterator
    dyn_one_ring_face_at_point::find_face_in_linked_ptr(
        dyn_one_ring_face_at_point::linked_face_ptr_types & lfpt,
        const dyn_one_ring_face_at_point::face_type & face)const
    {
      assert(face.size() == 3);
      assert(face[1] - face[0] <= face[1]); // check it sorted
      assert(face[2] - face[1] <= face[2]);

      for(linked_face_iterator lfit = lfpt.begin();
          lfit != lfpt.end(); ++lfit){
          if(face == *(*lfit)) return lfit;
        }

      return lfpt.end();
    }

    int dyn_one_ring_face_at_point::convert_face_to_matrix(
        matrixst & faces) const
    {
      faces.resize(3, faces_.size());
      size_t face_idx = 0;
      for(face_list_type::const_iterator fcit = faces_.begin();
          fcit != faces_.end(); ++fcit, ++face_idx){
          const face_type &one_face = *fcit;
          copy(one_face.begin(), one_face.end(),
               faces(zjucad::matrix::colon(),face_idx).begin());
        }
      return 0;
    }

    int euler_number::add_face(const size_t * face, const size_t &vertex_num){
      assert(face);
      if(vertex_num != 3 && vertex_num != 4) {
          std::cerr << "# [error] sorry, I can only handle triangle face or "
                       "quad face"  << std::endl;
          return __LINE__;
        }
      if(vertex_num == 3) return add_tri_face(face);
      if(vertex_num == 4) return add_quad_face(face);
    }

    int euler_number::add_tri_face(const size_t * face){
      assert(face);
      std::set<size_t> one_face(face, face+3);
      vertex_.insert(one_face.begin(), one_face.end());
      for(size_t i = 0; i < 3; ++i){
          if(face[i] > face[(i+1)%3])
            edge_.insert(std::make_pair(face[(i+1)%3], face[i]));
          else
            edge_.insert(std::make_pair(face[i], face[(i+1)%3]));
        }
      face_.insert(one_face);
      return 0;
    }

    int euler_number::add_quad_face(const size_t * face){
      assert(face);
      std::set<size_t> one_face(face, face+4);
      vertex_.insert(one_face.begin(), one_face.end());
      for(size_t i = 0; i < 4; ++i){
          if(face[i] > face[(i+1)%4])
            edge_.insert(std::make_pair(face[(i+1)%4], face[i]));
          else
            edge_.insert(std::make_pair(face[i], face[(i+1)%4]));
        }
      face_.insert(one_face);
      return 0;
    }


    int find_common_edge(const matrixst & face0,
                         const matrixst & face1,
                         size_t * edge)
    {
      vector<size_t> edge_p;
      for(size_t t = 0; t < face0.size(); ++t){
          if(std::find(face1.begin(),face1.end(),face0[t]) != face1.end())
            edge_p.push_back(face0[t]);
        }
      if(edge_p.size() != 2) return 1; // two face share more than one edge
      copy(edge_p.begin(), edge_p.end(), edge);
      return 0;
    }

    one_ring_point_at_point * one_ring_point_at_point::create(
        const matrixst &faces)
    {
      unique_ptr<one_ring_point_at_point> orpap(new one_ring_point_at_point);
      if(orpap->init(faces))
        return 0;
      return orpap.release();
    }

    void one_ring_point_at_point::sort_into_loop(
        const matrixst & face,
        const edge2cell_adjacent & ea)
    {
      using namespace zjucad::matrix;
      for(p2p_type::iterator it = p2p_.begin(); it != p2p_.end(); ++it){
          vector<size_t> & around_points = it->second;

          vector<pair<size_t,size_t> > edges;
          for(size_t p0 = 0; p0 < around_points.size(); ++p0){
              for(size_t p1 = p0 + 1; p1 < around_points.size(); ++p1){
                  const size_t edge_idx = ea.get_edge_idx(around_points[p0],around_points[p1]);
                  if(edge_idx != -1){
                      const pair<size_t,size_t> & face_pair = ea.edge2cell_[edge_idx];
                      if(!ea.is_boundary_edge(face_pair)){
                          size_t  oppo_vertex[2] = {
                            std::accumulate(face(colon(),face_pair.first).begin(),
                            face(colon(),face_pair.first).end(),0) - around_points[p0] - around_points[p1],
                            std::accumulate(face(colon(),face_pair.second).begin(),
                            face(colon(),face_pair.second).end(),0) - around_points[p0] - around_points[p1],
                          };
                          if(oppo_vertex[0] != it->first && oppo_vertex[1] != it->first)
                            continue;
                        }
                            edges.push_back(make_pair(around_points[p0], around_points[p1]));
                    }
                }
            }

          vector<deque<pair<size_t,size_t> > > chains;
          jtf::util::extract_chain_from_edges(edges, chains);
          if(chains.size() != 1){
              cerr << "# [error] arounding point chain should be only one "
                   << "for a manifold mesh." << endl;
              return ;
            }
          around_points.clear();
          for(size_t ei = 0; ei < chains.front().size(); ++ei){
              around_points.push_back(chains[0][ei].first);
            }
          around_points.push_back(chains[0].back().second);
        }
    }

    int one_ring_point_at_point::init(const matrixst & faces)
    {
      for(size_t fi = 0; fi < faces.size(2); ++fi){
          for(size_t pi = 0; pi < faces.size(1); ++pi){
          	vector<size_t> & ff = p2p_[faces(pi,fi)];
              for(size_t lp = 0; lp < faces.size(1) -1 ; ++lp)
                ff.push_back(faces((pi+lp+1)%faces.size(1),fi));
            }
        }

      for(p2p_type::iterator it = p2p_.begin(); it != p2p_.end(); ++it){
          vector<size_t> & around_points = it->second;
          std::sort(around_points.begin(), around_points.end());
          auto itend = std::unique(around_points.begin(), around_points.end());
          around_points.resize(std::distance(around_points.begin(),itend));
        }
      return 0;
    }

    void one_ring_face_at_edge::build(
        const one_ring_tet_at_edge & ortae,
        const jtf::mesh::face2tet_adjacent & fa,
        const matrixst & tet)
    {
      set<size_t> one_ring_points;
      for(one_ring_tet_at_edge::e2tet_type::const_iterator cit = ortae.e2t_.begin();
          cit != ortae.e2t_.end(); ++cit){
          one_ring_points.clear();
          const vector<size_t> & loop = cit->second;
          for(size_t i = 0; i < loop.size(); ++i){
              if(loop[i] == -1) continue;
              for(size_t j = 0; j < tet.size(1); ++j){
                  if(tet(j, loop[i]) != cit->first.first &&
                     tet(j, loop[i]) != cit->first.second)
                    one_ring_points.insert(tet(j, loop[i]));
                }
            }
            vector<size_t> & e2f_cit = e2f_[cit->first];
          for(set<size_t>::const_iterator scit = one_ring_points.begin();
              scit != one_ring_points.end(); ++scit){
              const size_t face_idx =
                  fa.get_face_idx(cit->first.first,cit->first.second, *scit);
              if(face_idx == -1){
                  cerr << "# [error] strange can not find face " << face_idx << endl;
                  return ;
                }
              e2f_cit.push_back(face_idx);
            }
        }
    }

    void one_ring_tet_at_edge::add_one_tet(
        const matrixst & one_tet,
        const jtf::mesh::face2tet_adjacent &fa)
    {
      for(int a = 0; a < 4; ++a)  {
          for(int b = a + 1; b < 4; ++b)
            {
              pair<size_t, size_t> edge = make_edge(one_tet[a],one_tet[b]);
              if(edge.first > edge.second)
                swap(edge.first, edge.second);
              // for edge (a, b)
              size_t cd[2]; // the other two nodes
              for(int i = 0, k = 0; i < 4; ++i) {
                  if(i == a || i == b) continue;
                  cd[k] =one_tet[i];
                  ++k;
                }
              const pair<size_t,size_t> &tet_pair_i =
                  fa.face2tet_[fa.get_face_idx(edge.first, edge.second, cd[0])];

              assert(tet_pair_i.first != -1 || tet_pair_i.second != -1);

              const pair<size_t,size_t> &tet_pair_j =
                  fa.face2tet_[fa.get_face_idx(edge.first, edge.second, cd[1])];

              assert(tet_pair_j.first != -1 || tet_pair_j.second != -1);

              size_t tet_a = tet_pair_i.first;
              size_t tet_b = tet_pair_i.second;

              if(tet_a > tet_b) swap(tet_a,tet_b);

              e2t_[edge].push_back(tet_a);
              e2t_[edge].push_back(tet_b);

              tet_a = tet_pair_j.first;
              tet_b = tet_pair_j.second;

              if(tet_a > tet_b) swap(tet_a,tet_b);

              e2t_[edge].push_back(tet_a);
              e2t_[edge].push_back(tet_b);
            }
        }
    }

    void one_ring_tet_at_edge::add_tets(
        const matrixst & tet,
        const jtf::mesh::face2tet_adjacent & fa)
    {
      for(size_t t = 0; t < tet.size(2); ++t){
          add_one_tet(tet(colon(),t),fa);
        }
    }

    int one_ring_tet_at_edge::sort_into_loop(
        const matrixst &tet,
        const matrixd &node,
        const int dir)
    {
      for( e2tet_type::iterator eti = e2t_.begin(); eti != e2t_.end(); ++eti) { // for each edge
          vector<size_t> &tet_couple = eti->second;
          boost::unordered_set<pair<size_t, size_t> > segments;
          for(size_t i = 0; i < tet_couple.size()/2; ++i)
            segments.insert(make_pair(tet_couple[i*2], tet_couple[i*2+1]));

          deque<size_t> loop;
          {// get the first pair into loop
            boost::unordered_set<pair<size_t,size_t> >::iterator first_sp = segments.begin();
            loop.push_back(first_sp->first);
            loop.push_back(first_sp->second);
            segments.erase(first_sp);
          }

          while(!segments.empty()) {
              boost::unordered_set<pair<size_t, size_t> >::iterator i;
              for(i = segments.begin(); i != segments.end(); ++i) {
                  if(loop.front() != -1 && (loop.front() == i->first || loop.front() == i->second)) {
                      loop.push_front(i->first+i->second-loop.front());
                      break;
                    }
                  if(loop.back() != -1 && (loop.back() == i->first || loop.back() == i->second)) {
                      loop.push_back(i->first+i->second-loop.back());
                      break;
                    }
                }
              if(i == segments.end()) {
                  cerr << "# strange: can not extract arounding tet loop  for edge: <"
                       << eti->first.first << "," << eti->first.second << ">" << endl;
                  copy(tet_couple.begin(), tet_couple.end(), ostream_iterator<size_t>(cerr, " "));
                  cerr << endl;
                  return __LINE__;
                }
              segments.erase(i);
            }

          if(loop.front() != loop.back())
            {
              if(loop.size() == 2)
                {
                  assert(loop.front() == -1 || loop.back() == -1);
                  if(loop.front() == -1)
                    loop.push_back(-1);
                  else  if(loop.back() == -1)
                    loop.push_front(-1);
                  //cerr << "# edge (" << eti->first.first << "," << eti->first.second << ") is a surface edge" << endl;
                }
              else{
                  cerr << "# strange wrong loop" << endl;
                  copy(tet_couple.begin(), tet_couple.end(), ostream_iterator<size_t>(cerr, " "));
                  cerr << endl << "## loop order" << endl;
                  copy(loop.begin(), loop.end(), ostream_iterator<size_t>(cerr, " "));
                  cerr << endl;
                  return __LINE__;
                }
            }

          assert(loop.back() == loop.front());

          if(dir == 0) {
              tet_couple.resize(loop.size());
              copy(loop.begin(), loop.end(), tet_couple.begin());
              continue;// do no need to rerange the tet loop
            }
          // a sharp edge which contain only one tet can not be reranged
          if(loop.front() == -1 && loop.size() == 3){
              tet_couple.resize(loop.size());
              copy(loop.begin(), loop.end(), tet_couple.begin());
              continue;
            }

          {// adjust the loop order, right hand order around the   v_first --> v_second
            size_t t1,t2;
            size_t v0,v1,v2;
            vector<size_t> vt1,vt2;
            if(loop.front() == -1)
              {
                t1 = loop[1];
                t2 = loop[2];
              }else
              {
                t1 = loop[0];
                t2 = loop[1];
              }
            const pair<size_t,size_t> &pair_ti = eti->first;

            for(size_t t = 0; t < 4; ++t) {
                if(tet(t,t1) != pair_ti.first && tet(t,t1) != pair_ti.second)
                  vt1.push_back(tet(t,t1));
                if(tet(t,t2) != pair_ti.first && tet(t,t2) != pair_ti.second)
                  vt2.push_back(tet(t,t2));
              }
            if(find(vt2.begin(),vt2.end(),vt1[0]) != vt2.end()) {
                v0 = vt1[0];
                v1 = vt1[1];
                v2 = vt2[0] + vt2[1] - v0;
              }else if(find(vt2.begin(),vt2.end(), vt1[1]) != vt2.end()){
                v0 = vt1[1];
                v1 = vt1[0];
                v2 = vt2[0] + vt2[1] - v0;
              }else{
                cerr << "# strange adjacent tet" << endl;
                return __LINE__;
              }

            matrixd edge = node(colon(),pair_ti.second) - node(colon(),pair_ti.first);
            matrixd t1e1 = node(colon(),v0) - node(colon(),v1);
            matrixd t1e2 = node(colon(),pair_ti.second) - node(colon(),v0);
            if(dot(cross(t1e1,t1e2 ),edge) < -1e-8 ) {
                reverse(loop.begin(),loop.end());
              }
          }
          tet_couple.resize(loop.size());
          if(dir == 1) {
              copy(loop.begin(), loop.end(), tet_couple.begin());
            }// right hand order
          else if(dir == -1) {// left hand order
              reverse(loop.begin(), loop.end());
              copy(loop.begin(), loop.end(), tet_couple.begin());
            }
          assert(tet_couple.front() == tet_couple.back());
        }
      return 0;
    }

    int one_ring_tet_at_edge::get_one_ring_points_of_edge(
        const matrixst & tet,
        const std::pair<size_t,size_t> & edge,
        std::vector<size_t> & points) const
    {
      e2tet_type::const_iterator it = e2t_.find(edge);
      if(it == e2t_.end())
        it = e2t_.find(make_pair(edge.second, edge.first));
      if(it == e2t_.end())
        return __LINE__;
      const vector<size_t> & loop = it->second;
      set<size_t> point_set;
      for(size_t t = 0; t < loop.size(); ++t){
          if(loop[t] != -1){
              for(size_t pi = 0; pi < tet.size(1); ++pi){
                  if(tet(pi, loop[t]) != edge.first &&
                     tet(pi, loop[t]) != edge.second)
                    point_set.insert(tet(pi, loop[t]));
                }
            }
        }

      points.resize(point_set.size());
      copy(point_set.begin(), point_set.end(), points.begin());

      return 0;
    }

    bool one_ring_tet_at_edge::is_chain_in_right_hand_order(
        const std::vector<size_t> & loop,
        const std::pair<size_t,size_t> & edge,
        const matrixst & tet,
        const matrixd & node)const
    {
      // beginning check
      vector<size_t> outside_tet_idx;
      for(size_t t = 0; t < loop.size(); ++t){
          if(loop[t] == -1) outside_tet_idx.push_back(t);
        }
      if(!outside_tet_idx.empty()){
          if(outside_tet_idx.size() > 2) {
              cerr << "# [error] wrong tet loop, more than three outside tet." << endl;
              return false;
            }
          for(size_t i = 0; i < outside_tet_idx.size(); ++i){
              if(outside_tet_idx[i] != 0 && outside_tet_idx[i] != loop.size() - 1)
                {
                  cerr << "# [error] outside tet not located on loop beginning or end." << endl;
                  return false;
                }
            }
        }

      vector<vector<size_t> > tet_vertex(2);
      vector<size_t> first_two_tets;
      for(size_t t = 0; t < loop.size(); ++t){
          if(loop[t]  != -1) {
              first_two_tets.push_back(loop[t]);
            }
          if(first_two_tets.size() == 2)
            break;
        }
      assert(first_two_tets.size() == 2);

      for(size_t t = 0; t < first_two_tets.size(); ++t) {
          for(size_t i = 0; i < tet.size(1); ++i)
            if(tet(i,first_two_tets[t]) != edge.first &&
               tet(i,first_two_tets[t]) != edge.second){
                tet_vertex[t].push_back(tet(i,first_two_tets[t]));
              }
          assert(tet_vertex[t].size() == 2);
        }
      size_t shared_vertex = -1;
      if(tet_vertex[0].front() != tet_vertex[1].front() &&
         tet_vertex[0].front() != tet_vertex[1].back()){
          assert(tet_vertex[0].back() == tet_vertex[1].front() ||
              tet_vertex[0].back() == tet_vertex[1].back());
          shared_vertex = tet_vertex[0].back();
        }else{
          shared_vertex = tet_vertex[0].front();
        }
      size_t begin_vertex = tet_vertex[0].front() + tet_vertex[0].back()
          - shared_vertex;
      const matrixd edge_v = node(colon(),edge.second) -
          node(colon(),edge.first);
      const matrixd round_v =
          cross((node(colon(),begin_vertex) - node(colon(),edge.first)),
                (node(colon(),shared_vertex) - node(colon(),begin_vertex)));
      assert(norm(edge_v) > 0);
      if(dot(edge_v,round_v) > 0)
        return true;
      else
        return false;
    }


    /////////////////////////////////////////////////////////////
    ////  hex

    face2hex_adjacent* face2hex_adjacent::create(const matrixi & hex)
    {
      unique_ptr<face2hex_adjacent> fa(new face2hex_adjacent);
      if(fa->init(hex))
        return 0;
      else
        return fa.release();
    }

    int face2hex_adjacent::init(const matrixi &hex)
    {
      typedef std::map<std::vector<size_t>,tuple<size_t,size_t,vector<size_t> > > map_type;
      map_type adj_map;
      assert(hex.size(1) == 8);
      //matrixi<size_t> face(4,6);
      const size_t face[] = {0,2,3,1,
                             7,6,4,5,
                             0,1,5,4,
                             7,3,2,6,
                             0,4,6,2,
                             7,5,1,3};
      itr_matrix<const size_t *> face0(4,6,&face[0]);

      vector<size_t> original_quad(4),sort_quad(4);

      for(size_t hex_idx = 0; hex_idx < hex.size(2); ++hex_idx){
          for(size_t quad_idx = 0; quad_idx < face0.size(2); ++quad_idx){
              for(size_t v = 0; v < face0.size(1); ++v)
                original_quad[v] = hex(face0(v,quad_idx),hex_idx);

              sort_quad = original_quad;
              sort(sort_quad.begin(),sort_quad.end());
              map_type::iterator i = adj_map.find(sort_quad);
              if(i == adj_map.end()){ // can not find such a quad face
                  adj_map[sort_quad] = make_tuple(hex_idx,-1,original_quad);
                }else{
                  if(get<0>(i->second) == -1){
                      // adj_map[sort_quad]
                      size_t & idx = get<0>(i->second);
                      idx = hex_idx;
                    }else if(get<1>(i->second) == -1){
                      size_t & idx = get<1>(i->second);
                      idx = hex_idx;
                    }else{
                      std::cerr << "# [error] this face meets at least three hex: " << get<0>(i->second)
                                << " " << get<1>(i->second) << " " << hex_idx << std::endl;
                      std::cerr << "# [error] face point: " ;
                      const vector<size_t> &common_face = i->first;

                      for(size_t t = 0; t < common_face.size(); ++t) std::cerr << common_face[t] << " ";
                      std::cerr << std::endl;

                      std::cerr << "# [error] hex_id: " << get<0>(i->second) << hex(colon(),get<0>(i->second));
                      std::cerr << "# [error] hex_id: " << get<1>(i->second) << hex(colon(),get<1>(i->second));
                      std::cerr << "# [error] hex_id: " << hex_idx << hex(colon(),hex_idx);

                      return __LINE__;
                    }

                }
            }
        }

      map_type::const_iterator cit = adj_map.begin();
      while(cit != adj_map.end()){
          faces_.push_back(get<2>(cit->second));
          faces_sort_.push_back(cit->first);
          face2hex_.push_back(make_pair(get<0>(cit->second),get<1>(cit->second)));
          ++cit;
        }

      return 0;
    }

    size_t face2hex_adjacent::get_face_idx(size_t a, size_t b,size_t c, size_t d) const{
      vector<size_t> quad(4);
      quad[0] = a; quad[1] = b; quad[2] = c; quad[3] = d;
      sort(quad.begin(),quad.end());
      std::vector<std::vector<size_t> >::const_iterator i
          = lower_bound(faces_sort_.begin(),faces_sort_.end(),quad);

      if(i == faces_sort_.end()) return -1; // can not find such a face
      else{
          const std::vector<size_t> &face = *i;
          for(size_t j = 0; j < face.size(); ++j) if(face[j] != quad[j]) return -1;
        }
      return static_cast<size_t> (i - faces_sort_.begin());
    }

    size_t face2hex_adjacent::get_face_idx(const size_t *abcd) const{
      assert(abcd);
      return get_face_idx(abcd[0],abcd[1],abcd[2],abcd[3]);
    }

    pair<size_t,size_t> face2hex_adjacent::query(size_t a, size_t b, size_t c, size_t d) const{
      size_t idx = get_face_idx(a,b,c,d);
      if(idx == -1) return make_pair(-1,-1);
      return face2hex_[idx];
    }

    pair<size_t,size_t> face2hex_adjacent::query(const size_t * abcd) const{
      assert(abcd);
      return query(abcd[0],abcd[1],abcd[2],abcd[3]);
    }

    size_t face2hex_adjacent::get_common_face_idx(const matrixst &hex,
                                                  const std::pair<size_t,size_t> &two_hex) const
    {
      if(two_hex.first == -1 || two_hex.second == -1) return -1;

      matrixst hex_a = hex(colon(),two_hex.first);
      matrixst hex_b = hex(colon(),two_hex.second);
      sort(hex_a.begin(),hex_a.end());
      sort(hex_b.begin(),hex_b.end());
      vector<size_t> int_(hex.size(1) * 2);
      vector<size_t>::const_iterator it;
      it = set_intersection(hex_a.begin(),hex_a.end(),hex_b.begin(),hex_b.end(),int_.begin());
      if(it == int_.begin() || it != int_.begin() + 4) return -1; // not a valid face
      return get_face_idx(int_[0],int_[1],int_[2],int_[3]);
    }

    int get_outside_face(const face2hex_adjacent &fa,
                         matrixst &face)
    {
      vector<size_t> face0;
      for(size_t t = 0; t < fa.faces_.size(); ++t){
          const pair<size_t,size_t> &two_hex = fa.face2hex_[t];
          if(!fa.is_outside_face(two_hex)) continue;
          const vector<size_t> & quad = fa.faces_[t];
          face0.insert(face0.end(),quad.begin(),quad.end());
        }

      face.resize(4,face0.size()/4);
      copy(face0.begin(),face0.end(),face.begin());
      return 0;
    }

    int get_outside_face_idx(const face2hex_adjacent &fa,
                             matrixst &face_idx)
    {
      vector<size_t> face0;
      for(size_t t = 0; t < fa.faces_.size(); ++t){
          const pair<size_t,size_t> &two_hex = fa.face2hex_[t];
          if(fa.is_outside_face(two_hex)) face0.push_back(t);
        }

      face_idx.resize(face0.size());
      copy(face0.begin(),face0.end(),face_idx.begin());
      return 0;
    }

    void one_ring_hex_at_edge::add_one_hex(const matrixst & hex,
                                           const matrixd &node,
                                           const face2hex_adjacent &fa)
    {
      assert(hex.size() == 8); // only one hex
      matrixst faces;
      get_faces_for_one_hex(hex,faces);

      typedef map<pair<size_t,size_t>,vector<size_t> >::const_iterator mcit;
      size_t edge[2];
      map<pair<size_t,size_t>,vector<size_t> > edge2face;
      for(size_t t = 0; t < faces.size(2); ++t){
          for(size_t i = 0; i < faces.size(1); ++i){
              edge[0] = faces(i,t);
              edge[1] = faces((i+1)%4,t);
              if(edge[0] > edge[1]) swap(edge[0],edge[1]);
              edge2face[make_pair(edge[0],edge[1])].push_back(t);
            }
        }

#if 1
      for(mcit it = edge2face.begin(); it != edge2face.end(); ++it){
          if(it->second.size() != 2){
              cerr << "# [error] degenerated hex, edge<" << it->first.first << "," << it->first.second << ">"
                   << " belongs to " << it->second.size() << "faces." << endl;
              return;
            }
        }
#endif

      for(mcit it = edge2face.begin(); it != edge2face.end(); ++it){
          const pair<size_t,size_t> & edge_ = it->first;
          for(size_t t = 0; t < it->second.size(); ++t){
              const pair<size_t,size_t> & both_hex = fa.face2hex_[fa.get_face_idx(&faces(0,it->second[t]))];
              e2h_[make_pair(edge_.first,edge_.second)].push_back(both_hex.first);
              e2h_[make_pair(edge_.first,edge_.second)].push_back(both_hex.second);
            }
        }
    }

    int  one_ring_hex_at_edge::sort_into_loop(const matrixst &hex,
                                              const matrixd &node,
                                              const face2hex_adjacent &fa)
    {
      size_t idx = 0;
      for(e2hex_type::iterator it = e2h_.begin(); it != e2h_.end(); ++it){
          ++idx;
          vector<size_t> &loop = it->second;
          set<pair<size_t,size_t> > segments;
          for(size_t t = 0; t < loop.size(); t+=2){
              if(loop[t] > loop[t+1])
                segments.insert(make_pair(loop[t+1],loop[t]));
              else
                segments.insert(make_pair(loop[t],loop[t+1]));
            }

          deque<size_t> new_loop;

          while(!segments.empty()){
              set<pair<size_t,size_t> >::iterator sit;
              for(sit = segments.begin(); sit != segments.end(); ++sit){
                  if(new_loop.empty()) {
                      new_loop.push_back(sit->first);
                      new_loop.push_back(sit->second);
                      break;
                    }
                  if(new_loop.back() != -1 && (new_loop.back() == sit->first || new_loop.back() == sit->second)){
                      new_loop.push_back(sit->first + sit->second - new_loop.back());
                      break;
                    }
                  else if (new_loop.front() != -1 && (new_loop.front() == sit->first || new_loop.front() == sit->second)) {
                      new_loop.push_front(sit->first + sit->second - new_loop.front());
                      break;
                    }
                } //end for

              if(sit == segments.end()){
                  cerr << "# idx = " << idx << endl;
                  cerr << "# strange error for edge: <" << it->first.first << "," << it->first.second << ">" << endl;
                  copy(loop.begin(), loop.end(), ostream_iterator<size_t>(cerr, " "));
                  cerr << endl;
                  return __LINE__;
                }
              segments.erase(sit);
            }// end while

          // a sharp edge, we add -1 to make it an open loop
          if(new_loop.size() == 2 && (new_loop[0] == -1)) new_loop.push_back(-1);
          if(new_loop.size() == 2 && (new_loop[1] == -1)) new_loop.push_front(-1);

          loop.resize(new_loop.size());
          copy(new_loop.begin(),new_loop.end(),loop.begin());
          if(loop.front() == -1 && loop.back() == -1 && loop.size() == 3) continue;

#if 1 // check to adjust the order of loop around edge into right-hand mode

          // adjust the loop order to make sure it fit right hand
          //      b1
          //    / | \
          //  b --|--b0
          //  |   |  |
          //  |   a1 |
          //  | /   \|
          //  a------a0   edge<a,b> == a-->b
          // to determin whether loop around edge<a,b> is right hand
          // we can check dot(edge<a,b>, cross(edge<a,a0>,edge<a0,a1>)) > 0

          if(loop.front() == -1){ // handle the open loop
              //      \-1/
              //     fa|fb\...
              matrixst  faces;
              const pair<size_t,size_t> & edge_ab = it->first;
              matrixst edge_a0b0(2),edge_a1b1(2);
              matrixst face_ab = fa.get_face(fa.get_common_face_idx(hex,make_pair(loop[1],loop[2])));

              assert(find(face_ab.begin(),face_ab.end(),edge_ab.first) != face_ab.end()
                  && find(face_ab.begin(),face_ab.end(),edge_ab.second) != face_ab.end());

              const matrixd ab = node(colon(),edge_ab.second) - node(colon(),edge_ab.first);

              { // handle edge<a1,b1>
                for(size_t i = 0, j = 0; i < 4; ++i){
                    if(face_ab[i] != edge_ab.first && face_ab[i] != edge_ab.second){
                        edge_a1b1[j] = face_ab[i];
                        ++j;
                        if(j == 2) break;
                      }
                  }
                matrixd a1b1 = node(colon(),edge_a1b1[1]) - node(colon(),edge_a1b1[0]);
                if(dot(a1b1,ab) < 0) swap(edge_a1b1[1],edge_a1b1[0]);
              }
              { // handle edge<a0,b0>
                vector<size_t> face_contain_edge;
                get_faces_for_one_hex(hex(colon(),loop[1]),faces);
                for(size_t t = 0; t < faces.size(2); ++t){
                    if(find(faces(colon(),t).begin(),faces(colon(),t).end(),edge_ab.first) != faces(colon(),t).end()
                       &&find(faces(colon(),t).begin(),faces(colon(),t).end(),edge_ab.second) != faces(colon(),t).end()){
                        face_contain_edge.push_back(t);
                        if(face_contain_edge.size() == 2) break;
                      }
                  }
                size_t edge_ab_[2] = {edge_ab.first,edge_ab.second};
                for(size_t t = 0; t < face_contain_edge.size(); ++t){
                    const pair<size_t,size_t> & hex_pair = fa.face2hex_[fa.get_face_idx(&faces(0,face_contain_edge[t]))];
                    if(fa.is_outside_face(hex_pair)){ // yes! this face is the one we need
                        vector<size_t> other_points(4);
                        vector<size_t>::iterator it = find_difference_set(faces(colon(),face_contain_edge[t]).begin(),
                                                                          faces(colon(),face_contain_edge[t]).end(),
                                                                          edge_ab_,edge_ab_ + 2,
                                                                          other_points.begin());
                        assert(it == other_points.begin() + 2);
                        copy(other_points.begin(),it,edge_a0b0.begin());
                        matrixd a0b0 = node(colon(),edge_a0b0[1]) - node(colon(),edge_a0b0[0]);
                        if(dot(a0b0,ab) < 0) swap(edge_a0b0[1],edge_a0b0[0]);
                        break;
                      }
                  }
              }
              { // to check the right hand order
                matrixd aa0 = node(colon(),edge_a0b0[0]) - node(colon(),edge_ab.first);
                matrixd a0a1 = node(colon(),edge_a1b1[0]) - node(colon(),edge_a0b0[0]);
                if(dot(cross(aa0,a0a1),ab) < 0){ // need to reverse the loop
                    reverse(loop.begin(),loop.end());
                  }
              }

            }else{
              //      \c /
              //    fa|fb \..
              assert(loop.size() > 3);
              matrixst two_face_idx(2);
              vector<size_t> other_points(4);
              const pair<size_t,size_t> & edge_ab = it->first;
              const matrixd ab = node(colon(),edge_ab.second) - node(colon(),edge_ab.first);
              // edges(colon(),0) == edge<a0,b0>
              // edges(colon(),1) == edge<a1,b1>
              matrixst edges(2,2);
              size_t edge_ab_[2] = {edge_ab.first,edge_ab.second};
              two_face_idx[0] = fa.get_common_face_idx(hex,make_pair(loop[0],loop[1]));
              two_face_idx[1] = fa.get_common_face_idx(hex,make_pair(loop[1],loop[2]));
              for(size_t t = 0; t < 2; ++t){
                  matrixst face = fa.get_face(two_face_idx[t]);
                  //sort(face.begin(),face.end())
                  vector<size_t>::iterator it = find_difference_set(face.begin(),face.end(),edge_ab_, edge_ab_ + 2, other_points.begin());
                  assert(it == other_points.begin() + 2);
                  copy(other_points.begin(),it,edges(colon(),t).begin());
                }
              for(size_t t = 0; t < 2; ++t){
                  matrixd edge_ = node(colon(),edges(1,t)) - node(colon(),edges(0,t));
                  if(dot(edge_,ab) < 0) swap(edges(1,t),edges(0,t));
                }
              { // to check the right hand order
                matrixd aa0 = node(colon(),edges(0,0)) - node(colon(),edge_ab.first);
                matrixd a0a1 = node(colon(),edges(0,1)) - node(colon(),edges(0,0));
                if(dot(cross(aa0,a0a1),ab) < 0){ // need to reverse the loop
                    reverse(loop.begin(),loop.end());
                  }
              }
            }
        }
#endif
      return 0;
    }

    //##        hex vertex order
    //##          2---0
    //##         /|  /|      w
    //##        3---1 |      | v
    //##        | 6-|-4      |/
    //##        |/  |/       0-->u
    //##        7---5
    // the faces are stored with pair
    int get_faces_for_one_hex(const matrixst &one_hex,
                              matrixst &faces)
    {
      if(one_hex.size() != 8) return __LINE__;

      faces.resize(4,6);

      const size_t four_faces[] = {0,2,3,1,
                                   7,6,4,5,

                                   0,1,5,4,
                                   7,3,2,6,

                                   0,4,6,2,
                                   7,5,1,3};
      for(size_t t = 0; t < 24; ++t){
          faces[t] = one_hex[four_faces[t]];
        }
      return 0;
    }

    int get_edges_for_one_hex(const matrixst &one_hex,
                              matrixst & edges)
    {
      if(one_hex.size() != 8) return __LINE__;
      edges.resize(2, 12);
      const size_t edge_raw [] = {
        0,1, 1,3, 3,2, 2,0,
        4,5, 5,7, 7,6, 6,4,
        0,4, 1,5, 2,6, 3,7
      };

      for(size_t pi = 0; pi < 2 * 12; ++pi)
        edges[pi] = one_hex[edge_raw[pi]];
      return 0;
    }

  }
}
