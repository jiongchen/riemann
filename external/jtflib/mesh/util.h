#ifndef JTFLIB_TRIMESH_UTIL_H
#define JTFLIB_TRIMESH_UTIL_H

#include <vector>
#include <stack>
#include <map>
#include "mesh.h"
#include <zjucad/matrix/matrix.h>

#define _EXPORTING

namespace jtf{
  namespace mesh{

    typedef zjucad::matrix::matrix<double> matrixd;
    typedef zjucad::matrix::matrix<size_t> matrixst;

    /// @brief calculate angle of one face (triangle/quad)
    /// @param mesh: input one face 3*1 or 4*1
    /// @param node: all nodes
    /// @param angle: output angle degree
	int JTF_MESH_API cal_face_angle(const matrixst & mesh,
                       const matrixd & node,
                       std::vector<double> & angle);

    ///
    /// @brief cal_min_angle_of_one_face
    /// @param mesh
    /// @param node
    /// @return minimal angle in degree
    ///
	double JTF_MESH_API cal_min_angle_of_one_face(const matrixst & mesh,
                                     const matrixd & node);

    ///
    /// @brief cal_min_angle_of_one_face
    /// @param mesh
    /// @param node
    /// @return minimal angle in degree, and which one is the minimal angle
    ///
	double JTF_MESH_API cal_min_angle_of_one_face(const matrixst & mesh,
                                     const matrixd & node,
                                     size_t &min_idx);

    ///
    /// @brief cal area of one face
    /// @param one_face
    /// @param node
    /// @return area
    ///
	double JTF_MESH_API cal_face_area(const matrixst & one_face,
                         const matrixd & node);

    ///
    /// @brief cal area of one face
    /// @param node of face point
    /// @return area
    ///
	double JTF_MESH_API cal_face_area(const matrixd &node);

    ///
    /// @brief cal_face_area
    /// @param one_face face array
    /// @param face_points point numbers
    /// @param node
    /// @return area
    ///
	double JTF_MESH_API cal_face_area(const size_t * one_face,
                         const size_t & face_points,
                         const matrixd & node);

    /// @brief cal triangle are by three edge length
	inline double JTF_MESH_API cal_face_area(double a, double b, double c)
    {
      assert(a > 0 && b > 0 & c > 0);
      double p = (a+b+c)/2.0;
      return std::sqrt(p*(p-a)*(p-b)*(p-c));
    }
    ///
    /// @brief cal_face_normal
    /// @param mesh
    /// @param node
    /// @param normal output normal
    /// @param is_normalized
    ///
	void JTF_MESH_API cal_face_normal(const matrixst & mesh,
                         const matrixd & node,
                         matrixd & normal,
                         bool is_normalized = true,
                         double eps = 1e-8);

    ///
    /// @brief cal_face_normal
    /// @param node node of one face
    /// @param normal output normal
    /// @param is_normalized
    ///
	void JTF_MESH_API cal_face_normal(const matrixd & node,
                         matrixd & normal,
                         bool is_normalized = true,
                         double eps = 1e-8 );

    ///
    /// @brief cal_point_normal, each point normal is area weight face normal
    ///        around one point
    /// @param mesh
    /// @param node
    /// @param normal
    ///
	void JTF_MESH_API cal_point_normal(const matrixst & mesh,
                          const matrixd & node,
                          matrixd & normal,
                          double eps = 1e-8);

    ///
    /// @brief transfer a planar 3d mesh to 2d mesh, notice that the mesh should
    ///        be planar, or the result can not be revert anymore.
    /// @param mesh input mesh
    /// @param node_3d input mesh node
    /// @param axes output axes of 2d plane, and the original point
    /// @param node_2d output project of each node to 2d plane
    ///
	void JTF_MESH_API mesh_to_2d(const matrixst & mesh,
                    const matrixd & node_3d,
                    matrixd & axes,
                    matrixd & node_2d);

    ///
    /// @brief mesh_from_2d revert a planar mesh from 2d coordniates
    /// @param mesh
    /// @param node_2d
    /// @param axes
    /// @param node_3d
    ///
	void JTF_MESH_API mesh_from_2d(const matrixst & mesh,
                      const matrixd & node_2d,
                      const matrixd & axes,
                      matrixd & node_3d);

    ///
    /// @brief cal_average_edge
    /// @param mesh
    /// @param node
    /// @return average edge length
    ///
	double JTF_MESH_API cal_average_edge(const matrixst & mesh,
                            const matrixd &node);


	class JTF_MESH_API curvature
    {
    public:
      curvature(const matrixst & mesh, const matrixd & node);
      void generate();
      const matrixd & get_gauss_curvature()const{
        return kG_;
      }
      // return kmax_kin: 2*N
      const matrixd & get_principle_curvature()const{
        return k1k2_;
      }
      const matrixd& get_mean_curvature()const{
        return kH_;
      }
      const matrixd& get_principle_curvature_direction()const{
        return d1d2_;
      }
    private:
      void kG();
      void kH();
      void k1k2();
      void d1d2();
      double cal_xi_angle(const size_t xi, const matrixst & one_face,
                          const matrixd & node) const;
      void get_common_edge(const matrixst & f0, const matrixst & f1,
                           size_t * edge)const;
      double get_beta_angle(const zjucad::matrix::matrix<double> & d0,
                            const zjucad::matrix::matrix<double> & d1,
                            const zjucad::matrix::matrix<double> & axis)const;
      void sort_eig(zjucad::matrix::matrix<double> &C, zjucad::matrix::matrix<double> &EV,
                    const zjucad::matrix::matrix<double> &n) const;
    private:
      const matrixst & mesh_;
      const matrixd & node_;
      matrixd Amixed_;
      matrixd angle_defect_;
      matrixd face_area_;
      matrixd face_normal_;
      matrixd kG_,kH_;
      matrixd k1k2_;
      matrixd d1d2_;
      jtf::mesh::one_ring_face_at_point orfap_;
      std::shared_ptr<jtf::mesh::edge2cell_adjacent> ea_;
    };

    ///
    /// @brief This function is used to extract mesh feature line
    /// @param face input mesh
    /// @param node input mesh node
    /// @param feature_line output feature line, each colon is an edge
    /// @param cos_angle_threshold default is 0.83 (34 degree)
    /// @return return 0 if works fine, or return non-zeros
    ///
	int JTF_MESH_API extract_mesh_feature_line(
        const matrixst & face,
        const matrixd & node,
        matrixst & feature_line,
        const double cos_angle_threshold = 0.83);

    ///
    /// \brief extract_hex_singularity_lines
    /// \param ortae    input one_ring_hex_at_edge data
    /// \param singularity_edges ouptut singularity edges
    /// \return
    ///
	int JTF_MESH_API extract_hex_singularity_lines(
        const one_ring_hex_at_edge &ortae,
        std::vector<std::pair<size_t,size_t> > & singularity_edges);

    ///
    /// @brief return cal_tet_vol
    ///
    template <typename E>
    inline typename E::value_type cal_tet_vol(const zjucad::matrix::matrix_expression<E> &tet)
    {
      using namespace zjucad::matrix;
      assert(tet().size(1) == 3);
      if(tet().size(2) == 4) {
          matrix<typename E::value_type> edges
              = tet()(colon(), colon(1, 3)) - tet()(colon(), 0)*
              ones<typename E::value_type>(1, 3);
          return cal_tet_vol(edges);
        }
      if(tet().size(2) == 3) // fix by jtf, here tet means three edges (3*3)
        return dot(cross(tet()(colon(), 0), tet()(colon(), 1)), tet()(colon(), 2))/6.0;
    }

    template <typename val_type>
    bool is_acute_triange(const zjucad::matrix::matrix<size_t> &one_tri,
                          const zjucad::matrix::matrix<val_type> & node)
    {
      assert(one_tri.size() == 3);
      std::vector<val_type> edge_len = {
        zjucad::matrix::norm(node(zjucad::matrix::colon(), one_tri[0]) - node(zjucad::matrix::colon(), one_tri[1])),
        zjucad::matrix::norm(node(zjucad::matrix::colon(), one_tri[1]) - node(zjucad::matrix::colon(), one_tri[2])),
        zjucad::matrix::norm(node(zjucad::matrix::colon(), one_tri[2]) - node(zjucad::matrix::colon(), one_tri[0]))
      };
      for(auto & len : edge_len) len = len*len;
      std::sort(edge_len.begin(), edge_len.end());
      if(edge_len[0]+edge_len[1]>edge_len[2]) return true;
      return false;
    }

    template <typename T1, typename T2>
    double get_triangle_circumcenter(const zjucad::matrix::matrix_expression<T1> &tri_node,
                                     zjucad::matrix::matrix_expression<T2> &center)
    {
      // using namespace zjucad::matrix;
      // zjucad::matrix::matrix<double> M(3,2),MT,MTM;
      // M(colon(),0) = tri_node()(colon(),1) - tri_node()(colon(),0);
      // M(colon(),1) = tri_node()(colon(),2) - tri_node()(colon(),0);
      // MT = trans(M);
      // MTM = MT*M;
      // zjucad::matrix::matrix<double> b = MT*tri_node()(colon(),0);
      // b[0] -= 0.5*std::pow(norm(tri_node()(colon(),1)),2) - 0.5*std::pow(norm(tri_node()(colon(),0)),2);
      // b[1] -= 0.5*std::pow(norm(tri_node()(colon(),2)),2) - 0.5*std::pow(norm(tri_node()(colon(),0)),2);
      // zjucad::matrix::matrix<double> inv_MTM = MTM;
      // jtf::math::invert_2dmatrix(inv_MTM);
      // zjucad::matrix::matrix<double> lambda = inv_MTM * b;
      // center() = M*lambda+tri_node()(colon(),0);
      // return norm(tri_node()(colon(),0) - center());
      return 0;
    }


    //     class JTF_MESH_API patch_separater
    // {
    // public:
    //   patch_separater(const zjucad::matrix::matrix<size_t> & faces):faces_(faces){
    //     ea_.reset(jtf::mesh::edge2cell_adjacent::create(faces_));
    //   }
    //   patch_separater(const zjucad::matrix::matrix<size_t> &faces,
    //                   std::shared_ptr<const jtf::mesh::edge2cell_adjacent> ea)
    //     :faces_(faces), ea_(ea){}

    //   void separater(const zjucad::matrix::matrix<size_t> & face_type,
    //                  std::vector<std::deque<std::pair<size_t,size_t> > > & chains){
    //     if(face_type.size() != faces_.size(2))
    //       throw std::invalid_argument("wrong face number.");

    //     std::vector<std::pair<size_t,size_t> > patch_boundary_edges;
    //     for(size_t ei = 0; ei < ea_->edge2cell_.size(); ++ei){
    //         const std::pair<size_t,size_t> & tri_pair = ea_->edge2cell_[ei];
    //         if(ea_->is_boundary_edge(tri_pair)) {
    //             patch_boundary_edges.push_back(ea_->edges_[ei]);
    //             continue;
    //           }
    //         if(face_type[tri_pair.first] != face_type[tri_pair.second]){
    //             patch_boundary_edges.push_back(ea_->edges_[ei]);
    //           }
    //       }

    //     jtf::util::extract_chain_from_edges(patch_boundary_edges, chains);
    //     separater(chains);
    //   }

    //   void separater(const std::vector<std::pair<size_t,size_t> > & edges){
    //     patches_.clear();
    //     edge2patch_idx_.clear();
    //     face2patches_.clear();

    //     std::vector<bool> is_face_visited(faces_.size(2), false);
    //     std::stack<size_t> face_stack;

    //     std::set<std::pair<size_t,size_t> > edge_set;
    //     for(const auto & one_edge : edges){
    //         if(one_edge.first > one_edge.second) edge_set.insert(std::make_pair(one_edge.second, one_edge.first));
    //         else edge_set.insert(one_edge);
    //       }

    //     auto cit = find(is_face_visited.begin(), is_face_visited.end(), false);
    //     while(cit != is_face_visited.end()){
    //         std::set<size_t> one_patch;
    //         face_stack.push(cit-is_face_visited.begin());
    //         while(!face_stack.empty()){
    //             const size_t f_idx = face_stack.top();
    //             face_stack.pop();
    //             is_face_visited[f_idx] = true;
    //             one_patch.insert(f_idx);

    //             for(size_t pi = 0; pi < faces_.size(1); ++pi){
    //                 std::pair<size_t,size_t> one_edge(faces_(pi,f_idx),faces_((pi+1)%faces_.size(1),f_idx));
    //                 if(one_edge.first > one_edge.second) std::swap(one_edge.first, one_edge.second);
    //                 if(edge_set.find(one_edge) != edge_set.end()) continue;

    //                 const size_t edge_idx = ea_->get_edge_idx(one_edge.first, one_edge.second);
    //                 if(edge_idx == -1){
    //                     std::cerr << "# [error] strange can not find edge idx of "
    //                               << one_edge.first << " " << one_edge.second << std::endl;
    //                     return ;
    //                   }
    //                 const std::pair<size_t,size_t> & tri_pair = ea_->edge2cell_[edge_idx];
    //                 if(ea_->is_boundary_edge(tri_pair)) continue;
    //                 if(f_idx != tri_pair.first && f_idx != tri_pair.second){
    //                     std::cerr << "# [error] strange face pair of edge "
    //                               << one_edge.first << " " << one_edge.second << " is "
    //                               << tri_pair.first << " " << tri_pair.second << " without "
    //                               << f_idx << std::endl;
    //                     return ;
    //                   }
    //                 const size_t other_face_idx = tri_pair.first + tri_pair.second - f_idx;
    //                 if(is_face_visited[other_face_idx]) continue;
    //                 face_stack.push(other_face_idx);
    //               }
    //           }
    //         std::vector<size_t> one_patch_vec(one_patch.size());
    //         std::copy(one_patch.begin(), one_patch.end(), one_patch_vec.begin());
    //         patches_.push_back(one_patch_vec);
    //         cit = std::find(is_face_visited.begin(), is_face_visited.end(), false);
    //       }
    //     face2patches_.resize(faces_.size(2));
    //     for(size_t pi = 0; pi < patches_.size(); ++pi){
    //         const std::vector<size_t> & one_patch = patches_[pi];
    //         for(const auto & one_face : one_patch){
    //             face2patches_[one_face] = pi;
    //           }
    //       }
    //     for(size_t ei = 0; ei < ea_->edge2cell_.size(); ++ei){
    //         std::pair<size_t,size_t> one_edge = ea_->edges_[ei];
    //         if(one_edge.first > one_edge.second) std::swap(one_edge.first, one_edge.second);
    //         const std::pair<size_t,size_t> & one_face_pair = ea_->edge2cell_[ei];
    //         if(one_face_pair.first != -1){
    //             edge2patch_idx_[one_edge].insert(face2patches_[one_face_pair.first]);
    //             patch2edge_[face2patches_[one_face_pair.first]].push_back(one_edge);
    //           }
    //         if(one_face_pair.second != -1){
    //             edge2patch_idx_[one_edge].insert(face2patches_[one_face_pair.second]);
    //             patch2edge_[face2patches_[one_face_pair.second]].push_back(one_edge);
    //           }
    //       }
    //   }

    //   void separater(const std::vector<std::deque<std::pair<std::size_t,std::size_t> > > &chain){
    //     chain2patch_idx_.clear();
    //     std::vector<std::pair<size_t,size_t> > edges;
    //     for(size_t ci = 0; ci < chain.size(); ++ci){
    //         const std::deque<std::pair<size_t,size_t> > & one_chain = chain[ci];
    //         for(const auto & one_edge : one_chain) edges.push_back(one_edge);
    //       }
    //     separater(edges);
    //     for(size_t ci = 0; ci < chain.size(); ++ci){
    //         const std::deque<std::pair<size_t,size_t> > & one_chain = chain[ci];
    //         for(const auto & one_edge : one_chain) {
    //             std::set<size_t> edge_adj_patch = get_edge_adj_patches(one_edge);
    //             chain2patch_idx_[ci].insert(edge_adj_patch.begin(), edge_adj_patch.end());
    //             for(const auto & one_p : edge_adj_patch)
    //               patch2chain_[one_p].insert(ci);
    //           }
    //       }
    //   }

    //   void separater(){
    //     zjucad::matrix::matrix<size_t> boundary_edges;
    //     get_boundary_edge(*ea_, boundary_edges);
    //     std::vector<std::pair<size_t,size_t> > boundary_edges_vec;
    //     for(size_t ei = 0 ; ei < boundary_edges.size(2); ++ei){
    //         boundary_edges_vec.push_back(std::make_pair(boundary_edges(0,ei), boundary_edges(1,ei)));
    //       }
    //     separater(boundary_edges_vec);
    //   }
    //   const std::vector<std::vector<size_t> > & get_all_patches() const{return patches_;}
    //   const std::vector<size_t> & get_patch(const size_t i ) const {return patches_.at(i);}
    //   const size_t get_face2patch(const size_t i) const {return face2patches_.at(i);}
    //   const std::vector<size_t> & get_all_face2patch()const{return face2patches_;}
    //   const std::set<size_t> get_edge_adj_patches(const std::pair<size_t,size_t> & edge) const{
    //     std::map<std::pair<size_t,size_t>, std::set<size_t> >::const_iterator it;
    //     if(edge.first > edge.second)
    //       it = edge2patch_idx_.find(std::make_pair(edge.second, edge.first));
    //     else it = edge2patch_idx_.find(edge);
    //     if(it == edge2patch_idx_.end()) return std::set<size_t>();
    //     return it->second;
    //   }
    //   const std::set<size_t>& get_chain_adj_patches(const size_t ci) const{
    //     return chain2patch_idx_.at(ci);
    //   }
    //   const std::vector<std::pair<size_t,size_t> >& get_edges_of_patch(const size_t pi)const{
    //     return patch2edge_.at(pi);
    //   }
    //   const std::set<size_t>& get_chain_of_patch(const size_t pi) const{
    //     return patch2chain_.at(pi);
    //   }
    // private:
    //   const zjucad::matrix::matrix<size_t> &faces_;
    //   std::shared_ptr<const jtf::mesh::edge2cell_adjacent> ea_;
    //   std::vector<std::vector<size_t> > patches_;
    //   std::vector<size_t> face2patches_;
    //   std::map<std::pair<size_t,size_t>, std::set<size_t> > edge2patch_idx_;
    //   std::map<size_t, std::vector<std::pair<size_t,size_t> > > patch2edge_;
    //   std::map<size_t, std::set<size_t> > chain2patch_idx_;
    //   std::map<size_t, std::set<size_t> > patch2chain_;
    // };

  }
}

#endif // UTIL_H
