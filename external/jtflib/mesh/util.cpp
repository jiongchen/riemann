#include <iostream>
#include <boost/math/constants/constants.hpp>
//#include <hjlib/math/blas_lapack.h>
//#include <zjucad/matrix/lapack.h>
#include <zjucad/matrix/io.h>
#include "util.h"

using namespace std;
using namespace zjucad::matrix;

namespace jtf{
  namespace mesh{

    template<typename T>
    T cos2deg(const T & cos_)
    {
      return acos(std::min(1.0, std::max(cos_, -1.0))) * 180.0
          / boost::math::constants::pi<double>();;
    }

    int cal_face_angle(const matrixst & mesh,
                       const matrixd & node,
                       std::vector<double> & angle)
    {
      assert(mesh.size() == 3 || mesh.size() == 4); // one face

      angle.resize(mesh.size());

      matrixd two_edge(node.size(1),2);

      for(size_t ai = 0; ai < mesh.size(); ++ai){
          for(size_t ei = 0; ei < 2; ++ei){
              two_edge(colon(), ei) =
                  node(colon(), mesh[(ai + mesh.size() + (ei==0?-1:1)) % mesh.size()])
                  - node(colon(), mesh[ai]);
              const double len = norm(two_edge(colon(), ei));
              if( len > 1e-6)
                two_edge(colon(), ei) /= len;
            }
          const double cos_ = dot(two_edge(colon(),0),  two_edge(colon(),1));
          angle[ai] = cos2deg(cos_);
        }

      return 0;
    }


    double cal_face_area(const size_t * one_face,
                         const size_t & face_points,
                         const matrixd & node)
    {
      if(face_points == 3){
          matrixd edges[2] = {
            node(colon(), one_face[1])-node(colon(), one_face[0]),
            node(colon(), one_face[2])-node(colon(), one_face[0])
          };
          return norm(cross(edges[0], edges[1]))/2.0;
        }else if(face_points == 4){
          matrixst face = zeros<size_t>(3,2);
          face(0,0) = one_face[0];
          face(1,0) = one_face[1];
          face(2,0) = one_face[2];

          face(0,1) = one_face[2];
          face(1,1) = one_face[3];
          face(2,1) = one_face[0];

          return cal_face_area(&face(0,0), 3, node) +
              cal_face_area(&face(0,1), 3, node);
        }else {
          std::cerr << "# [error] unsupported face type." << endl;
          return 0;
        }
    }

    double cal_face_area(const matrixst & one_face,
                         const matrixd & node)
    {
      return cal_face_area(&one_face[0], one_face.size(), node);
    }

    double cal_face_area(const matrixd &node)
    {
      static vector<size_t> face(node.size(2));
      for(size_t i = 0; i < face.size(); ++i) face[i] = i;

      return cal_face_area(&face[0], face.size(), node);
    }

    void cal_face_normal(const matrixst & mesh,
                         const matrixd & node,
                         matrixd & normal,
                         bool is_normalized,
                         double eps)
    {
      if(normal.size(1) != 3 || normal.size(2) != mesh.size(2));
      normal.resize(3, mesh.size(2));

      matrixd edges[2] ;

      for(size_t fi = 0; fi < mesh.size(2); ++fi){
          edges[0] = node(colon(), mesh(1,fi))-node(colon(), mesh(0,fi));
          edges[1] = node(colon(), mesh(2,fi))-node(colon(), mesh(0,fi));

          normal(colon(),fi) = cross(edges[0], edges[1]);
          if(is_normalized){
              const double len = norm(normal(colon(),fi));
              if(len > eps)
                normal(colon(),fi) /= len;
            }
        }
    }

    void cal_face_normal(const matrixd & node,
                         matrixd & normal,
                         bool is_normalized,
                         double eps)
    {
      if(normal.size() != 3)
        normal.resize(3,1);

      matrixd edges[2] ;
      edges[0] = node(colon(), 1)-node(colon(), 0);
      edges[1] = node(colon(), 2)-node(colon(), 0);

      normal = cross(edges[0], edges[1]);
      if(is_normalized){
          const double len = norm(normal);
          if(len > eps)
            normal /= len;
        }
    }

    void cal_point_normal(const matrixst & mesh,
                          const matrixd & node,
                          matrixd & normal,
                          double eps)
    {
      normal = zeros<double>(3, node.size(2));

      matrixd face_normal, face_area(mesh.size(2),1);
      cal_face_normal(mesh, node, face_normal);
      for(size_t fi = 0; fi < mesh.size(2); ++fi){
          face_area[fi] = cal_face_area(mesh(colon(),fi), node);
        }

      vector<vector<size_t> > one_ring_face_of_p(node.size(2));
      for(size_t fi = 0; fi < mesh.size(2); ++fi){
          for(size_t di = 0; di < mesh.size(1); ++di){
              one_ring_face_of_p[mesh(di,fi)].push_back(fi);
            }
        }

      for(size_t ni = 0; ni < node.size(2); ++ni){
          double total_area = 0;
          if(one_ring_face_of_p[ni].empty()) continue;
          for(size_t ai = 0; ai < one_ring_face_of_p[ni].size(); ++ai){
              normal(colon(),ni) += face_normal(colon(), one_ring_face_of_p[ni][ai])
                  * face_area[one_ring_face_of_p[ni][ai]];
              total_area += face_area[one_ring_face_of_p[ni][ai]];
            }

          normal(colon(),ni) /= total_area;
          const double len = norm(normal(colon(),ni));
          if(len > eps)
            normal(colon(),ni) /= len;
        }
    }

    void mesh_to_2d(const matrixst & mesh,
                    const matrixd & node_3d,
                    matrixd & axes_original, // axes with original_point
                    matrixd & node_2d)
    {
      node_2d = zeros<double>(2, node_3d.size(2));
      axes_original = zeros<double>(3,3);
      axes_original(colon(),2) = node_3d(colon(),0);
      matrixd edge[] = {
        node_3d(colon(), mesh[1]) - node_3d(colon(),mesh[0]),
        node_3d(colon(),mesh[2]) - node_3d(colon(), mesh[0])
      };

      const double len[] = { norm(edge[0]),   norm(edge[1])};
      for(size_t ei = 0; ei < 2; ++ei){
          if(len[ei] > 1e-6)
            edge[ei] /= len[ei];
        }

      axes_original(colon(),0) = edge[0];
      matrixd normal = cross(edge[0], edge[1]);

      axes_original(colon(),1) = cross(axes_original(colon(),0), normal);
      const double len_a = norm(axes_original(colon(),1));
      if(len_a > 1e-6)
        axes_original(colon(),1) /= len_a;

      for(size_t ni = 0; ni < node_3d.size(2); ++ni){
          node_2d(0,ni) = dot(node_3d(colon(),ni) - axes_original(colon(),2),
                              axes_original(colon(),0));
          node_2d(1,ni) = dot(node_3d(colon(),ni) - axes_original(colon(),2),
                              axes_original(colon(),1));
        }
    }

    void mesh_from_2d(const matrixst & mesh,
                      const matrixd & node_2d,
                      const matrixd & axes_original,
                      matrixd & node_3d)
    {
      node_3d = zeros<double>(3, node_2d.size(2));
      for(size_t ni = 0; ni < node_2d.size(2); ++ni){
          node_3d(colon(), ni) =
              axes_original(colon(),2) +
              node_2d(0,ni) * axes_original(colon(),0) +
              node_2d(1,ni) * axes_original(colon(),1);
        }
    }


    int extract_mesh_feature_line(
        const matrixst &face,
        const matrixd &node,
        matrixst &feature_line,
        const double cosin_theta)
    {
      boost::shared_ptr<edge2cell_adjacent> edge_mesh(
		  edge2cell_adjacent::create(face));
      if(edge_mesh == NULL)
        return __LINE__;
      size_t face1_id,face2_id;
      matrixd face1_normal, face2_normal;
      vector<size_t> feature_vec;
      size_t feature_line_num;
      for(size_t i = 0; i < edge_mesh -> edges_.size(); ++i) {
          if(! edge_mesh -> is_boundary_edge(edge_mesh -> edge2cell_[i])) {
              face1_id = edge_mesh -> edge2cell_[i].first;
              face2_id = edge_mesh -> edge2cell_[i].second;
              face1_normal = cross( (node(colon(), face(0,face1_id)) -
                                     node(colon(), face(1,face1_id))),
                                    (node(colon(), face(1,face1_id)) -
                                     node(colon(), face(2,face1_id))) );

              face2_normal = cross( (node(colon(), face(0,face2_id)) -
                                     node(colon(), face(1,face2_id))),
                                    (node(colon(), face(1,face2_id)) -
                                     node(colon(), face(2,face2_id))) );
              if(norm(face1_normal) < 1e-8 || norm(face2_normal) < 1e-8) {
                  cout << "#[info] tri: the length of the feature line is zero" <<endl;
                  continue;
                  // return 1;
                }
              face1_normal = face1_normal / norm(face1_normal);
              face2_normal = face2_normal / norm(face2_normal);
              if(fabs( dot(face1_normal, face2_normal)) <= cosin_theta) {
                  feature_vec.push_back(edge_mesh -> edges_[i].first);
                  feature_vec.push_back(edge_mesh -> edges_[i].second);
                }
            }
        }
      feature_line_num = feature_vec.size() / 2;
      feature_line.resize(2,feature_line_num);
      for(size_t j = 0; j < feature_vec.size(); ++j)
        feature_line[j] = feature_vec[j];
      //cerr << "# [error] empty function." << endl;
      return 0;
    }

    double cal_average_edge(const matrixst & mesh,
                            const matrixd &node)
    {
      assert(mesh.size());
      double avg = 0;
      for(size_t ci = 0; ci < mesh.size(2); ++ci){
          for(size_t pi = 0; pi < mesh.size(1); ++pi){
              avg += norm(node(colon(), mesh(pi,ci)) - node(colon(), mesh((pi+1)%mesh.size(1), ci )));
            }
        }
      avg /= mesh.size();
      return avg;
    }

    double cal_min_angle_of_one_face(const matrixst & mesh,
                                     const matrixd & node) {
      vector<double> angle;
      cal_face_angle(mesh, node, angle);
      return *min_element(angle.begin(), angle.end());
    }

    double cal_min_angle_of_one_face(const matrixst & mesh,
                                     const matrixd & node,
                                     size_t &min_idx){
      vector<double> angle;
      cal_face_angle(mesh, node, angle);
      min_idx = min_element(angle.begin(), angle.end()) - angle.begin();
      return angle[min_idx];
    }

    int extract_hex_singularity_lines(
        const one_ring_hex_at_edge &ortae,
        std::vector<std::pair<size_t,size_t> > & singularity_edges)
    {
      for(one_ring_hex_at_edge::e2hex_type::const_iterator
          eit = ortae.e2h_.begin(); eit != ortae.e2h_.end(); ++eit){
          const vector<size_t> & loop = eit->second;
          if(find(loop.begin(), loop.end(),-1) != loop.end()) continue;
          const size_t adj_hex_num = loop.size() - count(loop.begin(), loop.end(), -1);
          if(adj_hex_num != 5){
              singularity_edges.push_back(eit->first);
            }
        }
	  return 0;
    }

#if 0
    curvature::curvature(const matrixst & mesh, const matrixd & node)
      :mesh_(mesh), node_(node){
      ea_.reset(jtf::mesh::edge2cell_adjacent::create(mesh));
      if(!ea_.get()) throw std::invalid_argument("wrong mesh.");
      face_area_.resize(mesh_.size(2),1);
      for(size_t fi = 0; fi < mesh_.size(2); ++fi)
        face_area_[fi] = jtf::mesh::cal_face_area(mesh_(zjucad::matrix::colon(), fi), node_);

      jtf::mesh::cal_face_normal(mesh_, node_, face_normal_);

      orfap_.add_all_faces(mesh_, *ea_);
      orfap_.sort_int_loop_with_normal_info(mesh_, node_, *ea_, face_normal_);
    }

    void curvature::generate(){
      kG();
      kH();
      k1k2();
      d1d2();
    }

    // void curvature::kG()
    // {
      // Amixed_.resize(node_.size(2),1);
      // Amixed_ *= 0;
      // angle_defect_ = ones<double>(node_.size(2),1)*2*jtf::math::My_PI();
      // kG_.resize(node_.size(2),1);
      // for(const auto & one_p : orfap_.p2f_){
      //     const vector<size_t> & one_ring = one_p.second;
      //     if(one_ring.front() == -1 || one_ring.back() == -1){
      //         Amixed_[one_p.first]= 1; continue;
      //       }

      //     for(size_t fi = 0; fi < one_ring.size()-1; ++fi){
      //         double angle = cal_xi_angle(one_p.first, mesh_(colon(), one_ring[fi]), node_);
      //         angle_defect_[one_p.first] -= angle;
      //         Amixed_[one_p.first] += face_area_[one_ring[fi]]/3.0;
      //       }
      //     kG_[one_p.first] = angle_defect_[one_p.first] /Amixed_[one_p.first];
      //   }
  //    }

    // template <typename val_type>
    // val_type cot(const val_type &v){
    //   return tan(jtf::math::My_PI()*0.5-v);
    // }

    // void curvature::kH()
    // {
    //   kH_.resize(node_.size(2),1);
    //   kH_ *= 0;
    //   matrixst edge(2,1);
    //   for(const auto & one_p : orfap_.p2f_){
    //       const vector<size_t> & one_ring = one_p.second;
    //       if(one_ring.front() == -1 || one_ring.back() == -1) continue;

    //       for(size_t fi = 0; fi < one_ring.size()-1; ++fi){
    //           get_common_edge(mesh_(colon(), one_ring[fi]),
    //                           mesh_(colon(), one_ring[fi+1]),
    //               &edge[0]);
    //           const size_t other_p = edge[0] + edge[1] - one_p.first;
    //           const double len = norm(node_(colon(), other_p) - node_(colon(), one_p.first));
    //           double normal_diff_angle = jtf::math::safe_acos(
    //                 dot(face_normal_(colon(), one_ring[fi]), face_normal_(colon(), one_ring[fi+1])));
    //           if(dot(cross(face_normal_(colon(), one_ring[fi]), face_normal_(colon(), one_ring[fi+1])),
    //                  node_(colon(), one_p.first) - node_(colon(), other_p)) < 0)
    //             normal_diff_angle *= -1;
    //           kH_[one_p.first] += normal_diff_angle*len*0.25;
    //         }
    //       kH_[one_p.first] /= Amixed_[one_p.first];
    //     }
    // }

    // void curvature::k1k2()
    // {
    //   k1k2_.resize(2, node_.size(2));
    //   for(size_t i = 0; i < node_.size(2); ++i){
    //       double delta = kH_[i]*kH_[i]-kG_[i];
    //       if(delta < 0) delta = 0;
    //       k1k2_(0,i) = kH_[i] + sqrt(delta);
    //       k1k2_(1,i) = kH_[i] - sqrt(delta);
    //     }
    // }
    // double curvature::get_beta_angle(const zjucad::matrix::matrix<double> & d0,
    //                                  const zjucad::matrix::matrix<double> & d1,
    //                                  const zjucad::matrix::matrix<double> &axis)const
    // {
    //   matrix<double> d0_unit = d0/norm(d0);
    //   matrix<double> d1_unit = d1/norm(d1);

    //   double angle = jtf::math::safe_acos(dot(d0_unit,d1_unit));

    //   if(dot(cross(d0_unit,d1_unit),axis) < 0) angle*=-1;
    //   return angle;
    // }

    // void curvature::d1d2()
    // {
    //   cerr << "# [warning] principal direction is experimental." << endl;
    //   d1d2_.resize(6,node_.size(2));
    //   matrix<double> C(3,3);
    //   matrixst common_edge(2,1);
    //   matrixd e(3,1);
    //   matrixd EV(3,1);
    //   matrixd point_normal;
    //   jtf::mesh::cal_point_normal(mesh_, node_, point_normal);
    //   for(const auto & one_p : orfap_.p2f_){
    //       const vector<size_t> & one_ring_f = one_p.second;
    //       if(one_ring_f.front() == -1 || one_ring_f.back() == -1) continue;
    //       C *= 0;
    //       for(size_t i = 0; i != one_ring_f.size()-1; ++i){
    //           get_common_edge(mesh_(colon(),one_ring_f[i]),
    //                           mesh_(colon(),one_ring_f[i+1]), &common_edge[0]);
    //           e = node_(colon(), one_p.first) - node_(colon(),common_edge[0]+common_edge[1]-one_p.first);
    //           const double len = norm(e);
    //           e /= len;
    //           const double beta = get_beta_angle(face_normal_(colon(), one_ring_f[i]),
    //                                              face_normal_(colon(), one_ring_f[i+1]),
    //               e);
    //           C += beta * 0.5*len*(e*trans(e));
    //         }
    //      C /= Amixed_[one_p.first];
    //       eig(C,EV);
    //       sort_eig(C,EV, point_normal(colon(), one_p.first));

    //       d1d2_(colon(0,2),one_p.first) = C(colon(),2);
    //       d1d2_(colon(3,5),one_p.first) = C(colon(),1);
    //     }
    // }
    void curvature::sort_eig(zjucad::matrix::matrix<double> &C, zjucad::matrix::matrix<double> &EV,
                             const zjucad::matrix::matrix<double> & n) const
    {
      vector<pair<double,size_t> > eig_fabs(3);
      for(size_t i = 0; i < 3; ++i) eig_fabs[i] = make_pair(fabs(EV[i]),i);
      sort(eig_fabs.begin(), eig_fabs.end());
      if(fabs(dot(C(colon(),eig_fabs[2].second),n)) > 0.7)
        swap(eig_fabs.front(), eig_fabs[2]);
      else if(fabs(dot(C(colon(),eig_fabs[1].second),n)) > 0.7)
        swap(eig_fabs.front(), eig_fabs[1]);
      matrix<double> old_C = C, old_EV = EV;
      for(size_t i = 0; i < 3; ++i){
          C(colon(), i) = old_C(colon(), eig_fabs[i].second);
          EV[i] = old_EV[eig_fabs[i].second];
        }
    }
    double curvature::cal_xi_angle(
        const size_t xi, const matrixst & one_face, const matrixd & node) const
    {
      auto it = std::find(one_face.begin(), one_face.end(), xi);
      if(it == one_face.end())
        throw std::logic_error("can not find xi in this face.");
      size_t idx = static_cast<size_t>(it-one_face.begin());
      size_t next = one_face[(idx+1)%one_face.size()];
      size_t prev = one_face[(idx+one_face.size()-1)%one_face.size()];
      matrixd dirs(3,2);
      dirs(colon(), 0) = node(colon(), next) - node(colon(),xi);
      dirs(colon(), 1) = node(colon(), prev) - node(colon(),xi);
      return jtf::math::safe_acos(dot(dirs(colon(),0), dirs(colon(),1))/(norm(dirs(colon(), 0)) * norm(dirs(colon(),1))));
    }

    void curvature::get_common_edge(const matrixst & f0, const matrixst & f1,
                                    size_t * edge)const
    {
      size_t k = 0;
      for(size_t i = 0; i < f0.size(); ++i)
        for(size_t j = 0; j < f1.size(); ++j){
            if(f0[i] == f1[j]) {edge[k++] = f0[i];break;}
          }
      if(k != 2) throw std::logic_error("this two faces do not share one common edge.");
    }
#endif

  }
}
