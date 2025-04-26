#ifndef JTF_MESH_H
#define JTF_MESH_H

#include <zjucad/matrix/matrix.h>
#include <zjucad/matrix/itr_matrix.h>
#include <utility>
#include <vector>
#include <map>
#include <list>
#include <set>
#include <iostream>
#include <boost/unordered_set.hpp>
#include <boost/unordered_map.hpp>
#include "config.h"

namespace jtf
{
  namespace mesh
  {

    typedef zjucad::matrix::matrix<double> matrixd;
    typedef zjucad::matrix::matrix<size_t> matrixst;

    //! @brief static mesh contains node/mesh matrix
	class JTF_MESH_API meshes
    {
    public:
      matrixd node_;
      matrixst mesh_;
    };

    ////////////////////////////////////////////////////////////////////////
    /// Mesh Topology Utils
    ////////////////////////////////////////////////////////////////////////


    ////////////////////////////////////////////////////////////////////////
    ///  2D MESH UTILS
    ///////////////////////////////////////////////////////////////////////
    //! @brief this class is used to build edge2cell relationship
    //!        it assume the mesh to be a manifold surface mesh
	class JTF_MESH_API edge2cell_adjacent
    {
    public:
      static edge2cell_adjacent *create(const matrixst & mesh,
                                        const bool &show_debug_info = true);

      static edge2cell_adjacent *create_on_dirty_mesh(const matrixst & mesh);

      size_t get_edge_idx(size_t vi, size_t vj) const;

      size_t get_edge_idx(const size_t *v) const;

      std::pair<size_t,size_t> query(size_t vi, size_t vj) const;

      std::pair<size_t,size_t> query(const size_t *v) const;

      static inline bool is_boundary_edge(const std::pair<size_t,size_t> &nb_tri_id) {
        return (nb_tri_id.first == -1) ^ (nb_tri_id.second == -1);
      }

      matrixst get_edge(size_t id) const {
        matrixst rtn(2);
        rtn[0] = edges_[id].first;
        rtn[1] = edges_[id].second;
        return rtn;
      }
      std::vector<std::pair<size_t,size_t> > edges_;
      std::vector<std::pair<size_t,size_t> > edge2cell_;
    private:
      edge2cell_adjacent(){}

      ///
      /// @brief init initialize the edge2cell relationship
      /// @param mesh input manifold and orientable surface mesh
      /// @param show_debug_info a switcher to show information for debugging
      /// @return return 0 if ok, or non-zeros
      ///
      int init(const matrixst &mesh, const bool & show_debug_info);

      ///
      /// @brief init initialize the edge2cell relationship
      /// @param mesh input manifold surface mesh
      /// @return return 0 if ok, or non-zeros
      ///
      int init_on_dirty_mesh(const matrixst &mesh);
    };


    //! @brief this class is used to build edge2cell relationship
    //!        it does not assume the mesh to be a manifold surface mesh
	class JTF_MESH_API edge2cell_adjacent_general
    {
    public:
      static edge2cell_adjacent_general *create(const matrixst &);

      size_t get_edge_idx(size_t vi, size_t vj) const;

      size_t get_edge_idx(const size_t *v) const;

      std::vector<size_t> query(size_t vi, size_t vj) const;

      std::vector<size_t> query(const size_t *v) const;

      static inline bool is_boundary_edge(const std::vector<size_t>& adj_faces) {
        return (adj_faces.size() > 1? false:true);
      }

      matrixst get_edge(size_t id) const {
        matrixst rtn(2);
        rtn[0] = edges_[id].first;
        rtn[1] = edges_[id].second;
        return rtn;
      }
      std::vector<std::pair<size_t,size_t> > edges_;
      std::vector<std::vector<size_t> > edge2cell_;
    private:

      edge2cell_adjacent_general(){}

      ///
      /// @brief init initialize the edge2cell relationship
      /// @param mesh input manifold surface mesh
      /// @return return 0 if ok, or non-zeros
      ///
      int init(const matrixst &mesh);
    };

    //! @brief This class is used to calculate the euler number of a surface mesh
	class JTF_MESH_API euler_number
    {
    public:

      euler_number(){}

      inline size_t operator()() const{
        return vertex_.size() + face_.size() - edge_.size();
      }

      ///
      /// @brief add_face add face, and automatically cound the face number,
      ///        vertex number, edge num
      /// @param face store the face points
      /// @param vertex_num store the number of points in a face
      /// @return return 0 if succeed, or return non-zeros
      ///
      int add_face(const size_t * face, const size_t &vertex_num);

    protected:

      int add_tri_face(const size_t * face);

      int add_quad_face(const size_t * face);

    private:
      std::set<size_t> vertex_;
      std::set<std::pair<size_t,size_t> > edge_;
      std::set<std::set<size_t> > face_;
    };

    //! @breif store N ring face around a point, each ring will be stored
    //!        in a vector, so for each point there will be N vectors in map.
	class JTF_MESH_API N_ring_face_at_point{
    public:
      void add_all_faces(const matrixst & face,
                         const edge2cell_adjacent & ea,
                         const size_t N);

      typedef boost::unordered_map<size_t,std::vector<std::vector<size_t> > > p2f_type;

      p2f_type p2f_;
    };


    //! @breif store N ring face around a point, each ring will be stored
    //!        in a vector, so for each point there will be N vectors in map.
	class JTF_MESH_API one_ring_face_at_point
    {
    public:
      void add_all_faces(const matrixst & faces,
                         const edge2cell_adjacent & ea)
      {
        for(size_t t = 0; t < faces.size(2); ++t)
          add_face(faces(zjucad::matrix::colon(),t),ea);
      }

      void add_face(const matrixst & one_face,
                    const edge2cell_adjacent & ea);

      ///
      /// @brief sort_int_loop it will sort faces around one point into a loop
      /// @param face input face
      /// @param node input node
      /// @param return 0 for success, others error
      int sort_int_loop(const matrixst &face,
                         const matrixd & node);

      void sort_int_loop_with_normal_info(
          const matrixst &face,
          const matrixd & node,
          const edge2cell_adjacent &ea,
          const matrixd & normal);

      typedef std::map<size_t,std::vector<size_t> > p2f_type;
      std::map<size_t,std::vector<size_t> > p2f_;
    private:
      void fix_boundary_loop();
    };

    //! @brief store one ring point of a point
	class JTF_MESH_API one_ring_point_at_point
    {
    public:
      static one_ring_point_at_point * create(const matrixst & faces);
      void sort_into_loop(const matrixst & face, const edge2cell_adjacent & ea);
      typedef boost::unordered_map<size_t, std::vector<size_t> > p2p_type;
      p2p_type p2p_;
    private:
      int init(const matrixst & faces);
      one_ring_point_at_point(){}
    };


    //! @brief  dynamically maintain the one ring face information of one point
    //! WARNING: this class will be removed later
	class JTF_MESH_API dyn_one_ring_face_at_point
    {
    public:

      typedef std::vector<size_t> face_type; // face should be sorted
      typedef std::list<face_type> face_list_type;
      typedef std::list<face_list_type::iterator>  linked_face_ptr_types;
      typedef std::list<face_list_type::iterator>::const_iterator  linked_face_const_iterator;
      typedef std::list<face_list_type::iterator>::iterator  linked_face_iterator;
      typedef boost::unordered_map<size_t,linked_face_ptr_types> point2face_type;
      typedef point2face_type::const_iterator pcit;
      typedef point2face_type::iterator pit;

      void add_all_faces(const matrixst & faces);

      void add_face(const matrixst & one_face);

      ///
      /// @brief get_one_ring_face_ptr
      /// @param point
      /// @return linked_face_ptr_types  if point can be found, return its one ring
      ///                            around linking faces, or return empty container
      ///
      linked_face_ptr_types get_one_ring_face_ptr(const size_t & point) const{
        point2face_type::const_iterator pcit = p2f_.find(point);
        if(pcit == p2f_.end()){
            linked_face_ptr_types lfpt;
            return lfpt;
          }else
          return pcit->second;
      }

      ///
      /// @brief is_face_inside: check whether this face is inside or not
      /// @param sorted_face: WARNING!!! input face should be sorted
      /// @return return true if can find this face, or return false
      ///
      bool is_face_inside(const std::vector<size_t> & sorted_face) const ;


      ///
      /// @brief To find a face in linked_face_ptr_types
      /// @param lfpt input linked_face_ptr_types
      /// @param face input face
      /// @return if found this face, return linked_face_iterator, or return lfpt.end()
      ///
      linked_face_iterator find_face_in_linked_ptr(
          linked_face_ptr_types & lfpt,
          const face_type & face)const;

      ///
      /// @brief remove_face
      /// @param face  input on face
      /// @param point_num input face point number
      /// @return return 0 if succeed, or return non-zeros
      ///
      int remove_face(const size_t *face, const size_t point_num);

      ///
      /// @brief glue_two_faces_at_shared_edge:
      ///       face0:edge and point_a; face1:edge and point_b
      ///       and update its one ring face info
      /// @param edge
      /// @param point_from
      /// @param point_to
      /// @return  return 0 if succeed, or return non-zeros
      ///
      int glue_two_faces_at_shared_edge(
          const std::pair<size_t,size_t> & edge,
          const size_t & point_from,const size_t & point_to);

      ///
      /// @brief convert_face_to_matrix
      /// @param faces
      /// @return  return 0 if succeed, or return non-zeros
      ///
      int convert_face_to_matrix(matrixst & faces) const;

      face_list_type faces_; /** stores all faces */
      point2face_type p2f_; /** stores all point's linking info */
    };



    ///
    /// @brief check whether the surface is manifold, by now I use edge adjacent faces
    ///        to check, if there is an edge associated with more than 2 faces, it's
    ///        non-manifold, or it's manifold
    /// @param surface
    /// @param dorfap_ptr used to check whether the point is non-manifold,
    ///         if dorfap_ptr is given, this function assume that
    ///         it corresponds to surface
    /// @return bool return true if it's manifold, or return false
    ///
    bool is_surface_mesh_manifold(
        const matrixst & surface,
        const dyn_one_ring_face_at_point * dorfap_ptr = 0);

    ///
    /// @brief extract all boundary edges from ea.
    /// @param ea input edge2cell adjacent relationship
    /// @param edge output edges
    /// @return return 0 if succeed, or non-zeros
    ///
    int JTF_MESH_API get_boundary_edge(const edge2cell_adjacent& ea,
                          matrixst &edge);


    ///
    /// @brief get_boundary_edge_idx, index is recorded in ea
    /// @param ea
    /// @param edge_idx
    /// @return
    ///
    int get_boundary_edge_idx(const edge2cell_adjacent& ea,
                              matrixst &edge_idx);

    ////
    /// @brief find one common edge between two faces
    /// @param face0
    /// @param face1
    /// @param edge
    /// @return 0 if find one edge, or return non-zeros
    ///
	int JTF_MESH_API find_common_edge(const matrixst & face0,
                         const matrixst & face1,
                         size_t * edge);

    ///
    /// @brief reorder_face points to make each face normal is unified
    /// @param face: input triangle face
    /// @return 0 if ok, or non-zeros
    ///
	int JTF_MESH_API reorder_face(matrixst & face, const matrixd * node = nullptr);


    /////////////////////////////////////////////////////////////////////////
    /// 3D MESH UTILS
    /////////////////////////////////////////////////////////////////////////

	class JTF_MESH_API face2tet_adjacent
    {
    public:
      static face2tet_adjacent *create(
          const matrixst &tet,
          const std::string & strategy = "topology");

      std::string get_strategy() const {return strategy_;}
      size_t get_face_idx(size_t a, size_t b, size_t c) const;
      size_t get_face_idx(const size_t *abc) const;

      std::pair<size_t, size_t> query(size_t a, size_t b, size_t c) const;
      std::pair<size_t, size_t> query(const size_t *abc) const;
      static inline bool is_outside_face(const std::pair<size_t, size_t> &nb_tet_id) {
        return (nb_tet_id.first  == -1) ^ (nb_tet_id.second == -1);
      }

      matrixst get_face(size_t id) const {
        matrixst rtn(3);
        std::copy(faces_[id].begin(), faces_[id].end(), rtn.begin());
        return rtn;
      }

      // vector of length 3 vector
      std::vector<std::vector<size_t> > faces_;
      std::vector<std::pair<size_t, size_t> > face2tet_;
    private:
      std::map<std::vector<size_t>,size_t> face2idx_;

      face2tet_adjacent(){}
      int init_geometry(const matrixst &tet);
      int init_topology(const matrixst &tet);
      std::string strategy_;
    };

	class JTF_MESH_API one_ring_tet_at_point
    {
    public:
      void add_tets(const matrixst & tets){
        for(size_t ti = 0; ti < tets.size(2); ++ti){
            for(size_t pi = 0; pi < tets.size(1); ++pi){
                p2t_[tets(pi,ti)].push_back(ti);
              }
          }
      }
      std::map<size_t,std::vector<size_t> > p2t_;
    };

	class JTF_MESH_API one_ring_tet_at_edge
    {
    public:
      ///
      /// @brief add_all_tets,This function is used to add all tets into edge map
      /// @param tets, input all tets
      /// @param fa, input face2tet_adjacent
      ///
      void add_tets(const matrixst & tets,
                    const jtf::mesh::face2tet_adjacent & fa);

      ///
      /// @brief This function is used to sort the tets  arounding edge
      ///        into a loop, and direction will be controlled by dir
      /// @param tet, input tet
      /// @param node, input node
      /// @param dir, input direction, default is 1 which stands for right
      ///        hand direction, wile -1 stands for left hand direction,
      ///        and 0 stands for no-direction
      /// @return int return 0 if works fine, or retun non-zeros
      ///
      int sort_into_loop( const matrixst &tet,
                          const matrixd &node, const int dir = 1);

      ///
      /// @brief get_one_ring_points_of_edge,
      /// @param tet, input tet
      /// @param edge, input edge
      /// @param points, output one ring points of edge
      /// @return 0 if ok, or non-zeros if error.
      ///
      int get_one_ring_points_of_edge(
          const matrixst & tet,
          const std::pair<size_t,size_t> & edge,
          std::vector<size_t> & points) const;

      ///
      /// @brief is_inner_edge, test whether a loop imply an inner edge
      /// @param loop input loop tets
      /// @return true if is inner edge, or return false
      ///
      bool is_inner_edge(const std::vector<size_t> & loop) const
      {
        return (loop.front() == -1?false:true);
      }

      ///
      /// @brief is_chain_in_right_hand_order
      /// @param loop, input tets loop
      /// @param edge, input one edge
      /// @param tet, input a tet
      /// @param node, input node
      /// @return true if is in right hand order, or return false
      ///
      bool is_chain_in_right_hand_order(const std::vector<size_t> & loop,
                                        const std::pair<size_t,size_t> & edge,
                                        const matrixst & tet,
                                        const matrixd & node)const;

      // map: stores the begin and end vertex() of edge and the arounding tets;
      typedef boost::unordered_map<std::pair<size_t,size_t>, std::vector<size_t> > e2tet_type;
      e2tet_type e2t_;
    private:
      ///
      /// @brief add_tet, This function is used to add each tet into the edge map
      /// @param one_tet input one_tet
      /// @param fa input face2tet_adjacent
      ///
      void add_one_tet(const matrixst & one_tet,
                       const jtf::mesh::face2tet_adjacent &fa);

    };


	class JTF_MESH_API one_ring_face_at_edge
    {
    public:
      void build(const one_ring_tet_at_edge & ortae,
                 const jtf::mesh::face2tet_adjacent & fa,
                 const matrixst & tet);
      typedef boost::unordered_map<std::pair<size_t,size_t>,std::vector<size_t> >
      edge2_face_type;
      boost::unordered_map<std::pair<size_t,size_t>,std::vector<size_t> > e2f_;
    };

    ///
    /// @brief get_outside_face
    /// @param fa
    /// @param face
    /// @param check_order whethet to order point of one face
    ///
    void get_outside_face(const face2tet_adjacent &fa,
                          matrixst &face,
                          bool check_order = false,
                          const matrixd * node = nullptr);

    ///
    /// @brief get_outside_face_idx
    /// @param fa
    /// @param face_idx
    ///
    void get_outside_face_idx(const face2tet_adjacent &fa,
                              matrixst &face_idx);

    ///
    /// @brief find_common_face between two tets
    /// @param tet_i input tet_i
    /// @param tet_j input tet_j
    /// @param face output face
    /// @return 0 if find a face, non-zeros if no common face
    ///
    template <typename E1,typename E2, typename E3>
    int find_common_face(const zjucad::matrix::matrix_expression<E1> &tet_i,
                         const zjucad::matrix::matrix_expression<E2> &tet_j,
                         zjucad::matrix::matrix_expression<E3> &face )
    {
      std::map<size_t,size_t> tmp;
      for(size_t t = 0; t < 4; ++t)
        {
          if(tmp.find(tet_i()[t]) == tmp.end())
            tmp[tet_i()[t]] = 1;
          else
            tmp[tet_i()[t]] += 1;
          if(tmp.find(tet_j()[t]) == tmp.end())
            tmp[tet_j()[t]] = 1;
          else
            tmp[tet_j()[t]] += 1;
        }
      if(tmp.size() != 5) return __LINE__; // have no common face
      size_t idx = 0;
      face().resize(3,1);
      for(std::map<size_t,size_t>::const_iterator mci = tmp.begin();
          mci != tmp.end(); ++mci){
          if(mci->second == 2){
              face()[idx++] = mci->first;
            }
        }
      return 0;
    }

    ///
    /// @brief find_common_face between two tets, a wrapper
    /// @param tet_i input tet_i
    /// @param tet_j input tet_j
    /// @param face output face
    /// @return 0 if find a face, non-zeros if no common face
    ///
    template <typename E1, typename E2>
    int find_common_face(const zjucad::matrix::matrix_expression<E1> &tet_i,
                         const zjucad::matrix::matrix_expression<E2> &tet_j,
                         size_t *face )
    {
      zjucad::matrix::itr_matrix<size_t*> face0(3,1,face);
      return find_common_face(tet_i, tet_j, face0);
    }


    ////////////////////////////////////////////////////////////////////////////
    ///  face 2 hex
    ///
	class JTF_MESH_API face2hex_adjacent{
    public:
      typedef matrixst matrixi;
      static face2hex_adjacent *create(const matrixi &hex);

      size_t get_face_idx(size_t a, size_t b, size_t c, size_t d) const;
      size_t get_face_idx(const size_t *abcd) const;
      size_t get_common_face_idx(const matrixst &hex,
                                 const std::pair<size_t,size_t> &two_hex) const ;
      std::pair<size_t,size_t> query(size_t a, size_t b, size_t c, size_t d) const;
      std::pair<size_t,size_t> query(const size_t *abcd) const;

      static inline bool is_outside_face(const std::pair<size_t,size_t> &nb_hex_idx){
        return (nb_hex_idx.first == -1) ^ (nb_hex_idx.second == -1);
      }

      matrixi get_face(size_t id) const{
        matrixi rtn(4);
        copy(faces_[id].begin(),faces_[id].end(),rtn.begin());
        return rtn;
      }

      // vector of length 4 vector
      std::vector<std::vector<size_t> > faces_;
      std::vector<std::pair<size_t,size_t> > face2hex_;
    private:
      face2hex_adjacent(){}
      int init(const matrixi &hex);
      std::vector<std::vector<size_t> > faces_sort_;
    };


	class JTF_MESH_API one_ring_hex_at_edge
    {
    public:
      void add_hex(const zjucad::matrix::matrix<size_t> & hex,
                   const zjucad::matrix::matrix<double> &node,
                   const face2hex_adjacent &fa){
        for(size_t hi = 0; hi < hex.size(2); ++hi)
          add_one_hex(hex(zjucad::matrix::colon(),hi), node, fa);
      }

      int sort_into_loop(const zjucad::matrix::matrix<size_t> &hex,
                         const zjucad::matrix::matrix<double> &node,
                         const face2hex_adjacent &fa);

      typedef std::map<std::pair<size_t, size_t>, std::vector<size_t> > e2hex_type;
      e2hex_type e2h_;
    private:
      void add_one_hex(const zjucad::matrix::matrix<size_t> & one_hex,
                       const zjucad::matrix::matrix<double> &node,
                       const face2hex_adjacent &fa);
    };

    //! face is 4 * n
	int JTF_MESH_API get_outside_face(const face2hex_adjacent &fa,
                         matrixst &face);
    //! face idx size is n
	int JTF_MESH_API get_outside_face_idx(const face2hex_adjacent &fa,
                             matrixst &face_idx);

	int JTF_MESH_API get_faces_for_one_hex(const matrixst &one_hex,
                              matrixst &faces);

	int JTF_MESH_API get_edges_for_one_hex(const matrixst &one_hex,
                              matrixst & edges);

  }
}

#endif // MESH_H
