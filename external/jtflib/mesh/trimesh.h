#ifndef TRIMESH_TRIMESH_H
#define TRIMESH_TRIMESH_H

#include <jtflib/mesh/mesh.h>
#include <jtflib/mesh/io.h>
#include <jtflib/mesh/util.h>
#include "config.h"

namespace jtf{
namespace mesh{
	class JTF_MESH_API tri_mesh
  {
  public:
    tri_mesh(const char *file){init(file);}
    tri_mesh(const jtf::mesh::meshes & trm):trimesh_(trm){init(trimesh_);}
    tri_mesh(){}

    //////////////////////////////////////////////////////////////////////////////
    jtf::mesh::meshes trimesh_;
    zjucad::matrix::matrix<double> face_normal_;
    zjucad::matrix::matrix<double> face_area_;
    std::shared_ptr<jtf::mesh::edge2cell_adjacent> ea_;

  private:
    void init(const char * file){
      using namespace zjucad::matrix;
      if(jtf::mesh::load_obj(file, trimesh_.mesh_, trimesh_.node_))
        throw std::invalid_argument("invalide tetmesh.");

      init(trimesh_);
    }
    void init(const jtf::mesh::meshes &trm){
        using namespace zjucad::matrix;
        jtf::mesh::cal_face_normal(trm.mesh_,trm.node_,face_normal_, true);

        face_area_.resize(trm.mesh_.size(2),1);
        for(size_t fi = 0; fi < trm.mesh_.size(2); ++fi){
            face_area_[fi] = jtf::mesh::cal_face_area(trm.mesh_(colon(),fi),trm.node_);
          }
        ea_.reset(jtf::mesh::edge2cell_adjacent::create(trm.mesh_));
        if(!ea_.get()){
            throw std::invalid_argument("invalide outside_face.");
          }
    }
  };
				}
}

#endif // TRIMESH_H
