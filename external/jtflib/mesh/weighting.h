#ifndef WEIGHTING_H
#define WEIGHTING_H

#include <iostream>
#include <numeric>
#include <zjucad/matrix/matrix.h>

#include "util.h"

enum PARAM {
  // CELL_TYPE
  TRI,
  QUAD,

  //OBJ_TYPE
  NODE,
  EDGE,
  FACE,

  //STRATEGY
  AVG,
  CON
};

template<PARAM CELL_TYPE, PARAM OBJ_TYPE, PARAM STRATEGY>
class mesh_weight
{
public:
  mesh_weight(const zjucad::matrix::matrix<size_t>& mesh,
              const zjucad::matrix::matrix<double> & node){}
  double get_w(const size_t idx);
};


template <>
class mesh_weight<TRI, NODE, AVG>
{
public:
  mesh_weight(const zjucad::matrix::matrix<size_t>& mesh,
              const zjucad::matrix::matrix<double> & node)
    : mesh_(mesh), node_(node){init();}
  double get_w(const size_t idx){
    if(idx < node_w_.size())
      return node_w_[idx];
    else{
        std::cerr << "# [error] can not locate node " << idx << std::endl;
        return 0;
      }
  }
protected:
  void init(){
    using namespace zjucad::matrix;
    node_w_ = zeros<double>(node_.size(2),1);
    for(size_t fi = 0; fi < mesh_.size(2); ++fi){
        const double area = jtf::mesh::cal_face_area(mesh_(colon(), fi), node_);
        for(size_t di = 0; di < 3; ++di){
            node_w_[mesh_(di,fi)] += area/3.0;
          }
      }
    const double total_area = std::accumulate(node_w_.begin(), node_w_.end(),0.0);
    node_w_ /= total_area;
  }
private:
  const zjucad::matrix::matrix<size_t> & mesh_;
  const zjucad::matrix::matrix<double> & node_;
  zjucad::matrix::matrix<double> node_w_;
};
#endif // WEIGHTING_H
