#ifndef JTF_MESH_FEATURE_H
#define JTF_MESH_FEATURE_H

#include <string>
#include <vector>
#include <zjucad/matrix/matrix.h>

namespace jtf{

namespace mesh{

class feature
{
public:
  virtual std::string get_type() {return "";}
  virtual ~feature(){}
};


class feature_line: public feature
{
public:
  explicit feature_line(const zjucad::matrix::matrix<size_t> & line_idx,
                        const zjucad::matrix::matrix<double> & direction)
    :line_idx_(line_idx), direction_(direction){
  }
  virtual ~feature_line(){}
  virtual std::string get_type() {return "line";}
  const zjucad::matrix::matrix<size_t> & get_line_idx() const {
    return line_idx_;
  }
  zjucad::matrix::matrix<size_t> & get_line_idx() {
    return line_idx_;
  }
  const zjucad::matrix::matrix<double> & get_direction() const {
    return direction_;
  }
private:
  zjucad::matrix::matrix<size_t> line_idx_;
  zjucad::matrix::matrix<double> direction_;
};

class feature_point: public feature
{
public:
  explicit feature_point(const size_t & point_idx,
                         const zjucad::matrix::matrix<double> & node)
    :point_idx_(point_idx), position_(node){
  }
  virtual ~feature_point(){}
  virtual std::string get_type(){return "point";}
  size_t get_node_idx() const {return point_idx_;}
  const zjucad::matrix::matrix<double> & get_position() const { return position_;}
private:
  size_t point_idx_;
  const zjucad::matrix::matrix<double> position_;
};

}

}
#endif // JTF_MESH
