#ifndef GREEN_COORD_DEFORM_H
#define GREEN_COORD_DEFORM_H

#include <Eigen/Sparse>
#include <zjucad/matrix/matrix.h>

namespace geom_deform {

class green_deform
{
public:
  green_deform() {}
  virtual ~green_deform() {}
  virtual int load_sample_points(const char *file) = 0;
  virtual int load_cage(const char *file) = 0;
  virtual int calc_green_coords() = 0;
  virtual int move_cage(const size_t id, const double *dx, bool disp) = 0;
  virtual int deform() = 0;
  virtual int dump(const char *file) = 0;
protected:
  virtual int calc_outward_normal() = 0;
};

class green_deform_2d : public green_deform
{
public:
  green_deform_2d();
  int load_sample_points(const char *file);
  int load_cage(const char *file);
  int calc_green_coords();
  int move_cage(const size_t id, const double *dx, bool disp);
  int deform();
  int dump(const char *file);
  int dump_normal(const char *file);
  int dump_cage(const char *file);
protected:
  int calc_outward_normal();
  int calc_cage_edge_length(Eigen::VectorXd &len);
private:
  Eigen::MatrixXi cell_;
  Eigen::MatrixXd nods_;

  Eigen::MatrixXi cage_cell_;
  Eigen::MatrixXd cage_nods_;
  Eigen::MatrixXd cage_normal_;

  Eigen::VectorXd rest_len_;
  Eigen::VectorXd curr_len_;
  Eigen::MatrixXd phi_;
  Eigen::MatrixXd psi_;
};

class green_deform_3d : public green_deform
{
public:
  typedef zjucad::matrix::matrix<size_t> mati_t;
  typedef zjucad::matrix::matrix<double> matd_t;
  green_deform_3d();
  int load_sample_points(const char *file);
  int load_cage(const char *file);
  int calc_green_coords();
  int move_cage(const size_t id, const double *dx, bool disp);
  int deform();
  int dump(const char *file);
  int dump_normal(const char *file);
  int dump_cage(const char *file);
protected:
  int calc_outward_normal();
  int calc_stretch_ratio(matd_t &s);
  int calc_cage_edge(matd_t &uv);
  double GCTriInt(const matd_t &p, const matd_t &v1,
                  const matd_t &v2, const matd_t &eta);
private:
  mati_t cell_;
  matd_t nods_;

  mati_t cage_cell_;
  matd_t cage_nods_;
  matd_t cage_normal_;
  matd_t cage_rest_uv_;
  matd_t cage_curr_uv_;

  matd_t rest_area_;
  matd_t phi_;
  matd_t psi_;
};

}

#endif
