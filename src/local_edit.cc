#include "local_edit.h"

#include <jtflib/mesh/mesh.h>

#include "def.h"
#include "util.h"

using namespace std;
using namespace Eigen;
using namespace zjucad::matrix;

namespace riemann {

class handle_energy : public Functional<double>
{
public:
  handle_energy(const matd_t &nods, unordered_map<size_t, Vector3d> &hans, const double w)
    : dim_(nods.size()), hans_(hans), w_(w) {}
  size_t Nx() const {
    return dim_;
  }
  int Val(const double *x, double *val) const {
    Map<const MatrixXd> d(x, 3, dim_/3);
    for (auto &elem : hans_) {
      *val += 0.5*w_*(d.col(elem.first)-elem.second).squaredNorm();
    }
    return 0;
  }
  int Gra(const double *x, double *gra) const {
    Map<const MatrixXd> d(x, 3, dim_/3);
    Map<MatrixXd> G(gra, 3, dim_/3);
    for (auto &elem : hans_) {
      G.col(elem.first) += w_*(d.col(elem.first)-elem.second);
    }
    return 0;
  }
  int Hes(const double *x, vector<Triplet<double>> *hes) const {
    for (auto &elem : hans_) {
      hes->push_back(Triplet<double>(3*elem.first+0, 3*elem.first+0, w_));
      hes->push_back(Triplet<double>(3*elem.first+1, 3*elem.first+1, w_));
      hes->push_back(Triplet<double>(3*elem.first+2, 3*elem.first+2, w_));
    }
    return 0;
  }
private:
  const size_t dim_;
  double w_;
  unordered_map<size_t, Vector3d> &hans_;
};

class sparsity_energy : public Functional<double>
{
public:
  sparsity_energy(const matd_t &nods, unordered_map<size_t, Vector3d> &hans, const double w)
    : dim_(nods.size()), hans_(hans), w_(w) {}
  size_t Nx() const {
    return dim_;
  }
  int Val(const double *x, double *val) const {
    Map<const MatrixXd> d(x, 3, dim_/3);
    for (size_t i = 0; i < d.cols(); ++i) {
      if ( hans_.find(i) == hans_.end() )
        *val += 0.5*w_*d.col(i).norm();
    }
    return 0;
  }
  int Gra(const double *x, double *gra) const {
    Map<const MatrixXd> d(x, 3, dim_/3);
    for (size_t i = 0; i < d.cols(); ++i) {
      if ( hans_.find(i) == hans_.end() ) {

      }
    }
    return 0;
  }
  int Hes(const double *x, vector<Triplet<double>> *hes) const {
    return __LINE__;
  }
private:
  const size_t dim_;
  double w_;
  unordered_map<size_t, Vector3d> &hans_;
};

class fairness_energy : public Functional<double>
{
public:
  fairness_energy(const mati_t &quad, const matd_t &nods, const double w);
  size_t Nx() const {
    return dim_;
  }
  int Val(const double *x, double *val) const {
    return 0;
  }
  int Gra(const double *x, double *gra) const {
    return 0;
  }
  int Hes(const double *x, vector<Triplet<double>> *hes) const {
    return 0;
  }
private:
  const size_t dim_;
  double w_;
};

class subsitiution_planar_cons : public Constraint<double>
{

};

//==============================================================================
constrained_mesh_editor::constrained_mesh_editor(const mati_t &tris, const matd_t &nods) {

}

int constrained_mesh_editor::set_handles(const Json::Value &json) {
  if ( json.empty() ) {
    cerr << "[INFO] no handles\n";
    return __LINE__;
  }
  for (int i = 0; i < json.size(); ++i) {
    size_t vid = json[i]["vid"].asInt();
    double dx = json[i]["disp"][0].asDouble();
    double dy = json[i]["disp"][1].asDouble();
    double dz = json[i]["disp"][2].asDouble();
    handle_.insert(make_pair(vid, Vector3d(dx, dy, dz)));
  }
  return 0;
}

int constrained_mesh_editor::deform(double *d) const {
  return 0;
}

}
