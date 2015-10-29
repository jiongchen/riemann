#include "shell.h"

#include <iostream>
#include <unordered_set>
#include <zjucad/matrix/itr_matrix.h>
#include <jtflib/mesh/mesh.h>

#include "def.h"

using namespace std;
using namespace zjucad::matrix;
using namespace Eigen;
using jtf::mesh::edge2cell_adjacent;

namespace riemann {

extern "C" {

void calc_edge_length_(double *val, const double *x);
void calc_edge_length_jac_(double *jac, const double *x);

void calc_dihedral_angle_(double *val, const double *x);
void calc_dihedral_angle_jac_(double *jac, const double *x);

void calc_volume_(double *val, const double *x);
void calc_volume_jac_(double *jac, const double *x);

}

static void get_edge_elem(const mati_t &tris, mati_t &edge) {
  edge2cell_adjacent *e2c = edge2cell_adjacent::create(tris, false);
  edge.resize(2, e2c->edges_.size());
  for (size_t i = 0; i < edge.size(2); ++i) {
    edge(0, i) = e2c->edges_[i].first;
    edge(1, i) = e2c->edges_[i].second;
  }
  delete e2c;
}

static void get_diam_elem(const mati_t &tris, mati_t &diam) {
  edge2cell_adjacent *ea = edge2cell_adjacent::create(tris, false);
  mati_t bd_ed_id;
  get_boundary_edge_idx(*ea, bd_ed_id);
  diam.resize(4, ea->edges_.size()-bd_ed_id.size());
  for(size_t ei = 0, di = 0; ei < ea->edges_.size(); ++ei) {
    pair<size_t, size_t> nb_tr_id = ea->edge2cell_[ei];
    if( ea->is_boundary_edge(nb_tr_id) ) continue;
    diam(colon(1, 2), di) = ea->get_edge(ei);
    // orient
    bool need_swap = true;
    for(size_t k = 0; k < 3; ++k) {
      if( diam(1, di) == tris(k, nb_tr_id.first) ) {
        if( diam(2, di) != tris((k+1)%3, nb_tr_id.first) )
          need_swap = false;
      }
    }
    if( need_swap )
      swap(diam(1, di), diam(2, di));
    diam(0, di) = zjucad::matrix::sum(tris(colon(), nb_tr_id.first))
        - zjucad::matrix::sum(diam(colon(1, 2), di));
    diam(3, di) = zjucad::matrix::sum(tris(colon(), nb_tr_id.second))
        - zjucad::matrix::sum(diam(colon(1, 2), di));
    ++di;
  }
  delete ea;
}

class stretch_constraint : public Constraint<double>
{
public:
  stretch_constraint(const mati_t &edge, const matd_t &nods, const double w)
    : dim_(nods.size()), edges_(edge), w_(sqrt(w)) {
    len_.resize(edges_.size(2), 1);
#pragma omp parallel for
    for (size_t i = 0; i < edges_.size(2); ++i) {
      matd_t vert = nods(colon(), edges_(colon(), i));
      calc_edge_length_(&len_[i], &vert[0]);
    }
  }
  size_t Nx() const {
    return dim_;
  }
  size_t Nf() const {
    return edges_.size(2);
  }
  int Val(const double *x, double *val) const {
    itr_matrix<const double *> X(3, dim_/3, x);
#pragma omp parallel for
    for (size_t i = 0; i < edges_.size(2); ++i) {
      matd_t vert = X(colon(), edges_(colon(), i));
      double curr_len = 0;
      calc_edge_length_(&curr_len, &vert[0]);
      val[i] += w_/len_[i]*(curr_len-len_[i]);
    }
    return 0;
  }
  int Jac(const double *x, const size_t off, vector<Triplet<double>> *jac) const {
    itr_matrix<const double *> X(3, dim_/3, x);
    for (size_t i = 0; i < edges_.size(2); ++i) {
      matd_t vert = X(colon(), edges_(colon(), i));
      matd_t grad = zeros<double>(6, 1);
      calc_edge_length_jac_(&grad[0], &vert[0]);
      for (size_t j = 0; j < 6; ++j) {
        if ( grad[j] != 0.0 )
          jac->push_back(Triplet<double>(off+i, 3*edges_(j/3, i)+j%3, w_/len_[i]));
      }
    }
    return 0;
  }
private:
  const size_t dim_;
  const double w_;
  const mati_t &edges_;
  matd_t len_;
};

class bending_constraint : public Constraint<double>
{
public:
  bending_constraint(const mati_t &diamond, const matd_t &nods, const double w)
    : dim_(nods.size()), w_(sqrt(w)), dias_(diamond) {
    len_.resize(dias_.size(2), 1);
    area_.resize(dias_.size(2), 1);
#pragma omp parallel for
    for (size_t i = 0; i < dias_.size(2); ++i) {
      len_[i] = 0;
      area_[i] = 1;
    }
  }
  size_t Nx() const {
    return dim_;
  }
  size_t Nf() const {
    return dias_.size(2);
  }
  int Val(const double *x, double *val) const {
    itr_matrix<const double *> X(3, dim_/3, x);
#pragma omp parallel for
    for (size_t i = 0; i < dias_.size(2); ++i) {
      matd_t vert = X(colon(), dias_(colon(), i));
    }
    return 0;
  }
  int Jac(const double *x, const size_t off, vector<Triplet<double>> *jac) const {
    itr_matrix<const double *> X(3, dim_/3, x);
    for (size_t i = 0; i < dias_.size(2); ++i) {
      matd_t vert = X(colon(), dias_(colon(), i));
      matd_t grad = zeros<double>(12, 1);
      calc_dihedral_angle_jac_(&grad[0], &vert[0]);
    }
    return 0;
  }
private:
  const size_t dim_;
  const double w_;
  const mati_t &dias_;
  matd_t len_, area_;
};

class volume_constraint : public Constraint<double>
{
public:
};

class position_constraint : public Constraint<double>
{
public:
  position_constraint(const matd_t &nods, const double w)
    : dim_(nods.size()), w_(sqrt(w)) {}
  size_t Nx() const {
    return dim_;
  }
  size_t Nf() const {
    return 3*idx_.size();
  }
  int Val(const double *x, double *val) const {
    Map<const VectorXd> X(x, Nx());
    Map<VectorXd> fx(val, Nf());
    size_t cnt = 0;
    for (auto &id : idx_) {
      fx.segment<3>(3*cnt) = w_*(X.segment<3>(id)-po_.segment<3>(id));
      ++cnt;
    }
    return 0;
  }
  int Jac(const double *x, const size_t off, vector<Triplet<double>> *jac) const {
    size_t cnt = 0;
    for (auto &id : idx_) {
      jac->push_back(Triplet<double>(off+cnt++, 3*id+0, w_));
      jac->push_back(Triplet<double>(off+cnt++, 3*id+1, w_));
      jac->push_back(Triplet<double>(off+cnt++, 3*id+2, w_));
    }
    return 0;
  }
  void add(const size_t idx, const double *coords) {
    idx_.insert(idx);
    po_.segment<3>(3*idx) = Vector3d(coords);
  }
  void del(const size_t idx) {
    auto it = idx_.find(idx);
    if ( it != idx_.end() )
      idx_.erase(it);
  }
private:
  const size_t dim_;
  double w_;
  unordered_set<size_t> idx_;
  VectorXd po_;
};
//==============================================================================
void shell_deformer::temp_test() const {
  const double x[6] = {0,0,0, 3,4,5};
  double len = 0;
  calc_edge_length_(&len, x);
  cout << "len^2: " << len*len << endl;

  const double y[12] = {1,0,0, 0,0,0, 0,1,0, 0,0,-11};
  double angle = 0;
  calc_dihedral_angle_(&angle, y);
  cout << "angle: " << 180-angle/M_PI*180 << endl;
}

}
