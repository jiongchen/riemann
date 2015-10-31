#include "shell.h"

#include <iostream>
#include <unordered_set>
#include <zjucad/matrix/itr_matrix.h>
#include <zjucad/matrix/io.h>
#include <jtflib/mesh/mesh.h>

#include "def.h"
#include "config.h"

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
#pragma omp parallel
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

static inline double calc_tri_area(const matd_t &vert) {
  matd_t edge = vert(colon(), colon(1, 2))-vert(colon(), colon(0, 1));
  return 0.5*norm(cross(edge(colon(), 0), edge(colon(), 1)));
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
    itr_matrix<const double *> X(3, Nx()/3, x);
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
    itr_matrix<const double *> X(3, Nx()/3, x);
    for (size_t i = 0; i < edges_.size(2); ++i) {
      matd_t vert = X(colon(), edges_(colon(), i));
      matd_t grad = zeros<double>(6, 1);
      calc_edge_length_jac_(&grad[0], &vert[0]);
      grad *= w_/len_[i];
      for (size_t j = 0; j < 6; ++j) {
        if ( grad[j]*grad[j] != 0.0 )
          jac->push_back(Triplet<double>(off+i, 3*edges_(j/3, i)+j%3, grad[j]));
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
    angle_.resize(dias_.size(2), 1);
    len_.resize(dias_.size(2), 1);
    area_.resize(dias_.size(2), 1);
#pragma omp parallel for
    for (size_t i = 0; i < dias_.size(2); ++i) {
      matd_t vert = nods(colon(), dias_(colon(), i));
      calc_dihedral_angle_(&angle_[i], &vert[0]);
      calc_edge_length_(&len_[i], &vert(0, 1));
      area_[i] = calc_tri_area(vert(colon(), colon(0, 2)))
          + calc_tri_area(vert(colon(), colon(1, 3)));
    }
  }
  size_t Nx() const {
    return dim_;
  }
  size_t Nf() const {
    return dias_.size(2);
  }
  int Val(const double *x, double *val) const {
    itr_matrix<const double *> X(3, Nx()/3, x);
#pragma omp parallel for
    for (size_t i = 0; i < dias_.size(2); ++i) {
      matd_t vert = X(colon(), dias_(colon(), i));
      double curr = 0;
      calc_dihedral_angle_(&curr, &vert[0]);
      val[i] += w_*len_[i]/sqrt(area_[i])*(curr-angle_[i]);
    }
    return 0;
  }
  int Jac(const double *x, const size_t off, vector<Triplet<double>> *jac) const {
    itr_matrix<const double *> X(3, Nx()/3, x);
    for (size_t i = 0; i < dias_.size(2); ++i) {
      matd_t vert = X(colon(), dias_(colon(), i));
      matd_t grad = zeros<double>(12, 1);
      calc_dihedral_angle_jac_(&grad[0], &vert[0]);
      grad *= w_*len_[i]/sqrt(area_[i]);
      for (size_t j = 0; j < 12; ++j) {
        if ( grad[j]*grad[j] != 0.0 )
          jac->push_back(Triplet<double>(off+i, 3*dias_(j/3, i)+j%3, grad[j]));
      }
    }
    return 0;
  }
private:
  const size_t dim_;
  const double w_;
  const mati_t &dias_;
  matd_t angle_, len_, area_;
};

class volume_constraint : public Constraint<double>
{
public:
private:
};

class position_constraint : public Constraint<double>
{
public:
  position_constraint(const matd_t &nods, const double w)
    : dim_(nods.size()), w_(sqrt(w)) {
    po_.setZero(dim_);
  }
  size_t Nx() const {
    return dim_;
  }
  size_t Nf() const {
    return 3*pid_.size();
  }
  int Val(const double *x, double *val) const {
    if ( Nf() == 0 )
      return 1;
    Map<const VectorXd> X(x, Nx());
    Map<VectorXd> fx(val, Nf());
    size_t cnt = 0;
    for (auto &id : pid_) {
      fx.segment<3>(3*cnt) += w_*(X.segment<3>(3*id)-po_.segment<3>(3*id));
      ++cnt;
    }
    return 0;
  }
  int Jac(const double *x, const size_t off, vector<Triplet<double>> *jac) const {
    if ( Nf() == 0 )
      return 1;
    size_t cnt = 0;
    for (auto &id : pid_) {
      jac->push_back(Triplet<double>(off+cnt++, 3*id+0, w_));
      jac->push_back(Triplet<double>(off+cnt++, 3*id+1, w_));
      jac->push_back(Triplet<double>(off+cnt++, 3*id+2, w_));
    }
    return 0;
  }
  void add(const size_t id, const double *coords) {
    pid_.insert(id);
    po_.segment<3>(3*id) = Vector3d(coords);
  }
  void remove(const size_t id) {
    auto it = pid_.find(id);
    if ( it != pid_.end() )
      pid_.erase(it);
  }
private:
  const size_t dim_;
  double w_;
  unordered_set<size_t> pid_;
  VectorXd po_;
};
//==============================================================================
shell_deformer::shell_deformer(const mati_t &tris, const matd_t &nods, shell_args &args)
  : args_(args) {
  get_edge_elem(tris, edges_);
  get_diam_elem(tris, diams_);
  cbf_.resize(3);
  cbf_[0] = make_shared<position_constraint>(nods, args_.wp);
  cbf_[1] = make_shared<stretch_constraint>(edges_, nods, args_.ws);
  cbf_[2] = make_shared<bending_constraint>(diams_, nods, args_.wb);
}

void shell_deformer::fix_vert(const size_t id, const double *coords) {
  dynamic_pointer_cast<position_constraint>(cbf_[0])->add(id, coords);
}

void shell_deformer::free_vert(const size_t id) {
  dynamic_pointer_cast<position_constraint>(cbf_[0])->remove(id);
}

int shell_deformer::prepare() {
  try {
    constraint_ = make_shared<constraint_t<double>>(cbf_);
  } catch ( exception &e ){
    cerr << "[exception] " << e.what() << endl;
    exit(EXIT_FAILURE);
  }
  return 0;
}

int shell_deformer::solve(double *x) {
  Map<VectorXd> X(x, constraint_->Nx());
  VectorXd xstar = X;
  // gauss-newton
  for (size_t iter = 0; iter < args_.max_iter; ++iter) {
    VectorXd fc = VectorXd::Zero(constraint_->Nf()); {
      constraint_->Val(&xstar[0], &fc[0]);
      if ( iter % 10 == 0 ) {
        cout << "\t@iter " << iter << " error: " << fc.squaredNorm() << endl;
      }
    }
    SparseMatrix<double> J(constraint_->Nf(), constraint_->Nx()); {
      vector<Triplet<double>> trips;
      constraint_->Jac(&xstar[0], 0, &trips);
      J.reserve(trips.size());
      J.setFromTriplets(trips.begin(), trips.end());
    }
    solver_.compute(J.transpose()*J);
    ASSERT(solver_.info() == Success);
    VectorXd dx = solver_.solve(-J.transpose()*fc);
    ASSERT(solver_.info() == Success);
    const double xnorm = xstar.norm();
    // line search
    double h = 2.0;
    VectorXd xnew(constraint_->Nx()), fnew(constraint_->Nf());
    do {
      h *= 0.5;
      xnew = xstar+h*dx;
      fnew.setZero();
      constraint_->Val(&xnew[0], &fnew[0]);
    } while ( fnew.squaredNorm() >= fc.squaredNorm() && h > 1e-10 );
    xstar += h*dx;
    if ( h*dx.norm() <= args_.tolerance*xnorm ) {
      cout << "\t@converged, relative error is below " << args_.tolerance << "\n";
      break;
    }
  }
  X = xstar;
  return 0;
}

void shell_deformer::unit_test() const {
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
