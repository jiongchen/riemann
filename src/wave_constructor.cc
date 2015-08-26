#include "wave_constructor.h"

#include <iostream>
#include <zjucad/matrix/itr_matrix.h>
#include <zjucad/matrix/io.h>
#include <jtflib/mesh/io.h>
#include <Eigen/Geometry>

#include "def.h"
#include "write_vtk.h"
#include "config.h"
#include "util.h"
#include "vtk.h"

using namespace std;
using namespace zjucad::matrix;
using namespace Eigen;
using namespace jtf::mesh;
using mati_t = zjucad::matrix::matrix<size_t>;
using matd_t = zjucad::matrix::matrix<double>;

namespace riemann {

extern "C" {

void trans_condition_(double *val, const double *f, const double *Cij, const double *Cji);
void trans_condition_jac_(double *jac, const double *f, const double *Cij, const double *Cji);
void trans_condition_hes_(double *hes, const double *f, const double *Cij, const double *Cji);

void modulus_condition_(double *val, const double *f);
void modulus_condition_jac_(double *jac, const double *f);
void modulus_condition_hes_(double *hes, const double *f);

void phase_condition_(double *val, const double *f);
void phase_condition_jac_(double *jac, const double *f);
void phase_condition_hes_(double *hes, const double *f);

}

//class wave_value_energy : public Functional<double>
//{
//public:
//  wave_value_energy(const mati_t &edge, const matd_t &f, const matd_t &c, const double w=1.0)
//    : edge_(edge), dim_(f.size()), c_(c), w_(w) {}
//  size_t Nx() const {
//    return dim_;
//  }
//  int Val(const double *x, double *val) const {
//    itr_matrix<const double *> X(4, dim_/4, x);
//    for (size_t i = 0; i < edge_.size(2); ++i) {
//      matd_t vert = X(colon(), edge_(colon(), i));
//      double value = 0;
//      wave_value_condition_(&value, &vert[0], &c_(0, i));
//      *val += w_*value;
//    }
//    return 0;
//  }
//  int Gra(const double *x, double *gra) const {
//    itr_matrix<const double *> X(4, dim_/4, x);
//    itr_matrix<double *> G(4, dim_/4, gra);
//    for (size_t i = 0; i <edge_.size(2); ++i) {
//      matd_t vert = X(colon(), edge_(colon(), i));
//      matd_t g = zeros<double>(4, 2);
//      wave_value_condition_jac_(&g[0], &vert[0], &c_(0, i));
//      G(colon(), edge_(colon(), i)) += w_*g;
//    }
//    return 0;
//  }
//  int Hes(const double *x, vector<Triplet<double>> *hes) const {
//    itr_matrix<const double *> X(4, dim_/4, x);
//    for (size_t i = 0; i < edge_.size(2); ++i) {
//      matd_t vert = X(colon(), edge_(colon(), i));
//      matd_t H = zeros<double>(8, 8);
//      wave_value_condition_hes_(&H[0], &vert[0], &c_(0, i));
//      for (size_t q = 0; q < H.size(2); ++q) {
//        for (size_t p = 0; p < H.size(1); ++p) {
//          if ( H(p, q) != 0.0 ) {
//            size_t I = 4*edge_(p/4, i)+p%4;
//            size_t J = 4*edge_(q/4, i)+q%4;
//            hes->push_back(Triplet<double>(I, J, w_*H(p, q)));
//          }
//        }
//      }
//    }
//    return 0;
//  }
//private:
//  const mati_t &edge_;
//  const matd_t &c_;
//  const size_t dim_;
//  const double w_;
//};

//class modulus_energy : public Functional<double>
//{
//public:
//  modulus_energy(const matd_t &f, const double w=0.15)
//    : dim_(f.size()), w_(w) {}
//  size_t Nx() const {
//    return dim_;
//  }
//  int Val(const double *x, double *val) const {
//    itr_matrix<const double *> X(4, dim_/4, x);
//    for (size_t i = 0; i < X.size(2); ++i) {
//      double value = 0;
//      modulus_condition_(&value, &X(0, i));
//      *val += w_*value;
//    }
//    return 0;
//  }
//  int Gra(const double *x, double *gra) const {
//    itr_matrix<const double *> X(4, dim_/4, x);
//    itr_matrix<double *> G(4, dim_/4, gra);
//#pragma omp parallel for
//    for (size_t i = 0; i < X.size(2); ++i) {
//      matd_t g  = zeros<double>(4, 1);
//      modulus_condition_jac_(&g[0], &X(0, i));
//      G(colon(), i) += w_*g;
//    }
//    return 0;
//  }
//  int Hes(const double *x, vector<Triplet<double>> *hes) const {
//    itr_matrix<const double *> X(4, dim_/4, x);
//    for (size_t i = 0; i < X.size(2); ++i) {
//      matd_t H = zeros<double>(4, 4);
//      modulus_condition_hes_(&H[0], &X(0, i));
//      for (size_t q = 0; q < H.size(2); ++q) {
//        for (size_t p = 0; p < H.size(1); ++p)  {
//          if ( H(p, q) != 0.0 )
//            hes->push_back(Triplet<double>(4*i+p, 4*i+q, w_*H(p, q)));
//        }
//      }
//    }
//    return 0;
//  }
//private:
//  const size_t dim_;
//  const double w_;
//};

//class phase_energy : public Functional<double>
//{
//public:
//  phase_energy(const matd_t &f, const double w=0.15)
//    : dim_(f.size()), w_(w) {}
//  size_t Nx() const {
//    return dim_;
//  }
//  int Val(const double *x, double *val) const {
//    itr_matrix<const double *> X(4, dim_/4, x);
//    for (size_t i = 0; i < X.size(2); ++i) {
//      double value = 0;
//      phase_condition_(&value, &X(0, i));
//      *val += w_*value;
//    }
//    return 0;
//  }
//  int Gra(const double *x, double *gra) const {
//    itr_matrix<const double *> X(4, dim_/4, x);
//    itr_matrix<double *> G(4, dim_/4, gra);
//#pragma omp parallel for
//    for (size_t i = 0; i < X.size(2); ++i) {
//      matd_t g  = zeros<double>(4, 1);
//      phase_condition_jac_(&g[0], &X(0, i));
//      G(colon(), i) += w_*g;
//    }
//    return 0;
//  }
//  int Hes(const double *x, vector<Triplet<double>> *hes) const {
//    itr_matrix<const double *> X(4, dim_/4, x);
//    for (size_t i = 0; i < X.size(2); ++i) {
//      matd_t H = zeros<double>(4, 4);
//      phase_condition_hes_(&H[0], &X(0, i));
//      for (size_t q = 0; q < H.size(2); ++q) {
//        for (size_t p = 0; p < H.size(1); ++p)  {
//          if ( H(p, q) != 0.0 )
//            hes->push_back(Triplet<double>(4*i+p, 4*i+q, w_*H(p, q)));
//        }
//      }
//    }
//    return 0;
//  }
//private:
//  const size_t dim_;
//  const double w_;
//};

class wave_value_cons : public Constraint<double>
{
public:
  wave_value_cons(const mati_t &edge, const matd_t &f, const matd_t &cIJ, const matd_t &cJI, const double w=1.0)
    : edge_(edge), cIJ_(cIJ), cJI_(cJI), nx_(f.size()), nf_(2*edge.size(2)), w_(sqrt(w)) {}
  size_t Nx() const {
    return nx_;
  }
  size_t Nf() const {
    return nf_;
  }
  int Val(const double *x, double *val) const {
    itr_matrix<const double *> X(4, nx_/4, x);
#pragma omp parallel for
    for (size_t i = 0; i < edge_.size(2); ++i) {
      matd_t vert = X(colon(), edge_(colon(), i));
      val[2*i+0] += w_*(vert(0, 1)-dot(cIJ_(colon(), i), vert(colon(), 0)));
      val[2*i+1] += w_*(vert(0, 0)-dot(cJI_(colon(), i), vert(colon(), 1)));
    }
    return 0;
  }
  int Jac(const double *x, const size_t off, vector<Triplet<double>> *jac) const {
    for (size_t i = 0; i < edge_.size(2); ++i) {
      // I->J
      jac->push_back(Triplet<double>(off+2*i+0, 4*edge_(1, i)+0, w_*1.0));
      jac->push_back(Triplet<double>(off+2*i+0, 4*edge_(0, i)+0, -w_*cIJ_(0, i)));
      jac->push_back(Triplet<double>(off+2*i+0, 4*edge_(0, i)+1, -w_*cIJ_(1, i)));
      jac->push_back(Triplet<double>(off+2*i+0, 4*edge_(0, i)+2, -w_*cIJ_(2, i)));
      jac->push_back(Triplet<double>(off+2*i+0, 4*edge_(0, i)+3, -w_*cIJ_(3, i)));
      // J->I
      jac->push_back(Triplet<double>(off+2*i+1, 4*edge_(0, i)+0, w_*1.0));
      jac->push_back(Triplet<double>(off+2*i+1, 4*edge_(1, i)+0, -w_*cJI_(0, i)));
      jac->push_back(Triplet<double>(off+2*i+1, 4*edge_(1, i)+1, -w_*cJI_(1, i)));
      jac->push_back(Triplet<double>(off+2*i+1, 4*edge_(1, i)+2, -w_*cJI_(2, i)));
      jac->push_back(Triplet<double>(off+2*i+1, 4*edge_(1, i)+3, -w_*cJI_(3, i)));
    }
    return 0;
  }
private:
  const mati_t &edge_;
  const matd_t &cIJ_, &cJI_;
  const size_t nx_, nf_;
  const double w_;
};

class modulus_cons : public Constraint<double>
{
public:
  modulus_cons(const matd_t &f, const double w)
    : nx_(f.size()), nf_(f.size(2)), w_(sqrt(w)) {}
  size_t Nx() const {
    return nx_;
  }
  size_t Nf() const {
    return nf_;
  }
  int Val(const double *x, double *val) const {
    itr_matrix<const double *> X(4, nx_/4, x);
#pragma omp parallel for
    for (size_t i = 0; i < X.size(2); ++i) {
      matd_t vert = X(colon(), i);
      val[i] += w_*(dot(vert, vert)-1.0);
    }
    return 0;
  }
  int Jac(const double *x, const size_t off, vector<Triplet<double>> *jac) const {
    itr_matrix<const double *> X(4, nx_/4, x);
    for (size_t i = 0; i < X.size(2); ++i) {
      jac->push_back(Triplet<double>(off+i, 4*i+0, 2*w_*X(0, i)));
      jac->push_back(Triplet<double>(off+i, 4*i+1, 2*w_*X(1, i)));
      jac->push_back(Triplet<double>(off+i, 4*i+2, 2*w_*X(2, i)));
      jac->push_back(Triplet<double>(off+i, 4*i+3, 2*w_*X(3, i)));
    }
    return 0;
  }
private:
  const size_t nx_, nf_;
  const double w_;
};

class phase_cons : public Constraint<double>
{
public:
  phase_cons(const matd_t &f, const double w)
    : nx_(f.size()), nf_(f.size(2)), w_(sqrt(w)) {}
  size_t Nx() const {
    return nx_;
  }
  size_t Nf() const {
    return nf_;
  }
  int Val(const double *x, double *val) const {
    itr_matrix<const double *> X(4, nx_/4, x);
#pragma omp parallel for
    for (size_t i = 0; i < X.size(2); ++i) {
      val[i] += w_*(X(0, i)*X(3, i)-X(1, i)*X(2, i));
    }
    return 0;
  }
  int Jac(const double *x, const size_t off, vector<Triplet<double> > *jac) const {
    itr_matrix<const double *> X(4, nx_/4, x);
    for (size_t i = 0; i < X.size(2); ++i) {
      jac->push_back(Triplet<double>(off+i, 4*i+0, w_*X(3, i)));
      jac->push_back(Triplet<double>(off+i, 4*i+1, -w_*X(2, i)));
      jac->push_back(Triplet<double>(off+i, 4*i+2, -w_*X(1, i)));
      jac->push_back(Triplet<double>(off+i, 4*i+3, w_*X(0, i)));
    }
    return 0;
  }
private:
  const size_t nx_, nf_;
  const double w_;
};

class feature_cons : public Constraint<double>
{
public:
  feature_cons(const matd_t &f, const map<size_t, size_t> &vert_on_line, const double w=1.0)
    : nx_(f.size()), nf_(0), vert_on_line_(vert_on_line), w_(sqrt(w)) {
    for (auto &it : vert_on_line_) {
      if ( it.second == 1 )
        nf_ += 1;
      else if ( it.second == 2 )
        nf_ += 3;
    }
  }
  size_t Nx() const {
    return nx_;
  }
  size_t Nf() const {
    return nf_;
  }
  int Val(const double *x, double *val) const {
    itr_matrix<const double *> X(4, nx_/4, x);
    size_t cnt = 0;
    for (auto &it : vert_on_line_) {
      switch ( it.second ) {
        case 1:
          val[cnt++] += w_*X(3, it.first);
          break;
        case 2:
          val[cnt++] += w_*X(1, it.first);
          val[cnt++] += w_*X(2, it.first);
          val[cnt++] += w_*X(3, it.first);
          break;
        default:
          break;
      }
    }
    return 0;
  }
  int Jac(const double *x, const size_t off, vector<Triplet<double>> *jac) const {
    size_t cnt = 0;
    for (auto &it : vert_on_line_) {
      switch ( it.second ) {
        case 1:
          jac->push_back(Triplet<double>(off+cnt++, 4*it.first+3, w_));
          break;
        case 2:
          jac->push_back(Triplet<double>(off+cnt++, 4*it.first+1, w_));
          jac->push_back(Triplet<double>(off+cnt++, 4*it.first+2, w_));
          jac->push_back(Triplet<double>(off+cnt++, 4*it.first+3, w_));
          break;
        default:
          break;
      }
    }
    return 0;
  }
private:
  const size_t nx_;
  size_t nf_;
  const double w_;
  const map<size_t, size_t> &vert_on_line_;
};

//==============================================================================

void rotate_frame(const double *oriX, const double *oriY,
                  const double *refX, const double *refY,
                  double       *rotX, double       *rotY) {
  Map<const Vector3d> ori_x(oriX), ori_y(oriY), ref_x(refX), ref_y(refY);
  Map<Vector3d> rot_x(rotX), rot_y(rotY);
  Vector3d ori_nor = ori_x.cross(ori_y), ref_nor = ref_x.cross(ref_y);
  ori_nor /= ori_nor.norm();
  ref_nor /= ref_nor.norm();
  Matrix3d R;
  if ( std::fabs(ori_nor.dot(ref_nor)-1.0) < 1e-10 ) {
    R.setIdentity();
  } else if ( std::fabs(ori_nor.dot(ref_nor)+1.0) < 1e-10 ) {
    R = -Matrix3d::Identity();
  } else {
    Vector3d axis = ori_nor.cross(ref_nor);
    axis /= axis.norm();
    double angle = safe_acos(ori_nor.dot(ref_nor));
    R = AngleAxisd(angle, axis);
  }
  rot_x = R*ori_x;
  rot_y = R*ori_y;
}

/// @brief oriX, oriY are coplanar with refX, refY
int align_vector_field(const double *oriX, const double *oriY,
                       const double *refX, const double *refY,
                       double     *alignX, double     *alignY) {
  Map<const Vector3d> ori_x(oriX), ori_y(oriY), ref_x(refX), ref_y(refY);
  Map<Vector3d> align_x(alignX), align_y(alignY);
  VectorXd nor = ori_x.cross(ori_y);
  nor /= nor.norm();
  Matrix3d R;
  R = AngleAxisd(M_PI/2, nor);
  vector<double> buff(4);
  buff[0] = ref_x.dot(ori_x);
  buff[1] = ref_x.dot(R*ori_x);
  buff[2] = ref_x.dot(R*R*ori_x);
  buff[3] = ref_x.dot(R*R*R*ori_x);

  size_t kappa = std::max_element(buff.begin(), buff.end())-buff.begin();
  ASSERT(kappa >= 0 && kappa < 4);
  switch ( kappa ) {
    case 0:
      align_x = ori_x;
      align_y = ori_y;
      break;
    case 1:
      align_x = ori_y;
      align_y = -ori_x;
      break;
    case 2:
      align_x = -ori_x;
      align_y = -ori_y;
      break;
    case 3:
      align_x = -ori_y;
      align_y = ori_x;
      break;
    default: break;
  }
  return 0;
}

//==============================================================================
int wave_constructor::load_model_from_obj(const char *filename) {
  jtf::mesh::load_obj(filename, tris_, nods_);
  extract_edges();
  return 0;
}

int wave_constructor::extract_edges() {
  e2c_.reset(edge2cell_adjacent::create(tris_, false));
  if ( !e2c_.get() )
    return __LINE__;
  edges_.resize(2, e2c_->edges_.size());
#pragma omp parallel for
  for (size_t i = 0; i < e2c_->edges_.size(); ++i) {
    edges_(0, i) = e2c_->edges_[i].first;
    edges_(1, i) = e2c_->edges_[i].second;
  }
  return 0;
}

int wave_constructor::load_frame_field(const char *filename) {
  if ( !e2c_.get() )
    return __LINE__;
  ifstream ifs(filename);
  if ( ifs.fail() ) {
    cerr << "[error] can not open " << filename << endl;
    return __LINE__;
  }
  size_t edge_num;
  ifs >> edge_num;
  ASSERT(edge_num == e2c_->edges_.size());
  edge_frm_.resize(6, edge_num);
  size_t p, q, eid;
  for (size_t i = 0; i < edge_num; ++i) {
    ifs >> p >> q;
    eid = e2c_->get_edge_idx(p, q);
    ifs >> edge_frm_(0, eid) >> edge_frm_(1, eid) >> edge_frm_(2, eid)
        >> edge_frm_(3, eid) >> edge_frm_(4, eid) >> edge_frm_(5, eid);
  }
  e2c_.release();
  ifs.close();
  return 0;
}

int wave_constructor::vis_edge_frame_field(const char *file_x, const char *file_y, const double scale) const {
  if ( edges_.size() == 0 ) {
    cerr << "[error] edge is not initialized\n";
    return __LINE__;
  }
  matd_t X = scale*edge_frm_(colon(0, 2), colon());
  draw_edge_direct_field(file_x, &nods_[0], nods_.size(2), &edges_[0], edges_.size(2), &X[0]);
  matd_t Y = scale*edge_frm_(colon(3, 5), colon());
  draw_edge_direct_field(file_y, &nods_[0], nods_.size(2), &edges_[0], edges_.size(2), &Y[0]);
  return 0;
}

int wave_constructor::load_feature_line(const char *filename) {
  ifstream inf(filename);
  if ( inf.fail() ) {
    cerr << "[error] can not open file: " << filename << endl;
    return __LINE__;
  }

  size_t feature_num;
  inf >> feature_num;
  cout << "[info] feature chains num: " << feature_num << endl;
  if ( feature_num == 0 )
    return 0;

  std::vector<size_t> feature_edges;
  for (size_t fi = 0; fi < feature_num; ++fi) {
    size_t vnum, prev;
    inf >> vnum;
    if (vnum > 0) inf >> prev;
    for (size_t vi=1; vi < vnum; ++vi) {
      size_t curv;
      inf >> curv;
      feature_edges.push_back(prev);
      feature_edges.push_back(curv);
      prev = curv;
    }
  }
  lines_.resize(2, feature_edges.size()/2);
  std::copy(feature_edges.begin(), feature_edges.end(), lines_.begin());

  inf.close();
  return 0;
}

int wave_constructor::save_feature_to_vtk(const char *filename) const {
  ofstream os(filename);
  if ( os.fail() ) {
    cerr << "[error] can not open " << filename << endl;
    return __LINE__;
  }
  line2vtk(os, &nods_[0], nods_.size(2), &lines_[0], lines_.size(2));
  os.close();
  return 0;
}

int wave_constructor::save_model_to_obj(const char *filename) const {
  return jtf::mesh::save_obj(filename, tris_, nods_);
}

int wave_constructor::save_wave_to_vtk(const char *filename) const {
  if ( f_.size(1) != 4 || f_.size(2) != nods_.size(2) )
    return __LINE__;
  matd_t data = f_(0, colon());
  return draw_vert_value_to_vtk(filename, &nods_[0], nods_.size(2), &tris_[0], tris_.size(2), &data[0]);
}

int wave_constructor::scale_frame_field(const double scale) {
  edge_frm_ *= scale;
  return 0;
}

int wave_constructor::build_frame_on_vert() {
  vector<bool> vis(nods_.size(2));
  std::fill(vis.begin(), vis.end(), false);
  vert_frm_.resize(6, nods_.size(2));
  for (size_t i = 0; i < edges_.size(2); ++i) {
    for (size_t j = 0; j < 2; ++j) {
      if ( !vis[edges_(j, i)] ) {
        vis[edges_(j, i)] = true;
        vert_frm_(colon(), edges_(j, i)) = edge_frm_(colon(), i);
      }
    }
  }
  return 0;
}

int wave_constructor::solve_phase_transition() {
  build_frame_on_vert();
  matd_t alpha(2, edges_.size(2));
  matd_t beta(2, edges_.size(2));
#pragma omp parallel for
  for (size_t i = 0; i < edges_.size(2); ++i) {
    for (size_t j = 0; j < 2; ++j) {
      const size_t I = edges_(j, i);
      const size_t J = edges_(1-j, i);
      matd_t rotx(3, 1), roty(3, 1);
      rotate_frame(&edge_frm_(0, i), &edge_frm_(3, i), &vert_frm_(0, I), &vert_frm_(3, I), &rotx[0], &roty[0]);
      matd_t alignx(3, 1), aligny(3, 1);
      align_vector_field(&rotx[0], &roty[0], &vert_frm_(0, I), &vert_frm_(3, I), &alignx[0], &aligny[0]);
      alpha(j, i) = M_PI*dot(nods_(colon(), J)-nods_(colon(), I), alignx)/dot(alignx, alignx);
      beta(j, i) = M_PI*dot(nods_(colon(), J)-nods_(colon(), I), aligny)/dot(aligny, aligny);
    }
  }
  cIJ_.resize(4, edges_.size(2));
  cJI_.resize(4, edges_.size(2));
#pragma omp parallel for
  for (size_t i = 0; i < edges_.size(2); ++i) {
    double a = alpha(0, i), b = beta(0, i);
    cIJ_(0, i) = cos(a)*cos(b);
    cIJ_(1, i) = -cos(a)*sin(b);
    cIJ_(2, i) = -sin(a)*cos(b);
    cIJ_(3, i) = sin(a)*sin(b);
    a = alpha(1, i); b = beta(1, i);
    cJI_(0, i) = cos(a)*cos(b);
    cJI_(1, i) = -cos(a)*sin(b);
    cJI_(2, i) = -sin(a)*cos(b);
    cJI_(3, i) = sin(a)*sin(b);
  }
  return 0;
}

inline bool is_inflexion(const matd_t &p, const matd_t &front, const matd_t &back) {
  return dot(front-p, back-p)/(norm(front-p)*norm(back-p)) > -0.1;
}

void wave_constructor::count_vert_show_on_feature() {
  vert_count_.clear();
  std::map<size_t, std::set<size_t>> tmp_vert_count;
  for(size_t i = 0; i < lines_.size(2); ++i) {
    tmp_vert_count[lines_(0, i)].insert(lines_(1, i));
    tmp_vert_count[lines_(1, i)].insert(lines_(0, i));
  }
  for(auto i=tmp_vert_count.begin(); i != tmp_vert_count.end(); ++i) {
    const std::set<size_t> &adj_vids = i->second;
    switch (adj_vids.size()) {
      case 0:
        break;
      case 1:
        vert_count_[i->first] = 1;
        break;
      case 2:
        if ( is_inflexion(nods_(colon(), i->first),
                          nods_(colon(), *adj_vids.begin()),
                          nods_(colon(), *adj_vids.rbegin())) )
          vert_count_[i->first] = 2;
        else
          vert_count_[i->first] = 1;
        break;
      default:
        vert_count_[i->first] = 2;
        break;
    }
  }
}

int wave_constructor::prepare() {
#define SOFT_FEATURE
  f_.resize(4, nods_.size(2));
  buff_.resize(4);
  buff_[0] = std::make_shared<wave_value_cons>(edges_, f_, cIJ_, cJI_, 1.0);
  buff_[1] = std::make_shared<modulus_cons>(f_, 0.15);
  buff_[2] = std::make_shared<phase_cons>(f_, 0.15);
  count_vert_show_on_feature();
  if ( !vert_count_.empty() ) {
#ifdef SOFT_FEATURE
//    buff_[3] = std::make_shared<feature_cons>(f_, vert_count_, 1e5);
#else
    feature_cons_ = std::make_shared<feature_cons>(f_, vert_count_);
#endif
  }
  try {
    constraint_ = std::make_shared<constraint_t<double>>(buff_);
  } catch ( exception &e ) {
    cerr << "[exception] " << e.what() << endl;
    exit(EXIT_FAILURE);
  }
  cout << "[info] number of energy term: " << constraint_->Nf() << endl;
  return 0;
}

int wave_constructor::give_an_initial_value(const size_t idx) {
  if ( idx >= f_.size(2) ) {
    cerr << "[error] point ID is out of range\n";
    return __LINE__;
  }
  f_(0, idx) = 1.0;
  return 0;
}

int wave_constructor::solve_wave() {
#ifdef SOFT_FEATURE
  return solve_wave_soft_feature();
#else
  return solve_wave_hard_feature();
#endif
}

int wave_constructor::solve_wave_soft_feature() {
  // Gauss Newton
  const size_t xdim = constraint_->Nx();
  const size_t fdim = constraint_->Nf();
  Map<VectorXd> X(&f_[0], xdim);

  for (size_t iter = 0; iter < 50; ++iter) {
    VectorXd cv = VectorXd::Zero(fdim); {
      constraint_->Val(&X[0], cv.data());
      if ( iter % 1 == 0 ) {
        cout << "\t@energy value: " << cv.squaredNorm() << endl;
      }
      if ( cv.lpNorm<Infinity>() < 1e-8 ) {
        cout << "[info] converged after " << iter << " iteration\n";
        break;
      }
    }
    SparseMatrix<double> J(fdim, xdim); {
      vector<Triplet<double>> trips;
      constraint_->Jac(&X[0], 0, &trips);
      J.reserve(trips.size());
      J.setFromTriplets(trips.begin(), trips.end());
    }
    SparseMatrix<double> LHS = J.transpose()*J;
    VectorXd rhs = -J.transpose()*cv;
    ltl_solver_.compute(LHS);
    ASSERT(ltl_solver_.info() == Success);
    VectorXd dx = ltl_solver_.solve(rhs);
    ASSERT(ltl_solver_.info() == Success);
    X += dx;
  }
  return 0;
}

int wave_constructor::solve_wave_hard_feature() {
  const size_t xdim = constraint_->Nx();
  const size_t fdim = constraint_->Nf();
  const size_t cdim = feature_cons_->Nf();
  Map<VectorXd> X(&f_[0], xdim);
  VectorXd unkown = VectorXd::Zero(xdim+cdim);
  unkown.head(xdim) = X;

  // solve KKT
  cout.precision(15);
  for (size_t iter = 0; iter < 100; ++iter) {
    VectorXd cv = VectorXd::Zero(fdim); {
      constraint_->Val(&unkown[0], cv.data());
      if ( iter % 1 == 0 ) {
        cout << "\t@iter " << iter << endl;
        cout << "\t@energy value: " << cv.squaredNorm() << endl;
      }
    }
    VectorXd fv = VectorXd::Zero(cdim); {
      feature_cons_->Val(&unkown[0], fv.data());
      if ( iter % 1 == 0 )
        cout << "\t@feature constraint: " << fv.lpNorm<Infinity>() << endl << endl;
    }
    SparseMatrix<double> J(fdim, xdim); {
      vector<Triplet<double>> trips;
      constraint_->Jac(&X[0], 0, &trips);
      J.reserve(trips.size());
      J.setFromTriplets(trips.begin(), trips.end());
    }
    SparseMatrix<double> LHS(xdim+cdim, xdim+cdim); {
      vector<Triplet<double>> trips;
      SparseMatrix<double> JTJ = J.transpose()*J;
      extract_triplets_from_spmat<double, ColMajor>(JTJ, &trips);
      vector<Triplet<double>> cons_trips;
      feature_cons_->Jac(nullptr, xdim, &cons_trips);
      for (auto &it : cons_trips) {
        trips.push_back(Triplet<double>(it.row(), it.col(), it.value()));
        trips.push_back(Triplet<double>(it.col(), it.row(), it.value()));
      }
      LHS.reserve(trips.size());
      LHS.setFromTriplets(trips.begin(), trips.end());
    }
    VectorXd rhs = VectorXd::Zero(xdim+cdim); {
      rhs.head(xdim) = -J.transpose()*cv;
      rhs.tail(cdim) = -fv;
    }
    lu_solver_.compute(LHS);
    ASSERT(lu_solver_.info() == Success);
    VectorXd dx = lu_solver_.solve(rhs);
    ASSERT(lu_solver_.info() == Success);
    unkown += dx;
  }
  X = unkown.head(xdim);
  return 0;
}

//==============================================================================
//void wave_constructor::test_wave_conditions() const {
//  const double f1[4] = {1, 2, 3, 4};
//  const double f2[4] = {1, -1, 1, -1};
//  const double c[4] = {1, 1, 1, 1};
// {
//    cout << "#TEST 1:\n";
//    double val = 0;
//    modulus_condition_(&val, f1);
//    cout << "modulus: " << val << endl;
//    phase_condition_(&val, f1);
//    cout << "phase: " << val << endl << endl;
//  }
//  {
//    cout << "#TEST 2:\n";
//    double val = 0;
//    modulus_condition_(&val, f2);
//    cout << "modulus: " << val << endl;
//    phase_condition_(&val, f2);
//    cout << "phase: " << val << endl << endl;
//  }
//  {
//    cout << "#TEST 3:\n";
//    const double f[8] = {1, 2, 3, 4, 5, 6, 7, 9};
//    double val = 0;
//    wave_value_condition_(&val, f, c);
//    cout << "transition: " << val << endl << endl;
//  }
//}

}
