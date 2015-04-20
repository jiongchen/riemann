#include "green_coord_deform.h"

#include <zjucad/matrix/matrix.h>
#include <jtflib/mesh/io.h>
#include <jtflib/mesh/util.h>
#include <zjucad/matrix/io.h>
#include <boost/math/special_functions/sign.hpp>

#include "config.h"
#include "vtk.h"

using namespace std;
using namespace zjucad::matrix;
using namespace Eigen;
using boost::math::sign;

namespace geom_deform {

green_deform_2d::green_deform_2d() {}

int green_deform_2d::load_sample_points(const char *file) {
    matrix<size_t> cell;
    matrix<double> nods;
    if ( jtf::mesh::load_obj(file, cell, nods) ) {
        return __LINE__;
    }
    cell_.resize(cell.size(1), cell.size(2));
    nods_.resize(2, nods.size(2));
    std::copy(cell.begin(), cell.end(), cell_.data());
#pragma omp parallel for
    for (size_t i = 0; i < nods.size(2); ++i) {
        nods_(0, i) = nods(0, i);
        nods_(1, i) = nods(2, i);
    }
    return 0;
}

int green_deform_2d::load_cage(const char *file) {
    // parse manually
    ifstream is(file);
    if ( is.fail() ) {
        cerr << "# can not open " << file << "\n";
        return __LINE__;
    }
    size_t ele_num;
    string CELL;
    is >> CELL >> ele_num;
    cage_cell_.resize(2, ele_num);
    for (size_t i = 0; i < ele_num; ++i)
        is >> cage_cell_(0, i) >> cage_cell_(1, i);

    size_t nods_num;
    string NODES;
    is >> NODES >> nods_num;
    cage_nods_.resize(2, nods_num);
    for (size_t i = 0; i < nods_num; ++i)
        is >> cage_nods_(0, i) >> cage_nods_(1, i);

    // get rest segment element length
    rest_len_.resize(cage_cell_.cols());
    curr_len_.resize(cage_cell_.cols());
    calc_cage_edge_length(rest_len_);

    // calc initial normal
    cage_normal_.resize(2, cage_cell_.cols());
    calc_outward_normal();
    return 0;
}

int green_deform_2d::calc_cage_edge_length(VectorXd &len) {
#pragma omp parallel for
    for (size_t i = 0; i < cage_cell_.cols(); ++i)
        len[i] = (cage_nods_.col(cage_cell_(0, i))-cage_nods_.col(cage_cell_(1, i))).norm();
    return 0;
}

int green_deform_2d::calc_outward_normal() {
#pragma omp parallel for
    for (size_t i = 0; i < cage_cell_.cols(); ++i) {
        Vector2d dir = cage_nods_.col(cage_cell_(1, i)) - cage_nods_.col(cage_cell_(0, i));
        dir /= dir.norm();
        cage_normal_(0, i) = -dir[1];
        cage_normal_(1, i) = dir[0];
    }
    return 0;
}

int green_deform_2d::calc_green_coords() {
    phi_ = MatrixXd::Zero(cage_nods_.cols(), nods_.cols());
    psi_ = MatrixXd::Zero(cage_normal_.cols(), nods_.cols());
#pragma omp parallel for
    for (size_t pid = 0; pid < nods_.cols(); ++pid) {
        for (size_t i = 0; i < cage_cell_.cols(); ++i) {
            Vector2d a = cage_nods_.col(cage_cell_(1, i))-cage_nods_.col(cage_cell_(0, i));
            Vector2d b = cage_nods_.col(cage_cell_(0, i))-nods_.col(pid);
            double Q, S, R, BA, SRT, L0, L1, A0, A1, A10, L10, psi_entry, phi_entry1, phi_entry0;
            Q = a.dot(a);
            S = b.dot(b);
            R = 2*a.dot(b);
            BA = b.dot(a.norm()*cage_normal_.col(i));
            SRT = sqrt(4*S*Q-R*R);
            L0 = log(S);
            L1 = log(S+Q+R);
            A0 = atan2(R, SRT)/SRT;
            A1 = atan2(2*Q+R, SRT)/SRT;
            A10 = A1 - A0;
            L10 = L1 - L0;
            psi_entry = -a.norm()/(4*M_PI)*((4*S-R*R/Q)*A10 + R/(2*Q)*L10 + L1 - 2);
            phi_entry1 = -BA/(2*M_PI)*(L10/(2*Q) - A10*R/Q);
            phi_entry0 = +BA/(2*M_PI)*(L10/(2*Q) - A10*(2+R/Q));
#pragma omp critical
            {
                psi_(i, pid) += psi_entry;
                phi_(cage_cell_(1, i), pid) += phi_entry1;
                phi_(cage_cell_(0, i), pid) += phi_entry0;
            }
        }
    }
    ///> =_=#
    for (size_t j = 0; j < phi_.cols(); ++j) {
        double sum = phi_.col(j).sum();
        phi_.col(j) /= sum;
    }
    return 0;
}

int green_deform_2d::move_cage(const size_t id, const double *dx, bool disp) {
    Vector2d X(dx[0], dx[1]);
    if ( disp )
        cage_nods_.col(id) += X;
    else
        cage_nods_.col(id) = X;
    return 0;
}

int green_deform_2d::deform() {
    calc_outward_normal();
    calc_cage_edge_length(curr_len_);
    VectorXd ratio = curr_len_.cwiseQuotient(rest_len_);
    ASSERT(cage_normal_.cols() == ratio.rows());
    nods_ = cage_nods_*phi_ + (cage_normal_*ratio.asDiagonal())*psi_;
    return 0;
}
//==============================================================================
int green_deform_2d::dump(const char *file) {
    MatrixXd nods_3d(3, nods_.cols());
#pragma omp parallel for
    for (size_t i = 0; i < nods_.cols(); ++i) {
        nods_3d(0, i) = nods_(0, i);
        nods_3d(1, i) = 0.0;
        nods_3d(2, i) = nods_(1, i);
    }
    ofstream os(file);
    tri2vtk(os, nods_3d.data(), nods_3d.cols(), cell_.data(), cell_.cols());
    return 0;
}

int green_deform_2d::dump_cage(const char *file) {
    MatrixXd cage_nods_3d(3, cage_nods_.cols());
#pragma omp parallel for
    for (size_t i = 0; i < cage_nods_.cols(); ++i) {
        cage_nods_3d(0, i) = cage_nods_(0, i);
        cage_nods_3d(1, i) = 0.0;
        cage_nods_3d(2, i) = cage_nods_(1, i);
    }
    ofstream os(file);
    line2vtk(os, cage_nods_3d.data(), cage_nods_3d.cols(), cage_cell_.data(), cage_cell_.cols());
    return 0;
}

int green_deform_2d::dump_normal(const char *file) {
    matrix<size_t> normal_cell(2, cage_cell_.cols());
    for (size_t i = 0; i < normal_cell.size(); ++i)
        normal_cell[i] = i;
    MatrixXd normal_nods(3, 2 * normal_cell.size(2));
    normal_nods.setZero();
#pragma omp parallel for
    for (size_t i = 0; i < normal_nods.cols()/2; ++i) {
        Vector2d mid = 0.5 * (cage_nods_.col(cage_cell_(0, i)) + cage_nods_.col(cage_cell_(1, i)));
        Vector2d end = mid + cage_normal_.col(i);
        normal_nods(0, 2*i+0) = mid[0];
        normal_nods(2, 2*i+0) = mid[1];
        normal_nods(0, 2*i+1) = end[0];
        normal_nods(2, 2*i+1) = end[1];
    }
    ofstream os(file);
    line2vtk(os, normal_nods.data(), normal_nods.cols(), normal_cell.begin(), normal_cell.size(2));
    return 0;
}
//==============================================================================
green_deform_3d::green_deform_3d() {}

int green_deform_3d::load_sample_points(const char *file) {
    if ( jtf::mesh::load_obj(file, cell_, nods_) )
        return __LINE__;
    return 0;
}

int green_deform_3d::load_cage(const char *file) {
    if ( jtf::mesh::load_obj(file, cage_cell_, cage_nods_) )
        return __LINE__;
    // compute initial outward normal
    cage_normal_.resize(3, cage_cell_.size(2));
    calc_outward_normal();
    // compute initial face area
    rest_area_.resize(1, cage_cell_.size(2));
#pragma omp parallel for
    for (size_t i = 0; i < cage_cell_.size(2); ++i) {
        matd_t vt = cage_nods_(colon(), cage_cell_(colon(), i));
        rest_area_[i] = jtf::mesh::cal_face_area(vt);
    }
    // compute cage edge
    cage_rest_uv_.resize(3, 2*cage_cell_.size(2));
    cage_curr_uv_.resize(3, 2*cage_cell_.size(2));
    calc_cage_edge(cage_rest_uv_);
    return 0;
}

int green_deform_3d::calc_cage_edge(matd_t &uv) {
    if ( uv.size(1) != 3 || uv.size(2) != 2*cage_cell_.size(2) )
        uv.resize(3, 2*cage_cell_.size(2));
#pragma omp parallel for
    for (size_t i = 0; i < cage_cell_.size(2); ++i) {
        matd_t VJ = cage_nods_(colon(), cage_cell_(colon(), i));
        uv(colon(), colon(2*i, 2*i+1)) = VJ(colon(), colon(1, 2))-VJ(colon(), colon(0, 1));
    }
    return 0;
}

int green_deform_3d::calc_stretch_ratio(matd_t &s) {
    s.resize(cage_normal_.size(2), 1);
    calc_cage_edge(cage_curr_uv_);
    ASSERT(cage_curr_uv_.size(2) == 2*s.size());
#pragma omp parallel for
    for (size_t i = 0; i < s.size(); ++i) {
        matd_t u = cage_rest_uv_(colon(), 2*i+0);
        matd_t v = cage_rest_uv_(colon(), 2*i+1);
        matd_t ucap = cage_curr_uv_(colon(), 2*i+0);
        matd_t vcap = cage_curr_uv_(colon(), 2*i+1);
        s[i] = sqrt(dot(ucap, ucap)*dot(v, v)-2*dot(ucap, vcap)*dot(u,v)+dot(vcap, vcap)*dot(u, u))
                / (sqrt(8)*rest_area_[i]);
    }
    return 0;
}

int green_deform_3d::calc_outward_normal() {
    jtf::mesh::cal_face_normal(cage_cell_, cage_nods_, cage_normal_, true);
    return 0;
}

int green_deform_3d::dump(const char *file) {
    ofstream os(file);
    tri2vtk(os, &nods_[0], nods_.size(2), &cell_[0], cell_.size(2));
    return 0;
}

int green_deform_3d::dump_cage(const char *file) {
    ofstream os(file);
    tri2vtk(os, &cage_nods_[0], cage_nods_.size(2), &cage_cell_[0], cage_cell_.size(2));
    return 0;
}

int green_deform_3d::dump_normal(const char *file) {
    mati_t normal_cell(2, cage_cell_.size(2));
    matd_t normal_nods(3, 2*cage_cell_.size(2));
#pragma omp parallel for
    for (size_t i = 0; i < normal_cell.size(2); ++i) {
        normal_cell(0, i) = 2*i+0;
        normal_cell(1, i) = 2*i+1;
        normal_nods(colon(), normal_cell(0, i)) = 1.0/3*cage_nods_(colon(), cage_cell_(colon(), i))*ones<double>(3, 1);
        normal_nods(colon(), normal_cell(1, i)) = normal_nods(colon(), normal_cell(0, i))+cage_normal_(colon(), i);
    }
    ofstream os(file);
    line2vtk(os, &normal_nods[0], normal_nods.size(2), &normal_cell[0], normal_cell.size(2));
    return 0;
}
//==============================================================================
double green_deform_3d::GCTriInt(const matd_t &p, const matd_t &v1,
                                 const matd_t &v2, const matd_t &eta) {
    double alpha = acos(std::min(1.0, std::max(-1.0, dot(v2-v1, p-v1)/norm(v2-v1)/norm(p-v1))));
    double beta = acos(std::min(1.0, std::max(-1.0, dot(v1-p, v2-p)/norm(v1-p)/norm(v2-p))));
    double lambda = dot(p-v1, p-v1)*sin(alpha)*sin(alpha);
    double c = dot(p-eta, p-eta);

    double theta[2] = {M_PI-alpha, M_PI-alpha-beta};
    double I[2];
    double sqrt_c = sqrt(c);
    double sqrt_lambda = sqrt(lambda);
    for (int i = 0; i < 2; ++i) {
        double S = sin(theta[i]);
        double C = cos(theta[i]);
        I[i] = -sign(S)/2.0*(2.0*sqrt_c*atan2(sqrt_c*C, sqrt(lambda+S*S*c))
                            + sqrt_lambda*log(2.0*sqrt_lambda*S*S/((1-C)*(1-C))
                                              *(1.0-2*c*C/(c+c*C+lambda+sqrt(lambda*lambda+lambda*c*S*S)))));
    }
    return -1.0/(4*M_PI)*fabs(I[0]-I[1]-sqrt_c*beta);
}

int green_deform_3d::calc_green_coords() {
    matd_t zero3d = zeros<double>(3, 1);
    phi_ = zeros<double>(cage_nods_.size(2), nods_.size(2));
    psi_ = zeros<double>(cage_normal_.size(2), nods_.size(2));
    for (size_t pid = 0; pid < nods_.size(2); ++pid) {
        for (size_t i = 0; i < cage_cell_.size(2); ++i) {
            matd_t VJ = cage_nods_(colon(), cage_cell_(colon(), i))-nods_(colon(), pid)*ones<double>(1, 3);
            matd_t p = dot(VJ(colon(), 0), cage_normal_(colon(), i))*cage_normal_(colon(), i);
            matd_t s(3, 1), I(3, 1), II(3, 1), N(3, 3);
            for (size_t j = 0; j < 3; ++j) {
                s[j] = sign(dot(cross(VJ(colon(), j)-p, VJ(colon(), (j+1)%3)-p), cage_normal_(colon(), i)));
                I[j] = GCTriInt(p, VJ(colon(), j), VJ(colon(), (j+1)%3), zero3d);
                II[j] = GCTriInt(zero3d, VJ(colon(), (j+1)%3), VJ(colon(), j), zero3d);
                matd_t q = cross(VJ(colon(), (j+1)%3), VJ(colon(), j));
                N(colon(), j) = q/norm(q);
            }
            double I_ = -fabs(dot(s, I));
            psi_(i, pid) += -I_;
            matd_t w = I_*cage_normal_(colon(), i) + N * II;
            if ( norm(w) > 1e-8 ) {
                for (size_t j = 0; j < 3; ++j)
                    phi_(cage_cell_(j, i), pid) += dot(N(colon(), (j+1)%3), w)/dot(N(colon(), (j+1)%3), VJ(colon(), j));
            }
        }
    }
    ///> for translation invariance
#pragma omp parallel for
    for (size_t j = 0; j < phi_.size(2); ++j) {
        double col_sum = sum(phi_(colon(), j));
        phi_(colon(), j) /= col_sum;
    }
    return 0;
}

int green_deform_3d::move_cage(const size_t id, const double *dx, bool disp) {
    itr_matrix<const double *> X(3, 1, dx);
    if ( disp )
        cage_nods_(colon(), id) += X;
    else
        cage_nods_(colon(), id) = X;
    return 0;
}

int green_deform_3d::deform() {
    calc_outward_normal();
    matd_t s;
    calc_stretch_ratio(s);
    ASSERT(s.size() == cage_normal_.size(2));
    matd_t stretch_psi = psi_;
    for (size_t row = 0; row < stretch_psi.size(1); ++row)
        stretch_psi(row, colon()) *= s[row];
    nods_ = cage_nods_*phi_ + cage_normal_*psi_;
    return 0;
}

}
