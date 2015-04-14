#include "energy.h"

#include <jtflib/mesh/mesh.h>

#include "config.h"
#include "cotmatrix.h"

using namespace std;
using namespace Eigen;
using namespace jtf::mesh;

namespace surfparam {

dirichlet_energy::dirichlet_energy(const mati_t &tris, const matd_t &nods, const double w)
    : dim_(2*nods.size(2)), w_(w) {
    cotmatrix(tris, nods, 2, &L_);
    //
}

size_t dirichlet_energy::Nx() const {
    return dim_;
}

int dirichlet_energy::Val(const double *x, double *val) const {
    Map<const VectorXd> X(x, dim_);
    *val += 0.5 * w_ * X.dot(L_ * X);
    return dim_;
}

int dirichlet_energy::Gra(const double *x, double *gra) const {
    Map<const VectorXd> X(x, dim_);
    Map<VectorXd> grad(gra, dim_);
    grad += w_ * L_ * X;
    return 0;
}

int dirichlet_energy::Hes(const double *x, vector<Triplet<double>> *hes) const {
    for (size_t j = 0; j < L_.outerSize(); ++j)
        for (SparseMatrix<double>::InnerIterator it(L_, j); it; ++it)
            hes->push_back(Triplet<double>(it.row(), it.col(), it.value()));
    return 0;
}
//==============================================================================
param_area::param_area(const mati_t &tris, const matd_t &nods, const double w)
    : dim_(2*nods.size(2)), w_(w) {
    mati_t bnd_edge;
    int no_boundary = GetBoundaryEdge(tris, bnd_edge);
    ASSERT(no_boundary == 0);   // for disk topology
    vector<Triplet<double>> trips;
    for (size_t i = 0; i < bnd_edge.size(2); ++i) {
        size_t pi = bnd_edge(0, i);
        size_t qi = bnd_edge(1, i);
        trips.push_back(Triplet<double>(2*pi+0, 2*qi+1, 0.5));
        trips.push_back(Triplet<double>(2*pi+1, 2*qi+0, -0.5));
        trips.push_back(Triplet<double>(2*qi+0, 2*pi+1, -0.5));
        trips.push_back(Triplet<double>(2*qi+1, 2*pi+0, 0.5));
    }
    A_.resize(dim_, dim_);
    A_.reserve(trips.size());
    A_.setFromTriplets(trips.begin(), trips.end());
}

size_t param_area::Nx() const {
    return dim_;
}

int param_area::Val(const double *x, double *val) const {
    Map<const VectorXd> X(x, dim_);
    *val += 0.5 * w_ * X.dot(A_ * X);
    return 0;
}

int param_area::Gra(const double *x, double *gra) const {
    Map<const VectorXd> X(x, dim_);
    Map<VectorXd> grad(gra, dim_);
    grad += w_ * A_ * X;
    return 0;
}

int param_area::Hes(const double *x, vector<Triplet<double>> *hes) const {
    for (size_t j = 0; j < A_.outerSize(); ++j)
        for (SparseMatrix<double>::InnerIterator it(A_, j); it; ++it)
            hes->push_back(Triplet<double>(it.row(), it.col(), w_*it.value()));
    return 0;
}

int param_area::GetBoundaryEdge(const mati_t &tris, mati_t &bnd_edge) {
    vector<size_t> edge_buffer;
    shared_ptr<edge2cell_adjacent> e2c(edge2cell_adjacent::create(tris, false));
    for (size_t i = 0; i < tris.size(2); ++i) {
        for (size_t j = 0; j < 3; ++j) {
            size_t pi = tris(j, i);
            size_t qi = tris((j+1)%3, i);
            pair<size_t, size_t> fa = e2c->query(pi, qi);
            if ( e2c->is_boundary_edge(fa) ) {
                edge_buffer.push_back(pi);
                edge_buffer.push_back(qi);
            }
        }
    }
    if ( edge_buffer.size() == 0 )
        return __LINE__;
    bnd_edge.resize(2, edge_buffer.size()/2);
    copy(edge_buffer.begin(), edge_buffer.end(), bnd_edge.begin());
    return 0;
}

}
