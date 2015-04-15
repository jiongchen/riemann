#ifndef PARAM_ENERGY_H
#define PARAM_ENERGY_H

#include <zjucad/matrix/matrix.h>
#include "def.h"

namespace surfparam {

class dirichlet_energy : public Functional<double>
{
public:
    typedef zjucad::matrix::matrix<size_t> mati_t;
    typedef zjucad::matrix::matrix<double> matd_t;
    dirichlet_energy(const mati_t &tris, const matd_t &nods, const double w=1.0);
    size_t Nx() const;
    int Val(const double *x, double *val) const;
    int Gra(const double *x, double *gra) const;
    int Hes(const double *x, std::vector<Eigen::Triplet<double>> *hes) const;
private:
    const size_t dim_;
    const double w_;
    Eigen::SparseMatrix<double> L_;
};

class param_area : public Functional<double>
{
public:
    typedef zjucad::matrix::matrix<size_t> mati_t;
    typedef zjucad::matrix::matrix<double> matd_t;
    param_area(const mati_t &tris, const matd_t &nods, const double w=1.0);
    size_t Nx() const;
    int Val(const double *x, double *val) const;
    int Gra(const double *x, double *gra) const;
    int Hes(const double *x, std::vector<Eigen::Triplet<double>> *hes) const;
public:
    int GetBoundaryEdge(const mati_t &tris, mati_t &bnd_edge);
    const size_t dim_;
    const double w_;
    Eigen::SparseMatrix<double> A_;
};

}

#endif
