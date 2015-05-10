#ifndef GRADIENT_BASED_DEFORM_H
#define GRADIENT_BASED_DEFORM_H

#include <Eigen/Sparse>
#include <zjucad/matrix/matrix.h>

namespace geom_deform {

/// Lp = div w;
class gradient_field_deform
{
public:
    typedef zjucad::matrix::matrix<size_t> mati_t;
    typedef zjucad::matrix::matrix<double> matd_t;
    gradient_field_deform(const mati_t &tris, const matd_t &nods);
private:
    int compute_divergence(double *div);
    const mati_t tris_;
    const matd_t nods_;
    Eigen::SparseMatrix<double> L_;
};

}
#endif
