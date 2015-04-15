#ifndef PARAM_LSCM_H
#define PARAM_LSCM_H

#include <zjucad/matrix/matrix.h>
#include <unordered_set>

#include "def.h"

namespace surfparam {

class lscm_param
{
public:
    typedef zjucad::matrix::matrix<size_t> mati_t;
    typedef zjucad::matrix::matrix<double> matd_t;
    lscm_param(const mati_t &tris, const matd_t &nods);
    void set_fixed_bnd_vert(const size_t id, const double *x);
    int apply();
    int get_param_mesh(mati_t *param_tris, matd_t *param_nods);
private:
    std::vector<std::shared_ptr<Functional<double>>> buff_;
    std::shared_ptr<Functional<double>> conformal_energy_;
    const mati_t &tris_;
    std::unordered_set<size_t> fixed_dofs_;
    std::vector<size_t> g2l_;
    Eigen::VectorXd uv_;
};

}
#endif
