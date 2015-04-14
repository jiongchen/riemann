#ifndef PARAM_LSCM_H
#define PARAM_LSCM_H

#include <zjucad/matrix/matrix.h>
#include "def.h"

namespace surfparam {

class lscm_param
{
public:
    typedef zjucad::matrix::matrix<size_t> mati_t;
    typedef zjucad::matrix::matrix<double> matd_t;
    lscm_param(const mati_t &tris, const matd_t &nods);
    int apply();
    int get_param_mesh(mati_t *param_tris, matd_t *param_nods);
private:
    std::vector<Functional<double>> buff_;
    std::shared_ptr<Functional<double>> energy_;
    const mati_t &tris_;
    Eigen::VectorXd uv_;
};

}
#endif
