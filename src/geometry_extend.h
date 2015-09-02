#ifndef GEOMETRY_EXTEND_H
#define GEOMETRY_EXTEND_H

#include <zjucad/matrix/matrix.h>

namespace riemann {

void gen_rand_orth_vec_3d(const double *x, double *xT);

int calc_vert_local_frame(const zjucad::matrix::matrix<size_t> &tris,
                          const zjucad::matrix::matrix<double> &nods,
                          zjucad::matrix::matrix<double> &frame);

void project_point_on_plane(const double *x, const double *origin,
                            const double *tan0, const double *tan1,
                            double *proj_x);
}

#endif
