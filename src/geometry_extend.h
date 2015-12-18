#ifndef GEOMETRY_EXTEND_H
#define GEOMETRY_EXTEND_H

#include <zjucad/matrix/matrix.h>

namespace riemann {

void gen_rand_orth_vec_3d(const double *x, double *xT);

int calc_vert_local_frame(const zjucad::matrix::matrix<size_t> &tris,
                          const zjucad::matrix::matrix<double> &nods,
                          zjucad::matrix::matrix<double> &frame);

void project_vector_on_subspace(const double *vect, const size_t dim,
                                const double *basis, const size_t sub_dim,
                                double *proj_vect);

void project_point_on_subspace(const double *x, const size_t dim,
                               const double *origin, const double *basis, const size_t sub_dim,
                               double *proj_x);

void project_point_on_plane(const double *x, const double *origin,
                            const double *tan0, const double *tan1,
                            double *proj_x);

void rotate_and_scale_triangle(const double *x, const double *quaternion, const double s, const double *xnew);

int calc_one_ring_face(const zjucad::matrix::matrix<size_t> &tris,
                       std::vector<std::vector<size_t>> &p2f);
}

#endif
