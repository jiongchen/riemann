#ifndef VERT_LOCAL_FRAME_H
#define VERT_LOCAL_FRAME_H

#include <zjucad/matrix/matrix.h>

namespace riemann {

int calc_vert_local_frame(const zjucad::matrix::matrix<size_t> &tris,
                          const zjucad::matrix::matrix<double> &nods,
                          zjucad::matrix::matrix<double> &frame);

}

#endif
