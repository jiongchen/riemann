#include "vert_local_frame.h"

#include <jtflib/mesh/util.h>

using namespace std;
using namespace zjucad::matrix;
using mati_t=zjucad::matrix::matrix<size_t>;
using matd_t=zjucad::matrix::matrix<double>;

namespace riemann {

int calc_vert_local_frame(const mati_t &tris, const matd_t &nods, matd_t &frame) {
  frame.resize(9, nods.size(2));
  {
    matd_t normal;
    jtf::mesh::cal_point_normal(tris, nods, normal);
    frame(colon(6, 8), colon()) = normal;
  }
#pragma omp parallel for
  for (size_t i = 0; i < nods.size(2); ++i) {
    do {
      matd_t rand_vec = rand<double>(3, 1);
      frame(colon(0, 2), i) = cross(frame(colon(6, 8), i), rand_vec);
    } while ( norm(frame(colon(0, 2), i)) < 1e-8 );
    frame(colon(3, 5), i) = cross(frame(colon(6, 8), i), frame(colon(0, 2), i));
    frame(colon(0, 2), i) /= norm(frame(colon(0, 2), i));
    frame(colon(3, 5), i) /= norm(frame(colon(3, 5), i));
    frame(colon(6, 8), i) /= norm(frame(colon(6, 8), i));
  }
  return 0;
}

}
