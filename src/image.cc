#include "image.h"

#include <algorithm>
#include <cmath>

using namespace std;

namespace riemann {

/// ---------------- 2D grid --------------------
/// (0,h-1)->(1,h-1)->(2,h-1)-> ... -> (w-1,h-1)-
///  ...   ->  ...  ->  ...  -> ... ->   ...    -
/// (0,2)  -> (1,2) -> (2,2) -> ... -> (w-1,2)  -
/// (0,1)  -> (1,1) -> (2,1) -> ... -> (w-1,1)  -
/// (0,0)  -> (1,0) -> (2,0) -> ... -> (w-1,0)  -
/// ---------------------------------------------

struct tga_header {

};

image_t::~image_t() {
  if ( pixels_ ) free(pixels_);
}

int image_t::read(const char *filename) {
  return 0;
}

double& image_t::operator ()(const int x, const int y) {
  return pixels_[y*w_+x];
}

double image_t::operator ()(const int x, const int y) const {
  return pixels_[y*w_+x];
}

double image_t::sample(const double x, const double y) const {
  // bilinear interpolation
  double ax = x-std::floor(x), ay = y-std::floor(y);
  double bx = 1.0-ax, by = 1.0-ay;
  int x0 = (int)std::floor(x), y0 = (int)std::floor(y);
  int x1 = x0+1, y1 = y0+1;
  clamp(x0, y0);
  clamp(x1, y1);
  return by*(bx*(*this)(x0, y0)+ax*(*this)(x0, y1))+
         ay*(bx*(*this)(x0, y1)+ax*(*this)(x1, y1));
}

void image_t::clamp(int &x, int &y) const {
  x = std::max(0, std::min(w_-1, x));
  y = std::max(0, std::min(h_-1, y));
}

}
