#ifndef IMAGE_H
#define IMAGE_H

namespace riemann {

class image_t
{
public:
  ~image_t();
  int read(const char *filename);
  double operator()(const int x, const int y) const;
  double& operator()(const int x, const int y);
  double sample(const double x, const double y) const;
  int width() const { return w_; }
  int height() const { return h_; }
private:
  void clamp(int &x, int &y) const;
  double *pixels_;
  int w_, h_;
};

}

#endif
