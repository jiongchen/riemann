#ifndef NUMERIC_DEF_H
#define NUMERIC_DEF_H

#include <Eigen/Sparse>
#include <memory>

namespace riemann {

template <typename T>
class Functional
{
public:
  virtual ~Functional() {}
  virtual size_t Nx() const = 0;
  virtual int Val(const T *x, T *val) const = 0;
  virtual int Gra(const T *x, T *gra) const = 0;
  virtual int Hes(const T *x, std::vector<Eigen::Triplet<T>> *hes) const = 0;
  virtual void ResetWeight(const double w) {}
  virtual int operator ()(const T *x, T *val, T *gra, const T step, bool graON) { // for LBFGS
    this->Val(x, val);
    if ( graON )
      this->Gra(x, gra);
    return 0;
  }
};

template <typename T>
class Constraint
{
public:
  virtual ~Constraint() {}
  virtual size_t Nx() const = 0;
  virtual size_t Nf() const = 0;
  virtual int Val(const T *x, T *val) const = 0;
  virtual int Jac(const T *x, const size_t off, std::vector<Eigen::Triplet<T>> *jac) const = 0;
  virtual int Hes(const T *x, const size_t off, std::vector<std::vector<Eigen::Triplet<T>>> *hes) const {
    return __LINE__;
  }
};

template <typename T>
class energy_t : public Functional<T>
{
public:
  class null_input_exception : public std::exception {
  public :
    const char* what() const throw() {
      return "null input exception";
    }
  };
  class compatibility_exception : public std::exception {
  public :
    const char* what() const throw() {
      return "compatibility exception";
    }
  };
  energy_t(const std::vector<std::shared_ptr<Functional<T>>> &buffer)
    : buffer_(buffer), dim_(-1) {
    for (auto &e : buffer_) {
      if ( e.get() ) {
        dim_ = e->Nx();
        break;
      }
    }
    if ( dim_ == -1 ) {
      throw null_input_exception();
    }
    for (auto &e : buffer_) {
      if ( e.get() && e->Nx() != dim_ ) {
        throw compatibility_exception();
      }
    }
  }
  size_t Nx() const {
    return dim_;
  }
  int Val(const T *x, T *val) const {
    for (auto &e : buffer_) {
      if ( e.get() ) {
        e->Val(x, val);
      }
    }
    return 0;
  }
  int Gra(const T *x, T *gra) const {
    for (auto &e : buffer_) {
      if ( e.get() ) {
        e->Gra(x, gra);
      }
    }
    return 0;
  }
  int Hes(const T *x, std::vector<Eigen::Triplet<T>> *hes) const {
    for (auto &e : buffer_) {
      if ( e.get() ) {
        e->Hes(x, hes);
      }
    }
    return 0;
  }
protected:
  const std::vector<std::shared_ptr<Functional<T>>> &buffer_;
  size_t dim_;
};

template <typename T>
class constraint_t : public Constraint<T>
{
public:
  class null_input_exception : std::exception {
  public :
    const char* what() const throw() {
      return "null input exception";
    }
  };
  class compatibility_exception : std::exception {
  public :
    const char* what() const throw() {
      return "compatibility exception";
    }
  };
  constraint_t(const std::vector<std::shared_ptr<Constraint<T>>> &buffer)
    : buffer_(buffer), xdim_(-1), fdim_(0) {
    for (auto &e : buffer_) {
      if ( e.get() ) {
        xdim_ = e->Nx();
        break;
      }
    }
    if ( xdim_ == -1 )
      throw null_input_exception();
    bool compatible = true;
    for (auto &c : buffer_) {
      if ( c.get() ) {
        fdim_ += c->Nf();
        if ( c->Nx() != xdim_ )
          compatible = false;
      }
    }
    if ( !compatible )
      throw compatibility_exception();
  }
  size_t Nx() const {
    return xdim_;
  }
  size_t Nf() const {
    size_t fdim = 0;
    for (auto &c : buffer_) {
      if ( c.get() )
        fdim += c->Nf();
    }
    return fdim;
  }
  int Val(const T *x, T *val) const {
    Eigen::Map<Eigen::Matrix<T, -1, 1>> v(val, Nf());
    size_t offset = 0;
    for (auto &c : buffer_) {
      if ( c.get() ) {
        const size_t nf = c->Nf();
        Eigen::Matrix<T, -1, 1> value(nf);
        value.setZero();
        c->Val(x, value.data());
        v.segment(offset, nf) += value;
        offset += nf;
      }
    }
    return 0;
  }
  int Jac(const T *x, const size_t off, std::vector<Eigen::Triplet<T>> *jac) const {
    size_t offset = off;
    for (auto &c : buffer_) {
      if ( c.get() ) {
        c->Jac(x, offset, jac);
        offset += c->Nf();
      }
    }
    return 0;
  }
  int Hes(const T *x, const size_t off, std::vector<std::vector<Eigen::Triplet<T>>> *hes) const {
    if ( hes->size() != fdim_ )
      hes->resize(fdim_);
    size_t offset = 0;
    for (auto &c : buffer_) {
      if ( c.get() ) {
        c->Hes(x, offset, hes);
        offset += c->Nf();
      }
    }
    return 0;
  }
protected :
  const std::vector<std::shared_ptr<Constraint<T>>> &buffer_;
  size_t xdim_, fdim_;
};

}

#endif // NUMERIC_DEF_H
