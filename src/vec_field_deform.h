#ifndef VEC_FIELD_DEFORM_H
#define VEC_FIELD_DEFORM_H

#include <memory>
#include <Eigen/Dense>
#include <Eigen/Sparse>

namespace geom_deform {

static double b(const double r, const double ri, const double ro) {
    return 4*std::pow(r-ri, 3)*(1-(r-ri)/(ro-ri))/std::pow(ro-ri, 3) + std::pow(r-ri, 4)/std::pow(ro-ri, 4);
}

static double  db_dr(const double r, const double ri, const double ro) {
    return 12*std::pow(r-ri, 2)*(1-(r-ri)/(ro-ri))/std::pow(ro-ri, 3);
}

typedef Eigen::Vector3d Vec3;

class scalar_field
{
public:
    virtual ~scalar_field() {}
    virtual int eval_val(const double *x, double *val) const {}
    virtual int eval_gra(const double *x, double *gra) const {}
    virtual void set_direct(const double *u) {}
    virtual void set_center(const double *c) {}
};

class sphere_distace_field : public scalar_field
{
public:
    int eval_val(const double *x, double *val) const {

    }
    int eval_gra(const double *x, double *gra) const {

    }
    void set_center(const double *c) { c_ = Vec3(c[0], c[1], c[2]); }
private:
    Vec3 c_;
};

class linear_scalar_field : public scalar_field
{
public:
    int eval_val(const double *x, double *val) const {
        Vec3 X(x[0], x[1], x[2]);
        *val = u_.dot(X-c_);
        return 0;
    }
    int eval_gra(const double *x, double *gra) const {
        Vec3 X(x[0], x[1], x[2]);
        Eigen::Map<Eigen::VectorXd>(gra, 3) = u_;
        return 0;
    }
    void set_direct(const double *u) { u_ = Vec3(u[0], u[1], u[2]); }
    void set_center(const double *c) { c_ = Vec3(c[0], c[1], c[2]); }
private:
    Vec3 u_, c_;
};

class quadratic_scalar_field : public scalar_field
{
public:
    int eval_val(const double *x, double *val) const {
        Vec3 X(x[0], x[1], x[2]);
        *val = (u_.cross(X-c_)).squaredNorm();
        return 0;
    }
    int eval_gra(const double *x, double *gra) const {

    }
    void set_direct(const double *u) { u_ = Vec3(u[0], u[1], u[2]); }
    void set_center(const double *c) { c_ = Vec3(c[0], c[1], c[2]); }
private:
    Vec3 u_, c_;
};

class vector_field {
public:
    vector_field() {

    }
    virtual void set_trans_direc(const double *v) {
        Vec3 u, w;
        // uv = uw = wv = 0;
        e_->set_direct();
        f_->set_direct();
    }
    virtual void set_center(const double *x) {
        c_ = Vec3(x[0], x[1], x[2]);
        r_->set_center();
    }
    virtual void set_range(const double ri, const double ro) { ri_ = ri; ro_ = ro; }
    virtual int eval_val(const double *x, double *v) const {
        double dist = 0;
        r_->eval_val(x, &dist);
        Vec3 nabla_e, nabla_f;
        if ( dist < ri_ ) {

        } else if ( ri_ <= dist && dist < ro_ ) {

        } else {

        }
        return 0;
    }
private:
    double ri_, ro_;
    Vec3 c_;
    /// for translation
    std::shared_ptr<scalar_field> r_, e_, f_;
};

}
#endif
