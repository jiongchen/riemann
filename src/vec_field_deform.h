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

static void get_orthonormal_basis(const double *v, double *u, double *w) {

}

extern "C" {

void quad_scalar_field_(double *val, const double *x, const double *a, const double *c);
void quad_scalar_field_jac_(double *jac, const double *x, const double *a, const double *c);

}

typedef Eigen::Vector3d Vec3;

class scalar_field
{
public:
    virtual ~scalar_field() {}
    scalar_field(const Vec3 &a, const Vec3 &c) : a_(a), c_(c) {}
    virtual int eval_val(const double *x, double *val) const {}
    virtual int eval_gra(const double *x, double *gra) const {}
    virtual int set_axis(const double *a) {}
    virtual void set_center(const double *c) {}
protected:
    Vec3 a_, c_;
};

class linear_scalar_field : public scalar_field
{
public:
//    linear_scalar_field(const Vec3 &a, const Vec3 &c);
//    int eval_val(const double *x, double *val) const {
//        Vec3 X(x[0], x[1], x[2]);
//        *val = u_.dot(X-c_);
//        return 0;
//    }
//    int eval_gra(const double *x, double *gra) const {
//        Vec3 X(x[0], x[1], x[2]);
//        Eigen::Map<Eigen::VectorXd>(gra, 3) = u_;
//        return 0;
//    }
//    void set_direct(const double *u) {
//        u_ = Vec3(u[0], u[1], u[2]);
//    }
//    void set_center(const double *c) {
//        c_ = Vec3(c[0], c[1], c[2]);
//    }
};

class quadratic_scalar_field : public scalar_field
{
public:
    quadratic_scalar_field(const Vec3 &a, const Vec3 &c);
    int eval_val(const double *x, double *val) const {
        quad_scalar_field_(val, x, &a_[0], &c_[0]);
        return 0;
    }
    int eval_gra(const double *x, double *gra) const {
        quad_scalar_field_jac_(gra, x, &a_[0], &c_[0]);
        return 0;
    }
//    void set_direct(const double *u) {
//        u_ = Vec3(u[0], u[1], u[2]);
//    }
//    void set_center(const double *c) {
//        c_ = Vec3(c[0], c[1], c[2]);
//    }
};

class implicit_tool
{
public:
    virtual ~implicit_tool() {}
    implicit_tool(const Vec3 &c, const double ri, const double ro)
        : c_(c), ri_(ri), ro_(ro) {}
    virtual int eval_val(const double *x, double *val) const {}
    virtual int eval_gra(const double *x, double *gra) const {}
    virtual bool inner(const double *x) const {}
    virtual bool outer(const double *x) const {}
    virtual bool inter(const double *x) const {}
    virtual int set_center(const double *c) {}
    virtual int set_range(const double ri, const double ro) {}
protected:
    Vec3 c_;
    double ri_, ro_;
};

class sphere_tool : public implicit_tool
{
public:
    sphere_tool(const Vec3 &c, const double ri, const double ro)
        : implicit_tool(c, ri, ro) {}
    int eval_val(const double *x, double *val) const;
    int eval_gra(const double *x, double *gra) const;
    virtual bool inner(const double *x) const;
    virtual bool outer(const double *x) const;
    virtual bool inter(const double *x) const;
    virtual int set_center(const double *c);
    virtual int set_range(const double ri, const double ro);
};


class vector_field
{
public:
//    vector_field() {
//        r_ = std::make_shared<sphere_distace_field>();
//        e_ = std::make_shared<linear_scalar_field>();
//        f_ = std::make_shared<linear_scalar_field>();
//    }
//    virtual void set_center(const double *x) {
//        c_ = Vec3(x[0], x[1], x[2]);
//        e_->set_center(&c_[0]);
//        f_->set_center(&c_[0]);
//        r_->set_center(&c_[0]);
//    }
//    virtual void set_range(const double ri, const double ro) {
//        ri_ = ri;
//        ro_ = ro;
//    }
//    virtual void set_trans_direct(const double *v) {
//        Vec3 u, w;
//        get_orthonormal_basis(v, &u[0], &w[0]);
//        e_->set_direct(&u[0]);
//        f_->set_direct(&w[0]);
//    }
//    virtual int eval_val(const double *x, double *v) const {
//        double dist = 0;
//        r_->eval_val(x, &dist);
//        Vec3 nabla_e, nabla_f;
//        if ( dist < ri_ ) {

//        } else if ( ri_ <= dist && dist < ro_ ) {

//        } else {

//        }
//        return 0;
//    }
private:
    std::shared_ptr<implicit_tool> tools_;
    std::shared_ptr<scalar_field> e_, f_;
};

class deformer
{
public:
    deformer();
    int load_model(const char *file);
    int set_translate_tools();
    int advect();
    int save_model(const char *file);
private:
    std::vector<std::shared_ptr<vector_field>> vf_;
};

}
#endif
