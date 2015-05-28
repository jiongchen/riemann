#ifndef VEC_FIELD_DEFORM_H
#define VEC_FIELD_DEFORM_H

#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <zjucad/matrix/matrix.h>

namespace geom_deform {

double eval_blend_val(const double r, const double ri, const double ro);
double eval_blend_jac(const double r, const double ri, const double ro);
int get_ortn_basis(const double *v, double *u, double *w);

extern "C" {
void quad_scalar_field_(double *val, const double *x, const double *a, const double *c);
void quad_scalar_field_jac_(double *jac, const double *x, const double *a, const double *c);
}

typedef Eigen::Vector3d Vec3;

//==============================================================================
class scalar_field
{
public:
    virtual ~scalar_field() {}
    scalar_field(const Vec3 &a, const Vec3 &c) : a_(a), c_(c) {}
    virtual int eval_val(const double *x, double *val) const = 0;
    virtual int eval_gra(const double *x, double *gra) const = 0;
    virtual void set_axis(const double *a) {
        std::copy(a, a+3, a_.data());
    }
    virtual void set_center(const double *c) {
        std::copy(c, c+3, c_.data());
    }
protected:
    Vec3 a_, c_;
};

class linear_scalar_field : public scalar_field
{
public:
    linear_scalar_field(const Vec3 &a, const Vec3 &c) : scalar_field(a, c) {}
    virtual int eval_val(const double *x, double *val) const {
        Vec3 X(x[0], x[1], x[2]);
        *val = a_.dot(X-c_);
        return 0;
    }
    virtual int eval_gra(const double *x, double *gra) const {
        std::copy(a_.data(), a_.data()+a_.size(), gra);
        return 0;
    }
    void set_axis(const double *a) {
        scalar_field::set_axis(a);
    }
    void set_center(const double *c) {
        scalar_field::set_center(c);
    }
};

class cross_scalar_field : public scalar_field
{
public:
    cross_scalar_field(const Vec3 &a, const Vec3 &c) : scalar_field(a, c) {}
    int eval_val(const double *x, double *val) const {
        quad_scalar_field_(val, x, &a_[0], &c_[0]);
        return 0;
    }
    int eval_gra(const double *x, double *gra) const {
        quad_scalar_field_jac_(gra, x, &a_[0], &c_[0]);
        return 0;
    }
    void set_axis(const double *a) {
        scalar_field::set_axis(a);
    }
    void set_center(const double *c) {
        scalar_field::set_center(c);
    }
};

class dot_scalar_field : public scalar_field
{
public:
    dot_scalar_field(const Vec3 &a, const Vec3 &c) : scalar_field(a, c) {}
    int eval_val(const double *x, double *val) const {
        Vec3 X(x[0], x[1], x[2]);
        *val = std::pow(a_.dot(X-c_), 2);
        return 0;
    }
    int eval_gra(const double *x, double *gra) const {
        Vec3 X(x[0], x[1], x[2]);
        Vec3 g = 2*a_.dot(X-c_) * a_;
        std::copy(g.data(), g.data()+3, gra);
        return 0;
    }
    void set_axis(const double *a) {
        scalar_field::set_axis(a);
    }
    void set_center(const double *c) {
        scalar_field::set_center(c);
    }
};

//==============================================================================
class metaphor
{
public:
    virtual ~metaphor() {}
    metaphor(const Vec3 &c, const double ri, const double ro)
        : c_(c), ri_(ri), ro_(ro) {}
    virtual int eval_val(const double *x, double *val) const = 0;
    virtual int eval_gra(const double *x, double *gra) const = 0;
    virtual bool inner(const double *x) const = 0;
    virtual bool outer(const double *x) const = 0;
    virtual bool intermediate(const double *x) const = 0;
    virtual void set_center(const double *c) {
        std::copy(c, c+3, c_.data());
    }
    virtual void set_range(const double ri, const double ro) {
        ri_ = ri;
        ro_ = ro;
    }
    virtual double get_ri() const { return ri_; }
    virtual double get_ro() const { return ro_; }
    virtual Vec3 get_center() const { return c_; }
protected:
    Vec3 c_;
    double ri_, ro_;
};

class sphere_tool : public metaphor
{
public:
    sphere_tool(const Vec3 &c, const double ri, const double ro)
        : metaphor(c, ri, ro) {}
    int eval_val(const double *x, double *val) const {
        Vec3 X(x[0], x[1], x[2]);
        *val = (X-c_).norm();
        return 0;
    }
    int eval_gra(const double *x, double *gra) const {
        Vec3 X(x[0], x[1], x[2]);
        const double dist = (X-c_).norm();
        if ( dist <= 1e-16 ) {
            std::fill(gra, gra+3, 0);
        } else {
            gra[0] = (X[0]-c_[0]) / dist;
            gra[1] = (X[1]-c_[1]) / dist;
            gra[2] = (X[2]-c_[2]) / dist;
        }
        return 0;
    }
    bool inner(const double *x) const {
        Vec3 X(x[0], x[1], x[2]);
        const double dist = (X-c_).squaredNorm();
        return (dist < ri_*ri_);
    }
    bool intermediate(const double *x) const {
        Vec3 X(x[0], x[1], x[2]);
        const double dist = (X-c_).squaredNorm();
        return (ri_*ri_ <= dist && dist < ro_*ro_);
    }
    bool outer(const double *x) const {
        Vec3 X(x[0], x[1], x[2]);
        const double dist = (X-c_).squaredNorm();
        return (dist >= ro_*ro_);
    }
    void set_center(const double *c) {
        metaphor::set_center(c);
    }
    void set_range(const double ri, const double ro) {
        metaphor::set_range(ri, ro);
    }
    double get_ri() const {
        return metaphor::get_ri();
    }
    double get_ro() const {
        return metaphor::get_ro();
    }
    Vec3 get_center() const {
        return metaphor::get_center();
    }
};

class twist_cylinder : public metaphor
{
public:
    twist_cylinder(const Vec3 &c, const double ri, const double ro, const Vec3 &n)
        : metaphor(c, ri, ro), n_(n) {}
    int eval_val(const double *x, double *val) const {
        Vec3 X(x[0], x[1], x[2]);
        *val = n_.dot(X-c_);
        return 0;
    }
    int eval_gra(const double *x, double *gra) const {
        std::copy(n_.data(), n_.data()+3, gra);
        return 0;
    }
    bool inner(const double *x) const {
        double dist = 0;
        eval_val(x, &dist);
        return dist < ri_;
    }
    bool intermediate(const double *x) const {
        double dist = 0;
        eval_val(x, &dist);
        return (ri_ <= dist && dist < ro_);
    }
    bool outer(const double *x) const {
        double dist = 0;
        eval_val(x, &dist);
        return dist >= ro_;
    }
    void set_center(const double *c) {
        metaphor::set_center(c);
    }
    void set_range(const double ri, const double ro) {
        metaphor::set_range(ri, ro);
    }
    double get_ri() const {
        return metaphor::get_ri();
    }
    double get_ro() const {
        return metaphor::get_ro();
    }
    Vec3 get_center() const {
        return metaphor::get_center();
    }
private:
    const Vec3 n_;
};

//==============================================================================
class vector_field
{
public:
    vector_field(const Vec3 &center, const double ri, const double ro, const Vec3 &v, const std::string &type) {
        if ( type == "translate" ) {
            Vec3 u, w;
            get_ortn_basis(&v[0], &u[0], &w[0]);
            e_ = std::make_shared<linear_scalar_field>(u, center);
            f_ = std::make_shared<linear_scalar_field>(w, center);
            tools_ = std::make_shared<sphere_tool>(center, ri, ro);
        } else if ( type == "twist" ) {
            e_ = std::make_shared<dot_scalar_field>(v, center);
            f_ = std::make_shared<cross_scalar_field>(v, center);
            tools_ = std::make_shared<twist_cylinder>(center, ri, ro, v);
        }
    }
    virtual Vec3 operator ()(const Vec3 &x) const {
        Vec3 rtn;
        Vec3 e_gra = Vec3::Zero();
        Vec3 f_gra = Vec3::Zero();
        Vec3 r_gra = Vec3::Zero();

        if ( tools_->inner(&x[0]) ) {
            e_->eval_gra(&x[0], &e_gra[0]);
            f_->eval_gra(&x[0], &f_gra[0]);
            rtn = e_gra.cross(f_gra);
        } else if ( tools_->intermediate(&x[0]) ) {
            double e_val = 0, f_val = 0, r_val = 0;
            double ri = tools_->get_ri(), ro = tools_->get_ro();

            e_->eval_val(&x[0], &e_val);
            e_->eval_gra(&x[0], &e_gra[0]);
            f_->eval_val(&x[0], &f_val);
            f_->eval_gra(&x[0], &f_gra[0]);
            tools_->eval_val(&x[0], &r_val);
            tools_->eval_gra(&x[0], &r_gra[0]);

            double blend_val = eval_blend_val(r_val, ri, ro);
            double blend_jac = eval_blend_jac(r_val, ri, ro);
            Vec3 p_gra = (1.0-blend_val)*e_gra - blend_jac*e_val*r_gra;
            Vec3 q_gra = (1.0-blend_val)*f_gra - blend_jac*f_val*r_gra;
            rtn = p_gra.cross(q_gra);
        } else if ( tools_->outer(&x[0]) ) {
            rtn.setZero();
        }
        return rtn;
    }
private:
    std::shared_ptr<metaphor> tools_;
    std::shared_ptr<scalar_field> e_, f_;
};

class advector
{
public:
    typedef Eigen::Vector3d Vec3;
    advector(const std::vector<std::shared_ptr<vector_field>> &vfs)
        : vfs_(vfs) { }
    Vec3 operator()(const Vec3 &pos) const {
        const double h = 0.001;
        if ( vfs_.size() < 3 ) {
            const Vec3 vel = (*vfs_[vfs_.size()-1])(pos);
            return vel*h;
        } else {
            const Vec3 k1 = (*vfs_[vfs_.size()-3])(pos);
            const Vec3 k2 = (*vfs_[vfs_.size()-2])(pos+h*k1);
            const Vec3 k3 = (*vfs_[vfs_.size()-1])(pos+2*h*(-k1+2*k2));
            return (k1+4*k2+k3)/6*h;
        }
        std::cerr << "# error in advector." << std::endl;
        exit(0);
    }
private:
    const std::vector<std::shared_ptr<vector_field>> &vfs_;
};

class vel_field_deform
{
public:
    vel_field_deform();
    int load_model(const char *file);
    int translate_deform(const Vec3 &src, const Vec3 &des, const double ri, const double ro);
    int twist_deform(const Vec3 &center, const double ri, const double ro, const Vec3 &n, const size_t times);
    int bend_deform(const Vec3 &center, const double ri, const double ro);
    int save_model(const char *file);
private:
    Eigen::MatrixXi cell_;
    Eigen::MatrixXd nods_;
    std::vector<std::shared_ptr<vector_field>> vf_;
    std::shared_ptr<advector> advect_;
};

}
#endif
