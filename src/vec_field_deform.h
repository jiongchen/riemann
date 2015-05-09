#ifndef VEC_FIELD_DEFORM_H
#define VEC_FIELD_DEFORM_H

#include <memory>
#include <Eigen/Sparse>

namespace geom_deform {

static double b(const double r, const double ri, const double ro) {
    return 4*std::pow(r-ri, 3)*(1-(r-ri)/(ro-ri))/std::pow(ro-ri, 3) + std::pow(r-ri, 4)/std::pow(ro-ri, 4);
}

static double  db_dr(const double r, const double ri, const double ro) {
    return 12*std::pow(r-ri, 2)*(1-(r-ri)/(ro-ri))/std::pow(ro-ri, 3);
}

class scalar_field {
public:
    virtual int value(const double *x, double *val) const = 0;
    virtual int gradient(const double *x, double *gra) const = 0;
};

class constant_scalar_field : public scalar_field {

};

class linear_scalar_field : public scalar_field {

};

class quadratic_scalar_field : public scalar_field {

};

class piecewise_scalar_field {
public :
    piecewise_scalar_field(const std::shared_ptr<scalar_field> &e,
                           const std::shared_ptr<scalar_field> &r);
    int value(const double *x, const double ri, const double ro, double *val) const {
        double d = 0;
        r_->value(x, &d);
        if ( d < ri ) {
            e_->value(x, val);
        } else if ( ri <= d && d < ro ) {
            e_->value(x, val);
            *val *= 1.0 - b(d, ri, ro);
        } else {
            *val = 0;
        }
        return 0;
    }
    int gradient(const double *x, const double ri, const double ro, double *gra) const {
        double d = 0;
        r_->value(x, &d);
        if ( d < ri ) {
            e_->gradient(x, gra);
        } else if ( ri <= d && d< ro ) {

        } else {
            gra[0] = gra[1] = gra[2] = 0;
        }
        return 0;
    }
private:
    const std::shared_ptr<scalar_field> r_;
    const std::shared_ptr<scalar_field> e_;
};

class vector_field {
public:
    vector_field();
private:
    std::shared_ptr<piecewise_scalar_field> p_, q_;
    std::shared_ptr<scalar_field> r_, e_, f_;
};

}
#endif
