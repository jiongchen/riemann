#ifndef GREEN_COORD_DEFORM_H
#define GREEN_COORD_DEFORM_H

#include <Eigen/Sparse>

namespace geom_deform {

class green_deform
{
public:
    green_deform() {}
    virtual ~green_deform() {}
    virtual int load_mesh(const char *file) = 0;
    virtual int load_cage(const char *file) = 0;
    virtual int calc_green_coords() = 0;
    virtual int move_cage(const size_t id, const double *dx) = 0;
    virtual int deform() = 0;
    virtual int dump(const char *file) = 0;
protected:
    virtual int calc_outward_normal() = 0;
};

class green_deform_2d : public green_deform
{
public:
    green_deform_2d();
    int load_mesh(const char *file);
    int load_cage(const char *file);
    int calc_green_coords();
    int move_cage(const size_t id, const double *dx);
    int deform();
    int dump(const char *file);
    int dump_normal(const char *file);
    int dump_cage(const char *file);
protected:
    int calc_outward_normal();
private:
    Eigen::MatrixXi cell_;
    Eigen::VectorXd X_;

    Eigen::MatrixXi cage_cell_;
    Eigen::VectorXd Xcage_;
    Eigen::VectorXd N_;

    Eigen::VectorXd d_;
    Eigen::VectorXd s_;
    Eigen::SparseMatrix<double> phi_;
    Eigen::SparseMatrix<double> psi_;
};

class green_deform_3d : public green_deform
{
public:
};

}

#endif
