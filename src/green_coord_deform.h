#ifndef GREEN_COORD_DEFORM_H
#define GREEN_COORD_DEFORM_H

#include <Eigen/Sparse>

namespace geom_deform {

class green_deform
{
public:
    green_deform() {}
    virtual ~green_deform() {}
    virtual int load_sample_points(const char *file) = 0;
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
    int load_sample_points(const char *file);
    int load_cage(const char *file);
    int calc_green_coords();
    int move_cage(const size_t id, const double *dx);
    int deform();
    int dump(const char *file);
    int dump_normal(const char *file);
    int dump_cage(const char *file);
protected:
    int calc_outward_normal();
    int update_cage_edge_length();
private:
    Eigen::MatrixXi cell_;
    Eigen::MatrixXd nods_;

    Eigen::MatrixXi cage_cell_;
    Eigen::MatrixXd cage_nods_;
    Eigen::MatrixXd cage_normal_;

    Eigen::VectorXd rest_len_;
    Eigen::VectorXd curr_len_;
    Eigen::SparseMatrix<double> phi_;
    Eigen::SparseMatrix<double> psi_;
};

class green_deform_3d : public green_deform
{
public:
};

}

#endif
