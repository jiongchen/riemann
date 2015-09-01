#ifndef WRITE_VTK_H
#define WRITE_VTK_H

#include <vector>
#include <string>

namespace riemann {

int draw_face_value_to_vtk(const char *filename,
                           const double *vert, const size_t vert_num,
                           const size_t *face, const size_t face_num,
                           const double *data, const size_t type_num,
                           const std::vector<std::string> &data_name);

int draw_vert_value_to_vtk(const char *filename,
                           const double *vert, const size_t vert_num,
                           const size_t *face, const size_t face_num,
                           const double *data);

int draw_edge_value_to_vtk(const char *filename,
                           const double *vert, const size_t vert_num,
                           const size_t *edge, const size_t edge_num,
                           const double *data);

int draw_face_direct_field(const char *filename,
                           const double *vert, const size_t vert_num,
                           const size_t *face, const size_t face_num,
                           const double *field);

int draw_edge_direct_field(const char *filename,
                           const double *vert, const size_t vert_num,
                           const size_t *edge, const size_t edge_num,
                           const double *field);

int draw_vert_direct_field(const char *filename,
                           const double *vert, const size_t vert_num,
                           const double *field);
}

#endif
