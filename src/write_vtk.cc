#include "write_vtk.h"

#include <iostream>
#include <fstream>
#include <zjucad/matrix/matrix.h>
#include <zjucad/matrix/itr_matrix.h>

#include "vtk.h"

using namespace std;
using namespace zjucad::matrix;

namespace riemann {

int draw_face_value_to_vtk(const char *filename,
                           const double *vert, const size_t vert_num,
                           const size_t *face, const size_t face_num,
                           const double *data, const size_t type_num,
                           const vector<string> &data_name) {
  ofstream os(filename);
  if ( os.fail() )
    return __LINE__;
  os.precision(15);
  tri2vtk(os, vert, vert_num, face, face_num);
  if ( type_num == 0 )
    return 0;
  cell_data(os, data, face_num, data_name[0].c_str(), data_name[0].c_str());
  for (size_t i = 1; i < type_num; ++i) {
    vtk_data(os, data+i*face_num, face_num, data_name[i].c_str(), data_name[i].c_str());
  }
  os.close();
  return 0;
}

int draw_vert_value_to_vtk(const char *filename,
                           const double *vert, const size_t vert_num,
                           const size_t *face, const size_t face_num,
                           const double *data) {
  ofstream os(filename);
  if ( os.fail() )
    return __LINE__;
  os.precision(15);
  tri2vtk(os, vert, vert_num, face, face_num);
  point_data(os, data, vert_num, "vert_value", "vert_value");
  os.close();
  return 0;
}

int draw_edge_value_to_vtk(const char *filename,
                           const double *vert, const size_t vert_num,
                           const size_t *edge, const size_t edge_num,
                           const double *data) {
  ofstream os(filename);
  if ( os.fail() )
    return __LINE__;
  os.precision(15);
  line2vtk(os, vert, vert_num, edge, edge_num);
  cell_data(os, data, edge_num, "edge_value", "edge_value");
  os.close();
  return 0;
}

int draw_face_direct_field(const char *filename,
                           const double *vert, const size_t vert_num,
                           const size_t *face, const size_t face_num,
                           const double *field) {
  ofstream os(filename);
  if ( os.fail() )
    return __LINE__;

  itr_matrix<const double *> nods(3, vert_num, vert);
  itr_matrix<const size_t *> cell(3, face_num, face);
  itr_matrix<const double *> df(3, face_num, field);

  matrix<size_t> line(2, face_num);
  line(0, colon()) = colon(0, face_num-1);
  line(1, colon()) = colon(face_num, 2*face_num-1);
  matrix<double> pts(3, 2*face_num);
#pragma omp parallel for
  for (size_t i = 0; i < face_num; ++i) {
    pts(colon(), i) = nods(colon(), cell(colon(), i))*ones<double>(3, 1)/3.0;
    pts(colon(), i+face_num) = pts(colon(), i)+df(colon(), i);
  }
  line2vtk(os, pts.begin(), pts.size(2), line.begin(), line.size(2));
  os.close();
  return 0;
}

int draw_edge_direct_field(const char *filename,
                           const double *vert, const size_t vert_num,
                           const size_t *edge, const size_t edge_num,
                           const double *field) {
  ofstream os(filename);
  if ( os.fail() )
    return __LINE__;

  itr_matrix<const double *> nods(3, vert_num, vert);
  itr_matrix<const size_t *> cell(2, edge_num, edge);
  itr_matrix<const double *> df(3, edge_num, field);

  matrix<size_t> line(2, edge_num);
  line(0, colon()) = colon(0, edge_num-1);
  line(1, colon()) = colon(edge_num, 2*edge_num-1);
  matrix<double> pts(3, 2*edge_num);
#pragma omp parallel for
  for (size_t i = 0; i < edge_num; ++i) {
    pts(colon(), i) = nods(colon(), cell(colon(), i))*ones<double>(2, 1)/2.0;
    pts(colon(), i+edge_num) = pts(colon(), i)+df(colon(), i);
  }
  line2vtk(os, pts.begin(), pts.size(2), line.begin(), line.size(2));
  os.close();
  return 0;
}

int draw_vert_direct_field(const char *filename,
                           const double *vert, const size_t vert_num,
                           const double *field) {
  ofstream os(filename);
  if ( os.fail() )
    return __LINE__;

  itr_matrix<const double *> nods(3, vert_num, vert);
  itr_matrix<const double *> df(3, vert_num, field);

  matrix<size_t> line(2, vert_num);
  line(0, colon()) = colon(0, vert_num-1);
  line(1, colon()) = colon(vert_num, 2*vert_num-1);
  matrix<double> pts(3, 2*vert_num);
  pts(colon(), colon(0, vert_num-1)) = nods-df;
  pts(colon(), colon(vert_num, 2*vert_num-1)) = nods+df;
  line2vtk(os, pts.begin(), pts.size(2), line.begin(), line.size(2));
  os.close();
  return 0;
}

}
