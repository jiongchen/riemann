#ifndef JTF_MESH_IO_H
#define JTF_MESH_IO_H

#include <iostream>
#include <fstream>
#include <deque>
#include "config.h"
#include <zjucad/matrix/matrix.h>

namespace jtf{
  namespace mesh{

    typedef zjucad::matrix::matrix<double> matrixd;
    typedef zjucad::matrix::matrix<size_t> matrixst;

    struct material {
      double color_[3][3], shi_;	// amb, dif, spe
      std::string name_, texture_;
    };

    ///
    /// @brief load_obj, it only support pure triangle/quad mesh
    /// @param filename
    /// @param mesh output mesh matrix
    /// @param node output node matrix
    /// @return return 0 if fine or non-zeros
    ///
	int JTF_MESH_API load_obj(const char * filename, matrixst & mesh, matrixd & node);

    ///
    /// @brief save_obj
    /// @param filename output filename
    /// @param mesh input mesh matrix
    /// @param node input node matrix
    /// @return return 0 if fine or non-zeros
    ///
	int JTF_MESH_API save_obj(const char * filename, const matrixst & mesh, const matrixd & node);

    ///
    /// @brief load_feature_line
    /// @param filename
    /// @param feature_lines
    /// @return return 0 if fine or non-zeros
    ///
	int JTF_MESH_API load_feature_line(
        const char *filename,
        std::vector<std::vector<size_t> > & feature_lines);

    ///
    /// @brief load_feature_line
    /// @param feature_line_file
    /// @param point_mapping
    /// @param feature_line
    /// @return return 0 if ok
    ///
	int JTF_MESH_API load_feature_line(
        const char * feature_line_file,
        const char * point_mapping,
        std::vector<std::vector<size_t> > &feature_line);

    ///
    /// @brief load_feature_lines
    /// @param feature_line_file
    /// @param feature_lines
    /// @return return 0 of or
    ///
	int JTF_MESH_API load_feature_line(
        const char * feature_line_file,
        std::vector<std::deque<std::pair<size_t,size_t> > > & feature_lines);

    ///
    /// @brief save_feature_line
    /// @param filename
    /// @param feature_lines
    /// @return return 0 if fine or non-zeros
    ///
	int JTF_MESH_API save_feature_line(
        const char *filename,
        const std::vector<std::vector<size_t> > & feature_lines);

    ///
    /// @brief tet_mesh_read_from_zjumat, each tet has positive volume,
    ///        i.e. <(p1-p0)x(p2-p0),(p3-p0)> >= 0
    /// @param path
    /// @param node
    /// @param tet
    /// @param tri
    /// @return 0 if succeed, or non-zeros
    ///
	int JTF_MESH_API tet_mesh_read_from_zjumat(const char *path, matrixd *node = 0,
                                  matrixst *tet = 0, matrixst *tri = 0);


    ///
    /// @brief tet_mesh_write_to_zjumat, each tet has positive volume,
    ///        i.e. <(p1-p0)x(p2-p0),(p3-p0)> >= 0
    /// @param path
    /// @param node
    /// @param tet
    /// @param tri
    /// @return 0 if succeed, or non-zeros
    ///
	int JTF_MESH_API tet_mesh_write_to_zjumat(const char *path, const matrixd *node = 0,
                                 const matrixst *tet = 0, const matrixst *tri = 0);

    ///
    /// @brief tet_mesh_read_from_vtk
    /// @param path
    /// @param node
    /// @param tet
    /// @return 0 if succeed, or non-zeros
    ///
	int JTF_MESH_API tet_mesh_read_from_vtk(const char *path, matrixd *node = 0,
                               matrixst *tet = 0);

    ///////////////////////////////////////////////////////////////////
    ///  Low level reading: matrix
    //////////////////////////////////////////////////////////////////

    template <typename T>
    int read_matrix(std::istream &is, zjucad::matrix::matrix<T> &m)
    {
      int nrow, ncol;
      is.read((char *)&nrow, sizeof(int));
      is.read((char *)&ncol, sizeof(int));
      m.resize(nrow, ncol);
      is.read((char *)&m[0], sizeof(T)*m.size());
      return is.fail();
    }

    template <typename T>
    int read_matrix(std::istream &is, zjucad::matrix::matrix<T> &m, size_t nrow)
    {
      zjucad::matrix::matrix<T> tmp;
      read_matrix(is, tmp);
      if((tmp.size() % nrow) != 0)
        return -1;
      m.resize(nrow, tmp.size()/nrow);
      std::copy(tmp.begin(), tmp.end(), m.begin());
      return is.fail();
    }

    template <typename T>
    int write_matrix(std::ostream &os, const zjucad::matrix::matrix<T> &m)
    {
      int nrow = m.size(1), ncol = m.size(2);
      os.write((const char *)&nrow, sizeof(int));
      os.write((const char *)&ncol, sizeof(int));
      os.write((const char *)&m[0], sizeof(T)*m.size());
      return os.fail();
    }

    template <typename T>
    inline int read_matrix(const char *path, zjucad::matrix::matrix<T> &m)
    {
      using namespace std;
      ifstream ifs(path, ifstream::binary);
      if(ifs.fail()) {
          std::cerr << "open " << path << " for read fail." << endl;
          return __LINE__;
        }
      read_matrix(ifs, m);
      return ifs.fail();
    }

    template <typename T>
    inline int write_matrix(const char *path, const zjucad::matrix::matrix<T> &m)
    {
      using namespace std;
      ofstream ofs(path, ofstream::binary);
      if(ofs.fail()) {
          std::cerr << "open " << path << " for write fail." << endl;
          return __LINE__;
        }
      write_matrix(ofs, m);
      return ofs.fail();
    }
  }
}

#endif //TRIMESH_IO_H
