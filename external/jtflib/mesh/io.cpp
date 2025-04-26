#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <set>
#include <vector>

#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/classification.hpp>

#include "io.h"

using namespace std;
using namespace zjucad::matrix;

namespace jtf{
  namespace mesh{

    static int analysis_obj(ifstream & obj_ifs, size_t & vertex_num,
                            size_t & face_num, size_t & face_type) // 3:triangle, 4: quad
    {
      std::string strtemp;
      vertex_num = 0;
      face_num = 0;
      face_type = -1;

      vector<string> string_q;
      while(!obj_ifs.eof()) {
          std::getline(obj_ifs,strtemp);
          if (strtemp[0] == 'v' || strtemp[0] == 'V'){
              if(strtemp.size() >= 2 && strtemp[1] == ' ')
                ++vertex_num;
            }
          if (strtemp[0] == 'f' || strtemp[0] == 'F'){
              ++face_num;
              boost::split(string_q, strtemp, boost::is_any_of("\t, ,\r"));
	      while(!string_q.empty()){
		if(string_q.back().empty() || string_q.back() == "\r")
		  string_q.pop_back();
		else
		  break;
	      }
              if(face_type == -1)
                face_type = string_q.size() - 1;
              else{
                  if(face_type != string_q.size()-1){
                      cerr << "# [error] unsupported mesh." << endl;
                      return __LINE__;
                    }
                }
            }
          strtemp.clear();
        }

      obj_ifs.clear();
      obj_ifs.seekg(0,std::ios::beg);

      if(face_type == 3)
        cerr << "# [info] triangle mesh " << endl;
      else if(face_type == 4)
        cerr << "# [info] quad mesh " << endl;
      cerr << "# [info] vertex " << vertex_num  << " face " << face_num  << endl;

      return 0;
    }

    int load_obj(const char * filename, matrixst & faces, matrixd & nodes)
    {
      ifstream obj_ifs(filename);
      if(obj_ifs.fail()) {
          std::cerr << "# open " << filename << " fail. " << std::endl;
          return __LINE__;
        }

      std::string strtemp;
      size_t pointsnum = 0;
      size_t polygonsnum = 0;
      size_t face_type = -1;

      if(analysis_obj(obj_ifs, pointsnum, polygonsnum, face_type))
        return __LINE__;

      assert(face_type == 3 || face_type == 4); // only support triangle/quad mesh

      nodes.resize(3 , pointsnum);
      faces.resize(face_type, polygonsnum);
      size_t idxpoints = 0;
      size_t idxpolygons = 0;

      while(!obj_ifs.eof()) {
          obj_ifs >> strtemp;
          if(obj_ifs.eof()) break;
          if(strtemp.empty()) {
              break;
            } if(strtemp == "#" || strtemp[0] == '#') {
              std::getline(obj_ifs,strtemp);
              strtemp.clear();
              continue;
            } if (strtemp == "v" || strtemp == "V") {
              for(size_t i = 0; i < 3; ++i){
                  obj_ifs  >> nodes(i, idxpoints);
                }
              ++idxpoints;
            } if (strtemp == "f" || strtemp == "F") {
              std::getline(obj_ifs,strtemp);
              std::stringstream sss(strtemp);
              std::string temp_str;
              size_t vidx = 0;
              while(sss >> temp_str) {
                  if(vidx == face_type) {
                      std::cerr << "# not trimesh" << std::endl;
                      return __LINE__;
                    }
                  int num = count(temp_str.begin(),temp_str.end(),'/');
                  if (num == 0) {
                      faces(vidx,idxpolygons) = atoi(temp_str.c_str()) - 1;
                    }else if(num == 1 || num == 2) {
                      temp_str.erase(temp_str.begin()+temp_str.find('/',0),temp_str.end());
                      faces(vidx,idxpolygons) = atoi(temp_str.c_str()) - 1;
                      //faces[idxpolygons].push_back( atoi(temp_str.c_str())-1);
                    }
                  ++vidx;
                }
              ++idxpolygons;
            }
          strtemp.clear();
        }
      obj_ifs.close(); // finish reading obj
      return 0;
    }

    int save_obj(const char * filename, const matrixst & face, const matrixd & node)
    {
      ofstream ofs(filename);
      if(ofs.fail()){
          cerr << "# [error] can not open file " << filename << endl;
          return __LINE__;
        }
      for(size_t t = 0; t < node.size(2); ++t){
          ofs << "v ";
          for(size_t i = 0; i < node.size(1); ++i)
            ofs << node(i,t) << " " ;
          ofs << endl;
        }

      for(size_t t = 0; t < face.size(2); ++t){
          ofs << "f ";
          for(size_t i = 0; i < face.size(1); ++i){
              ofs << face(i,t) + 1<< " ";
            }
          ofs << endl;
        }
      return 0;
    }

    int load_feature_line(const char *filename,
                          std::vector<std::vector<size_t> > &feature_lines)
    {
      ifstream ifs(filename);

      if(ifs.fail()){
          cerr << "# [error] can not load feature line." << endl;
          return __LINE__;
        }

      size_t line_num = 0;
      ifs >> line_num;
      size_t point_num = 0;

      feature_lines.resize(line_num);
      for(size_t li = 0; li < line_num; ++li){
          ifs >> point_num;
          feature_lines[li].resize(point_num);
          for(size_t pi = 0; pi < point_num; ++pi){
              ifs >> feature_lines[li][pi];
            }
        }

      cerr << "# [info] load feature line succeed: " << line_num << " lines." << endl;
      return 0;
    }

    int load_feature_line(
        const char * feature_line_file,
        const char * s2v_file,
        vector<vector<size_t> > &feature_line)
    {
      ifstream ifs_f(feature_line_file), ifs_s2v(s2v_file);
      if(ifs_f.fail()){
          cerr << "# [error] can not open feature line file." << endl;
          return __LINE__;
        }

      if(ifs_s2v.fail()){
          cerr << "# [error] can not open s2v_file file." << endl;
          return __LINE__;
        }

      matrix<int32_t> surface_node_idx2vol_node_idx;
      if(jtf::mesh::read_matrix(ifs_s2v, surface_node_idx2vol_node_idx))
        return __LINE__;

      feature_line.clear();
      size_t feature_line_num;
      ifs_f >> feature_line_num;
      if(feature_line_num == 0) return 0;

      feature_line.resize(feature_line_num);
      size_t pnum = 0;
      for(size_t li = 0; li < feature_line_num; ++li){
          ifs_f >> pnum;
          feature_line[li].resize(pnum);
          for(size_t pi = 0; pi < pnum; ++pi){
              ifs_f >> feature_line[li][pi];
              feature_line[li][pi] = surface_node_idx2vol_node_idx[feature_line[li][pi]];
            }
        }
      return 0;
    }

    int load_feature_line(const char * feature_line_file,
                          vector<deque<pair<size_t,size_t> > > & feature_lines)
    {
      ifstream ifs(feature_line_file);
      if(ifs.fail()){
          cerr << "# [error] can not open feature line file." << endl;
          return __LINE__;
        }

      size_t feature_line_num;
      ifs >> feature_line_num;
      feature_lines.clear();
      feature_lines.reserve(feature_line_num);

      size_t point_num;
      vector<size_t> one_line;
      deque<pair<size_t,size_t> > one_chain;
      for(size_t li = 0; li < feature_line_num; ++li){
          ifs >> point_num;
          one_line.resize(point_num);
          for(size_t pi = 0; pi < point_num; ++pi){
              ifs >> one_line[pi];
            }
          if(point_num == 1) continue;
          one_chain.clear();
          for(size_t pi = 0; pi < point_num - 1; ++pi){
              one_chain.push_back(make_pair(one_line[pi], one_line[pi+1]));
            }
          feature_lines.push_back(one_chain);
        }

      return 0;
    }


    int save_feature_line(const char * filename,
                          const vector<vector<size_t> > & lines)
    {
      ofstream ofs(filename);

      if(ofs.fail()){
          cerr << "# [error] can not save feature line." << endl;
          return __LINE__;
        }

      ofs << lines.size() << endl;

      for(size_t li = 0 ; li < lines.size(); ++li){
          ofs << lines[li].size() << endl;
          for(size_t pi = 0; pi < lines[li].size(); ++pi)
            ofs << lines[li][pi] << " ";
          ofs << endl;
        }

      cerr << "# [info] save feature line succeed: " << lines.size() << " lines." << endl;
      return 0;
    }

    int  tet_mesh_read_from_zjumat(
        const char *path,
        matrixd *node,
        matrixst *tet,
        matrixst *tri)
    {
      ifstream ifs(path, ifstream::binary);
      if(ifs.fail()) {
          cerr << "[info] " << "can not open file" << path << endl;
          return __LINE__;
        }

      matrixd node0;
      matrix<int> tet1, tri1;
      if(!node) node = &node0;

      read_matrix(ifs, *node);
      read_matrix(ifs, tet1);
      read_matrix(ifs, tri1);
      if(tet)
        *tet = tet1;
      if(tri)
        *tri = tri1;

      if(max(*tet) >= node->size(2)){
          cerr << "# [error] tet index beyond node size " << endl;
          return __LINE__;
        }
      return 0;
    }

    int tet_mesh_write_to_zjumat(
        const char *path,
        const matrixd *node,
        const matrixst *tet,
        const matrixst *tri)
    {
      //cout << "[info]" << "path = " << path << endl;
      ofstream ofs(path, ofstream::binary);
      matrix<int> tet1, tri1;
      if(node)
        write_matrix(ofs, *node);
      if(tet){
          tet1 = *tet;
          write_matrix(ofs, tet1);
        }
      if(tri){
          tri1 = *tri;
          write_matrix(ofs, tri1);
        }

      return 0;
    }


    int tet_mesh_read_from_vtk( const char *path, matrixd *node, matrixst *tet )
    {
      ifstream ifs(path);
      if(ifs.fail()) {
          cerr << "[info] " << "can not open file" << path << endl;
          return __LINE__;
        }

      matrixd node0;
      matrix<int> tet1;

      string str;
      int point_num = 0,cell_num = 0;

      while(!ifs.eof()){
          ifs >> str;
          if(str == "POINTS"){
              ifs >> point_num >> str;
              node0.resize(3, point_num);
              for(size_t i = 0;i < point_num; ++i){
                  for(size_t j = 0;j < 3; ++j)
                    ifs >> node0(j, i);
                }
              continue;
            }
          if(str == "CELLS"){
              ifs >> cell_num >> str;
              int point_number_of_cell = 0;
              vector<size_t> tet_temp;
              for(size_t ci = 0; ci < cell_num; ++ci){
                  ifs >> point_number_of_cell;
                  if(point_number_of_cell != 4){
                      for(size_t i = 0; i < point_number_of_cell; ++i)
                        ifs >> str;
                    }else{
                      int p;
                      for(size_t i = 0; i < point_number_of_cell; ++i){
                          ifs >> p;
                          tet_temp.push_back(p);
                        }
                    }
                }
              tet1.resize(4, tet_temp.size()/4);
              copy(tet_temp.begin(), tet_temp.end(), tet1.begin());
            }
        }

      *node = node0;
      *tet = tet1;
      return 0;
    }
  }
}
