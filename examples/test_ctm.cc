#include <iostream>
#include <openctm.h>
#include <jtflib/mesh/io.h>

using namespace std;
using namespace zjucad::matrix;

int main(int argc, char *argv[])
{
  if ( argc != 3 ) {
    cerr << "#usage: ./test_ctm model.obj model.ctm\n";
    return __LINE__;
  }
  matrix<size_t> tris;
  matrix<double> nods;
  jtf::mesh::load_obj(argv[1], tris, nods);

  CTMcontext context;
  CTMuint vertCount, triCount, * indices;
  CTMfloat *vertices;

  // Create our mesh in memory
  vertCount = nods.size(2);
  triCount = tris.size(2);
  vertices = (CTMfloat *)malloc(3*sizeof(CTMfloat)*vertCount);
  indices = (CTMuint *)malloc(3*sizeof(CTMuint)*triCount);
  std::copy(tris.begin(), tris.end(), indices);
  std::copy(nods.begin(), nods.end(), vertices);
  cout << "[Info] total memory cost: " << 3*sizeof(CTMfloat)*vertCount+3*sizeof(CTMuint)*triCount << endl;

  // Create a new context
  context = ctmNewContext(CTM_EXPORT);

  // Config compression method
  ctmCompressionMethod(context, CTM_METHOD_MG2);
  ctmCompressionLevel(context, 1);

  // Define our mesh representation to OpenCTM (store references to it in
  // the context)
  ctmDefineMesh(context, vertices, vertCount, indices, triCount, NULL);

  // Save the OpenCTM file
  ctmSave(context, argv[2]);

  // Free the context
  ctmFreeContext(context);

  // Free our mesh
  free(indices);
  free(vertices);

  cout << "[Info] done\n";
  return 0;
}
