#ifndef HJ_ZJUCAD_MATRIX_EIGEN_H_
#define HJ_ZJUCAD_MATRIX_EIGEN_H_

namespace zjucad { namespace matrix {

template <typename E, typename Derived>
void assign(matrix_expression<E> &dst, const Eigen::MatrixBase<Derived> &src)
{
  dst().resize(src.rows(), src.cols());
  for(size_t c = 0; c < dst().size(2); ++c)
    for(size_t r = 0; r < dst().size(1); ++r)
      dst()(r, c) = src(r, c);
}

}}

#endif
