#ifndef HJ_ZJUCAD_MATRIX_CONVERSION_H_
#define HJ_ZJUCAD_MATRIX_CONVERSION_H_

namespace zjucad { namespace matrix {

template <typename Con>
class container_1d
{
public:
  container_1d(const Con &con)
    :con_(con){}
  size_t size() const { return con_.size(); }
  typename Con::const_iterator begin() const { return con_.begin(); }
  typename Con::const_iterator end() const { return con_.end(); }
private:
  const Con &con_;
};

template<typename Con>
container_1d<Con>
to_mat(const Con &c)
{
  return container_1d<Con>(c);
}

template <typename E, typename Con>
void assign(matrix_expression<E> &dst, const container_1d<Con> &src)
{
  dst().resize(src.size());
  copy(src.begin(), src.end(), dst().begin());
}

}}

#endif
