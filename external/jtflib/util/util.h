#ifndef COMMON_UTIL_H
#define COMMON_UTIL_H

#include <vector>
#include <deque>
#include <cstddef>
#include <limits>
#include <time.h>

namespace jtf{
  namespace util{
    ///
    /// @brief extract_chain_from_edges
    /// @param edges input edges
    /// @param chain output chain
    /// @param point_mapping
    ///
    void extract_chain_from_edges(
        const std::vector<std::pair<std::size_t,std::size_t> > & edges,
        std::vector<std::deque<std::pair<std::size_t,std::size_t> > > & chain,
		const std::vector<size_t> * point_mapping = 0);

    template <typename FUNC>
    double func_time(const FUNC &f)
    {
      clock_t start = clock();
      f();
      clock_t end = clock();
      double total = end - start;
      return total/CLOCKS_PER_SEC;
    }
  }
}
#endif // common_util_h
