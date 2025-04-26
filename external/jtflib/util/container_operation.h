#ifndef CONTAINER_OPERATION_H
#define CONTAINER_OPERATION_H

#include <algorithm>
#include <vector>
#include <list>
#include <deque>
#include <set>
#include <boost/unordered_set.hpp>

template <class InputIterator1, class InputIterator2, class OutputIterator>
OutputIterator find_difference_set_with_sorted_input(
    InputIterator1 first1, InputIterator1 last1,
    InputIterator2 first2, InputIterator2 last2,
    OutputIterator result )
{
  return set_difference(first1,last1,first2,last2,result);
}

template <class InputIterator1, class InputIterator2, class OutputIterator>
OutputIterator find_difference_set(
    InputIterator1 first1, InputIterator1 last1,
    InputIterator2 first2, InputIterator2 last2,
    OutputIterator result )
{
	std::sort(first1,last1);
	std::sort(first2,last2);
  return find_difference_set_with_sorted_input(first1,last1,
                                               first2,last2,result);
}

template <class InputIterator1, class InputIterator2, class OutputIterator>
OutputIterator find_intersection_set_with_sorted_input(
    InputIterator1 first1, InputIterator1 last1,
    InputIterator2 first2, InputIterator2 last2,
    OutputIterator result)
{
  return set_intersection(first1,last1,first2,last2,result);
}

template <class InputIterator1, class InputIterator2, class OutputIterator>
OutputIterator find_intersection_set(
    InputIterator1 first1, InputIterator1 last1,
    InputIterator2 first2, InputIterator2 last2,
    OutputIterator result)
{
	std::sort(first1,last1);
	std::sort(first2,last2);
  return find_intersection_set_with_sorted_input(first1,last1,
                                                 first2,last2,result);
}

template <class T>
std::vector<T>& operator << (std::vector<T> & container, const T & data){
  container.push_back(data);
  return container;
}

template <class T>
std::list<T>& operator << (std::list<T> & container, const T & data){
  container.push_back(data);
  return container;
}

template <class T>
std::deque<T>& operator << (std::deque<T> & container, const T & data){
  container.push_back(data);
  return container;
}

template <class T>
std::set<T>& operator << (std::set<T> & container, const T & data){
  container.insert(data);
  return container;
}

template <class T>
boost::unordered_set<T>& operator << (boost::unordered_set<T> container, const T & data){
  container.insert(data);
  return container;
}


template <class T>
std::list<T>& unite(std::list<T> & to, const std::list<T> & from){
  for(typename std::list<T>::const_iterator lcit = from.begin();
      lcit != from.end(); ++lcit){
    if(find(to.begin(), to.end(), *lcit) == to.end()){
      to.push_back(*lcit);
    }
  }
  return to;
}

template <typename T>
void generate_N_number(std::vector<T>& containter, const size_t & num)
{
  containter.reserve(num);
  T data = 0;
  for(size_t t = 0; t < num; ++t)
    containter << data++;
}

template <typename T>
void generate_N_number(std::set<T>& containter, const size_t & num)
{
  T data = 0;
  for(size_t t = 0; t < num; ++t)
    containter << data++;
}

template <typename T>
void generate_N_number(boost::unordered_set<T>& containter, const size_t & num)
{
  T data = 0;
  for(size_t t = 0; t < num; ++t)
    containter << data++;
}

template <typename T>
std::vector<T> make_vector()
{
  std::vector<T> empty;
  return empty;
}

template <typename T>
std::vector<T> make_vector(const T * v, const size_t &num)
{
  assert(v);
  std::vector<T> vec;
  for(size_t t = 0; t < num; ++t)
    vec << v[t];
  return vec;
}

template <typename T>
std::pair<T, T> make_edge(T a, T b) {
    assert(a != b);
    if(a > b)
        std::swap(a, b);
    return std::make_pair(a, b);
}

template <typename T>
double vec_norm(const std::vector<T> & a, const std::vector<T> & b)
{
  assert(a.size() == b.size());
  double r = 0;
  for(size_t i = 0; i < a.size(); ++i){
    r += (a[i] - b[i]) * (a[i] - b[i]);
  }
  return r;
}

#endif // CONTAINER_OPERATION_H
