/*
	State Key Lab of CAD&CG Zhejiang Unv.

	Author: Jin Huang (hj@cad.zju.edu.cn)

	Copyright (c) 2004-2011 <Jin Huang>
	All rights reserved.

	Redistribution and use in source and binary forms, with or without
	modification, are permitted provided that the following conditions
	are met:
	1. Redistributions of source code must retain the above copyright
	notice, this list of conditions and the following disclaimer.
	2. Redistributions in binary form must reproduce the above copyright
	notice, this list of conditions and the following disclaimer in the
	documentation and/or other materials provided with the distribution.
	3. The name of the author may not be used to endorse or promote products
	derived from this software without specific prior written permission.
*/

#ifndef _ZJUCAD_MATRIX_EXPRESSION_H_
#define _ZJUCAD_MATRIX_EXPRESSION_H_

#include "configure.h"
#include "iterator.h"
#include "functional.h"
#include <algorithm>
#include <cassert>
#include "colon.h" //fix compile bug for gcc4.3.3-5ubuntu4.

namespace zjucad {namespace matrix {

template <typename E>
struct matrix_expression;
class colon_void;

//2D //cv_ or _cv
template <typename E>
class idx_cv_i;
template <typename E1, typename E2>
class idx_cv_J;
template <typename E, typename R>
class idx_cv_range;
template <typename E, typename R>
class idx_cv_iv;

template <typename E>
class idx_i_cv;
template <typename E1, typename E2>
class idx_J_cv;
template <typename E, typename R>
class idx_range_cv;
template <typename E, typename R>
class idx_iv_cv;

//2D different mode with out cv
//template <typename E>
//class idx_i_i;
template <typename E1, typename E2>
class idx_J_i;
template <typename E, typename R>
class idx_range_i;
template <typename E, typename R>
class idx_iv_i;
template <typename E1, typename E2>
class idx_i_J;
template <typename E1, typename E2, typename E3>
class idx_J_J;
template <typename E1, typename R, typename E2>
class idx_range_J;
template <typename E1, typename R, typename E2>
class idx_iv_J;
template <typename E, typename R>
class idx_i_range;
template <typename E1, typename E2, typename R>
class idx_J_range;
template <typename E, typename R1, typename R2>
class idx_range_range;
template <typename E, typename R1, typename R2>
class idx_iv_range;
template <typename E, typename R>
class idx_i_iv;
template <typename E1, typename E2, typename R>
class idx_J_iv;
template <typename E, typename R1, typename R2>
class idx_range_iv;
template <typename E, typename R1, typename R2>
class idx_iv_iv;

//1D
template <typename E, typename E1>
class idx_I;
template <typename E>
class idx_cv;
template <typename E, typename T>
class idx_range;
template <typename E, typename T>
class idx_iv;

template <typename T>
struct range;

#define CONST_PROXY_ACCESS\
	\
	const expression_type &\
		operator()(colon_void, colon_void) const\
  { return *this; }                         \
	const idx_cv_i<const expression_type>\
		operator()(colon_void, size_type col) const\
		{ return idx_cv_i<const expression_type>(*this, col); }\
	template <typename E9>\
	const idx_cv_J<const expression_type, E9>\
		operator()(colon_void, const matrix_expression<E9> &col) const\
		{ return idx_cv_J<const expression_type, E9>(*this, col()); }\
	template <typename T8>\
	const idx_cv_range<const expression_type, const range<T8> >\
		operator()(colon_void, const range<T8> &range_col) const\
		{ return idx_cv_range<const expression_type, const range<T8> >(*this, range_col); }\
	template <typename T8>\
	const idx_cv_iv<const expression_type, const interleaved<T8> >\
		operator()(colon_void, const interleaved<T8> &interleaved_col) const\
		{ return idx_cv_iv<const expression_type, const interleaved<T8> >(*this, interleaved_col); }\
	\
	const idx_i_cv<const expression_type>\
		operator()(size_type row, colon_void) const\
		{ return idx_i_cv<const expression_type>(*this, row); }\
	template <typename E9>\
	const idx_J_cv<const expression_type, E9>\
		operator()(const matrix_expression<E9> &row, colon_void) const\
        { return idx_J_cv<const expression_type, E9>(*this, row()); }\
	template <typename T8>\
    const idx_range_cv<const expression_type, const range<T8> >\
		operator()(const range<T8> &range_row, colon_void) const\
		{ return idx_range_cv<const expression_type, const range<T8> >(*this, range_row); }\
	template <typename T8>\
	const idx_iv_cv<const expression_type, const interleaved<T8> >\
		operator()(const interleaved<T8> &interleaved_row, colon_void) const\
		{ return idx_iv_cv<const expression_type, const interleaved<T8> >(*this, interleaved_row); }\
	\
	template <typename E9>\
	const idx_J_i<const expression_type, E9>\
		operator()(const matrix_expression<E9> &row, size_type col) const\
		{ return idx_J_i<const expression_type, E9>(*this, row(), col); }\
	template <typename T8>\
	const idx_range_i<const expression_type, const range<T8> >\
		operator()(const range<T8> &range_row, size_type i) const\
		{ return idx_range_i<const expression_type, const range<T8> >(*this, range_row, i); }\
	template <typename T8>\
	const idx_iv_i<const expression_type, const interleaved<T8> >\
		operator()(const interleaved<T8> &interleaved_row, size_type i) const\
		{ return idx_iv_i<const expression_type, const interleaved<T8> >(*this, interleaved_row, i); }\
	template <typename E9>\
	const idx_i_J<const expression_type, E9>\
		operator()(size_type row, const matrix_expression<E9> &col) const\
		{ return idx_i_J<const expression_type, E9>(*this, row, col()); }\
	template <typename E8, typename E9>\
	const idx_J_J<const expression_type, E8, E9>\
		operator()(const matrix_expression<E8> &row, const matrix_expression<E9> &col) const\
		{ return idx_J_J<const expression_type, E8, E9>(*this, row(), col()); }\
	template <typename T8, typename E8>\
	const idx_range_J<const expression_type, const range<T8>, E8 >\
		operator()(const range<T8> &range_row, const matrix_expression<E8> &col) const\
		{ return idx_range_J<const expression_type, const range<T8>, E8 >(*this, range_row, col()); }\
	template <typename T8, typename E8>\
	const idx_iv_J<const expression_type, const interleaved<T8>, E8 >\
		operator()(const interleaved<T8> &interleaved_row, const matrix_expression<E8> &col) const\
		{ return idx_iv_J<const expression_type, const interleaved<T8>, E8 >(*this, interleaved_row, col()); }\
	template <typename T8>\
	const idx_i_range<const expression_type, const range<T8> >\
		operator()(size_type i, const range<T8> &range_col) const\
		{ return idx_i_range<const expression_type, const range<T8> >(*this, i, range_col); }\
	template <typename E8, typename T8>\
	const idx_J_range<const expression_type, E8, const range<T8> >\
		operator()(const matrix_expression<E8> &row, const range<T8> &range_col) const\
		{ return idx_J_range<const expression_type, E8, const range<T8> >(*this, row(), range_col); }\
	template <typename T8, typename T9>\
	const idx_range_range<const expression_type, const range<T8>, const range<T9> >\
		operator()(const range<T8> &range1, const range<T9> &range2) const\
		{ return idx_range_range<const expression_type, const range<T8>, const range<T9> >(*this, range1, range2); }\
	template <typename T8, typename T9>\
	const idx_iv_range<const expression_type, const interleaved<T8>, const range<T9> >\
		operator()(const interleaved<T8> &interleaved_row, const range<T9> &range_col) const\
		{ return idx_iv_range<const expression_type, const interleaved<T8>, const range<T9> >(*this, interleaved_row, range_col); }\
	template <typename T8>\
	const idx_i_iv<const expression_type, const interleaved<T8> >\
		operator()(size_type row, const interleaved<T8> &interleaved_col ) const\
		{ return idx_i_iv<const expression_type, const interleaved<T8> >(*this, row, interleaved_col); }\
	template <typename E8, typename T8>\
	const idx_J_iv<const expression_type, E8, const interleaved<T8> >\
		operator()(const matrix_expression<E8> &row, const interleaved<T8> &interleaved_col) const\
		{ return idx_J_iv<const expression_type, E8, const interleaved<T8> >(*this, row(), interleaved_col); }\
	template <typename T8, typename T9>\
	const idx_range_iv<const expression_type, const range<T8>, const interleaved<T9> >\
		operator()(const range<T8> &range_row, const interleaved<T9> &interleaved_col) const\
		{ return idx_range_iv<const expression_type, const range<T8>, const interleaved<T9> >(*this, range_row, interleaved_col); }\
	template <typename T8, typename T9>\
	const idx_iv_iv<const expression_type, const interleaved<T8>, const interleaved<T9> >\
		operator()(const interleaved<T8> &interleaved_row, const interleaved<T9> &interleaved_col) const\
		{ return idx_iv_iv<const expression_type, const interleaved<T8>, const interleaved<T9> >(*this, interleaved_row, interleaved_col); }\
	\
	template <typename E9>\
	const idx_I<const expression_type, E9>\
		operator()(const matrix_expression<E9> &I) const\
		{ return idx_I<const expression_type, E9>(*this, I()); }\
	const idx_cv<const expression_type>\
		operator()(colon_void) const\
		{ return idx_cv<const expression_type>(*this); }\
	template <typename T9>\
	const idx_range<const expression_type, T9>\
		operator()(const range<T9> &range) const\
		{ return idx_range<const expression_type, T9>(*this, range); }\
	template <typename T9>\
	const idx_iv<const expression_type, T9>\
		operator()(const interleaved<T9> &interleaved) const\
		{return idx_iv<const expression_type, T9>(*this, interleaved); }

#define NON_CONST_PROXY_ACCESS\
	\
	expression_type &\
		operator()(colon_void, colon_void)\
  { return *this; }                         \
	idx_cv_i<expression_type>\
		operator()(colon_void, size_type col)\
		{ return idx_cv_i<expression_type>(*this, col); }\
	template <typename E9>\
	idx_cv_J<expression_type, E9>\
		operator()(colon_void, const matrix_expression<E9> &col)\
		{ return idx_cv_J<expression_type, E9>(*this, col()); }\
	template <typename T8>\
	idx_cv_range<expression_type, const range<T8> >\
		operator()(colon_void, const range<T8> &range_col) \
		{ return idx_cv_range<expression_type, const range<T8> >(*this, range_col); }\
	template <typename T8>\
	idx_cv_iv<expression_type, const interleaved<T8> >\
		operator()(colon_void, const interleaved<T8> &interleaved_col) \
		{ return idx_cv_iv<expression_type, const interleaved<T8> >(*this, interleaved_col); }\
	\
	idx_i_cv<expression_type>\
		operator()(size_type row, colon_void)\
		{ return idx_i_cv<expression_type>(*this, row); }\
	template <typename E9>\
	idx_J_cv<expression_type, E9>\
		operator()(const matrix_expression<E9> &row, colon_void)\
		{ return idx_J_cv<expression_type, E9>(*this, row()); }\
	template <typename T8>\
	idx_range_cv<expression_type, const range<T8> >\
		operator()(const range<T8> &range_row ,colon_void) \
		{ return idx_range_cv<expression_type, const range<T8> >(*this, range_row); }\
	template <typename T8>\
	idx_iv_cv<expression_type, const interleaved<T8> >\
		operator()(const interleaved<T8> &inter_row ,colon_void) \
		{ return idx_iv_cv<expression_type, const interleaved<T8> >(*this, inter_row); }\
	\
	template <typename E9>\
	idx_J_i<expression_type, E9>\
		operator()(const matrix_expression<E9> &row, size_type col)\
		{ return idx_J_i<expression_type, E9>(*this, row(), col); }\
	template <typename T8>\
	idx_range_i<expression_type, const range<T8> >\
		operator()(const range<T8> &range_row, size_type i) \
		{ return idx_range_i<expression_type, const range<T8> >(*this, range_row, i); }\
	template <typename T8>\
	idx_iv_i<expression_type, const interleaved<T8> >\
		operator()(const interleaved<T8> &interleaved_row, size_type i)\
		{ return idx_iv_i<expression_type, const interleaved<T8> >(*this, interleaved_row, i); }\
	template <typename E9>\
	idx_i_J<expression_type, E9>\
		operator()(size_type row, const matrix_expression<E9> &col)\
		{ return idx_i_J<expression_type, E9>(*this, row, col()); }\
	template <typename E8, typename E9>\
	idx_J_J<expression_type, E8, E9>\
		operator()(const matrix_expression<E8> &row, const matrix_expression<E9> &col)\
		{ return idx_J_J<expression_type, E8, E9>(*this, row(), col()); }\
	template <typename T8, typename E8>\
	idx_range_J<expression_type, const range<T8>, E8 >\
		operator()(const range<T8> &range_row, const matrix_expression<E8> &col)\
		{ return idx_range_J<expression_type, const range<T8>, E8 >(*this, range_row, col()); }\
	template <typename T8, typename E8>\
	idx_iv_J<expression_type, const interleaved<T8>, E8 >\
		operator()(const interleaved<T8> &interleaved_row, const matrix_expression<E8> &col)\
		{ return idx_iv_J<expression_type, const interleaved<T8>, E8 >(*this, interleaved_row, col()); }\
	template <typename T8>\
	idx_i_range<expression_type, const range<T8> >\
		operator()(size_type i, const range<T8> &range_col)\
		{ return idx_i_range<expression_type, const range<T8> >(*this, i, range_col); }\
	template <typename E8,typename T8>\
	idx_J_range<expression_type, E8, const range<T8> >\
		operator()(const matrix_expression<E8> &row, const range<T8> &range_col)\
		{ return idx_J_range<expression_type, E8, const range<T8> >(*this, row(), range_col); }\
	template <typename T8, typename T9>\
	idx_range_range<expression_type, const range<T8>, const range<T9> >\
		operator()(const range<T8> &range1, const range<T9> &range2)\
		{ return idx_range_range<expression_type, const range<T8>, const range<T9> >(*this, range1, range2); }\
	template <typename T8, typename T9>\
	idx_iv_range<expression_type, const interleaved<T8>, const range<T9> >\
		operator()(const interleaved<T8> &interleaved_row, const range<T9> &range_col)\
		{ return idx_iv_range<expression_type, const interleaved<T8>, const range<T9> >(*this, interleaved_row, range_col); }\
	template <typename T8>\
	idx_i_iv<expression_type, const interleaved<T8> >\
		operator()(size_type row, const interleaved<T8> &interleaved_col )\
		{ return idx_i_iv<expression_type, const interleaved<T8> >(*this, row, interleaved_col); }\
	template <typename E8, typename T8>\
	idx_J_iv<expression_type, E8, const interleaved<T8> >\
		operator()(const matrix_expression<E8> &row, const interleaved<T8> &interleaved_col)\
		{ return idx_J_iv<expression_type, E8, const interleaved<T8> >(*this, row(), interleaved_col); }\
	template <typename T8, typename T9>\
	idx_range_iv<expression_type, const range<T8>, const interleaved<T9> >\
		operator()(const range<T8> &range_row, const interleaved<T9> &interleaved_col)\
		{ return idx_range_iv<expression_type, const range<T8>, const interleaved<T9> >(*this, range_row, interleaved_col); }\
	template <typename T8, typename T9>\
	idx_iv_iv<expression_type, const interleaved<T8>, const interleaved<T9> >\
		operator()(const interleaved<T8> &interleaved_row, const interleaved<T9> &interleaved_col)\
		{ return idx_iv_iv<expression_type, const interleaved<T8>, const interleaved<T9> >(*this, interleaved_row, interleaved_col); }\
	\
	template <typename E9>\
	idx_I<expression_type, E9>\
		operator()(const matrix_expression<E9> &I)\
		{ return idx_I<expression_type, E9>(*this, I()); }\
	idx_cv<expression_type>\
		operator()(colon_void) \
		{ return idx_cv<expression_type>(*this); }\
	template <typename T1>\
	idx_range<expression_type, T1>\
		operator()(const range<T1> &range)\
		{ return idx_range<expression_type, T1>(*this, range); }\
	template <typename T1>\
	idx_iv<expression_type, T1>\
		operator()(const interleaved<T1> &interleaved)\
		{return idx_iv<expression_type, T1>(*this, interleaved); }

#define PROXY_ACCESS\
	CONST_PROXY_ACCESS;\
	NON_CONST_PROXY_ACCESS;

#define MATRIX_SELF_SCALAR_OP(op) \
const expression_type& operator op (const value_type &a) {\
	iterator i1 = begin(), end1 = end();\
	for(; i1!=end1; ++i1) *i1 op a;\
	return *this;\
}

#define MATRIX_SELF_MATRIX_OP(op) \
template <typename EXX> \
const expression_type& operator op (const matrix_expression<EXX> &e) {\
	assert(size(1) == e().size(1) && size(2) == e().size(2));\
	iterator i1 = begin(), end1 = end();\
	typename EXX::const_iterator i2 = e().begin(), end2 = e().end();\
	for(; i1!=end1; ++i1, ++i2) *i1 op *i2;\
	return *this;\
}

#define MATRIX_SELF_OP\
	MATRIX_SELF_SCALAR_OP(+=);\
	MATRIX_SELF_SCALAR_OP(-=);\
	MATRIX_SELF_SCALAR_OP(*=);\
	MATRIX_SELF_SCALAR_OP(/=);\
	MATRIX_SELF_SCALAR_OP(%=);\
	MATRIX_SELF_SCALAR_OP(&=);\
	MATRIX_SELF_SCALAR_OP(^=);\
	MATRIX_SELF_SCALAR_OP(|=);\
	MATRIX_SELF_SCALAR_OP(<<=);\
	MATRIX_SELF_SCALAR_OP(>>=);\
\
	MATRIX_SELF_MATRIX_OP(+=);\
	MATRIX_SELF_MATRIX_OP(-=);\
	MATRIX_SELF_MATRIX_OP(*=);\
	MATRIX_SELF_MATRIX_OP(/=);\
	MATRIX_SELF_MATRIX_OP(%=);\
	MATRIX_SELF_MATRIX_OP(&=);\
	MATRIX_SELF_MATRIX_OP(^=);\
	MATRIX_SELF_MATRIX_OP(|=);\
	MATRIX_SELF_MATRIX_OP(<<=);\
	MATRIX_SELF_MATRIX_OP(>>=);

template <typename E>
struct matrix_expression {
	typedef E expression_type;
	
	expression_type &operator()(void) {return *static_cast<expression_type *>(this);}
	const expression_type &operator()(void) const {return *static_cast<const expression_type *>(this);}

	template <typename RHS>
	const expression_type& operator=(const RHS &rhs) {
    assign((*this)(), rhs);
		return (*this)();
	}
};

template <typename E1, typename E2>
void assign(matrix_expression<E1> &lhs, const matrix_expression<E2> &rhs) {
  lhs().resize(rhs().size(1), rhs().size(2));
  std::copy(rhs().begin(), rhs().end(), lhs().begin());
}

template <typename E1>
void assign(matrix_expression<E1> &lhs, const typename E1::value_type &T) {
  std::fill(lhs().begin(), lhs().end(), T);
}

template <typename ARG, typename F>
struct unary_funcion : public matrix_expression<F> {
	unary_funcion(ARG e):m_e(e){}
	ARG m_e;
};

template <typename ARG1, typename ARG2, typename F>
struct binary_function : public matrix_expression<F> {
	binary_function(ARG1 e1, ARG2 e2):m_e1(e1), m_e2(e2){}
	ARG1 m_e1;
	ARG2 m_e2;
};

template <typename E, typename F>
struct matrix_each_element : public unary_funcion<E&, matrix_each_element<E, F> > {
	typedef matrix_each_element<E, F> expression_type;
	typedef typename E::value_type value_type;
	typedef typename E::size_type size_type;
	typedef value_type const_reference;
	typedef void reference;
    typedef void iterator;

	matrix_each_element(E &e):unary_funcion<E&, matrix_each_element<E, F> >(e){}

	size_type size(void) const {return this->m_e.size();}
	//! @param dim 0 return nrows, 1 return ncols
	size_type size(int dim) const {return this->m_e.size(dim);}
	
	struct const_iterator : public default_matrix_iterator<value_type>  {
		const_iterator(){}
		const_iterator(typename E::const_iterator i):m_i(i){}
		value_type operator*(void) const {return f(*m_i);}
		const_iterator operator++(void) {++m_i; return *this;}
		const_iterator operator+=(size_type n) {m_i+=n; return *this;}
		const_iterator operator+(size_type n) const {const_iterator i(m_i); i+=n; return i;}
		bool operator != (const const_iterator &i) const {return m_i != i.m_i;}
		bool operator == (const const_iterator &i) const {return !this->operator!=(i);}
		typename E::const_iterator m_i;
		F f;
	};

	const_iterator begin(void) const {return const_iterator(this->m_e.begin());}
	const_iterator end(void) const {return const_iterator(this->m_e.end());}

	value_type operator()(idx_type i) const {return F()(this->m_e(i));}
	value_type operator[](idx_type i) const {return F()(this->m_e(i));}
	value_type operator()(idx_type row, idx_type col) const {return F()(this->m_e1(row, col));}

	CONST_PROXY_ACCESS;
};

template<class E1, class E2, class F>
struct matrix_matrix_each_element : public binary_function<const E1&, const E2&, matrix_matrix_each_element<E1, E2, F> > {
    typedef matrix_matrix_each_element expression_type;
	typedef typename E1::size_type size_type;
	typedef typename E1::value_type const_reference;
	typedef typename E1::value_type value_type;
	// for proxy access
	typedef void reference;
	typedef void iterator;

	value_type operator()(idx_type i) const {return F()(this->m_e1(i), this->m_e2(i));}
	value_type operator[](idx_type i) const {return F()(this->m_e1(i), this->m_e2(i));}
	value_type operator()(idx_type row, idx_type col) const {return F()(this->m_e1(row, col), this->m_e2(row, col));}

	struct const_iterator : public default_matrix_iterator<value_type>  {
		typedef typename E1::const_iterator const_iterator1;
		typedef typename E2::const_iterator const_iterator2;
		const_iterator(){}
		const_iterator(const const_iterator1 &i1, const const_iterator2 &i2):m_i1(i1), m_i2(i2){}
		const_iterator operator++(void) {++m_i1; ++m_i2; return *this;}
		const_iterator operator+=(size_type n) {m_i1+=n; m_i2+=n; return *this;}
		const_iterator operator+(size_type n) const {const_iterator i(m_i1, m_i2); i+=n; return i;}
		bool operator!=(const const_iterator &i) const {return m_i1!=i.m_i1 || m_i2!=m_i2;}
		bool operator==(const const_iterator &i) const {return !this->operator!=(i);}
		value_type operator*(void) const {return f(*m_i1, *m_i2);}

		const_iterator1 m_i1;
		const_iterator2 m_i2;
		F f;
	};
	matrix_matrix_each_element(const E1 &e1, const E2 &e2)
		:binary_function<const E1&, const E2&, matrix_matrix_each_element<E1, E2, F> >(e1, e2) {
		assert(e1.size(1) == e2.size(1) && e1.size(2) == e2.size(2));
	}

	size_type size(void) const {return this->m_e1.size();}
	size_type size(int dim) const {return this->m_e1.size(dim);}

	const_iterator begin(void) const {return const_iterator(this->m_e1.begin(), this->m_e2.begin());}
	const_iterator end(void) const {return const_iterator(this->m_e1.end(), this->m_e2.end());}

	CONST_PROXY_ACCESS;
};

template<class E, class F>
struct matrix_scalar_function : public binary_function<const E&, typename E::value_type, matrix_scalar_function<E, F> > {
    typedef matrix_scalar_function<E, F> expression_type;
	typedef typename E::size_type size_type;
	typedef typename E::value_type const_reference;
	typedef typename E::value_type value_type;
	// for proxy access
	typedef void reference;
	typedef void iterator;

	struct const_iterator : public default_matrix_iterator<value_type>  {
		typedef typename E::const_iterator const_iterator1;
		const_iterator(){}
		const_iterator(const const_iterator1 &i, const value_type &a):m_i(i), m_a(a){}
		const_iterator operator++(void) {++m_i; return *this;}
		const_iterator operator+=(size_type n) {m_i+=n; return *this;}
		const_iterator operator +(size_type n) const {const_iterator i(m_i, m_a); i+=n; return i;}
		bool operator!=(const const_iterator &i) const {return m_i!=i.m_i;}
		bool operator==(const const_iterator &i) const {return !this->operator!=(i);}
		value_type operator*(void) const {return f(*m_i, m_a);}

		const_iterator1 m_i;
		value_type m_a;
		F f;
	};
	matrix_scalar_function(const E &e, typename E::value_type &a)
		:binary_function<const E&, typename E::value_type, matrix_scalar_function<E, F> >(e, a) {}

	size_type size(void) const {return this->m_e1.size();}
	size_type size(int dim) const {return this->m_e1.size(dim);}

	const_iterator begin(void) const {
		return const_iterator(this->m_e1.begin(), this->m_e2);
	}
	const_iterator end(void) const {return const_iterator(this->m_e1.end(), this->m_e2);}

	value_type operator()(idx_type i) const {return F()(this->m_e1(i), this->m_e2);}
	value_type operator[](idx_type i) const {return F()(this->m_e1[i], this->m_e2);}
	value_type operator()(idx_type i, idx_type j) const {return F()(this->m_e1(i, j), this->m_e2);}

	CONST_PROXY_ACCESS;
};

template<class E, class F>
struct scalar_matrix_function : public binary_function<typename E::value_type, const E&, scalar_matrix_function<E, F> > {
    typedef scalar_matrix_function<E, F> expression_type;
	typedef typename E::size_type size_type;
	typedef typename E::value_type const_reference;
	typedef typename E::value_type value_type;
	// for proxy access
	typedef void reference;
	typedef void iterator;

	struct const_iterator : public default_matrix_iterator<value_type>  {
		typedef typename E::const_iterator const_iterator1;
		const_iterator(){}
		const_iterator(const value_type &a, const const_iterator1 &i):m_i(i), m_a(a){}
		const_iterator operator++(void) {++m_i; return *this;}
		const_iterator operator+=(size_type n) {m_i+=n; return *this;}
		const_iterator operator +(size_type n) const {const_iterator i(m_i, m_a); i+=n; return i;}
		bool operator!=(const const_iterator &i) const {return m_i!=i.m_i;}
		bool operator==(const const_iterator &i) const {return !this->operator!=(i);}
		value_type operator*(void) const {return f(m_a, *m_i);}

		const_iterator1 m_i;
		value_type m_a;
		F f;
	};
	scalar_matrix_function(typename E::value_type &a, const E &e)
		:binary_function<typename E::value_type, const E&, scalar_matrix_function<E, F> >(a, e) {}

	size_type size(void) const {return this->m_e2.size();}
	size_type size(int dim) const {return this->m_e2.size(dim);}

	const_iterator begin(void) const {
		return const_iterator(this->m_e1, this->m_e2.begin());
	}
	const_iterator end(void) const {return const_iterator(this->m_e1, this->m_e2.end());}

	value_type operator()(idx_type i) const {return F()(this->m_e1, this->m_e2(i));}
	value_type operator[](idx_type i) const {return F()(this->m_e1, this->m_e2[i]);}
	value_type operator()(idx_type i, idx_type j) const {return F()(this->m_e1, this->m_e2(i, j));}

	CONST_PROXY_ACCESS;
};

template<class E, class F>
struct matrix_scalar_logic_function : public binary_function<const E&,
	typename E::value_type, matrix_scalar_logic_function<E, F> > {
    typedef matrix_scalar_logic_function<E, F> expression_type;
	typedef typename E::size_type size_type;
	typedef bool const_reference;
	typedef bool value_type;
	// for proxy access
	typedef void reference;
	typedef void iterator;

	struct const_iterator : public default_matrix_iterator<value_type>  {
		typedef typename E::const_iterator const_iterator1;
		const_iterator(){}
		const_iterator(const const_iterator1 &i, const typename E::value_type &a):m_i(i), m_a(a){}
		const_iterator operator++(void) {++m_i; return *this;}
		const_iterator operator+=(size_type n) {m_i+=n; return *this;}
		bool operator!=(const const_iterator &i) const {return m_i!=i.m_i;}
		bool operator==(const const_iterator &i) const {return !this->operator!=(i);}
		bool operator*(void) const {return f(*m_i, m_a);}

		const_iterator1 m_i;
		typename E::value_type m_a;
		F f;
	};
	matrix_scalar_logic_function(const E &e, typename E::value_type &a)
		:binary_function<const E&, typename E::value_type, matrix_scalar_logic_function<E, F> >(e, a) {}

	size_type size(void) const {return this->m_e1.size();}
	size_type size(int dim) const {return this->m_e1.size(dim);}

	const_iterator begin(void) const {
		return const_iterator(this->m_e1.begin(), this->m_e2);
	}
	const_iterator end(void) const {return const_iterator(this->m_e1.end(), this->m_e2);}

	value_type operator()(idx_type i) const {return F()(this->m_e1(i), this->m_e2);}
	value_type operator[](idx_type i) const {return F()(this->m_e1[i], this->m_e2);}
	value_type operator()(idx_type i, idx_type j) const {return F()(this->m_e1(i, j), this->m_e2);}

	CONST_PROXY_ACCESS;
};

}	// namespace matrix
}	// namespace zjucad

#endif
