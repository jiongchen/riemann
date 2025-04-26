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

#pragma once

#include <cassert>

#include "storage.h"

#include "matrix_expression.h"
#include "matrix_proxy.h"

namespace zjucad {namespace matrix {

#define default_container unbounded_array
struct column_major;
typedef column_major default_format;

template <typename T, typename F, typename A>
struct matrix_traits {
	typedef typename F::template matrix_type<T, A> matrix_type;
};

template <typename T, typename F=default_format, typename A=default_container<T>, bool tmp = false>
struct matrix : public matrix_expression<matrix<T, F, A, tmp> > {
	typedef typename matrix_traits<T, F, A>::matrix_type matrix_type;
	typedef matrix<T, F, A, tmp> expression_type;
	typedef matrix<T, F, A, tmp> E;
	typedef typename matrix_type::raw_data_type raw_data_type;

	typedef typename matrix_type::value_type value_type;
	typedef typename matrix_type::reference reference;
	typedef typename matrix_type::const_reference const_reference;
	typedef typename matrix_type::const_iterator const_iterator;
	typedef typename matrix_type::iterator iterator;
	typedef typename matrix_type::size_type size_type;

	// constructor
	matrix(){}
  explicit
	matrix(size_type size){m_matrix.resize(size, 1);}
	matrix(size_type row, size_type col){m_matrix.resize(row, col);}

	template <typename E>
	matrix(const matrix_expression<E> &e):m_matrix(){*this=e();}
	template <typename T1, typename F1, typename A1, bool tmp1>
	matrix(const matrix<T1, F1, A1, tmp1> &m):m_matrix(){*this=m;}
	matrix(matrix<T, F, A, true> &m):m_matrix(){*this=m;}

	// size
	size_type size(void) const {return m_matrix.size();} 
	size_type size(int dim) const {return (dim==1)?m_matrix.size(1):m_matrix.size(2);}
	void resize(size_type size) {m_matrix.resize(size, 1);}
	void resize(size_type row, size_type col) {m_matrix.resize(row, col);}

	// element access
	const_reference operator[](idx_type i) const {return m_matrix(i);}
	const_reference operator()(idx_type i) const {return m_matrix(i);}
	const_reference operator()(idx_type row, idx_type col) const {return m_matrix(row, col);}
	reference operator[](idx_type i) {return m_matrix(i);}
	reference operator()(idx_type i) {return m_matrix(i);}
	reference operator()(idx_type row, idx_type col) {return m_matrix(row, col);}

	PROXY_ACCESS;

	// iterator access
	const_iterator begin(void) const {return m_matrix.begin();}
	const_iterator end(void) const {return m_matrix.end();}
	iterator begin(void) {return m_matrix.begin();}
	iterator end(void) {return m_matrix.end();}

	template <typename RHS>
	const matrix& operator=(const RHS &e){
    assign((*this), e);
    return *this;
	}
	const matrix& operator=(matrix<T, F, A, true> &m){m_matrix.swap(m.m_matrix); return *this;}
	const matrix& operator=(value_type v){this->resize(1); (*this)[0] = v; return *this;}

	MATRIX_SELF_OP;

	const raw_data_type &data(void) const {return m_matrix.data();}
	raw_data_type &data(void) {return m_matrix.data();}

	matrix_type m_matrix;
};

template <typename OS, typename T, typename F, typename A, bool tmp>
void write(OS &os, const matrix<T, F, A, tmp> &m) {
	typename matrix<T, F, A, tmp>::size_type row = m.size(1), col = m.size(2);
	os.write((const char *)&row, sizeof(typename matrix<T, F, A, tmp>::size_type));
	os.write((const char *)&col, sizeof(typename matrix<T, F, A, tmp>::size_type));
	write(os, m.data());
}

template <typename IS, typename T, typename F, typename A, bool tmp>
int read(IS &is, matrix<T, F, A, tmp> &m) {
	typename matrix<T, F, A, tmp>::size_type row, col;
	is.read((char *)&row, sizeof(typename matrix<T, F, A, tmp>::size_type));
	if(is.fail())    return 1;
	is.read((char *)&col, sizeof(typename matrix<T, F, A, tmp>::size_type));
	if(is.fail())    return 2;
	m.resize(row, col);
	return read(is, m.data());
}

struct column_major {
	template <typename T, typename A>
	struct matrix_type {
		typedef zjucad::matrix::size_type size_type;
		typedef typename A::value_type value_type;
		typedef typename A::reference reference;
		typedef typename A::const_reference const_reference;
		typedef typename A::const_iterator const_iterator;
		typedef typename A::iterator iterator;
		typedef A raw_data_type;

		matrix_type():m_row(0), m_col(0){}

		// size
		size_type size(void) const {return m_row*m_col;} 
		size_type size(int dim) const {return (dim==1)?m_row:m_col;} 

		// element access
		const_reference operator()(idx_type i) const {return m_data[i];}
		const_reference operator()(idx_type row, idx_type col) const {return m_data[row+col*m_row];}
		reference operator()(idx_type i) {return m_data[i];}
		reference operator()(idx_type row, idx_type col) {return m_data[row+col*m_row];}

		// iterator access
		const_iterator begin(void) const {return m_data.begin();}
		const_iterator end(void) const {return m_data.end();}
		iterator begin(void) {return m_data.begin();}
		iterator end(void) {return m_data.end();}

		// used by matrix
		void resize(size_type row, size_type col) {
			m_data.resize(row*col);
			m_row = row; m_col = col;
		}
		void swap(matrix_type<T, A> &m) {
			m_data.swap(m.m_data);
			std::swap(m_row, m.m_row);
			std::swap(m_col, m.m_col);
		}

		// raw data access
		const A &data() const {return m_data;}
		A &data() {return m_data;}

		// data
		A m_data;
		size_type m_row, m_col;
	};
};

struct row_major {
	template <typename T, typename A>
	struct matrix_type {
		typedef zjucad::matrix::size_type size_type;
		typedef typename A::value_type value_type;
		typedef typename A::reference reference;
		typedef typename A::const_reference const_reference;
		typedef A raw_data_type;
		typedef matrix_type<T, A> expression_type;

		matrix_type():m_row(0), m_col(0){}

		// size
		size_type size(void) const {return m_row*m_col;} 
		size_type size(int dim) const {return (1==dim)?m_row:m_col;} 

		// element access
		const_reference operator()(idx_type i) const {return operator()(i%m_row, i/m_row);}
		const_reference operator()(idx_type row, idx_type col) const {return m_data[row*m_col+col];}
		reference operator()(idx_type i) {return operator()(i%m_row, i/m_row);}
		reference operator()(idx_type row, idx_type col) {return m_data[row*m_col+col];}

		// iterator access
		template <typename M, typename I>
		struct iterator_base : public default_matrix_iterator<value_type> {
			iterator_base():m_e(0){}
			iterator_base(M *e, size_type row, size_type col):m_row(row), m_col(col), m_e(e)
			{ m_i = m_e->data().begin(); m_i+=row*m_e->size(2)+col; }
			const_reference operator*() const {return *m_i;}
			iterator_base operator++(void) {
				++m_row;
				m_i+=m_e->size(2);
				if(m_row == m_e->size(1)) {
					m_row = 0;
					++m_col;
					m_i = m_e->data().begin();
					m_i+=m_col;
				}
				return *this;
			}
			iterator_base operator +=(offset_type n) {
				if(!n) return *this;
				offset_type r, c, t = m_row+n;
				r = t%m_e->size(1);
				c =(t/m_e->size(1));
				m_i+=(c+(r-m_row)*m_e->size(2));
				m_row = r;
				m_col += c;
				return *this;
			}
			iterator_base operator+(offset_type n) const {
				iterator_base i=*this;
				i+=n;
				return i;
			}
			template <typename I1>
			bool operator != (const I1 &i) const {assert(m_e == i.m_e); return m_row!=i.m_row || m_col != i.m_col;}
			template <typename I1>
			bool operator == (const I1 &i) const { return !this->operator!=(i);}

			I m_i;
			size_type m_row, m_col;
			M *m_e;
		};
		typedef iterator_base<const expression_type, typename A::const_iterator> const_iterator;
		struct iterator : public iterator_base<expression_type, typename A::iterator> {
			iterator(){}
			iterator(expression_type *e, size_type row, size_type col)
				:iterator_base<expression_type, typename A::iterator>(e, row, col){}
			reference operator*(void) {return *this->m_i;}
			iterator operator+(size_type n) const {
				iterator i=*this;
				i+=n;
				return i;
			}
		};

		const_iterator begin(void) const {return const_iterator(this, 0, 0);}
		const_iterator end(void) const {return const_iterator(this, 0, size(2));}
		iterator begin(void) {return iterator(this, 0, 0);}
		iterator end(void) {return iterator(this, 0, size(2));}

		// used by matrix
		void resize(size_type row, size_type col) {
			m_data.resize(row*col);
			m_row = row; m_col = col;
		}
		void swap(matrix_type<T, A> &m) {
			m_data.swap(m.m_data);
			std::swap(m_row, m.m_row);
			std::swap(m_col, m.m_col);
		}

		// raw data access
		const A &data() const {return m_data;}
		A &data() {return m_data;}

		// data
		A m_data;
		size_type m_row, m_col;
	};
};

struct bmp_order {
	template <typename T, typename A>
	struct matrix_type {
		typedef zjucad::matrix::size_type size_type;
		typedef typename A::value_type value_type;
		typedef typename A::reference reference;
		typedef typename A::const_reference const_reference;
		typedef A raw_data_type;
		typedef matrix_type<T, A> expression_type;

		matrix_type():m_row(0), m_col(0){}

		// size
		size_type size(void) const {return m_row*m_col;} 
		size_type size(int dim) const {return (1==dim)?m_row:m_col;} 

		// element access
		const_reference operator()(idx_type i) const {return (*this)(i%m_row, i/m_row);}
		const_reference operator()(idx_type row, idx_type col) const {return m_data[(m_row-row-1)*m_col+col];}
		reference operator()(idx_type i) {return (*this)(i%m_row, i/m_row);}
		reference operator()(idx_type row, idx_type col) {return m_data[(m_row-row-1)*m_col+col];}

		// iterator access
		template <typename M, typename I>
		struct iterator_base : public default_matrix_iterator<value_type> {
			iterator_base(){}
			iterator_base(M *e, size_type row, size_type col):m_row(row), m_col(col), m_e(e)
			{ m_i = m_e->data().begin(); m_i+=(m_e->size(1)-row-1)*m_e->size(2)+col; }
			const_reference operator*() const {return *m_i;}
			iterator_base operator++(void) {
				++m_row;
				m_i-=m_e->size(2);
				if(m_row == m_e->size(1)) {
					m_row = 0;
					++m_col;
					m_i = m_e->data().begin();
					m_i+=(m_e->size(1)-1)*m_e->size(2)+m_col;
				}
				return *this;
			}
			iterator_base operator +=(offset_type n) {
				if(!n) return *this;
				offset_type r, c, t = m_row+n;
				r = t%m_e->size(1);
				c =(t/m_e->size(1));
				m_i+=(c-(r-m_row)*m_e->size(2));
				m_row = r;
				m_col += c;
				return *this;
			}
			iterator_base operator+(offset_type n) const {
				iterator_base i=*this;
				i+=n;
				return i;
			}
			template <typename I1>
			bool operator != (const I1 &i) const {assert(m_e == i.m_e); return m_row!=i.m_row || m_col != i.m_col;}
			template <typename I1>
			bool operator == (const I1 &i) const { return !this->operator!=(i);}
			I m_i;
			size_type m_row, m_col;
			M *m_e;
		};
		typedef iterator_base<const expression_type, typename A::const_iterator> const_iterator;
		struct iterator : public iterator_base<expression_type, typename A::iterator> {
			iterator(){}
			iterator(expression_type *e, size_type row, size_type col)
				:iterator_base<expression_type, typename A::iterator>(e, row, col){}
			reference operator*(void) {return *this->m_i;}
			iterator operator+(size_type n) const {
				iterator i=*this;
				i+=n;
				return i;
			}
			operator const_iterator() const { return const_iterator(this->m_e, m_row, m_col); }
		};

		const_iterator begin(void) const {return const_iterator(this, 0, 0);}
		const_iterator end(void) const {return const_iterator(this, 0, size(2));}
		iterator begin(void) {return iterator(this, 0, 0);}
		iterator end(void) {return iterator(this, 0, size(2));}

		// used by matrix
		void resize(size_type row, size_type col) {
			m_data.resize(row*col);
			m_row = row; m_col = col;
		}
		void swap(matrix_type<T, A> &m) {
			m_data.swap(m.m_data);
			std::swap(m_row, m.m_row);
			std::swap(m_col, m.m_col);
		}

		// raw data access
		const A &data() const {return m_data;}
		A &data() {return m_data;}

		// data
		A m_data;
		size_type m_row, m_col;
	};
};

template <typename T>
struct default_tmp_matrix {
	typedef matrix<T, default_format, default_container<T>, false> type;
};

template <typename E>
inline typename default_tmp_matrix<typename E::value_type>::type
temp(const matrix_expression<E> &e) {
	return typename default_tmp_matrix<typename E::value_type>::type(e());
}

template <typename T>
struct identity_matrix : matrix_expression<identity_matrix<T> > {
public:
	typedef identity_matrix expression_type;
	typedef T value_type;
	typedef T const_reference;
	typedef void reference;
	typedef zjucad::matrix::size_type size_type;

	identity_matrix(size_type size):m_size(size){}

	size_type size() const {return m_size*m_size;}
	size_type size(int dim) const { return m_size;}

	struct const_iterator : public default_matrix_iterator<value_type>  {
		const_iterator(){}
		const_iterator(size_type size, size_type row, size_type col):m_size(size), m_row(row), m_col(col){}
		const_iterator operator++(void) {
			++m_row;
			if(m_row == m_size) {
				m_row = 0;
				++m_col;
			}
			return *this;
		}
		const_iterator operator+=(offset_type n) {
			offset_type t = m_row+n;
			m_row = t%m_size;
			m_col+= t/m_size;
			return *this;
		}
		value_type operator*(void) const {return m_row == m_col;}
		bool operator != (const const_iterator &i) const 
		{assert(m_size == i.m_size); return m_row!=i.m_row || m_col != i.m_col;}
		bool operator == (const const_iterator &i) const 
		{ return !this->operator!=(i);}

		size_type m_size, m_row, m_col;
	};
	typedef const_iterator iterator;

	const_iterator begin(void) const {return const_iterator(m_size, 0, 0);}
	const_iterator end(void) const {return const_iterator(m_size, 0, m_size);}

	CONST_PROXY_ACCESS

	size_type m_size;
};

template <typename T>
inline identity_matrix<T> eye(size_type n) {return identity_matrix<T>(n);}

//inline identity_matrix<int> eye(size_type n) {return identity_matrix<int>(n);}

template <typename T, typename F>
struct value_matrix : public matrix_expression<value_matrix<T, F> > {
	typedef value_matrix expression_type;
	typedef T value_type;
	typedef T const_reference;
	typedef void reference;
	typedef zjucad::matrix::size_type size_type;

	value_matrix(size_type row, size_type col):m_row(row), m_col(col){}

	size_type size() const {return m_row*m_col;}
	size_type size(int dim) const {return (1==dim)?m_row:m_col;}

	struct const_iterator : public default_matrix_iterator<value_type>  {
		const_iterator(){}
		const_iterator(size_type i):m_i(i){}
		const_iterator operator++(void) {++m_i; return *this;}
		const_iterator operator+=(offset_type n) {m_i+=n; return *this;}
		value_type operator*(void) const {return static_cast<value_type>(F::value());}
		bool operator != (const const_iterator &i) const {return m_i != i.m_i;}
		bool operator == (const const_iterator &i) const { return !this->operator!=(i);}
		size_type m_i;
	};
	typedef const_iterator iterator;

	const_iterator begin(void) const {return const_iterator(0);}
	const_iterator end(void) const {return const_iterator(size());}

	CONST_PROXY_ACCESS

	size_type m_row, m_col;
};

struct zero_functor {
	static int value(void) {return 0;}
};

struct one_functor {
	static int value(void) {return 1;}
};

struct rand_functor {
	static double value(void) {return ::rand()/double(RAND_MAX);}
};

template <typename T>
inline value_matrix<T, zero_functor> zeros(size_type n) {return value_matrix<T, zero_functor>(n, n);}
template <typename T>
inline value_matrix<T, zero_functor> zeros(size_type row, size_type col) {return value_matrix<T, zero_functor>(row, col);}
inline value_matrix<int, zero_functor> zeros(size_type n) {return value_matrix<int, zero_functor>(n, n);}
inline value_matrix<int, zero_functor> zeros(size_type row, size_type col) {return value_matrix<int, zero_functor>(row, col);}

template <typename T>
inline value_matrix<T, one_functor> ones(size_type n) {return value_matrix<T, one_functor>(n, n);}
template <typename T>
inline value_matrix<T, one_functor> ones(size_type row, size_type col) {return value_matrix<T, one_functor>(row, col);}
inline value_matrix<int, one_functor> ones(size_type n) {return value_matrix<int, one_functor>(n, n);}
inline value_matrix<int, one_functor> ones(size_type row, size_type col) {return value_matrix<int, one_functor>(row, col);}

template <typename T>
inline value_matrix<T, rand_functor> rand(size_type n) {return value_matrix<T, rand_functor>(n, n);}
template <typename T>
inline value_matrix<T, rand_functor> rand(size_type row, size_type col) {return value_matrix<T, rand_functor>(row, col);}
inline value_matrix<double, rand_functor> rand(size_type n) {return value_matrix<double, rand_functor>(n, n);}
inline value_matrix<double, rand_functor> rand(size_type row, size_type col) {return value_matrix<double, rand_functor>(row, col);}

}	// namespace matrix
}	// namespace zjucad

#include "operation.h"
