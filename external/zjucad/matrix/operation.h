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
#include <cmath>
#include <functional>
#include <algorithm>
#include <vector>

#include <boost/shared_ptr.hpp>
#include <boost/type_traits/remove_const.hpp>

using std::plus;
using std::minus;
using std::multiplies;
using std::divides;
using std::logical_not;
using std::negate;
using std::modulus;

#include "matrix.h"

namespace zjucad { namespace matrix {

#define FUNC_EACH_ELEMENT(f)\
FUNC_SCALAR_FUNCTOR(f);\
template <typename E>\
matrix_each_element<const E, f##_scalar_functor<typename E::value_type> >\
f(const matrix_expression<E> &e) {\
	return matrix_each_element<const E, f##_scalar_functor<typename E::value_type> >(e());\
}

FUNC_EACH_ELEMENT(abs);
FUNC_EACH_ELEMENT(acos);
FUNC_EACH_ELEMENT(acosf);
FUNC_EACH_ELEMENT(asin);
FUNC_EACH_ELEMENT(asinf);
FUNC_EACH_ELEMENT(atan);
FUNC_EACH_ELEMENT(atanf);
FUNC_EACH_ELEMENT(ceil);
FUNC_EACH_ELEMENT(ceilf);
FUNC_EACH_ELEMENT(cos);
FUNC_EACH_ELEMENT(cosf);
FUNC_EACH_ELEMENT(cosh);
FUNC_EACH_ELEMENT(coshf);
FUNC_EACH_ELEMENT(exp);
FUNC_EACH_ELEMENT(expf);
FUNC_EACH_ELEMENT(fabs);
FUNC_EACH_ELEMENT(fabsf);
FUNC_EACH_ELEMENT(floor);
FUNC_EACH_ELEMENT(floorf);
FUNC_EACH_ELEMENT(labs);
FUNC_EACH_ELEMENT(log);
FUNC_EACH_ELEMENT(logf);
FUNC_EACH_ELEMENT(log10);
FUNC_EACH_ELEMENT(log10f);
FUNC_EACH_ELEMENT(sin);
FUNC_EACH_ELEMENT(sinf);
FUNC_EACH_ELEMENT(sinh);
FUNC_EACH_ELEMENT(sinhf);
FUNC_EACH_ELEMENT(sqrt);
FUNC_EACH_ELEMENT(tan);
FUNC_EACH_ELEMENT(tanf);
FUNC_EACH_ELEMENT(tanh);
FUNC_EACH_ELEMENT(tanhf);
#ifdef _MSC_VER
FUNC_EACH_ELEMENT(_chgsign);
FUNC_EACH_ELEMENT(_finite);
#endif

template <typename E>
matrix_each_element<const E, negate<typename E::value_type> >
operator-(const matrix_expression<E> &e) {
	return matrix_each_element<const E, negate<typename E::value_type> >(e());
}

template <typename E1, typename E2>
matrix_matrix_each_element<const E1, const E2, plus<typename E1::value_type> >
operator+(const matrix_expression<E1> &e1, const matrix_expression<E2> &e2) {
	return matrix_matrix_each_element<const E1, const E2, plus<typename E1::value_type> >(e1(), e2());
}

template <typename E1, typename E2>
matrix_matrix_each_element<const E1, const E2, minus<typename E1::value_type> >
operator-(const matrix_expression<E1> &e1, const matrix_expression<E2> &e2) {
	return matrix_matrix_each_element<const E1, const E2, minus<typename E1::value_type> >(e1(), e2());
}

template <typename E>
matrix_scalar_function<const E, plus<typename E::value_type> >
operator+(const matrix_expression<E> &e1, typename E::value_type a) {
	return matrix_scalar_function<const E, plus<typename E::value_type> >(e1(), a);
}

template <typename E>
scalar_matrix_function<const E, plus<typename E::value_type> >
operator+(typename E::value_type a, const matrix_expression<E> &e1) {
	return scalar_matrix_function<const E, plus<typename E::value_type> >(a, e1());
}

template <typename E>
matrix_scalar_function<const E, minus<typename E::value_type> >
operator-(const matrix_expression<E> &e1, typename E::value_type a) {
	return matrix_scalar_function<const E, minus<typename E::value_type> >(e1(), a);
}

template <typename E>
scalar_matrix_function<const E, minus<typename E::value_type> >
operator-(typename E::value_type a, const matrix_expression<E> &e1) {
	return scalar_matrix_function<const E, minus<typename E::value_type> >(a, e1());
}

template <typename E>
matrix_scalar_function<const E, multiplies<typename E::value_type> >
operator*(const matrix_expression<E> &e1, typename E::value_type a) {
	return matrix_scalar_function<const E, multiplies<typename E::value_type> >(e1(), a);
}

template <typename E>
scalar_matrix_function<const E, multiplies<typename E::value_type> >
operator*(typename E::value_type a, const matrix_expression<E> &e1) {
	return scalar_matrix_function<const E, multiplies<typename E::value_type> >(a, e1());
}

template <typename E>
matrix_scalar_function<const E, divides<typename E::value_type> >
operator/(const matrix_expression<E> &e1, typename matrix_scalar_type_traits<E>::arg_type2 a) {
	return matrix_scalar_function<const E, divides<typename E::value_type> >(e1(), a);
}

template <typename E>
scalar_matrix_function<const E, divides<typename E::value_type> >
operator/(typename matrix_scalar_type_traits<E>::arg_type2 a, const matrix_expression<E> &e1) {
	return scalar_matrix_function<const E, divides<typename E::value_type> >(a, e1());
}

template <typename E>
matrix_scalar_function<const E, modulus<typename E::value_type> >
operator%(const matrix_expression<E> &e1, typename matrix_scalar_type_traits<E>::arg_type2 a) {
	return matrix_scalar_function<const E, modulus<typename E::value_type> >(e1(), a);
}

template <typename E>
scalar_matrix_function<const E, modulus<typename E::value_type> >
operator%(typename matrix_scalar_type_traits<E>::arg_type2 a, const matrix_expression<E> &e1) {
	return scalar_matrix_function<const E, modulus<typename E::value_type> >(a, e1());
}

template <typename E>
matrix_scalar_function<const E, pow_functor<typename E::value_type, typename E::value_type> >
pow(const matrix_expression<E> &e1, typename E::value_type a) {
	return matrix_scalar_function<const E, pow_functor<typename E::value_type, typename E::value_type> >(e1(), a);
}

template <typename E>
matrix_scalar_logic_function<const E, std::equal_to<typename E::value_type> >
operator==(const matrix_expression<E> &e1, typename E::value_type a) {
	return matrix_scalar_logic_function<const E, std::equal_to<typename E::value_type> >(e1(), a);
}

template <typename E>
matrix_scalar_logic_function<const E, std::greater_equal<typename E::value_type> >
operator>=(const matrix_expression<E> &e1, typename E::value_type a) {
	return matrix_scalar_logic_function<const E, std::greater_equal<typename E::value_type> >(e1(), a);
}

template <typename E>
matrix_scalar_logic_function<const E, std::greater<typename E::value_type> >
operator>(const matrix_expression<E> &e1, typename E::value_type a) {
	return matrix_scalar_logic_function<const E, std::greater<typename E::value_type> >(e1(), a);
}

template <typename E>
matrix_scalar_logic_function<const E, std::less_equal<typename E::value_type> >
operator<=(const matrix_expression<E> &e1, typename E::value_type a) {
	return matrix_scalar_logic_function<const E, std::less_equal<typename E::value_type> >(e1(), a);
}

template <typename E>
matrix_scalar_logic_function<const E, std::less<typename E::value_type> >
operator<(const matrix_expression<E> &e1, typename E::value_type a) {
	return matrix_scalar_logic_function<const E, std::less<typename E::value_type> >(e1(), a);
}

template <typename E>
matrix_scalar_logic_function<const E, std::not_equal_to<typename E::value_type> >
operator!=(const matrix_expression<E> &e1, typename E::value_type a) {
	return matrix_scalar_logic_function<const E, std::not_equal_to<typename E::value_type> >(e1(), a);
}

template <typename E>
matrix_each_element<const E, logical_not<typename E::value_type> >
operator!(const matrix_expression<E> &e) {
	return matrix_each_element<const E, logical_not<typename E::value_type> >(e());
}

template <typename E1, typename E2>
struct matrix_matrix_multiplies : public binary_function<const E1&, const E2&, matrix_matrix_multiplies<E1, E2> > {
    typedef matrix_matrix_multiplies expression_type;
	typedef typename E1::size_type size_type;
	typedef typename E1::value_type const_reference;
	typedef void reference;
	typedef typename E1::value_type value_type;

	struct const_iterator : public default_matrix_iterator<value_type>  {
		typedef typename E1::const_iterator const_iterator1;
		typedef typename E2::const_iterator const_iterator2;
		const_iterator(){}
		const_iterator(size_type row, size_type col, const expression_type *m)
			:m_row(row), m_col(col), m_mult(m) {
			m_i1 = m->m_e1.begin();
			m_i1+=row;
			m_i2 = m->m_e2.begin();
			m_i2+=col*m->m_e2.size(1);
		}
		const_iterator operator++(void) {
			++m_row;
			++m_i1;
			if(m_row == m_mult->size(1)) {
				m_row = 0;
				m_i1 = m_mult->m_e1.begin();
				m_i2+=m_mult->m_e2.size(1);
				++m_col;
			}
			return *this;
		}
		const_iterator operator+=(offset_type n) {
			if(!n) return *this;
			if(n == m_mult->m_e1.size(1)) {
				m_i2+=m_mult->m_e2.size(1);
				++m_col;
			}
			else {
				offset_type t = m_row+n, r, c;
				r = t % m_mult->size(1);
				c = t / m_mult->size(1);
				m_i1 = m_mult->m_e1.begin();
				m_i1 += r;
				m_i2 += c*m_mult->m_e2.size(1);
				m_row = r;
				m_col +=c;
			}
			return *this;
		}
		value_type operator*(void) const {
			const_iterator1 i1 = m_i1;
			const_iterator2 i2 = m_i2;
			typename boost::remove_const<value_type>::type r = 0;
			for(size_type i = 0; i < m_mult->m_e1.size(2); ++i, i1+=m_mult->m_e1.size(1), ++i2)
				r += *i1 * *i2;
			return r;
		}
		bool operator != (const const_iterator &i) const {return m_mult->size() && (m_row!=i.m_row || m_col!=i.m_col);}
		bool operator == (const const_iterator &i) const {return !this->operator!=(i);}

		const_iterator1 m_i1;
		const_iterator2 m_i2;
		size_type m_row, m_col;
		const expression_type *m_mult;
	};
	matrix_matrix_multiplies(const E1 &e1, const E2 &e2)
		:binary_function<const E1&, const E2&, matrix_matrix_multiplies<E1, E2> >(e1, e2) {
		assert(e1.size(2) == e2.size(1));
	}

	size_type size(void) const {return size(1)*size(2);}
	size_type size(int dim) const {return (1==dim)?this->m_e1.size(1):this->m_e2.size(2);}

	const_iterator begin(void) const {return const_iterator(0, 0, this);}
	const_iterator end(void) const {return const_iterator(0, size(2), this);}
};

template <typename E1, typename E2>
matrix_matrix_multiplies<const E1, const E2>
operator*(const matrix_expression<E1> &e1, const matrix_expression<E2> &e2) {
	return matrix_matrix_multiplies<const E1, const E2>(e1(), e2());
}

template <typename E1, typename E2>
matrix_matrix_each_element<const E1, const E2, multiplies<typename E1::value_type> >
ele_prod(const matrix_expression<E1> &e1, const matrix_expression<E2> &e2) {
	return matrix_matrix_each_element<const E1, const E2, multiplies<typename E1::value_type> >(e1(), e2());
}

template <typename E>
struct matrix_trans : public matrix_expression<matrix_trans<E> > {
	typedef typename E::const_reference const_reference;
	typedef typename E::reference reference;
	typedef typename E::size_type size_type;
	typedef typename E::value_type value_type;
	typedef matrix_trans<E> expression_type;

	matrix_trans(E& e):m_e(e){}

	size_type size(void) const {return m_e.size();}
	size_type size(int dim) const {return (1==dim)?m_e.size(2):m_e.size(1);}

	// element access
	const_reference operator()(idx_type i) const {return m_e(i/m_e.size(1), i%m_e.size(1));}
	const_reference operator[](idx_type i) const {return m_e(i/m_e.size(1), i%m_e.size(1));}
	const_reference operator()(idx_type row, idx_type col) const {return m_e(col, row);}
	reference operator()(idx_type i) {return m_e(i/m_e.size(1), i%m_e.size(1));}
	reference operator[](idx_type i) {return m_e(i/m_e.size(1), i%m_e.size(1));}
	reference operator()(idx_type row, idx_type col) {return m_e(col, row);}

	template <typename M, typename I>
	struct iterator_base : public default_matrix_iterator<value_type> {
		iterator_base():m_e(0){}
		iterator_base(M *e, size_type row, size_type col):m_row(row), m_col(col), m_e(e) 
		{ m_i = m_e->begin(); m_i+=row*m_e->size(1)+col; }
		const_reference operator*() const {return *m_i;}
		iterator_base operator++(void) {
			++m_row;
			m_i+=m_e->size(1);
			if(m_row == m_e->size(2)) {
				m_row = 0;
				++m_col;
				m_i = m_e->begin();
				m_i+=m_col;
			}
			return *this;
		}
		iterator_base operator +=(offset_type n) {
			if(!n) return *this;
			offset_type r, c, t = m_row+n;
			r = t%m_e->size(2);
			c =(t/m_e->size(2));
			m_i+=(c+(r-m_row)*m_e->size(1));
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
		bool operator == (const I1 &i) const {return !this->operator!=(i);}
		I m_i;
		size_type m_row, m_col;
		M *m_e;
	};
	typedef iterator_base<const E, typename E::const_iterator> const_iterator;

	struct iterator : public iterator_base<E, typename E::iterator> {
		iterator(){}
		iterator(E *e, size_type row, size_type col):iterator_base<E, typename E::iterator>(e, row, col){}
		operator const_iterator() const {return const_iterator(this->m_e, this->m_row, this->m_col);}
		reference operator*(void) {return *this->m_i;}
		iterator operator+(size_type n) const {
			iterator i=*this;
			i+=n;
			return i;
		}
	};

	const_iterator begin(void) const {return const_iterator(&m_e, 0, 0);}
	const_iterator end(void) const {return const_iterator(&m_e, 0, m_e.size(1));}
	iterator begin(void) {return iterator(&m_e, 0, 0);}
	iterator end(void) {return iterator(&m_e, 0, m_e.size(1));}

	PROXY_ACCESS;

	E &m_e;
};

template <typename E>
const matrix_trans<const E> trans(const matrix_expression<E> &e) { return matrix_trans<const E>(e());}
template <typename E>
matrix_trans<E> trans(matrix_expression<E> &e) { return matrix_trans<E>(e());}

template <typename T, typename E>
struct matrix_convert : public unary_funcion<E&, matrix_convert<T, E> > {
	typedef matrix_convert<T, E> expression_type;
	typedef T value_type;
	typedef typename E::size_type size_type;
	typedef value_type const_reference;
	typedef void reference;
	typedef void iterator;

	matrix_convert(E &e):unary_funcion<E&, matrix_convert<T, E> >(e){}

	size_type size(void) const {return this->m_e.size();}
	size_type size(int dim) const {return this->m_e.size(dim);}
	
	struct const_iterator : public default_matrix_iterator<value_type>  {
		const_iterator(){}
		const_iterator(typename E::const_iterator i):m_i(i){}
		value_type operator*(void) const {return static_cast<T>(*m_i);}
		const_iterator operator++(void) {++m_i; return *this;}
		const_iterator operator+=(size_type n) {m_i+=n; return *this;}
		const_iterator operator+(size_type n) const {const_iterator i(m_i); i+=n; return i;}
		bool operator != (const const_iterator &i) const {return m_i != i.m_i;}
		bool operator == (const const_iterator &i) const {return !this->operator!=(i);}
		typename E::const_iterator m_i;
	};

	const_iterator begin(void) const {return const_iterator(this->m_e.begin());}
	const_iterator end(void) const {return const_iterator(this->m_e.end());}

	value_type operator[](idx_type i) const {return static_cast<T>(this->m_e[i]);}
	value_type operator()(idx_type i) const {return static_cast<T>(this->m_e(i));}
	value_type operator()(idx_type row, idx_type col) const {return static_cast<T>(this->m_e(row, col));}

	CONST_PROXY_ACCESS;
};

template <typename T, typename E>
const matrix_convert<T, const E> to(const matrix_expression<E> &e) { return matrix_convert<T, const E>(e());}

template <typename E>
typename E::value_type norm(const matrix_expression<E> &e) {
	typename E::const_iterator i = e().begin(), end = e().end();
	typename E::value_type r = 0;
	for(; i != end; ++i) r+=(*i)*(*i);
	return ::sqrt(r);
}

template <typename E>
typename E::value_type sqnorm(const matrix_expression<E> &e) {
	typename E::const_iterator i = e().begin(), end = e().end();
	typename E::value_type r = 0;
	for(; i != end; ++i) r+=(*i)*(*i);
	return r;
}

template <typename T, typename E>
T norm(const matrix_expression<E> &e) {
	typename E::const_iterator i = e().begin(), end = e().end();
	T r = 0;
	for(; i != end; ++i) r+=(*i)*(*i);
	return ::sqrt(r);
}

template <typename E1, typename E2>
typename E1::value_type dot(const matrix_expression<E1> &e1, const matrix_expression<E2> &e2) {
	assert(e1().size() == e2().size());
	typename E1::const_iterator i1 = e1().begin(), end1 = e1().end();
	typename E2::const_iterator i2 = e2().begin();
	typename E1::value_type r = 0;
	for(; i1 != end1; ++i1, ++i2) r+=(*i1)*(*i2);
	return r;
}

template <typename E, typename E1, typename E2>
matrix_expression<E> &
cross(matrix_expression<E> *r, const matrix_expression<E1> &e1, const matrix_expression<E2> &e2) {
	assert(r && e1().size() == 3 && e2().size() == 3);
	(*r)().resize(3);
	(*r)()(0) = e1()(1)*e2()(2)-e1()(2)*e2()(1);
	(*r)()(1) = e1()(2)*e2()(0)-e1()(0)*e2()(2);
	(*r)()(2) = e1()(0)*e2()(1)-e1()(1)*e2()(0);
	return *r;
}

template <typename E1, typename E2>
typename default_tmp_matrix<typename E1::value_type>::type
cross(const matrix_expression<E1> &e1, const matrix_expression<E2> &e2) {
	typename default_tmp_matrix<typename E1::value_type>::type r(3);
	cross(&r, e1, e2);
	return r;
}

template <typename T, typename E>
T sum(const matrix_expression<E> &e) {
	T r = 0;
	typename E::const_iterator i = e().begin(), end = e().end();
	for(; i!=end; ++i) r+=*i;
	return r;
}

template <typename E>
typename E::value_type sum(const matrix_expression<E> &e) {
	return sum<typename E::value_type>(e());
}
/*
//matlab arbitary argument length
template <typename VT, typename CT>
struct matlab_arbitary_argument_length
{
	template <typename E>
	CT& operator,(const matrix_expression<E> &e){
        vec_.push_back( boost::shared_ptr<op_abstract<VT> >(new op<E,typename E::value_type>( &(e()) )) );
		return *((CT*)this);
	}

	virtual typename default_tmp_matrix<VT>::type fetch(){
		assert(false);
		return default_tmp_matrix<VT>::type(0);
	}
	
protected:

	template <typename T>
	struct op_abstract{
		virtual size_type size() const =0;
		virtual size_type size(size_type dim) const =0;
		virtual T operator()(size_type i) const =0;
	};

	template <typename E, typename T>
	struct op : public op_abstract<T>{
		op(const E* proxy):proxy_(proxy){assert(proxy);}
		virtual size_type size() const{return proxy_->size();}
		virtual size_type size(size_type dim) const{return proxy_->size(dim);}
		virtual T operator()(size_type i) const{return proxy_->operator()(i);}
		const E *proxy_;
	};

	std::vector<boost::shared_ptr<op_abstract<VT> > > vec_;
	
};

//matlab blkdiag
template <typename VT>
struct blkdiag : public matlab_arbitary_argument_length<VT,blkdiag<VT> >
{
    typedef matlab_arbitary_argument_length<VT,blkdiag<VT> > base_class;
    typedef typename base_class::op_abstract<VT> OP_TYPE;

	virtual typename default_tmp_matrix<VT>::type fetch(){
        assert(!base_class::vec_.empty());
		
		size_type sz_row=0,sz_col=0;
        for(std::vector<boost::shared_ptr<base_class::op_abstract<VT> > >::iterator
            begin=base_class::vec_.begin(),end=base_class::vec_.end();begin != end;begin++){
			sz_row+=(*begin)->size(1);
			sz_col+=(*begin)->size(2);
		}

		default_tmp_matrix<VT>::type ret(sz_row,sz_col);

		size_type csize1=0,csize2=0;
        for(std::vector<boost::shared_ptr<OP_TYPE> >::iterator
            begin=base_class::vec_.begin(),end=base_class::vec_.end();begin != end;begin++){
			size_type size1=(*begin)->size(1),size2=(*begin)->size(2);
			
			size_type id=0;
			for(size_type i=0;i < size2;i++)
				for(size_type j=0;j < size1;j++)
					ret(csize1+j,csize2+i)=(*begin)->operator()(id++);

			csize1+=size1;
			csize2+=size2;
		}
		
        base_class::vec_.swap(std::vector<boost::shared_ptr<OP_TYPE> >());

		return ret;
	}
};

//matlab cat
template <typename VT>
struct cat : public matlab_arbitary_argument_length<VT,cat<VT> >
{
	typename default_tmp_matrix<VT>::type fetch1(){
		assert(!vec_.empty());

		size_type sz_row=0,sz_col=vec_[0]->size(2);
		for(std::vector<boost::shared_ptr<op_abstract<VT> > >::iterator 
			begin=vec_.begin(),end=vec_.end();begin != end;begin++){
			sz_row+=(*begin)->size(1);
			assert(sz_col == (*begin)->size(2));
		}

		default_tmp_matrix<VT>::type ret(sz_row,sz_col);

		size_type csize1=0;
		for(std::vector<boost::shared_ptr<op_abstract<VT> > >::iterator 
			begin=vec_.begin(),end=vec_.end();begin != end;begin++){
			size_type size1=(*begin)->size(1);
			
			size_type id=0;
			for(size_type i=0;i < sz_col;i++)
				for(size_type j=0;j < size1;j++)
					ret(csize1+j,i)=(*begin)->operator()(id++);

			csize1+=size1;
		}
		
		vec_.swap(std::vector<boost::shared_ptr<op_abstract<VT> > >());

		return ret;
	}

	typename default_tmp_matrix<VT>::type fetch2(){
		assert(!vec_.empty());
		
		size_type sz_row=vec_[0]->size(1),sz_col=0;
		for(std::vector<boost::shared_ptr<op_abstract<VT> > >::iterator 
			begin=vec_.begin(),end=vec_.end();begin != end;begin++){
			assert(sz_row == (*begin)->size(1));
			sz_col+=(*begin)->size(2);
		}

		default_tmp_matrix<VT>::type ret(sz_row,sz_col);

		size_type csize2=0;
		for(std::vector<boost::shared_ptr<op_abstract<VT> > >::iterator 
			begin=vec_.begin(),end=vec_.end();begin != end;begin++){
			size_type size1=(*begin)->size(1),size2=(*begin)->size(2);
			
			size_type id=0;
			for(size_type i=0;i < size2;i++)
				for(size_type j=0;j < size1;j++)
					ret(j,csize2+i)=(*begin)->operator()(id++);

			csize2+=size2;
		}
		
		vec_.swap(std::vector<boost::shared_ptr<op_abstract<VT> > >());

		return ret;
	}

};
//matlab circshift
template <typename E>
typename default_tmp_matrix<typename E::value_type>::type 
circshift(const matrix_expression<E> &e,size_type row_shift,size_type col_shift = 0)
{
	const size_type sz_row=e().size(1);
	const size_type sz_col=e().size(2);

	row_shift=sz_row-row_shift;
	if(row_shift < 0)row_shift=(row_shift%sz_row+sz_row)%sz_row;
	
	col_shift=sz_col-col_shift;
	if(col_shift < 0)col_shift=(col_shift%sz_col+sz_col)%sz_col;

	default_tmp_matrix<typename E::value_type>::type ret(sz_row,sz_col);
	
	for(size_type i=0;i < sz_col;i++)
		for(size_type j=0;j < sz_row;j++)
			ret(j,i)=e()((j+row_shift)%sz_row,(i+col_shift)%sz_col);

	return ret;
}

//matlab diag
template <typename E>
typename default_tmp_matrix<typename E::value_type>::type
diag_gen(const matrix_expression<E> &e, size_type k = 0){
	const size_type sz=e().size();
	const size_type sz_out=sz+std::abs(k);
	const size_type sz_off=sz_out+1;
	
	assert(sz_out);

	default_tmp_matrix<E::value_type>::type out(sz_out,sz_out);
	default_tmp_matrix<E::value_type>::type::iterator out_begin=out.begin()+(k>0 ? k*sz_out : -k);
	E::const_iterator in_begin=e().begin(),in_end=e().end();

	while(in_begin != in_end)
	{
		*out_begin=*in_begin;
		out_begin+=sz_off;++in_begin;
	}

	return out;
}

template <typename E>
typename default_tmp_matrix<typename E::value_type>::type
diag_get(const matrix_expression<E> &e, size_type k = 0){
	const size_type sz_row=e().size(1),sz_col=e().size(2);
	const size_type sz_out=k>0 ? std::min(sz_row,sz_col-k) : std::min(sz_row+k,sz_col);
	const size_type sz_off=sz_row+1;
	
	assert(sz_out > 0);

	default_tmp_matrix<E::value_type>::type out(sz_out,1);
	default_tmp_matrix<E::value_type>::type::iterator out_begin=out.begin(),out_end=out.end();
	E::const_iterator in_begin=e().begin()+(k>0 ? k*sz_row : -k);

	while(out_begin != out_end)
	{
		*out_begin=*in_begin;
		++out_begin;in_begin+=sz_off;
	}

	return out;
}

//matlab flipdim
template <typename E>
typename default_tmp_matrix<typename E::value_type>::type
flipdim(const matrix_expression<E> &e, size_type dim = 1)
{
	const size_type sz_row=e().size(1);
	const size_type sz_col=e().size(2);
	const size_type sz_row_1=sz_row-1;
	const size_type sz_col_1=sz_col-1;

	default_tmp_matrix<typename E::value_type>::type ret(sz_row,sz_col);
	
	if(dim == 1){
		for(size_type i=0;i < sz_col;i++)
			for(size_type j=0;j < sz_row;j++)
				ret(j,i)=e()(sz_row_1-j,i);
	}else{
		assert(dim == 2);
		for(size_type i=0;i < sz_col;i++)
			for(size_type j=0;j < sz_row;j++)
				ret(j,i)=e()(j,sz_col_1-i);
	}

	return ret;
}

//matlab fliplr
template <typename E>
typename default_tmp_matrix<typename E::value_type>::type
fliplr(const matrix_expression<E> &e){return flipdim(e,2);}

//matlab flipud
template <typename E>
typename default_tmp_matrix<typename E::value_type>::type
flipud(const matrix_expression<E> &e){return flipdim(e,1);}

//matlab horzcat
template <typename VT>
struct horzcat : public cat<VT>
{
	virtual typename default_tmp_matrix<VT>::type fetch(){
		return fetch2();
	}
};

//matlab repmat
template <typename E>
typename default_tmp_matrix<typename E::value_type>::type
repmat(const matrix_expression<E> &e,size_type m,size_type n){

	typedef idx_range_range<default_tmp_matrix<typename E::value_type>::type,const range<int>,const range<int> > ret_range_type;
	typedef ret_range_type::iterator ret_range_iterator;

	const size_type sz_row_e=e().size(1);
	const size_type sz_col_e=e().size(2);
	
	const size_type sz_row=m*sz_row_e;
	const size_type sz_col=n*sz_col_e;

	default_tmp_matrix<typename E::value_type>::type ret(sz_row,sz_col);

	for(int i=0,cid=0;i<n;i++,cid+=sz_col_e)
	{
		for(int j=0,rid=0;j<m;j++,rid+=sz_row_e)
		{
			ret_range_type sub=ret(colon(rid,rid+sz_row_e-1),
								   colon(cid,cid+sz_col_e-1));
			
			E::const_iterator begin_e=e().begin();
			
			for(ret_range_iterator begin=sub.begin(),end=sub.end();
				begin != end;++begin,++begin_e)
			{*begin=*begin_e;}
		}
	}

	return ret;
}

//matlab reshape
template <typename E>
typename default_tmp_matrix<typename E::value_type>::type 
reshape(const matrix_expression<E> &e, size_type row, size_type col){
	assert(e().size() == row*col);
	
	default_tmp_matrix<E::value_type>::type out(row,col);
	default_tmp_matrix<E::value_type>::type::iterator out_begin=out.begin(),out_end=out.end();
	E::const_iterator in_begin=e().begin(),in_end=e().end();

	while(in_begin != in_end)
	{
		*out_begin=*in_begin;
		++out_begin;++in_begin;
	}

	return out;
}

template <typename E>
typename default_tmp_matrix<typename E::value_type>::type 
reshape(const matrix_expression<E> &e, size_type row){
	return reshape(e(),row,e().size()/row);
}

//matlab rot90
template <typename E>
typename default_tmp_matrix<typename E::value_type>::type 
rot_90(const matrix_expression<E> &e)
{
	const size_type sz_row=e().size(1);
	const size_type sz_col=e().size(2);

	default_tmp_matrix<typename E::value_type>::type ret(sz_col,sz_row);

	for(int i=0;i < sz_col;i++)
		for(int j=0;j < sz_row;j++)
			ret(sz_col-i-1,j)=e()(j,i);
	
	return ret;
}

template <typename E>
typename default_tmp_matrix<typename E::value_type>::type 
rot_180(const matrix_expression<E> &e)
{
	const size_type sz_row=e().size(1);
	const size_type sz_col=e().size(2);

	default_tmp_matrix<typename E::value_type>::type ret(sz_row,sz_col);

	for(int i=0;i < sz_col;i++)
		for(int j=0;j < sz_row;j++)
			ret(sz_row-j-1,sz_col-i-1)=e()(j,i);
	
	return ret;
}

template <typename E>
typename default_tmp_matrix<typename E::value_type>::type
rot_270(const matrix_expression<E> &e)
{
	const size_type sz_row=e().size(1);
	const size_type sz_col=e().size(2);

	default_tmp_matrix<typename E::value_type>::type ret(sz_col,sz_row);

	for(int i=0;i < sz_col;i++)
		for(int j=0;j < sz_row;j++)
			ret(i,sz_row-j-1)=e()(j,i);

	return ret;
}

template <typename E>
typename default_tmp_matrix<typename E::value_type>::type
rot90(const matrix_expression<E> &e,size_type k=1)
{
	k=(k%4+4)%4;

	assert(k);

	switch(k)
	{
	case 1:return rot_90(e);
	case 2:return rot_180(e);
	case 3:return rot_270(e);
	default:assert(false);
	}
}

//matlab vertcat
template <typename VT>
struct vertcat : public cat<VT>
{
	virtual typename default_tmp_matrix<VT>::type fetch(){
		return fetch1();
	}
};
*/

// added by ZeyunSHI
template <typename T, typename E>
  T trace(const matrix_expression<E> &e) {
  const size_t N = std::min(e().size(1), e().size(2));
  T tmp = 0;
  for (size_t i = 0; i < N; ++i)
    tmp += e()(i,i);
  return tmp;
}

template <typename E>
typename E::value_type trace(const matrix_expression<E> &e) {
  return trace<typename E::value_type>(e());
}


#undef max
#undef min
template <typename E>
typename E::value_type max(const matrix_expression<E> &m) {
	return *std::max_element(m().begin(), m().end());
}

template <typename E>
typename E::value_type min(const matrix_expression<E> &m) {
	return *std::min_element(m().begin(), m().end());
}

}	// namespace matrix
}	// namespace zjucad
