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

namespace zjucad { namespace matrix {

class colon_void{};
inline colon_void colon(void) {return colon_void();}

template <typename E>
struct matrix_expression;

template <typename T>
struct range : public matrix_expression<range<T> > {
	typedef range<T> expression_type;
	typedef zjucad::matrix::size_type size_type;
	typedef T value_type;
	typedef const T const_reference;
	typedef T reference;

	range(){}
	range(const T front, const T last):m_front(front), m_end(size_type(last-front)+front+1){}

	size_type size() const {return m_end-m_front;}
	size_type size(int dim) const {return (1==dim)?m_end-m_front:1;}

	struct const_iterator : public default_matrix_iterator<value_type> {
		const_iterator(){}
		const_iterator(const T& value):m_value(value){}
		value_type operator*(void) const {return m_value;}
		const_iterator operator++(void) {++m_value; return *this;}	// need traits for ++
		const_iterator operator+=(size_type n) {m_value+=n; return *this;}
		bool operator != (const const_iterator &i) const {return m_value != i.m_value;}	// when m_value over end?
		bool operator == (const const_iterator &i) const {return !this->operator!=(i);}
		T m_value;
	};
	const_iterator begin(void) const {return const_iterator(m_front);}
	const_iterator end(void) const {return const_iterator(m_end);}

	value_type operator()(size_type i) const {return m_front+i;}
	value_type operator[](size_type i) const {return m_front+i;}
	const expression_type& operator=(const expression_type &r) {
		m_front = r.m_front; m_end= r.m_end;
		return *this;
	}

	T m_front, m_end;
};

template <typename T>
struct interleaved : public matrix_expression<interleaved<T> > {
	typedef interleaved<T> expression_type;
	typedef zjucad::matrix::size_type size_type;
	typedef T value_type;
	typedef const T const_reference;
	typedef T reference;

	interleaved(){}
	interleaved(const T front, const T last, const T step)
		:m_front(front), m_step(step)
	{
		m_nr=(last-front)/step;
		m_nr++;
	}

	size_type size() const {return m_nr;}
	size_type size(int dim) const {return (1==dim)?m_nr:1;}

	struct const_iterator : public default_matrix_iterator<value_type> {
		
		const_iterator(){}
		const_iterator(const T& value,const T& step):m_value(value),m_step(step){}
		
		value_type operator*(void) const {return m_value;}
		
		const_iterator operator++(void) {m_value+=m_step; return *this;}	// need traits for ++
		const_iterator operator+=(size_type n) {m_value+=n*m_step; return *this;}
		
		bool operator != (const const_iterator &i) const {return m_value != i.m_value;}	// when m_value over end?
		bool operator == (const const_iterator &i) const {return !this->operator!=(i);}
		
		T m_value,m_step;
	};

	const_iterator begin(void) const {return const_iterator(m_front,m_step);}
	const_iterator end(void) const {return const_iterator(m_front+m_step*m_nr,m_step);}

	value_type operator()(size_type i) const {return m_front+i*m_step;}
	value_type operator[](size_type i) const {return m_front+i*m_step;}
	const expression_type& operator=(const expression_type &r) {
        m_front = r.m_front; m_nr= r.m_nr; m_step = r.m_step;
		return *this;
	}

	T m_front, m_nr, m_step;
};

template <typename T1, typename T2>
inline range<T1> colon(const T1 &front, const T2 &last) { return range<T1>(front, last); }

template <typename T1, typename T2, typename T3>
inline interleaved<T1> colon(const T1 &front, const T2 &step, const T3& last) { return interleaved<T1>(front, last, step); }

}	// namespace matrix
}	// namespace zjucad

