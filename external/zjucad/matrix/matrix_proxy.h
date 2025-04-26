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

#include "matrix_expression.h"
#include "colon.h"
#include <algorithm>

namespace zjucad { namespace matrix {

#define KEEP_SIZE_RESIZE \
	void resize(size_type size) {assert(size == size());} \
	void resize(size_type row, size_type col) {assert(row == size(1) && col == size(2));}

#define SCALAR_FILL\
	const expression_type& operator=(const_reference a) {\
		std::fill(begin(), end(), a);\
		return *this;\
	}

#define PROXY_ASSIGN \
	template <typename E9> \
	const expression_type& operator=(const matrix_expression<E9> &e) { \
		assign(*this, e());                                        \
		return *this; \
	} \
	const expression_type& operator=(const expression_type &e) { \
		assign(*this, e);                                    \
		return *this; \
	}

#define PROXY_OPS \
	PROXY_ACCESS; \
	SCALAR_FILL; \
	MATRIX_SELF_OP; \
	KEEP_SIZE_RESIZE; \
	PROXY_ASSIGN;

//2D cv_ or _cv
template <typename E>
class idx_cv_i : public matrix_expression<idx_cv_i<E> > {
public:
    typedef idx_cv_i<E> expression_type;
	typedef typename E::size_type size_type;
	typedef typename E::const_reference const_reference;
	typedef typename base_type<E>::reference reference;
	typedef typename base_type<E>::value_type value_type;

	idx_cv_i(){}
	idx_cv_i(E &m, size_type col):m_m(m), m_col(col) {assert(m_col >= 0 && m_col < m_m.size(2));}

	// basic operation
	// matrix operation
	size_type size(int dim) const {return (1==dim)?m_m.size(1):1;}
	const_reference operator()(idx_type i, idx_type j) const {assert(j == 0); return m_m(i, m_col);}
	reference operator()(idx_type i, idx_type j) {assert(j == 0); return m_m(i, m_col);}

	// vector operation
	size_type size(void) const {return m_m.size(1);}

	const_reference operator()(idx_type i) const {return m_m(i, m_col);}
	reference operator()(idx_type i) {return m_m(i, m_col);}
	const_reference operator[](idx_type i) const {return m_m(i, m_col);}
	reference operator[](idx_type i) {return m_m(i, m_col);}

	typedef typename E::const_iterator const_iterator;
	typedef typename base_type<E>::iterator iterator;

	const_iterator begin(void) const {
		assert(m_col>=0 && m_col<m_m.size(2));
		return const_iterator(m_m.begin()+m_col*m_m.size(1));
	}
	const_iterator end(void) const {return const_iterator(m_m.begin()+(m_col+1)*m_m.size(1));}

	iterator begin(void) {
		assert(m_col>=0 && m_col<m_m.size(2));
		return iterator(m_m.begin()+m_col*m_m.size(1));
	}
	iterator end(void) {return iterator(m_m.begin()+(m_col+1)*m_m.size(1));}

	PROXY_OPS;

	E &m_m;
	size_type m_col;
};

template <typename E1, typename E2>
class idx_cv_J : public matrix_expression<idx_cv_J<E1, E2> > {
public:
    typedef idx_cv_J<E1, E2> expression_type;
	typedef typename E1::size_type size_type;
	typedef typename E1::const_reference const_reference;
	typedef typename base_type<E1>::reference reference;
	typedef typename base_type<E1>::value_type value_type;

	idx_cv_J(E1 &m, const E2 &J):m_m(m), m_J(J) {}

	// basic operation
	// matrix operation
	size_type size(int dim) const {return (1==dim)?m_m.size(1):m_J.size();}
	const_reference operator()(idx_type i, idx_type j) const {return m_m(i, m_J[j]);}
	reference operator()(idx_type i, idx_type j) {return m_m(i, m_J[j]);}

	// vector operation
	size_type size(void) const {return m_m.size(1)*m_J.size();}

	const_reference operator()(idx_type i) const {return m_m(i%m_m.size(1), m_J[i/m_m.size(1)]);}
	reference operator()(idx_type i) {return m_m(i%m_m.size(1), m_J[i/m_m.size(1)]);}
	const_reference operator[](idx_type i) const {return m_m(i%m_m.size(1), m_J[i/m_m.size(1)]);}
	reference operator[](idx_type i) {return m_m(i%m_m.size(1), m_J[i/m_m.size(1)]);}

	template <typename M, typename I>
	struct iterator_base : public default_matrix_iterator<value_type> {
		iterator_base(){}
		iterator_base(M *m, size_type row, size_type col):m_m(m), m_row(row), m_col(col){
			if(m_col < m_m->m_J.size())
				m_i = m_m->m_m.begin()+
					m_row+m_m->m_J[m_col]*m_m->m_m.size(1);
		}
		iterator_base operator++(void) {
			if(m_m->m_m.size(1)-1 == m_row) {
				++m_col;
				m_row = 0;
				if(m_col < m_m->m_J.size())
					m_i=m_m->m_m.begin()+m_m->m_m.size(1)*m_m->m_J[m_col];
			}
			else {
				++m_i; ++m_row;
			}
			return *this;
		}
		iterator_base operator+=(offset_type n) {
			const offset_type rn = m_m->m_m.size(1);
			const offset_type n2 = m_row+m_col*rn+n;
			m_row = n2%rn;
			m_col = n2/rn;
			if(m_col < m_m->m_J.size())
				m_i = m_m->m_m.begin()+
					m_row+m_m->m_J[m_col]*m_m->m_m.size(1);
			return *this;
		}
		iterator_base operator+(offset_type n) const {
			iterator_base i = *this;
			i+=n;
			return i;
		}
		const_reference operator*(void) const {return *m_i;}
		template <typename I1>
		bool operator != (const I1 &i) const {assert(m_m == i.m_m);return m_row !=i.m_row || m_col != i.m_col;}
		template <typename I1>
		bool operator == (const I1 &i) const {return !this->operator!=(i);}
		M *m_m;
		I m_i;
		size_type m_row, m_col;
	};
	typedef iterator_base<const expression_type, typename E1::const_iterator> const_iterator;
	const_iterator begin(void) const {return const_iterator(this, 0, 0);}
	const_iterator end(void) const {return const_iterator(this, 0, size(2));}

	struct iterator : public iterator_base<expression_type, typename base_type<E1>::iterator> {
		iterator(){}
		iterator(expression_type *m, size_type row, size_type col)
			:iterator_base<expression_type, typename base_type<E1>::iterator>(m, row, col){}
		operator const_iterator() const {return const_iterator(this->m_m, this->m_row, this->m_col);}
		reference operator*(void) {return *this->m_i;}
		iterator operator+(size_type n) const {
			iterator i = *this;
			i+=n;
			return i;
		}
	};
	iterator begin(void) {return iterator(this, 0, 0);}
	iterator end(void) {return iterator(this, 0, size(2));}

	PROXY_OPS;

	E1 &m_m;
	const E2 &m_J;
};

template <typename E, typename R>
class idx_cv_range : public matrix_expression<idx_cv_range<E, R> > {
public:
    typedef idx_cv_range<E, R> expression_type;
	typedef typename E::size_type size_type;
	typedef typename E::const_reference const_reference;
	typedef typename base_type<E>::reference reference;
	typedef typename base_type<E>::value_type value_type;

	idx_cv_range(E &m, R &range):m_m(m), m_range(range) {}

	// basic operation
	// vector operation
	size_type size(void) const {return m_m.size(1)*m_range.size();}
	const_reference operator()(size_type i) const {
		return m_m(i%m_m.size(1), i/m_m.size(1)+m_range.m_front);
	}
	reference operator()(size_type i) {
		return m_m(i%m_m.size(1), i/m_m.size(1)+m_range.m_front);
	}
	const_reference operator[](size_type i) const {
		return m_m(i%m_m.size(1), i/m_m.size(1)+m_range.m_front);
	}
	reference operator[](size_type i) {
		return m_m(i%m_m.size(1), i/m_m.size(1)+m_range.m_front);
	}

	// 2D matrix operation
	size_type size(int dim) const {return (1==dim)?m_m.size(1):m_range.size();}
	const_reference operator()(size_type i, size_type j) const {
		return m_m(i, m_range[j]);
	}
	reference operator()(size_type i, size_type j) {
		return m_m(i, m_range[j]);
	}

	typedef typename E::const_iterator const_iterator;
	typedef typename base_type<E>::iterator iterator;
	const_iterator begin(void) const {return m_m.begin()+m_range.m_front*m_m.size(1);}
	const_iterator end(void) const {return m_m.begin()+m_range.m_end*m_m.size(1);}

	iterator begin(void) {return m_m.begin()+m_range.m_front*m_m.size(1);}
	iterator end(void) {return m_m.begin()+m_range.m_end*m_m.size(1);}

	PROXY_OPS;

	E &m_m;
	R m_range;
};

template <typename E, typename R>
class idx_cv_iv : public matrix_expression<idx_cv_iv<E, R> > {
public:
    typedef idx_cv_iv<E, R> expression_type;
	typedef typename E::size_type size_type;
	typedef typename E::const_reference const_reference;
	typedef typename base_type<E>::reference reference;
	typedef typename base_type<E>::value_type value_type;

	idx_cv_iv(E &m, R &inter):m_m(m), m_inter(inter) {}

	// basic operation
	// vector operation
	size_type size(void) const {return m_m.size(1)*m_inter.size();}
	const_reference operator()(size_type i) const {
		return m_m(i%m_m.size(1), m_inter[i/m_m.size(1)]);
	}
	reference operator()(size_type i) {
		return m_m(i%m_m.size(1), m_inter[i/m_m.size(1)]);
	}
	const_reference operator[](size_type i) const {
		return m_m(i%m_m.size(1), m_inter[i/m_m.size(1)]);
	}
	reference operator[](size_type i) {
		return m_m(i%m_m.size(1), m_inter[i/m_m.size(1)]);
	}

	// 2D matrix operation
	size_type size(int dim) const {return (1==dim)?m_m.size(1):m_inter.size();}
	const_reference operator()(size_type i, size_type j) const {
		return m_m(i, m_inter[j]);
	}
	reference operator()(size_type i, size_type j) {
		return m_m(i, m_inter[j]);
	}

	template <typename M, typename I>
	struct iterator_base : public default_matrix_iterator<value_type> {
		iterator_base(){}
		iterator_base(M *m, size_type row, size_type col):m_m(m), m_row(row), m_col(col){
			if(m_col < m_m->m_inter.size())
				m_i = m_m->m_m.begin()+
					  m_row+
					  m_m->m_inter[m_col]*m_m->m_m.size(1);
		}
		iterator_base operator++(void) {
			if(m_m->m_m.size(1)-1 == m_row) {
				++m_col;
				m_row = 0;
				if(m_col < m_m->m_inter.size())
					m_i+=(m_m->m_inter.m_step-1)*m_m->size(1)+1;
			}
			else {
				++m_i; ++m_row;
			}
			return *this;
		}
		iterator_base operator+=(offset_type n) {
			const offset_type rn = m_m->size(1);
			const offset_type col_adv = m_row+n;
			m_row = col_adv%rn;
			m_col +=col_adv/rn;
			
			if(m_col < m_m->m_inter.size())
				m_i = m_m->m_m.begin()+
					  m_row+
					  m_m->m_inter[m_col]*m_m->m_m.size(1);
			
			return *this;
		}
		iterator_base operator+(offset_type n) const {
			iterator_base i = *this;
			i+=n;
			return i;
		}

		const_reference operator*(void) const {return *m_i;}
		template <typename I1>bool operator != (const I1 &i) const {assert(m_m == i.m_m);return m_row !=i.m_row || m_col != i.m_col;}
		template <typename I1>bool operator == (const I1 &i) const {return !this->operator!=(i);}
		
		M *m_m;
		I m_i;
		size_type m_row, m_col;
	};
	typedef iterator_base<const expression_type, typename E::const_iterator> const_iterator;
	const_iterator begin(void) const {return const_iterator(this, 0, 0);}
	const_iterator end(void) const {return const_iterator(this, 0, size(2));}

	struct iterator : public iterator_base<expression_type, typename base_type<E>::iterator> {
		iterator(){}
		iterator(expression_type *m, size_type row, size_type col)
			:iterator_base<expression_type, typename base_type<E>::iterator>(m, row, col){}
		operator const_iterator() const {return const_iterator(this->m_m, this->m_row, this->m_col);}
		reference operator*(void) {return *this->m_i;}
		iterator operator+(size_type n) const {
			iterator i = *this;
			i+=n;
			return i;
		}
	};

	iterator begin(void) {return iterator(this, 0,0);}
	iterator end(void) {return iterator(this, 0,size(2));}

	PROXY_OPS;

	E &m_m;
	R m_inter;
};

template <typename E>
class idx_i_cv : public matrix_expression<idx_i_cv<E> > {
public:
    typedef idx_i_cv<E> expression_type;
	typedef typename E::size_type size_type;
	typedef typename E::const_reference const_reference;
	typedef typename base_type<E>::reference reference;
	typedef typename base_type<E>::value_type value_type;

	idx_i_cv(E &m, size_type row):m_m(m), m_row(row) {assert(m_row >= 0 && m_row < m_m.size(1));}

	// basic operation
	// vector operation
	size_type size(void) const {return m_m.size(2);}

	// 2D matrix operation
	size_type size(int dim) const {return (1==dim)?1:m_m.size(2);}
	const_reference operator()(idx_type i, idx_type j) const {assert(i == 0); return m_m(this->m_row, j);}
	reference operator()(idx_type i, idx_type j) {assert(i == 0); return m_m(this->m_row, j);}

	const_reference operator()(idx_type i) const {return m_m(m_row, i);}
	reference operator()(idx_type i) {return m_m(m_row, i);}
	const_reference operator[](idx_type i) const {return m_m(m_row, i);}
	reference operator[](idx_type i) {return m_m(m_row, i);}

	template <typename I>
	struct iterator_base : public default_matrix_iterator<value_type> {
		iterator_base(){}
		iterator_base(I i, size_type rows):m_i(i), m_rows(rows) {}
		const_reference operator *(void) const {return *m_i;}
		iterator_base operator++(void) {m_i+=m_rows; return *this;}
		iterator_base operator+=(size_type n) {m_i+=n*m_rows; return *this;}
		iterator_base operator+(size_type n) const {iterator_base i=*this; i+=n; return i;}
		template <typename I1>
		bool operator != (const I1 &i) const {return m_i != i.m_i;}
		template <typename I1>
		bool operator == (const I1 &i) const {return !this->operator!=(i);}
		I m_i;
		size_type m_rows;
	};
	typedef iterator_base<typename E::const_iterator> const_iterator;
	struct iterator : public iterator_base<typename base_type<E>::iterator> {
		iterator(){}
		iterator(typename E::iterator i, size_type rows)
      :iterator_base<typename base_type<E>::iterator>(i, rows) {}
		operator const_iterator() const {return const_iterator(this->m_i, this->m_rows);}
		iterator operator+(size_type n) const {iterator i=*this; i+=n; return i;}
		reference operator *(void) {return *this->m_i;}
	};

	const_iterator begin(void) const {return const_iterator(m_m.begin()+m_row, m_m.size(1));}
	const_iterator end(void) const {return const_iterator(m_m.begin()+m_row+m_m.size(2)*m_m.size(1), m_m.size(1));}
	iterator begin(void) {return iterator(m_m.begin()+m_row, m_m.size(1));}
	iterator end(void) {return iterator(m_m.begin()+m_row+m_m.size(2)*m_m.size(1), m_m.size(1));}

	PROXY_OPS;

	E &m_m;
	size_type m_row;
};

template <typename E1, typename E2>
class idx_J_cv : public matrix_expression<idx_J_cv<E1, E2> > {
public:
    typedef idx_J_cv<E1, E2> expression_type;
	typedef typename E1::size_type size_type;
	typedef typename E1::const_reference const_reference;
	typedef typename base_type<E1>::reference reference;
	typedef typename base_type<E1>::value_type value_type;

	idx_J_cv(E1 &m, const E2 &J):m_m(m), m_J(J) {}

	// basic operation
	// matrix operation
	size_type size(int dim) const {return (1==dim)?m_J.size():m_m.size(2);}
	const_reference operator()(idx_type i, idx_type j) const {return m_m(m_J[i], j);}
	reference operator()(idx_type i, idx_type j) {return m_m(m_J[i], j);}

	// vector operation
	size_type size(void) const {return m_J.size()*m_m.size(2);}

	const_reference operator()(idx_type i) const {return m_m(m_J[i%m_J.size()], i/m_J.size());}
	reference operator()(idx_type i) {return m_m(m_J[i%m_J.size()], i/m_J.size());}
	const_reference operator[](idx_type i) const {return m_m(m_J[i%m_J.size()], i/m_J.size());}
	reference operator[](idx_type i) {return m_m(m_J[i%m_J.size()], i/m_J.size());}

	template <typename M, typename I>
	struct iterator_base : public default_matrix_iterator<value_type> {
		iterator_base(){}
		iterator_base(M *m, size_type row, size_type col):m_m(m), m_row(row), m_col(col){
			if(m_col < m_m->size(2))
				m_i = m_m->m_m.begin()+
					  m_m->m_J[m_row]+m_col*m_m->m_m.size(1);
		}
		iterator_base operator++(void) {
			if(m_m->m_J.size()-1 == m_row) {
				++m_col;
				m_row = 0;
				if(m_col < m_m->size(2))
					m_i = m_m->m_m.begin()+
						  m_m->m_J[m_row]+m_col*m_m->m_m.size(1);
			}
			else {
				m_i+=m_m->m_J[m_row+1]-m_m->m_J[m_row]; 
				++m_row;
			}
			return *this;
		}
		iterator_base operator+=(offset_type n) {
			const offset_type rn = m_m->size(1);
			const offset_type n2 = m_row+m_col*rn+n;
			m_row = n2%rn;
			m_col = n2/rn;

			if(m_col < m_m->size(2))
				m_i = m_m->m_m.begin()+
					  m_m->m_J[m_row]+m_col*m_m->m_m.size(1);
			
			return *this;
		}
		iterator_base operator+(offset_type n) const {
			iterator_base i = *this;
			i+=n;
			return i;
		}
		const_reference operator*(void) const {return *m_i;}
		template <typename I1>
		bool operator != (const I1 &i) const {assert(m_m == i.m_m);return m_row !=i.m_row || m_col != i.m_col;}
		template <typename I1>
		bool operator == (const I1 &i) const {return !this->operator!=(i);}
		M *m_m;
		I m_i;
		size_type m_row, m_col;
	};
	typedef iterator_base<const expression_type, typename E1::const_iterator> const_iterator;
	const_iterator begin(void) const {return const_iterator(this, 0, 0);}
	const_iterator end(void) const {return const_iterator(this, 0, size(2));}

	struct iterator : public iterator_base<expression_type, typename base_type<E1>::iterator> {
		iterator(){}
		iterator(expression_type *m, size_type row, size_type col)
			:iterator_base<expression_type, typename base_type<E1>::iterator>(m, row, col){}
		operator const_iterator() const {return const_iterator(this->m_m, this->m_row, this->m_col);}
		reference operator*(void) {return *this->m_i;}
		iterator operator+(size_type n) const {
			iterator i = *this;
			i+=n;
			return i;
		}
	};
	iterator begin(void) {return iterator(this, 0, 0);}
	iterator end(void) {return iterator(this, 0, size(2));}

	PROXY_OPS;

	E1 &m_m;
	const E2 &m_J;
};

template <typename E, typename R>
class idx_range_cv : public matrix_expression<idx_range_cv<E, R> > {
public:
    typedef idx_range_cv<E, R> expression_type;
	typedef typename E::size_type size_type;
	typedef typename E::const_reference const_reference;
	typedef typename base_type<E>::reference reference;
	typedef typename base_type<E>::value_type value_type;

	idx_range_cv(E &m, R &range):m_m(m), m_range(range) {}

	// basic operation
	// vector operation
	size_type size(void) const {return m_range.size()*m_m.size(2);}
	const_reference operator()(size_type i) const {return m_m(i%m_range.size()+m_range.m_front, i/m_range.size());}
	reference operator()(size_type i) {return m_m(i%m_range.size()+m_range.m_front, i/m_range.size());}
	const_reference operator[](size_type i) const {return m_m(i%m_range.size()+m_range.m_front, i/m_range.size());}
	reference operator[](size_type i) {return m_m(i%m_range.size()+m_range.m_front, i/m_range.size());}

	// 2D matrix operation
	size_type size(int dim) const {return (1==dim)?m_range.size():m_m.size(2);}
	const_reference operator()(size_type i, size_type j) const {return m_m(m_range[i] ,j);}
	reference operator()(size_type i, size_type j) {return m_m(m_range[i] ,j);}

	template <typename M,typename I>
	struct iterator_base : public default_matrix_iterator<value_type> {
		
		iterator_base(){}
		iterator_base(M *m, size_type row, size_type col)
			:m_m(m),m_row(row),m_col(col)
		{
			if(m_col < m_m->size(2)){
				m_i=m_m->m_m.begin()+
					(m_row+m_m->m_range.m_front)+
					(m_col*m_m->m_m.size(2));
			}
		}
		
		const_reference operator *(void) const {return *m_i;}

		iterator_base operator++(void){
			if(m_row == m_m->size(1)-1){
				m_row=0;
				m_col++;
				if(m_col < m_m->size(2))
					m_i+=m_m->m_m.size(1)-
						 (m_m->m_range.m_end-m_m->m_range.m_front-1);
			}else{
				++m_row;
				++m_i;
			}
			return *this;
		}
		
		iterator_base operator+=(size_type n){
			offset_type new_row = n+m_row;
			const offset_type fn = m_m->size(1);
			const offset_type col_adv = new_row/fn;

			new_row %= fn;
			
			if(m_col+col_adv < m_m->size(2))
				m_i+=col_adv*m_m->m_m.size(1)+
					 (new_row-m_row);

			m_row = new_row;
			m_col += col_adv;
			
			return *this;
		}
		
		iterator_base operator+(size_type n) const{
			iterator_base i=*this; 
			i+=n; 
			return i;
		}
		
		template <typename I1>bool operator != (const I1 &i) const {assert(m_m == i.m_m);return m_row != i.m_row || m_col != i.m_col;}
		template <typename I1>bool operator == (const I1 &i) const {return !this->operator!=(i);}
		
		M *m_m;
		I m_i;
		size_type m_row,m_col;
	};
	typedef iterator_base<const expression_type, typename E::const_iterator> const_iterator;
	
	struct iterator : public iterator_base<expression_type, typename base_type<E>::iterator> {
		iterator(){}
		iterator(expression_type *m, size_type row, size_type col)
			:iterator_base<expression_type, typename base_type<E>::iterator>(m,row,col) {}
		operator const_iterator() const {return const_iterator(this->m_m, this->m_row, this->m_col);}
		iterator operator+(size_type n) const {iterator i=*this; i+=n; return i;}
		reference operator *(void) {return *this->m_i;}
	};

	const_iterator begin(void) const {return const_iterator(this, 0,0);}
	const_iterator end(void) const {return const_iterator(this, 0,m_m.size(2));}
	iterator begin(void) {return iterator(this, 0,0);}
	iterator end(void) {return iterator(this, 0,m_m.size(2));}

	PROXY_OPS;

	E &m_m;
	R m_range;
};

template <typename E, typename R>
class idx_iv_cv : public matrix_expression<idx_iv_cv<E, R> > {
public:
    typedef idx_iv_cv<E, R> expression_type;
	typedef typename E::size_type size_type;
	typedef typename E::const_reference const_reference;
	typedef typename base_type<E>::reference reference;
	typedef typename base_type<E>::value_type value_type;

	idx_iv_cv(E &m, R &inter):m_m(m), m_inter(inter) {}

	// basic operation
	// vector operation
	size_type size(void) const {return m_inter.size()*m_m.size(2);}
	const_reference operator()(size_type i) const {return m_m(m_inter[i%m_inter.size()], i/m_inter.size());}
	reference operator()(size_type i) {return m_m(m_inter[i%m_inter.size()], i/m_inter.size());}
	const_reference operator[](size_type i) const {return m_m(m_inter[i%m_inter.size()], i/m_inter.size());}
	reference operator[](size_type i) {return m_m(m_inter[i%m_inter.size()], i/m_inter.size());}

	// 2D matrix operation
	size_type size(int dim) const {return (1==dim)?m_inter.size():m_m.size(2);}
	const_reference operator()(size_type i, size_type j) const {
		return m_m(m_inter[i], j);
	}
	reference operator()(size_type i, size_type j) {
		return m_m(m_inter[i], j);
	}

	template <typename M, typename I>
	struct iterator_base : public default_matrix_iterator<value_type> {
		iterator_base(){}
		iterator_base(M *m, size_type row, size_type col):m_m(m), m_row(row), m_col(col){
			if(m_col < m_m->size(2))
				m_i = m_m->m_m.begin()+
					  m_m->m_inter[m_row]+
					  m_col*m_m->m_m.size(1);
		}
		iterator_base operator++(void) {
			if(m_m->m_inter.size()-1 == m_row) {
				++m_col;
				m_row = 0;
				if(m_col < m_m->size(2))
					m_i = m_m->m_m.begin()+
						  m_m->m_inter[m_row]+
						  m_col*m_m->m_m.size(1);
			}
			else {
				m_i+=m_m->m_inter.m_step; 
				++m_row;
			}
			return *this;
		}
		iterator_base operator+=(offset_type n) {
			const offset_type rn = m_m->size(1);
			const offset_type col_adv = m_row+n;
			m_row = col_adv%rn;
			m_col +=col_adv/rn;
			
			if(m_col < m_m->size(2))
				m_i = m_m->m_m.begin()+
					  m_m->m_inter[m_row]+
					  m_col*m_m->m_m.size(1);
			
			return *this;
		}
		iterator_base operator+(offset_type n) const {
			iterator_base i = *this;
			i+=n;
			return i;
		}

		const_reference operator*(void) const {return *m_i;}
		template <typename I1>bool operator != (const I1 &i) const {assert(m_m == i.m_m);return m_row !=i.m_row || m_col != i.m_col;}
		template <typename I1>bool operator == (const I1 &i) const {return !this->operator!=(i);}
		
		M *m_m;
		I m_i;
		size_type m_row, m_col;
	};
	typedef iterator_base<const expression_type, typename E::const_iterator> const_iterator;
	const_iterator begin(void) const {return const_iterator(this, 0, 0);}
	const_iterator end(void) const {return const_iterator(this, 0, size(2));}

	struct iterator : public iterator_base<expression_type, typename base_type<E>::iterator> {
		iterator(){}
		iterator(expression_type *m, size_type row, size_type col)
			:iterator_base<expression_type, typename base_type<E>::iterator>(m, row, col){}
		operator const_iterator() const {return const_iterator(this->m_m, this->m_row, this->m_col);}
		reference operator*(void) {return *this->m_i;}
		iterator operator+(size_type n) const {
			iterator i = *this;
			i+=n;
			return i;
		}
	};

	iterator begin(void) {return iterator(this, 0,0);}
	iterator end(void) {return iterator(this, 0,size(2));}

	PROXY_OPS;

	E &m_m;
	R m_inter;
};

//2D different mode
template <typename I, typename Arg1, typename Arg2>
class vp2
{
public:
	vp2(Arg1 &arg1,Arg2 &arg2)
	:inner(arg1,arg2){}
	I inner;
};

template <typename I, typename Arg1, typename Arg2, typename Arg3>
class vp3
{
public:
	vp3(Arg1 &arg1,Arg2 &arg2,Arg3 &arg3)
	:inner(arg1,arg2,arg3){}
	I inner;
};

//template <typename E>
//class idx_i_i : public matrix_expression<idx_i_i<E> >, public vp2<idx_cv_i<E>, E, size_type > {
//public:
//    typedef idx_i_i<E> expression_type;
//	typedef typename E::size_type size_type;
//	typedef typename E::const_reference const_reference;
//	typedef typename base_type<E>::reference reference;
//	typedef typename base_type<E>::value_type value_type;
//
//	typedef idx_cv_i<E> inner_proxy_type;
//	typedef idx_i_cv<inner_proxy_type> proxy_type;
//	typedef typename proxy_type::iterator iterator;
//	typedef typename proxy_type::const_iterator const_iterator;
//
//	idx_i_i(E &m, size_type row, size_type col)
//	:vp2(m,col),proxy(inner,row){}
//
//	// basic operation
//	// vector operation
//	size_type size(void) const {return proxy.size();}
//	const_reference operator()(size_type i) const {return proxy(i);}
//	reference operator()(size_type i) {return proxy(i);}
//	const_reference operator[](size_type i) const {return proxy[i];}
//	reference operator[](size_type i) {return proxy[i];}
//
//	// 2D matrix operation
//	size_type size(int dim) const {return proxy.size(dim);}
//	const_reference operator()(size_type i, size_type j) const {return proxy(i,j);}
//	reference operator()(size_type i, size_type j) {return proxy(i,j);}
//
//	const_iterator begin(void) const {return proxy.begin();}
//	const_iterator end(void) const {return proxy.end();}
//	iterator begin(void) {return proxy.begin();}
//	iterator end(void) {return proxy.end();}
//
//	PROXY_OPS;
//
//	proxy_type proxy;
//};

template <typename E1, typename E2>
class idx_J_i : public vp2<idx_cv_i<E1>, E1, size_type >, public matrix_expression<idx_J_i<E1, E2> > {
public:
    typedef idx_J_i<E1, E2> expression_type;
	typedef typename E1::size_type size_type;
	typedef typename E1::const_reference const_reference;
	typedef typename base_type<E1>::reference reference;
	typedef typename base_type<E1>::value_type value_type;
    typedef vp2<idx_cv_i<E1>, E1, size_type > base_class;

	typedef idx_cv_i<E1> inner_proxy_type;
	typedef idx_J_cv<inner_proxy_type, E2> proxy_type;
	typedef typename proxy_type::iterator iterator;
	typedef typename proxy_type::const_iterator const_iterator;

	idx_J_i(E1 &m, const E2 &row, size_type col)
    :base_class(m,col),proxy(base_class::inner,row){}

	// basic operation
	// vector operation
	size_type size(void) const {return proxy.size();}
	const_reference operator()(size_type i) const {return proxy(i);}
	reference operator()(size_type i) {return proxy(i);}
	const_reference operator[](size_type i) const {return proxy[i];}
	reference operator[](size_type i) {return proxy[i];}

	// 2D matrix operation
	size_type size(int dim) const {return proxy.size(dim);}
	const_reference operator()(size_type i, size_type j) const {return proxy(i,j);}
	reference operator()(size_type i, size_type j) {return proxy(i,j);}

	const_iterator begin(void) const {return proxy.begin();}
	const_iterator end(void) const {return proxy.end();}
	iterator begin(void) {return proxy.begin();}
	iterator end(void) {return proxy.end();}

	PROXY_OPS;

	proxy_type proxy;
};

template <typename E, typename R>
class idx_range_i : public matrix_expression<idx_range_i<E, R> > {
public:
    typedef idx_range_i<E, R> expression_type;
	typedef typename E::size_type size_type;
	typedef typename E::const_reference const_reference;
	typedef typename base_type<E>::reference reference;
	typedef typename base_type<E>::value_type value_type;

	idx_range_i(E &m, R &range_row, size_type col):m_m(m), m_range(range_row), m_col(col) {}

	// basic operation
	// vector operation
	size_type size(void) const {return m_range.size();}
	const_reference operator()(size_type i) const {
		return m_m(m_range(i), m_col);
	}
	reference operator()(size_type i) {
		return m_m(m_range(i), m_col);
	}
	const_reference operator[](size_type i) const {
		return m_m(m_range(i), m_col);
	}
	reference operator[](size_type i) {
		return m_m(m_range(i), m_col);
	}

	// 2D matrix operation
	size_type size(int dim) const {return (1==dim)?m_range.size():1;}
	const_reference operator()(size_type i, size_type j) const {
		assert(j == 0);
		return m_m(m_range[i], m_col);
	}
	reference operator()(size_type i, size_type j) {
		assert(j == 0);
		return m_m(m_range[i], m_col);
	}

	typedef typename E::const_iterator const_iterator;
	typedef typename base_type<E>::iterator iterator;
	const_iterator begin(void) const {return m_m.begin()+m_range.m_front+m_col*m_m.size(1);}
	const_iterator end(void) const {return m_m.begin()+m_range.m_end+m_col*m_m.size(1);}

	iterator begin(void) {return m_m.begin()+m_range.m_front+m_col*m_m.size(1);}
	iterator end(void) {return m_m.begin()+m_range.m_end+m_col*m_m.size(1);}

	PROXY_OPS;

	E &m_m;
	R m_range;
	size_type m_col;
};

template <typename E, typename R>
class idx_iv_i : public matrix_expression<idx_iv_i<E, R> >, public vp2<idx_cv_i<E>, E, size_type > {
public:
    typedef idx_iv_i<E, R> expression_type;
	typedef typename E::size_type size_type;
	typedef typename E::const_reference const_reference;
	typedef typename base_type<E>::reference reference;
	typedef typename base_type<E>::value_type value_type;
    typedef vp2<idx_cv_i<E>, E, size_type > base_class;

	typedef idx_cv_i<E> inner_proxy_type;
	typedef idx_iv_cv<inner_proxy_type, R> proxy_type;
	typedef typename proxy_type::iterator iterator;
	typedef typename proxy_type::const_iterator const_iterator;

	idx_iv_i(E &m, R &inter_row, size_type col)
    :base_class(m,col),proxy(base_class::inner,inter_row){}

	// basic operation
	// vector operation
	size_type size(void) const {return proxy.size();}
	const_reference operator()(size_type i) const {return proxy(i);}
	reference operator()(size_type i) {return proxy(i);}
	const_reference operator[](size_type i) const {return proxy[i];}
	reference operator[](size_type i) {return proxy[i];}

	// 2D matrix operation
	size_type size(int dim) const {return proxy.size(dim);}
	const_reference operator()(size_type i, size_type j) const {return proxy(i,j);}
	reference operator()(size_type i, size_type j) {return proxy(i,j);}

	const_iterator begin(void) const {return proxy.begin();}
	const_iterator end(void) const {return proxy.end();}
	iterator begin(void) {return proxy.begin();}
	iterator end(void) {return proxy.end();}

	PROXY_OPS;

	proxy_type proxy;
};

template <typename E1, typename E2>
class idx_i_J : public matrix_expression<idx_i_J<E1, E2> >, public vp2<idx_cv_J<E1, E2>, E1, const E2 > {
public:
    typedef idx_i_J<E1, E2> expression_type;
	typedef typename E1::size_type size_type;
	typedef typename E1::const_reference const_reference;
	typedef typename base_type<E1>::reference reference;
	typedef typename base_type<E1>::value_type value_type;
    typedef vp2<idx_cv_J<E1, E2>, E1, const E2 > base_class;

	typedef idx_cv_J<E1, E2> inner_proxy_type;
	typedef idx_i_cv<inner_proxy_type> proxy_type;
	typedef typename proxy_type::iterator iterator;
	typedef typename proxy_type::const_iterator const_iterator;

	idx_i_J(E1 &m, size_type row, const E2 &col)
    :base_class(m,col),proxy(base_class::inner,row){}

	// basic operation
	// vector operation
	size_type size(void) const {return proxy.size();}
	const_reference operator()(size_type i) const {return proxy(i);}
	reference operator()(size_type i) {return proxy(i);}
	const_reference operator[](size_type i) const {return proxy[i];}
	reference operator[](size_type i) {return proxy[i];}

	// 2D matrix operation
	size_type size(int dim) const {return proxy.size(dim);}
	const_reference operator()(size_type i, size_type j) const {return proxy(i,j);}
	reference operator()(size_type i, size_type j) {return proxy(i,j);}

	const_iterator begin(void) const {return proxy.begin();}
	const_iterator end(void) const {return proxy.end();}
	iterator begin(void) {return proxy.begin();}
	iterator end(void) {return proxy.end();}

	PROXY_OPS;

	proxy_type proxy;
};

template <typename E1, typename E2, typename E3>
class idx_J_J : public matrix_expression<idx_J_J<E1, E2, E3> >, public vp2<idx_cv_J<E1, E3>, E1, const E3 > {
public:
    typedef idx_J_J<E1, E2, E3> expression_type;
	typedef typename E1::size_type size_type;
	typedef typename E1::const_reference const_reference;
	typedef typename base_type<E1>::reference reference;
	typedef typename base_type<E1>::value_type value_type;
    typedef vp2<idx_cv_J<E1, E3>, E1, const E3 > base_class;

	typedef idx_cv_J<E1, E3> inner_proxy_type;
	typedef idx_J_cv<inner_proxy_type, E2> proxy_type;
	typedef typename proxy_type::iterator iterator;
	typedef typename proxy_type::const_iterator const_iterator;

	idx_J_J(E1 &m, const E2 &row, const E3 &col)
    :base_class(m,col),proxy(base_class::inner,row){}

	// basic operation
	// vector operation
	size_type size(void) const {return proxy.size();}
	const_reference operator()(size_type i) const {return proxy(i);}
	reference operator()(size_type i) {return proxy(i);}
	const_reference operator[](size_type i) const {return proxy[i];}
	reference operator[](size_type i) {return proxy[i];}

	// 2D matrix operation
	size_type size(int dim) const {return proxy.size(dim);}
	const_reference operator()(size_type i, size_type j) const {return proxy(i,j);}
	reference operator()(size_type i, size_type j) {return proxy(i,j);}

	const_iterator begin(void) const {return proxy.begin();}
	const_iterator end(void) const {return proxy.end();}
	iterator begin(void) {return proxy.begin();}
	iterator end(void) {return proxy.end();}

	PROXY_OPS;

	proxy_type proxy;
};

template <typename E1, typename R, typename E2>
class idx_range_J : public matrix_expression<idx_range_J<E1, R, E2> >, public vp2<idx_cv_J<E1, E2>, E1, const E2 > {
public:
    typedef idx_range_J<E1, R, E2> expression_type;
	typedef typename E1::size_type size_type;
	typedef typename E1::const_reference const_reference;
	typedef typename base_type<E1>::reference reference;
	typedef typename base_type<E1>::value_type value_type;
    typedef vp2<idx_cv_J<E1, E2>, E1, const E2 > base_class;

	typedef idx_cv_J<E1, E2> inner_proxy_type;
	typedef idx_range_cv<inner_proxy_type, R> proxy_type;
	typedef typename proxy_type::iterator iterator;
	typedef typename proxy_type::const_iterator const_iterator;

	idx_range_J(E1 &m, R &row, const E2 &col)
    :base_class(m,col),proxy(base_class::inner,row){}

	// basic operation
	// vector operation
	size_type size(void) const {return proxy.size();}
	const_reference operator()(size_type i) const {return proxy(i);}
	reference operator()(size_type i) {return proxy(i);}
	const_reference operator[](size_type i) const {return proxy[i];}
	reference operator[](size_type i) {return proxy[i];}

	// 2D matrix operation
	size_type size(int dim) const {return proxy.size(dim);}
	const_reference operator()(size_type i, size_type j) const {return proxy(i,j);}
	reference operator()(size_type i, size_type j) {return proxy(i,j);}

	const_iterator begin(void) const {return proxy.begin();}
	const_iterator end(void) const {return proxy.end();}
	iterator begin(void) {return proxy.begin();}
	iterator end(void) {return proxy.end();}

	PROXY_OPS;

	proxy_type proxy;
};

template <typename E1, typename R, typename E2>
class idx_iv_J : public matrix_expression<idx_iv_J<E1, R, E2> >, public vp2<idx_cv_J<E1, E2>, E1, const E2 > {
public:
    typedef idx_iv_J<E1, R, E2> expression_type;
	typedef typename E1::size_type size_type;
	typedef typename E1::const_reference const_reference;
	typedef typename base_type<E1>::reference reference;
	typedef typename base_type<E1>::value_type value_type;
    typedef vp2<idx_cv_J<E1, E2>, E1, const E2 > base_class;

	typedef idx_cv_J<E1, E2> inner_proxy_type;
	typedef idx_iv_cv<inner_proxy_type, R> proxy_type;
	typedef typename proxy_type::iterator iterator;
	typedef typename proxy_type::const_iterator const_iterator;

	idx_iv_J(E1 &m, R &row, const E2 &col)
    :base_class(m,col),proxy(base_class::inner,row){}

	// basic operation
	// vector operation
	size_type size(void) const {return proxy.size();}
	const_reference operator()(size_type i) const {return proxy(i);}
	reference operator()(size_type i) {return proxy(i);}
	const_reference operator[](size_type i) const {return proxy[i];}
	reference operator[](size_type i) {return proxy[i];}

	// 2D matrix operation
	size_type size(int dim) const {return proxy.size(dim);}
	const_reference operator()(size_type i, size_type j) const {return proxy(i,j);}
	reference operator()(size_type i, size_type j) {return proxy(i,j);}

	const_iterator begin(void) const {return proxy.begin();}
	const_iterator end(void) const {return proxy.end();}
	iterator begin(void) {return proxy.begin();}
	iterator end(void) {return proxy.end();}

	PROXY_OPS;

	proxy_type proxy;
};

template <typename E, typename R>
class idx_i_range : public matrix_expression<idx_i_range<E, R> >, public vp2<idx_cv_range<E, R>, E, R > {
public:
    typedef idx_i_range<E, R> expression_type;
	typedef typename E::size_type size_type;
	typedef typename E::const_reference const_reference;
	typedef typename base_type<E>::reference reference;
	typedef typename base_type<E>::value_type value_type;
    typedef vp2<idx_cv_range<E, R>, E, R >  base_class;

	typedef idx_cv_range<E, R> inner_proxy_type;
	typedef idx_i_cv<inner_proxy_type> proxy_type;
	typedef typename proxy_type::iterator iterator;
	typedef typename proxy_type::const_iterator const_iterator;

	idx_i_range(E &m, size_type i, R &range_col)
    :base_class(m,range_col),proxy(base_class::inner,i){}

	// basic operation
	// vector operation
	size_type size(void) const {return proxy.size();}
	const_reference operator()(size_type i) const {return proxy(i);}
	reference operator()(size_type i) {return proxy(i);}
	const_reference operator[](size_type i) const {return proxy[i];}
	reference operator[](size_type i) {return proxy[i];}

	// 2D matrix operation
	size_type size(int dim) const {return proxy.size(dim);}
	const_reference operator()(size_type i, size_type j) const {return proxy(i,j);}
	reference operator()(size_type i, size_type j) {return proxy(i,j);}

	const_iterator begin(void) const {return proxy.begin();}
	const_iterator end(void) const {return proxy.end();}
	iterator begin(void) {return proxy.begin();}
	iterator end(void) {return proxy.end();}

	PROXY_OPS;

	proxy_type proxy;
};

template <typename E1, typename E2, typename R>
class idx_J_range : public matrix_expression<idx_J_range<E1, E2, R> >, public vp2<idx_cv_range<E1, R>, E1, R > {
public:
    typedef idx_J_range<E1, E2, R> expression_type;
	typedef typename E1::size_type size_type;
	typedef typename E1::const_reference const_reference;
	typedef typename base_type<E1>::reference reference;
	typedef typename base_type<E1>::value_type value_type;
    typedef vp2<idx_cv_range<E1, R>, E1, R > base_class;

	typedef idx_cv_range<E1, R> inner_proxy_type;
	typedef idx_J_cv<inner_proxy_type, E2> proxy_type;
	typedef typename proxy_type::iterator iterator;
	typedef typename proxy_type::const_iterator const_iterator;

	idx_J_range(E1 &m, const E2 &row, R &range_col)
    :base_class(m,range_col),proxy(base_class::inner,row){}

	// basic operation
	// vector operation
	size_type size(void) const {return proxy.size();}
	const_reference operator()(size_type i) const {return proxy(i);}
	reference operator()(size_type i) {return proxy(i);}
	const_reference operator[](size_type i) const {return proxy[i];}
	reference operator[](size_type i) {return proxy[i];}

	// 2D matrix operation
	size_type size(int dim) const {return proxy.size(dim);}
	const_reference operator()(size_type i, size_type j) const {return proxy(i,j);}
	reference operator()(size_type i, size_type j) {return proxy(i,j);}

	const_iterator begin(void) const {return proxy.begin();}
	const_iterator end(void) const {return proxy.end();}
	iterator begin(void) {return proxy.begin();}
	iterator end(void) {return proxy.end();}

	PROXY_OPS;

	proxy_type proxy;
};

template <typename E, typename R1, typename R2>
class idx_range_range : public matrix_expression<idx_range_range<E, R1, R2> > {
public:
    typedef idx_range_range<E, R1, R2> expression_type;
	typedef typename E::size_type size_type;
	typedef typename E::const_reference const_reference;
	typedef typename base_type<E>::reference reference;
	typedef typename base_type<E>::value_type value_type;

	idx_range_range(E &m, R1 &range1, R2 &range2):m_m(m), m_range1(range1), m_range2(range2) {}

	// basic operation
	// vector operation
	size_type size(void) const {return m_range1.size()*m_range2.size();}
	const_reference operator()(size_type i) const {
		return (*this)(i%m_range1.size(), i/m_range1.size());
	}
	reference operator()(size_type i) {
		return (*this)(i%m_range1.size(), i/m_range1.size());
	}
	const_reference operator[](size_type i) const {
		return (*this)(i%m_range1.size(), i/m_range1.size());
	}
	reference operator[](size_type i) {
		return (*this)(i%m_range1.size(), i/m_range1.size());
	}

	// 2D matrix operation
	size_type size(int dim) const {return (1==dim)?m_range1.size():m_range2.size();}
	const_reference operator()(size_type i, size_type j) const {
		return m_m(m_range1[i], m_range2[j]);
	}
	reference operator()(size_type i, size_type j) {
		return m_m(m_range1[i], m_range2[j]);
	}

	template <typename M, typename I>
	struct iterator_base : public default_matrix_iterator<value_type> {
		iterator_base(){}
		iterator_base(M *m, size_type row, size_type col):m_m(m), m_row(row){
			m_i = m_m->m_m.begin()+
				(m_m->m_range1.m_front+m_row)+
				(m_m->m_range2.m_front+col)*m_m->m_m.size(1);
		}
//		template <typename I1>
//		iterator_base(M *m, size_type row, I1 i):m_m(m), m_row(row), m_i(i){}
		iterator_base operator++(void) {
			if(m_m->m_range1.size()-1 == m_row) {
				m_i+=(m_m->m_m.size(1)-m_m->size(1)+1); m_row = 0;
			}
			else {
				++m_i; ++m_row;
			}
			return *this;
		}
		iterator_base operator+=(offset_type n) {
			if(!n) return *this;
			if(n == m_m->size(1)) {
				m_i+=m_m->m_m.size(1);
				return *this;
			}
			offset_type row_adv = m_row+n;
			offset_type col_adv = row_adv/m_m->size(1);
			row_adv = (row_adv%m_m->size(1))-m_row;
			m_i+=(col_adv*m_m->m_m.size(1)+row_adv);
			m_row += row_adv;
			return *this;
		}
		iterator_base operator+(offset_type n) const {
			iterator_base i = *this;
			i+=n;
			return i;
		}
		const_reference operator*(void) const {return *m_i;}
		template <typename I1>
		bool operator != (const I1 &i) const {assert(m_m == i.m_m);return m_i !=i.m_i;}
		template <typename I1>
		bool operator == (const I1 &i) const {return !this->operator!=(i);}
		M *m_m;
		I m_i;
		size_type m_row;
	};

	typedef iterator_base<const expression_type, typename E::const_iterator> const_iterator;
	const_iterator begin(void) const {return const_iterator(this, 0, 0);}
	const_iterator end(void) const {return const_iterator(this, 0, size(2));}

	struct iterator : public iterator_base<expression_type, typename base_type<E>::iterator> {
		iterator(){}
		iterator(expression_type *m, size_type row, size_type col)
			:iterator_base<expression_type, typename base_type<E>::iterator>(m, row, col){}
		operator const_iterator() const {return const_iterator(this->m_m, this->m_row, this->m_i);}
		reference operator*(void) {return *this->m_i;}
		iterator operator+(size_type n) const {
			iterator i = *this;
			i+=n;
			return i;
		}
	};

	iterator begin(void) {return iterator(this, 0, 0);}
	iterator end(void) {return iterator(this, 0, size(2));}

	PROXY_OPS;

	E &m_m;
	R1 m_range1;
	R2 m_range2;
};

template <typename E, typename R1, typename R2>
class idx_iv_range : public matrix_expression<idx_iv_range<E, R1, R2> >, public vp2<idx_cv_range<E, R2>, E, R2 > {
public:
    typedef idx_iv_range<E, R1, R2> expression_type;
	typedef typename E::size_type size_type;
	typedef typename E::const_reference const_reference;
	typedef typename base_type<E>::reference reference;
	typedef typename base_type<E>::value_type value_type;
    typedef vp2<idx_cv_range<E, R2>, E, R2 > base_class;

	typedef idx_cv_range<E, R2> inner_proxy_type;
	typedef idx_iv_cv<inner_proxy_type, R1> proxy_type;
	typedef typename proxy_type::iterator iterator;
	typedef typename proxy_type::const_iterator const_iterator;

	idx_iv_range(E &m, R1 &interleaved_row, R2 &range_col)
    :base_class(m,range_col),proxy(base_class::inner,interleaved_row){}

	// basic operation
	// vector operation
	size_type size(void) const {return proxy.size();}
	const_reference operator()(size_type i) const {return proxy(i);}
	reference operator()(size_type i) {return proxy(i);}
	const_reference operator[](size_type i) const {return proxy[i];}
	reference operator[](size_type i) {return proxy[i];}

	// 2D matrix operation
	size_type size(int dim) const {return proxy.size(dim);}
	const_reference operator()(size_type i, size_type j) const {return proxy(i,j);}
	reference operator()(size_type i, size_type j) {return proxy(i,j);}

	const_iterator begin(void) const {return proxy.begin();}
	const_iterator end(void) const {return proxy.end();}
	iterator begin(void) {return proxy.begin();}
	iterator end(void) {return proxy.end();}

	PROXY_OPS;

	proxy_type proxy;
};

template <typename E, typename R>
class idx_i_iv : public matrix_expression<idx_i_iv<E, R> >, public vp2<idx_cv_iv<E, R>, E, R > {
public:
    typedef idx_i_iv<E, R> expression_type;
	typedef typename E::size_type size_type;
	typedef typename E::const_reference const_reference;
	typedef typename base_type<E>::reference reference;
	typedef typename base_type<E>::value_type value_type;
    typedef vp2<idx_cv_iv<E, R>, E, R > base_class;

	typedef idx_cv_iv<E, R> inner_proxy_type;
	typedef idx_i_cv<inner_proxy_type> proxy_type;
	typedef typename proxy_type::iterator iterator;
	typedef typename proxy_type::const_iterator const_iterator;

	idx_i_iv(E &m, size_type row, R &interleaved_col)
    :base_class(m,interleaved_col),proxy(base_class::inner,row){}

	// basic operation
	// vector operation
	size_type size(void) const {return proxy.size();}
	const_reference operator()(size_type i) const {return proxy(i);}
	reference operator()(size_type i) {return proxy(i);}
	const_reference operator[](size_type i) const {return proxy[i];}
	reference operator[](size_type i) {return proxy[i];}

	// 2D matrix operation
	size_type size(int dim) const {return proxy.size(dim);}
	const_reference operator()(size_type i, size_type j) const {return proxy(i,j);}
	reference operator()(size_type i, size_type j) {return proxy(i,j);}

	const_iterator begin(void) const {return proxy.begin();}
	const_iterator end(void) const {return proxy.end();}
	iterator begin(void) {return proxy.begin();}
	iterator end(void) {return proxy.end();}

	PROXY_OPS;

	proxy_type proxy;
};

template <typename E1, typename E2, typename R>
class idx_J_iv : public matrix_expression<idx_J_iv<E1, E2, R> >, public vp2<idx_cv_iv<E1, R>, E1, R > {
public:
    typedef idx_J_iv<E1, E2, R> expression_type;
	typedef typename E1::size_type size_type;
	typedef typename E1::const_reference const_reference;
	typedef typename base_type<E1>::reference reference;
	typedef typename base_type<E1>::value_type value_type;
    typedef vp2<idx_cv_iv<E1, R>, E1, R > base_class;

	typedef idx_cv_iv<E1, R> inner_proxy_type;
	typedef idx_J_cv<inner_proxy_type,E2> proxy_type;
	typedef typename proxy_type::iterator iterator;
	typedef typename proxy_type::const_iterator const_iterator;

	idx_J_iv(E1 &m, const E2 &row, R &interleaved_col)
    :base_class(m,interleaved_col),proxy(base_class::inner,row){}

	// basic operation
	// vector operation
	size_type size(void) const {return proxy.size();}
	const_reference operator()(size_type i) const {return proxy(i);}
	reference operator()(size_type i) {return proxy(i);}
	const_reference operator[](size_type i) const {return proxy[i];}
	reference operator[](size_type i) {return proxy[i];}

	// 2D matrix operation
	size_type size(int dim) const {return proxy.size(dim);}
	const_reference operator()(size_type i, size_type j) const {return proxy(i,j);}
	reference operator()(size_type i, size_type j) {return proxy(i,j);}

	const_iterator begin(void) const {return proxy.begin();}
	const_iterator end(void) const {return proxy.end();}
	iterator begin(void) {return proxy.begin();}
	iterator end(void) {return proxy.end();}

	PROXY_OPS;

	proxy_type proxy;
};

template <typename E, typename R1, typename R2>
class idx_range_iv : public matrix_expression<idx_range_iv<E, R1, R2> >, public vp2<idx_cv_iv<E, R2>, E, R2 > {
public:
    typedef idx_range_iv<E, R1, R2> expression_type;
	typedef typename E::size_type size_type;
	typedef typename E::const_reference const_reference;
	typedef typename base_type<E>::reference reference;
	typedef typename base_type<E>::value_type value_type;
    typedef vp2<idx_cv_iv<E, R2>, E, R2 > base_class;

	typedef idx_cv_iv<E, R2> inner_proxy_type;
	typedef idx_range_cv<inner_proxy_type,R1> proxy_type;
	typedef typename proxy_type::iterator iterator;
	typedef typename proxy_type::const_iterator const_iterator;

	idx_range_iv(E &m, R1 &range_row, R2 &interleaved_col)
    :base_class(m,interleaved_col),proxy(base_class::inner,range_row){}

	// basic operation
	// vector operation
	size_type size(void) const {return proxy.size();}
	const_reference operator()(size_type i) const {return proxy(i);}
	reference operator()(size_type i) {return proxy(i);}
	const_reference operator[](size_type i) const {return proxy[i];}
	reference operator[](size_type i) {return proxy[i];}

	// 2D matrix operation
	size_type size(int dim) const {return proxy.size(dim);}
	const_reference operator()(size_type i, size_type j) const {return proxy(i,j);}
	reference operator()(size_type i, size_type j) {return proxy(i,j);}

	const_iterator begin(void) const {return proxy.begin();}
	const_iterator end(void) const {return proxy.end();}
	iterator begin(void) {return proxy.begin();}
	iterator end(void) {return proxy.end();}

	PROXY_OPS;

	proxy_type proxy;
};

template <typename E, typename R1, typename R2>
class idx_iv_iv : public matrix_expression<idx_iv_iv<E, R1, R2> >, public vp2<idx_cv_iv<E, R2>, E, R2 > {
public:
    typedef idx_iv_iv<E, R1, R2> expression_type;
	typedef typename E::size_type size_type;
	typedef typename E::const_reference const_reference;
	typedef typename base_type<E>::reference reference;
    typedef typename base_type<E>::value_type value_type;
    typedef  vp2<idx_cv_iv<E, R2>, E, R2 > base_class;

	typedef idx_cv_iv<E, R2> inner_proxy_type;
	typedef idx_iv_cv<inner_proxy_type,R1> proxy_type;
	typedef typename proxy_type::iterator iterator;
	typedef typename proxy_type::const_iterator const_iterator;

	idx_iv_iv(E &m, R1 &interleaved_row, R2 &interleaved_col)
    :base_class(m,interleaved_col),proxy(base_class::inner,interleaved_row){}

	// basic operation
	// vector operation
	size_type size(void) const {return proxy.size();}
	const_reference operator()(size_type i) const {return proxy(i);}
	reference operator()(size_type i) {return proxy(i);}
	const_reference operator[](size_type i) const {return proxy[i];}
	reference operator[](size_type i) {return proxy[i];}

	// 2D matrix operation
	size_type size(int dim) const {return proxy.size(dim);}
	const_reference operator()(size_type i, size_type j) const {return proxy(i,j);}
	reference operator()(size_type i, size_type j) {return proxy(i,j);}

	const_iterator begin(void) const {return proxy.begin();}
	const_iterator end(void) const {return proxy.end();}
	iterator begin(void) {return proxy.begin();}
	iterator end(void) {return proxy.end();}

	PROXY_OPS;

	proxy_type proxy;
};

//1D
template <typename E, typename E1>
class idx_I : public matrix_expression<idx_I<E, E1> > {
public:
    typedef idx_I<E, E1> expression_type;
	typedef typename E::size_type size_type;
	typedef typename E::const_reference const_reference;
	typedef typename base_type<E>::reference reference;
	typedef typename base_type<E>::value_type value_type;

	idx_I(E &m, const E1 &I):m_m(m), m_I(I) {}

	// basic operation
	// vector operation
	size_type size(void) const {return m_I.size();}
	const_reference operator[](idx_type i) const {return m_m(m_I(i));}
	const_reference operator()(idx_type i) const {return m_m(m_I(i));}
	reference operator[](idx_type i) {return m_m(m_I(i));}
	reference operator()(idx_type i) {return m_m(m_I(i));}

	// 2D matrix operation
	size_type size(int dim) const {return (1==dim)?m_I.size():1;}
	const_reference operator()(idx_type row, idx_type col) const {assert(col == 0);return m_m(m_I(row));}
	reference operator()(idx_type row, idx_type col) {assert(col == 0);return m_m(m_I(row));}

	template <typename E2>
	struct iterator_base : public default_matrix_iterator<value_type> {
		iterator_base(){}
		iterator_base(E2 *m, typename E1::const_iterator iter_I):m_m(m), m_iter_I(iter_I) {}
		const_reference operator *(void) const {return (*m_m)(*m_iter_I);}
		iterator_base operator++(void) {++m_iter_I; return *this;}
		iterator_base operator+=(size_type n) {m_iter_I+=n; return *this;}
		iterator_base operator+(size_type n) const {iterator_base i=*this; i+=n; return i;}
		template <typename I1>
		bool operator != (const I1 &i) const {assert(m_m == i.m_m); return m_iter_I != i.m_iter_I;}
		template <typename I1>
		bool operator == (const I1 &i) const {return !this->operator!=(i);}
		E2 *m_m;
		typename E1::const_iterator m_iter_I;
	};
	typedef iterator_base<const E> const_iterator;
	struct iterator : public iterator_base<E> {
		iterator(){}
		iterator(E *m, typename E1::const_iterator iter_I):iterator_base<E>(m, iter_I) {}
		operator const_iterator() const {return const_iterator(m_m, this->m_i);}
		iterator operator+(size_type n) const {iterator i=*this; i+=n; return i;}
		reference operator *(void) {return (*this->m_m)(*this->m_iter_I);}
	};

	const_iterator begin(void) const {return const_iterator(&m_m, m_I.begin());}
	const_iterator end(void) const {return const_iterator(&m_m, m_I.end());}
	iterator begin(void) {return iterator(&m_m, m_I.begin());}
	iterator end(void) {return iterator(&m_m, m_I.end());}

	PROXY_OPS;

	E &m_m;
	const E1 &m_I;
};

template <typename E>
class idx_cv : public matrix_expression<idx_cv<E> > {
public:
    typedef idx_cv<E> expression_type;
	typedef typename E::size_type size_type;
	typedef typename E::const_reference const_reference;
	typedef typename base_type<E>::reference reference;
	typedef typename base_type<E>::value_type value_type;

	idx_cv(E &m):m_m(m) {}

	// basic operation
	// vector operation
	size_type size(void) const {return m_m.size();}
	const_reference operator()(size_type i) const {return m_m(i);}
	reference operator()(size_type i) {return m_m(i);}
	const_reference operator[](idx_type i) const {return m_m[i];}
	reference operator[](idx_type i) {return m_m[i];}

	// 2D matrix operation
	size_type size(int dim) const {return (dim==1)?m_m.size():1;}
	const_reference operator()(size_type i, size_type j) const {
		assert(j == 0);
		return m_m(i);
	}
	reference operator()(size_type i, size_type j) {
		assert(j == 0);
		return m_m(i);
	}

	typedef typename E::const_iterator const_iterator;
	typedef typename base_type<E>::iterator iterator;

	const_iterator begin(void) const {return m_m.begin();}
	const_iterator end(void) const {return m_m.end();}
	iterator begin(void) {return m_m.begin();}
	iterator end(void) {return m_m.end();}

	PROXY_OPS;

	E &m_m;
};

template <typename E, typename T>
class idx_range : public matrix_expression<idx_range<E, T> > {
public:
    typedef idx_range<E, T> expression_type;
	typedef typename E::size_type size_type;
	typedef typename E::const_reference const_reference;
	typedef typename base_type<E>::reference reference;
	typedef typename base_type<E>::value_type value_type;

	idx_range(E &m, const range<T> &range):m_m(m), m_range(range) {}

	// basic operation
	// vector operation
	size_type size(void) const {return m_range.size();}
	const_reference operator()(size_type i) const {return m_m(m_range(i));}
	reference operator()(size_type i) {return m_m(m_range(i));}
	const_reference operator[](idx_type i) const {return m_m(m_range(i));}
	reference operator[](idx_type i) {return m_m(m_range(i));}

	// 2D matrix operation
	size_type size(int dim) const {return (dim==1)?m_range.size():1;}
	const_reference operator()(size_type i, size_type j) const {
		assert(j == 0);
		return m_m(m_range(i));
	}
	reference operator()(size_type i, size_type j) {
		assert(j == 0);
		return m_m(m_range(i));
	}

	typedef typename E::const_iterator const_iterator;
	typedef typename base_type<E>::iterator iterator;

	const_iterator begin(void) const {return m_m.begin()+m_range.m_front;}
	const_iterator end(void) const {return m_m.begin()+m_range.m_end;}
	iterator begin(void) {return m_m.begin()+m_range.m_front;}
	iterator end(void) {return m_m.begin()+m_range.m_end;}

	PROXY_OPS;

	E &m_m;
	range<T> m_range;
};

template <typename E, typename T>
class idx_iv : public matrix_expression<idx_iv<E, T> > {
public:
    typedef idx_iv<E, T> expression_type;
	typedef typename E::size_type size_type;
	typedef typename E::const_reference const_reference;
	typedef typename base_type<E>::reference reference;
	typedef typename base_type<E>::value_type value_type;

	idx_iv(E &m, const interleaved<T> &inter):m_m(m), m_inter(inter) {}

	// basic operation
	// vector operation
	size_type size(void) const {return m_inter.size();}
	const_reference operator()(size_type i) const {return m_m(m_inter(i));}
	reference operator()(size_type i) {return m_m(m_inter(i));}
	const_reference operator[](idx_type i) const {return m_m(m_inter(i));}
	reference operator[](idx_type i) {return m_m(m_inter(i));}

	// 2D matrix operation
	size_type size(int dim) const {return (dim==1)?m_inter.size():1;}
	const_reference operator()(size_type i, size_type j) const {
		assert(j == 0);
		return m_m(m_inter(i));
	}
	reference operator()(size_type i, size_type j) {
		assert(j == 0);
		return m_m(m_inter(i));
	}

	
    template <typename E2,typename I>
	struct iterator_base : public default_matrix_iterator<value_type> {
		
		iterator_base(){}
        iterator_base(E2 *m, size_type i):m_m(m),m_id(i) {
			m_i=m->m_m.begin()+m_m->m_inter[i];
		}
		
		const_reference operator *(void) const {return *m_i;}
		iterator_base operator++(void) {m_i+=m_m->m_inter.m_step;m_id++; return *this;}
		iterator_base operator+=(size_type n) {m_i+=m_m->m_inter.m_step*n;m_id+=n; return *this;}
		iterator_base operator+(size_type n) const {iterator_base i=*this; i+=n; return i;}
		
		template <typename I1>bool operator != (const I1 &i) const {assert(m_m == i.m_m); return m_i != i.m_i;}
		template <typename I1>bool operator == (const I1 &i) const {return !this->operator!=(i);}
		
        E2 *m_m;
		I m_i;
		size_type m_id;
	};

	typedef iterator_base<const expression_type,typename E::const_iterator > const_iterator;
	
	struct iterator : public iterator_base<expression_type, typename base_type<E>::iterator > {
        typedef iterator_base<expression_type, typename base_type<E>::iterator > base_class;
		iterator(){}
        iterator(expression_type *m, size_type i)
          :iterator_base<expression_type,typename base_type<E>::iterator>(m, i) {}
		operator const_iterator() const {return const_iterator(m_m, this->m_id);}
		iterator operator+(size_type n) const {iterator i=*this; i+=n; return i;}
        reference operator *(void) {return *base_class::m_i;}
	};

	const_iterator begin(void) const {return const_iterator(this, 0);}
	const_iterator end(void) const {return const_iterator(this, size());}
	iterator begin(void) {return iterator(this, 0);}
	iterator end(void) {return iterator(this, size());}

	PROXY_OPS;

	E &m_m;
	interleaved<T> m_inter;
};

}	// namespace matrix
}	// namespace zjucad
