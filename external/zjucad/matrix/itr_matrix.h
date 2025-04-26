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

#ifndef _ZJUCAD_MATRIX_ITR_MATRIX_H_
#define _ZJUCAD_MATRIX_ITR_MATRIX_H_

#include "matrix_expression.h"
#include "matrix_proxy.h"

namespace zjucad { namespace matrix {

// convert a random access iterator into a column major matrix

template <typename ITR>
class itr_matrix : public matrix_expression<itr_matrix<ITR> >
{
public:
	typedef itr_matrix<ITR> expression_type;
	typedef itr_matrix<ITR> E;

	typedef typename ITR::value_type value_type;
	typedef value_type& reference;
	typedef const value_type& const_reference;
	typedef ITR iterator;
	typedef ITR const_iterator;

	typedef zjucad::matrix::size_type size_type;

	// constructor
	itr_matrix():beg_(0), nrows_(0), ncols_(0){}
	itr_matrix(size_type size, ITR beg):beg_(beg), nrows_(size), ncols_(1){}
	itr_matrix(idx_type nrows, idx_type ncols, ITR beg):beg_(beg), nrows_(nrows), ncols_(ncols){}
	itr_matrix(const itr_matrix<ITR> &m) {
		*this = m;
	}

	// size
	size_type size(void) const {return nrows_*ncols_;} 
	size_type size(int dim) const {return (dim==1)?nrows_:ncols_;}

	void resize(size_type size) {
		nrows_ = size, ncols_ = 1;
	}
	void resize(size_type nrows, size_type ncols) {
		nrows_ = nrows, ncols_ = ncols;
	}

	// element access
	const_reference operator[](idx_type i) const {return beg_[i];}
	const_reference operator()(idx_type i) const {return beg_[i];}
	const_reference operator()(idx_type row, idx_type col) const {
		return beg_[row+col*nrows_];
	}
	reference operator[](idx_type i) {return beg_[i];}
	reference operator()(idx_type i) {return beg_[i];}
	reference operator()(idx_type row, idx_type col) {
		return beg_[row+col*nrows_];
	}

	// iterator access
	const_iterator begin(void) const {return beg_;}
	const_iterator end(void) const {return beg_+ncols_*nrows_;}
	iterator begin(void) {return beg_;}
	iterator end(void) {return beg_+ncols_*nrows_;}

	template <typename E>
	const itr_matrix<ITR>& operator=(const matrix_expression<E> &e){
		assert(size(1) == e().size(1) || size(2) == e().size(2));
		std::copy(e().begin(), e().end(), begin());
		return *this;
	}
	const itr_matrix<ITR>& operator=(const itr_matrix<ITR> &m){
		assert(size(1) == m.size(1) || size(2) == m.size(2));
		std::copy(m.begin(), m.end(), begin());
		return *this;
	}

	PROXY_ACCESS;
	MATRIX_SELF_OP;

	ITR beg_;	
	size_type nrows_, ncols_;
};


// e.g: itr_matrix<int *>, T=int
template <typename T>
class itr_matrix<T*> : public matrix_expression<itr_matrix<T*> >
{
public:
 	typedef itr_matrix<T*> expression_type;
 	typedef itr_matrix<T*> E;

 	typedef zjucad::matrix::size_type size_type;

	typedef T value_type;
	typedef T& reference;
	typedef const T& const_reference;
	typedef T* iterator;
	typedef T const * const_iterator;

	itr_matrix(size_type nrows, size_type ncols, T *beg)
		:nrows_(nrows), ncols_(ncols), beg_(beg)
		{}

	size_type size(void) const { return nrows_*ncols_; }
 	size_type size(int dim) const {return (dim==1)?nrows_:ncols_;}

	void resize(size_type size) {
		nrows_ = size, ncols_ = 1;
	}
	void resize(size_type nrows, size_type ncols) {
		nrows_ = nrows, ncols_ = ncols;
	}

	// element access
	const_reference operator[](idx_type i) const {return beg_[i];}
	const_reference operator()(idx_type i) const {return beg_[i];}
	const_reference operator()(idx_type row, idx_type col) const {
		return beg_[row+col*nrows_];
	}
	reference operator[](idx_type i) {return beg_[i];}
	reference operator()(idx_type i) {return beg_[i];}
	reference operator()(idx_type row, idx_type col) {
		return beg_[row+col*nrows_];
	}

	// iterator access
	const_iterator begin(void) const {return beg_;}
	const_iterator end(void) const {return beg_+size();}
	iterator begin(void) {return beg_;}
	iterator end(void) {return beg_+size();}

	template <typename E>
	const itr_matrix<T*>& operator=(const matrix_expression<E> &e){
		assert(size(1) == e().size(1) || size(2) == e().size(2));
		std::copy(e().begin(), e().end(), begin());
		return *this;
	}
	const itr_matrix<T*>& operator=(const itr_matrix<T*> &m){
		assert(size(1) == m.size(1) || size(2) == m.size(2));
		std::copy(m.begin(), m.end(), begin());
		return *this;
	}

 	PROXY_ACCESS;
 	MATRIX_SELF_OP;

protected:
	size_type nrows_, ncols_;
	T *beg_;
};

}}

#endif
