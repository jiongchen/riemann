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

#if defined __GNUC__
	#include <bits/allocator.h>
#elif defined _MSC_VER
	#include <xmemory>
#endif
using std::allocator;

#ifdef __GNUG__
#  define DEPRECATED __attribute__((deprecated))
#endif

#ifdef _MSC_VER
#  define DEPRECATED __declspec(deprecated)
#endif

namespace zjucad { namespace matrix {

// Unbounded array
template<typename T, typename A=allocator<T> >
class unbounded_array {
public:
	typedef zjucad::matrix::size_type size_type;
	typedef typename A::value_type value_type;
	typedef typename A::const_reference const_reference;
	typedef typename A::reference reference;
	typedef typename A::const_pointer const_pointer;
	typedef typename A::pointer pointer;

    // Construction and destruction
    unbounded_array ():size_ (0), data_ (0) {}
    unbounded_array (size_type size):size_ (0), data_(0){construct(size);}
    unbounded_array (const unbounded_array &a):
        size_ (0), data_ (0) {*this = a;}
    ~unbounded_array () {
		destroy();
	}

    // Resizing
    void resize (size_type size) {
        if (size != size_) {
			destroy();
			construct(size);
        }
    }

    size_type size () const {return size_;}

    // Element access
    const_reference operator [] (size_type i) const {
        assert(i >= 0 && i < size_);
        return data_ [i];
    }
    reference operator [] (size_type i) {
        assert(i >= 0 && i < size_);
        return data_ [i];
    }

    // Assignment
    unbounded_array &operator = (const unbounded_array &a) {
        if (this != &a) {
            resize (a.size_);
#if _MSC_VER // stupid vc stl checking
			if(data_)
#endif
            std::copy (a.data_, a.data_ + a.size_, data_);
        }
        return *this;
    }

    typedef const_pointer const_iterator;
    const_iterator begin () const {return data_;}
    const_iterator end () const {return data_ + size_;}
    typedef pointer iterator;
    iterator begin () {return data_;}
    iterator end () {return data_ + size_;}

	void swap (unbounded_array<T, A> &a) {
		if (this != &a) {
			std::swap (size_, a.size_);
			std::swap (data_, a.data_);
		}
	}

    size_type size_;
    pointer data_;
	A alloc_;
//protected: // cannot make friend template function in g++?
	void construct(size_type size) {
		assert(data_ == 0);
		data_ = alloc_.allocate(size);
		size_ = size;
		pointer data = data_;
		for(size_type i = 0; i < size_; ++i, ++data)
			alloc_.construct(data, value_type());
	}
	void destroy(void) {
		pointer data = data_;
		for(size_type i = 0; i < size_; ++i, ++data)
			alloc_.destroy(data);
		alloc_.deallocate(data_, size_);
		data_ = 0;
		size_ = 0;
	}
	//template <typename IS, typename T, typename A>
	//friend void read(IS &is, unbounded_array<T, A> &array);
};

template <typename OS, typename T, typename A>
void write(OS &os, const unbounded_array<T, A> &array);
template <typename IS, typename T, typename A>
int read(IS &is, unbounded_array<T, A> &array);

template <typename T1>
struct is_simple_obj {	// default is not simple_obj
	template <typename OS, typename T, typename A>
	static void write0(OS &os, const unbounded_array<T, A> &array) {
		os.write((const char *)&array.size_, sizeof(typename unbounded_array<T, A>::size_type));
		typename unbounded_array<T, A>::const_iterator iter = array.begin(), end = array.end();
		for(; iter != end; ++iter)
			write(os, *iter);
	}
	template <typename IS, typename T, typename A>
	static int read0(IS &is, unbounded_array<T, A> &array) {
		typename unbounded_array<T, A>::size_type size;
		is.read((char *)&size, sizeof(typename unbounded_array<T, A>::size_type));
		if(is.fail())    return 1;
		array.resize(size);
		typename unbounded_array<T, A>::iterator iter = array.begin(), end = array.end();
		for(; iter != end; ++iter) {
			read(is, *iter);
			if(is.fail())    return 1;
		}
		return 0;
	}
};

#define SIMPLE_OBJ(T1)\
template <>\
struct is_simple_obj<T1> {\
	template <typename OS, typename T, typename A>\
	static void write0(OS &os, const unbounded_array<T, A> &array) {\
		os.write((const char *)&array.size_, sizeof(typename unbounded_array<T, A>::size_type));\
		os.write((const char *)array.data_, array.size_*sizeof(typename unbounded_array<T, A>::value_type));\
	}\
	template <typename IS, typename T, typename A>\
	static int read0(IS &is, unbounded_array<T, A> &array) {\
		is.read((char *)&array.size_, sizeof(typename unbounded_array<T, A>::size_type));\
        if(is.fail())    return 1;\
		array.data_ = array.alloc_.allocate(array.size_);\
		is.read((char *)array.data_, array.size_*sizeof(typename unbounded_array<T, A>::value_type));\
        if(is.fail())    return 1;\
        return 0;\
	}\
};

SIMPLE_OBJ(char);
SIMPLE_OBJ(unsigned char);
SIMPLE_OBJ(signed char);
SIMPLE_OBJ(short);
SIMPLE_OBJ(unsigned short);
SIMPLE_OBJ(int);
SIMPLE_OBJ(unsigned int);
SIMPLE_OBJ(long);
SIMPLE_OBJ(unsigned long);
SIMPLE_OBJ(float);
SIMPLE_OBJ(double);
SIMPLE_OBJ(long double);

template <typename OS, typename T, typename A>
DEPRECATED void write(OS &os, const unbounded_array<T, A> &array) {
	is_simple_obj<T>::write0(os, array);
}

template <typename IS, typename T, typename A>
DEPRECATED int read(IS &is, unbounded_array<T, A> &array) {
	array.destroy();
	return is_simple_obj<T>::read0(is, array);
}

}	// namespace matrix
}	// namespace zjucad
