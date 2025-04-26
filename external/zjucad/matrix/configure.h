#ifndef _ZJUCAD_MATRIX_H_
#define _ZJUCAD_MATRIX_H_

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

#include <stddef.h>

namespace zjucad { namespace matrix {

//! use machine pointer size as the type
//! try to avoid using unsigned.

typedef ptrdiff_t size_type;
typedef ptrdiff_t idx_type;
typedef ptrdiff_t offset_type;

#define ZJUCAD_MATRIX_VERSION 1

}}

#endif
