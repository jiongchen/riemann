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

#ifndef ZJUCAD_MATRIX_BLAS_H_
#define ZJUCAD_MATRIX_BLAS_H_

#include <zjucad/matrix/matrix.h>

namespace zjucad { namespace matrix {

template <typename T>
struct value_type_of_pointer
{
};

template <typename T>
struct value_type_of_pointer<T *>
{
	typedef T value_type;
};

#ifndef F77_INT
typedef value_type_of_pointer<FINT>::value_type F77_INT;
#endif

static char __F77_TRANS[] = {'n', 't'};
// GEMV
void inline gemv(float alpha, bool trans, const matrix<float> &A,
	const matrix<float> &x, float beta, matrix<float> &Ax)
{
	Ax.resize(trans?A.size(2):A.size(1), 1);
	F77_INT MN1[] = {F77_INT(A.size(1)), F77_INT(A.size(2)), F77_INT(1)};
	F77_sgemv(&__F77_TRANS[trans], MN1, MN1+1,
			  &alpha, &A[0], MN1,
			  &x[0], MN1+2,
			  &beta, &Ax[0], MN1+2);
}

void inline gemv(double alpha, bool trans, const matrix<double> &A,
	const matrix<double> &x, double beta, matrix<double> &Ax)
{
	Ax.resize(trans?A.size(2):A.size(1), 1);
	F77_INT MN1[] = {F77_INT(A.size(1)), F77_INT(A.size(2)), 1};
	F77_dgemv(&__F77_TRANS[trans], MN1, MN1+1,
			  &alpha, &A[0], MN1,
			  &x[0], MN1+2,
			  &beta, &Ax[0], MN1+2);
}

template <typename T>
void inline gemv(bool trans, const matrix<T> &A,
	const matrix<T> &x, matrix<T> &Ax)
{
	gemv(1, trans, A, x, 0, Ax);
}

// GEMM
void inline gemm(float alpha, bool transA, const matrix<float> &A,
	bool transB, const matrix<float> &B, float beta, matrix<float> &C)
{
	F77_INT MNK[] = {
		F77_INT(transA?A.size(2):A.size(1)),
		F77_INT(transB?B.size(1):B.size(2)),
		F77_INT(transA?A.size(1):A.size(2))
	};
	C.resize(MNK[0], MNK[1]);
	F77_INT LD_ABC[] = {
		F77_INT(A.size(1)), F77_INT(B.size(1)), F77_INT(C.size(1))
	};
	F77_sgemm(&__F77_TRANS[transA], &__F77_TRANS[transB],
			  MNK, MNK+1, MNK+2,
			  &alpha, &A[0], LD_ABC+0,
			  &B[0], LD_ABC+1,
			  &beta, &C[0], LD_ABC+2);
}

void inline gemm(double alpha, bool transA, const matrix<double> &A,
	bool transB, const matrix<double> &B, double beta, matrix<double> &C)
{
	F77_INT MNK[] = {
		F77_INT(transA?A.size(2):A.size(1)),
		F77_INT(transB?B.size(1):B.size(2)),
		F77_INT(transA?A.size(1):A.size(2))
	};
	C.resize(MNK[0], MNK[1]);
	F77_INT LD_ABC[] = {
		F77_INT(A.size(1)), F77_INT(B.size(1)), F77_INT(C.size(1))
	};
	F77_dgemm(&__F77_TRANS[transA], &__F77_TRANS[transB],
			  MNK, MNK+1, MNK+2,
			  &alpha, &A[0], LD_ABC+0,
			  &B[0], LD_ABC+1,
			  &beta, &C[0], LD_ABC+2);
}

template <typename T>
void inline gemm(bool transA, const matrix<T> &A,
	bool transB, const matrix<T> &B, matrix<T> &C)
{
	gemm(1, transA, A, transB, B, 0, C);
}

}} // zjucad::matrix

#endif
