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

#include <zjucad/matrix/matrix.h>
#include <algorithm>

namespace zjucad { namespace matrix {

int inline lu(matrix<float> &A, matrix<int> &ipiv)
{
	int M = A.size(1), N = A.size(2), info;
	ipiv.resize(std::min(M, N));
        clapack::sgetrf_(&M, &N, A.data().begin(), &M, ipiv.data().begin(), &info);
	return info;
}

int inline lu(matrix<double> &A, matrix<int> &ipiv)
{
	int M = A.size(1), N = A.size(2), info;
	ipiv.resize(std::min(M, N));
	clapack::dgetrf_(&M, &N, A.data().begin(), &M, ipiv.data().begin(), &info);
	return info;
}

int inline inv(matrix<float> &A)
{
	matrix<int> ipiv;
	lu(A, ipiv);
	int N = A.size(1), info, lwork = std::max(int(1),N);
	matrix<float> work(lwork);
	clapack::sgetri_(&N, A.data().begin(), &N, ipiv.data().begin(), work.data().begin(), &lwork, &info);
	return info;
}

int inline inv(matrix<double> &A)
{
	matrix<int> ipiv;
	lu(A, ipiv);
	int N = A.size(1), info, lwork = std::max(int(1),N);
	matrix<double> work(lwork);
	clapack::dgetri_(&N, A.data().begin(), &N, ipiv.data().begin(), work.data().begin(), &lwork, &info);
	return info;
}

int inline svd(matrix<float> &A, matrix<float> &U, matrix<float> &S, matrix<float> &VT)
{
	int m = A.size(1), n = A.size(2), info = -1;
	int min_mn = std::min(m,n);
	int lwork = std::max(3*min_mn+std::max(m,n), 5*min_mn);
	matrix<float> w(lwork), diagS(m);
	U.resize(m, min_mn); VT.resize(min_mn, n);
  static char sS[] = "S";
	clapack::sgesvd_(sS, sS, &m, &n, A.data().begin(), &m, diagS.data().begin(),
		U.data().begin(), &m, VT.data().begin(), &min_mn, w.data().begin(), &lwork, &info);
	S = zeros<float>(min_mn);
	for(int i = 0; i < min_mn; ++i) S(i, i) = diagS[i];
	return info;
}

int inline svd(matrix<double> &A, matrix<double> &U, matrix<double> &S, matrix<double> &VT)
{
	int m = A.size(1), n = A.size(2), info = -1;
	int min_mn = std::min(m,n);
	int lwork = std::max(3*min_mn+std::max(m,n), 5*min_mn);
	matrix<double> w(lwork), diagS(m);
	U.resize(m, min_mn); VT.resize(min_mn, n);
  static char sS[] = "S";
	clapack::dgesvd_(sS, sS, &m, &n, A.data().begin(), &m, diagS.data().begin(),
		U.data().begin(), &m, VT.data().begin(), &min_mn, w.data().begin(), &lwork, &info);
	S = zeros<double>(min_mn);
	for(int i = 0; i < min_mn; ++i) S(i, i) = diagS[i];
	return info;
}

template <typename T>
T inline det(matrix<T> &A)
{
	assert(A.size(1) == A.size(2));
	matrix<int> ipiv;
	lu(A, ipiv);
	T rtn = 1.0f;
	zjucad::matrix::size_type i;
	for(i = 0; i < A.size(); i += A.size(1)+1)
		rtn *= A[i];
	for(i = 0; i < ipiv.size(); ++i)
		if(i != ipiv[i]-1)
			rtn *= -1.0f;
	return rtn;
};

int inline eig(matrix<double> &A, matrix<double> &e)
{
	int N = A.size(1), lwork = -1, info;
	matrix<double> work(1);
  static char V[]="V", U[]="U";
	clapack::dsyev_(V, U, &N, &A[0], &N, &e[0], &work[0], &lwork, &info);
	lwork = static_cast<int>(work[0]+0.5);
	work.resize(lwork);
	clapack::dsyev_(V, U, &N, &A[0], &N, &e[0], &work[0], &lwork, &info);
	return info;
}

int inline eig(matrix<float> &A, matrix<float> &e)
{
	int N = A.size(1), lwork = -1, info;
	matrix<float> work(1);
  static char V[]="V", U[]="U";
	clapack::ssyev_(V, U, &N, &A[0], &N, &e[0], &work[0], &lwork, &info);
	lwork = static_cast<int>(work[0]+0.5);
	work.resize(lwork);
	clapack::ssyev_(V, U, &N, &A[0], &N, &e[0], &work[0], &lwork, &info);
	return info;
}

int inline dgesv(matrix<double> &A, matrix<double> &bx)
{
	int n = A.size(1), nrhs = bx.size(2), info;
	matrix<int> lpiv(n);
	clapack::dgesv_(&n, &nrhs,
		   &A[0], &n,
		   &lpiv[0],
		   &bx[0], &n,
		   &info);
	return info;
}

}
}
