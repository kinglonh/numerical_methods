//This file contains some function to solve the linear equations directly.
//Author: Chen Jinlong
//Time: 19:08 2020/9/25

#ifndef LINEQADIRECT
#define LINEQADIRECT

#include "matrix.h"

using namespace matrix;

namespace lineqa
{
	Mat Gauss(const Mat& A, const Mat& b);
	Mat Gauss_MaxColumn(const Mat& A, const Mat& b);
	Mat Doolittle(const Mat& A);
	Mat Doolittle_Solve(const Mat& A, const Mat& b);
	Mat Cholesky(const Mat& A);
	Mat Cholesky_Solve(const Mat& A, const Mat& b);
	Mat LDL(const Mat& A);
	Mat LDL_Solve(const Mat& A, const Mat& b);
	void QR_Householder(const Mat& A, Mat& G, Mat& d, Mat& alpha);
	void QR_Householder_Show(const Mat& A);
	void QR_Householder_Show(const Mat& G, const Mat& d, const Mat& alpha);
	Mat QR_Householder_Solve(const Mat& A, const Mat& b);
}

#endif // !LINEQADIRECT
