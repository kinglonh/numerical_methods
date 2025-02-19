//This file contains some function to solve the linear equations by iterative method.
//Author: Chen Jinlong
//Time: 21:34 2020/10/16

#ifndef LINEQAITERATIVE
#define LINEQAITERATIVE

#include "matrix.h"

using namespace matrix;

namespace lineqa
{
	void Jacobi(const Mat& A, const Mat& b, const double& error, Mat& x);
	void GaussSeidel(const Mat& A, const Mat& b, const double& error, Mat& x);
	void SOR(const Mat& A, const Mat& b, const double& omega, const double& error, Mat& x);
	void Conjugate_Gradient(const Mat& A, const Mat& b, const double& error, Mat& x);
}
#endif // !LINEQAITERATIVE
