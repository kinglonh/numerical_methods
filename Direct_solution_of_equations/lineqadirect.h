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
}

#endif // !LINEQADIRECT
