//This file contains some function to solve function interpolation.
//Author: Chen Jinlong
//Time: 12:11 2021/1/4

#ifndef FUNCINTERPOLATION
#define FUNCINTERPOLATION

#include "matrix.h"

using namespace matrix;

namespace function
{
	//三次样条x(1:n):自变量,y(1,n):函数值,bond:插值条件，对应三种边界条件
	//left_bond,right_bond:左右边界值
	//默认情况：自然三次样条插值
	//返回值：插值函数系数，从左到右幂次从高到低
	Mat Cubic_Spline(Mat& x, Mat& y, int bond = 1, double left_bond = 0.0, double right_bond = 0.0);
}
#endif