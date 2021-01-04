//This file contains some function to solve function interpolation.
//Author: Chen Jinlong
//Time: 12:11 2021/1/4

#ifndef FUNCINTERPOLATION
#define FUNCINTERPOLATION

#include "matrix.h"

using namespace matrix;

namespace function
{
	//��������x(1:n):�Ա���,y(1,n):����ֵ,bond:��ֵ��������Ӧ���ֱ߽�����
	//left_bond,right_bond:���ұ߽�ֵ
	//Ĭ���������Ȼ����������ֵ
	//����ֵ����ֵ����ϵ�����������ݴδӸߵ���
	Mat Cubic_Spline(Mat& x, Mat& y, int bond = 1, double left_bond = 0.0, double right_bond = 0.0);
}
#endif