//This file contains some function and classes to store marix.
//Author: Chen Jinlong
//Time: 19:15 2020/9/25


#ifndef MATRIX
#define MATRIX

#include <iostream>

namespace matrix
{
	class Mat
	{
	private:
		int row_num;
		int column_num;
		double** p;
	public:
		Mat(int row = 0, int column = 0);
		Mat(const Mat& m);
		Mat(const Mat& m, const Mat& n);//增广矩阵 Augmented matrix
		~Mat();
		int row() const { return row_num; }
		int column() const { return column_num; }
		bool symmetric() const;//Judge if it is a  symmetric matrix.判断矩阵是否对称
		Mat T() const;//矩阵转置 Matrix transposition
		double norm_inf() const;//无穷范数
		Mat& operator=(const Mat& m);
		Mat& operator=(double* num);
		Mat operator+(const Mat& m) const;
		Mat operator-(const Mat& m) const;
		Mat operator*(const Mat& m) const;
		Mat operator*(const double& m) const;
		double* operator[](int n) const { return p[n]; }
		friend std::ostream& operator<<(std::ostream& os, const Mat& m);
		friend Mat operator*(const double& t, const Mat& m) { return m * t; }
	};

	Mat I(int n);//单位矩阵 Identity Matrix
}


#endif //MATRIX 

