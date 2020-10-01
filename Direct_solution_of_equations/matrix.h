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
		Mat(const Mat& m, const Mat& n);//Ôö¹ã¾ØÕó Augmented matrix
		~Mat();
		int row() const { return row_num; }
		int column() const { return column_num; }
		bool symmetric() const;//Judge if it is a  symmetric matrix.ÅĞ¶Ï¾ØÕóÊÇ·ñ¶Ô³Æ
		Mat T() const;//¾ØÕó×ªÖÃ Matrix transposition
		Mat& operator=(const Mat& m);
		Mat& operator=(double* num);
		double* operator[](int n) const { return p[n]; }
		friend std::ostream& operator<<(std::ostream& os, const Mat& m);
	};
}


#endif //MATRIX 

