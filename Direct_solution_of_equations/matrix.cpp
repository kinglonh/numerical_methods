//This file contains some function and classes to store marix.
//Author: Chen Jinlong
//Time: 19:36 2020/9/25

#include "matrix.h"

namespace matrix
{
	Mat::Mat(int row, int column)
	{
		row_num = row;
		column_num = column;
		p = new double* [row_num];
		int i, j;
		for (i = 0; i < row_num; i++)
		{
			p[i] = new double[column_num];
		}
		for (i = 0; i < row_num; i++)
			for (j = 0; j < column_num; j++)
				p[i][j] = 0;
	}

	Mat::Mat(const Mat& m)
	{
		row_num = m.row_num;
		column_num = m.column_num;
		p = new double* [m.row_num];
		int i, j;
		for (i = 0; i < m.row_num; i++)
		{
			p[i] = new double[m.column_num];
		}
		for (i = 0; i < m.row_num; i++)
			for (j = 0; j < m.column_num; j++)
				p[i][j] = m[i][j];
	}

	Mat::Mat(const Mat& m, const Mat& n)
	{
		if (m.row_num != n.row_num)
		{
			throw "m and n must have same rows.";
			exit(1);
		}
		else
		{
			int i, j;
			row_num = m.row_num;
			column_num = m.column_num + n.column_num;
			p = new double* [row_num];
			for (i = 0; i < m.row_num; i++)
			{
				p[i] = new double[column_num];
			}
			for (i = 0; i < m.row_num; i++)
			{
				for (j = 0; j < m.column_num; j++)
					p[i][j] = m[i][j];
				for (j = 0; j < n.column_num; j++)
					p[i][m.column_num + j] = n[i][j];
			}
		}
	}

	Mat::~Mat()
	{
		int i;
		for (i = 0; i < row_num; i++)
		{
			delete[] p[i];
		}
		delete[] p;
	}

	Mat& Mat::operator=(const Mat& m)
	{
		if (this == &m)
			return *this;
		else
		{
			int i, j;
			for (i = 0; i < row_num; i++)
			{
				delete[] p[i];
			}
			delete[] p;
			row_num = m.row_num;
			column_num = m.column_num;
			p = new double* [m.row_num];
			for (i = 0; i < m.row_num; i++)
			{
				p[i] = new double[m.column_num];
			}
			for (i = 0; i < m.row_num; i++)
				for (j = 0; j < m.column_num; j++)
					p[i][j] = m[i][j];
			return *this;
		}
	}

	Mat& Mat::operator=(double* num)
	{
		if (row_num == 0 || column_num == 0)
		{
			throw "The number of columns and columns in the matrix is not defined!";
			exit(1);
		}
		else
		{
			int i, j;
			for (i = 0; i < row_num; i++)
				for (j = 0; j < column_num; j++)
					p[i][j] = num[i*column_num + j];
		}
		return *this;
	}

	std::ostream& operator<<(std::ostream& os, const Mat& m)
	{
		int i, j;
		for (i = 0; i < m.row_num; i++)
		{
			os << std::endl;
			for (j = 0; j < m.column_num; j++)
				os << m[i][j] << "\t";
		}
		os << std::endl;
		return os;
	}
}