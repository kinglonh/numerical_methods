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

	bool Mat::symmetric() const
	{
		int i, j;
		if (row_num != column_num)
			return false;
		for (i = 0; i < row_num; i++)
		{
			for (j = i + 1; j < row_num; j++)
			{
				if (p[i][j] != p[j][i])
					return false;
			}
		}
		return true;
	}

	Mat Mat::T() const
	{
		Mat TA(column_num, row_num);
		int i, j;
		for (i = 0; i < row_num; i++)
		{
			for (j = 0; j < column_num; j++)
			{
				TA[j][i] = p[i][j];
			}
		}
		return TA;
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

	Mat Mat::operator+(const Mat& m) const
	{
		if ((row_num != m.row_num) || (column_num != m.column_num))
		{
			throw "Inconsistent matrix dimensions!";
			exit(1);
		}
		else
		{
			Mat result(row_num, column_num);
			int i, j;
			for (i = 0; i < row_num; i++)
				for (j = 0; j < column_num; j++)
					result[i][j] = p[i][j] + m[i][j];
			return result;
		}
	}

	Mat Mat::operator-(const Mat& m) const
	{
		if ((row_num != m.row_num) || (column_num != m.column_num))
		{
			throw "Inconsistent matrix dimensions!";
			exit(1);
		}
		else
		{
			Mat result(row_num, column_num);
			int i, j;
			for (i = 0; i < row_num; i++)
				for (j = 0; j < column_num; j++)
					result[i][j] = p[i][j] - m[i][j];
			return result;
		}
	}

	Mat Mat::operator*(const Mat& m) const
	{
		if (column_num != m.row_num)
		{
			throw "Inconsistent matrix dimensions!";
			exit(1);
		}
		else
		{
			Mat result(row_num, m.column_num);
			int i, j, k;
			for (i = 0; i < row_num; i++)
				for (j = 0; j < m.column_num; j++)
					for (k = 0; k < column_num; k++)
						result[i][j] += (p[i][k] * m[k][j]);
			return result;
		}
	}

	Mat Mat::operator*(const double& m) const
	{
		Mat result(row_num, column_num);
		int i, j;
		for (i = 0; i < row_num; i++)
			for (j = 0; j < column_num; j++)
				result[i][j] = p[i][j] * m;
		return result;
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

	Mat I(int n)
	{
		Mat I(n, n);
		for (int i = 0; i < n; i++)
			I[i][i] = 1;
		return I;
	}
}