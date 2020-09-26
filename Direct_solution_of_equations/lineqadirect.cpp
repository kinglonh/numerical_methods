//This file contains some function to solve the linear equations directly.
//Author: Chen Jinlong
//Time: 13:47 2020/9/26
#include <iostream>
#include <cmath>
#include "lineqadirect.h"

using namespace matrix;

namespace lineqa
{
	Mat Gauss(const Mat& A, const Mat& b)
	{
		if (A.row() == 0 || A.row() != A.column() || b.column() != 1 || A.row() != b.row())
		{
			throw "Fail to solve equations: Parameter error.";
			exit(1);
		}
		else
		{
			Mat x(b.row(), 1);
			Mat G = A, d = b;
			int i, j, k;
			int n = A.row();
			for (k = 1; k < n; k++) //step: 1 to n-1
			{
				if (abs(G[k - 1][k - 1]) < 1e-5)
				{
					throw "The value of A[k - 1][k - 1] is too small, please try Gauss Elimination with Maximal Column Pivoting.";
					exit(1);
				}
				else
				{
					for (i = k; i < n; i++)
					{
						double l = G[i][k - 1] / G[k - 1][k - 1];
						G[i][k - 1] = 0;
						for (j = k; j < n; j++)
						{
							G[i][j] = G[i][j] - l * G[k - 1][j];
						}
						d[i][0] = d[i][0] - l * d[k - 1][0];
					}
				}
			}
			
			for (k = n; k > 0; k--)
			{
				double s = d[k - 1][0];
				for (i = k; i < n; i++)
				{
					s = s - G[k - 1][i] * x[i][0];
				}
				x[k - 1][0] = s / G[k - 1][k - 1];
			}
			return x;
		}
	}

	Mat Gauss_MaxColumn(const Mat& A, const Mat& b)
	{
		if (A.row() == 0 || A.row() != A.column() || b.column() != 1 || A.row() != b.row())
		{
			throw "Fail to solve equations: Parameter error.";
			exit(1);
		}
		else
		{
			Mat x(b.row(), 1);
			Mat G = A, d = b;
			int i, j, k;
			int n = A.row();
			for (k = 1; k < n; k++) //step: 1 to n-1
			{
				double max_column = abs(G[k - 1][k - 1]);
				int m, max_column_num = k - 1;
				for (m = k; m < n; m++)
				{
					if (abs(G[k - 1][m]) > max_column)
					{
						max_column = abs(G[k - 1][m]);
						max_column_num = m;
					}
				}
				if (max_column_num != (k - 1))
				{
					double t;
					for (m = k - 1; m < n; m++)
					{
						t = G[k - 1][m];
						G[k - 1][m] = G[max_column_num][m];
						G[max_column_num][m] = t;
					}
					t = d[k - 1][0];
					d[k - 1][0] = d[max_column_num][0];
					d[max_column_num][0] = t;
				}
				if (abs(G[k - 1][k - 1]) < 1e-5)
				{
					throw "The value of A[k - 1][k - 1] is the largest number in the column. but it is too small. Please try another solution";
					exit(1);
				}
				else
				{
					for (i = k; i < n; i++)
					{
						double l = G[i][k - 1] / G[k - 1][k - 1];
						G[i][k - 1] = 0;
						for (j = k; j < n; j++)
						{
							G[i][j] = G[i][j] - l * G[k - 1][j];
						}
						d[i][0] = d[i][0] - l * d[k - 1][0];
					}
				}
			}

			for (k = n; k > 0; k--)
			{
				double s = d[k - 1][0];
				for (i = k; i < n; i++)
				{
					s = s - G[k - 1][i] * x[i][0];
				}
				x[k - 1][0] = s / G[k - 1][k - 1];
			}
			return x;
		}
	}
}