//This file contains some function to solve function interpolation.
//Author: Chen Jinlong
//Time: 12:20 2021/1/4

#include <iostream>
#include <cmath>
#include "lineqadirect.h"
#include "funcinterpolation.h"

using namespace matrix;
using namespace lineqa;

namespace function
{
	Mat Cubic_Spline(Mat& x, Mat& y, int bond, double left_bond, double right_bond)
	{
		if (x.row() != 1 || y.row() != 1 || x.column() < 3 || x.column() != y.column() || bond < 1 || bond >3)
		{
			throw "Fail to solve equations: Parameter error.";
			exit(1);
		}
		else
		{
			int n = x.column() - 1;
			Mat S(n, 4);
			Mat h(1, n);
			Mat miu(1, n);
			Mat lamda(1, n);
			Mat d(n + 1, 1);
			Mat M(n + 1, 1);
			int i;
			for (i = 0; i < n; i++)
			{
				h[0][i] = x[0][i + 1] - x[0][i];
			}
			for (i = 0; i < n - 1; i++)
			{
				miu[0][i] = h[0][i] / (h[0][i] + h[0][i + 1]);
				lamda[0][i] = 1 - miu[0][i];
			}
			for (i = 1; i < n; i++)
			{
				d[i][0] = 6 / (h[0][i - 1] + h[0][i]) * ((y[0][i + 1] - y[0][i]) / h[0][i] - (y[0][i] - y[0][i - 1]) / h[0][i - 1]);
			}
			if (bond == 1)
			{
				M[0][0] = left_bond;
				M[n][0] = right_bond;
				Mat A(1, n - 2);
				Mat B(1, n - 1);
				Mat C(1, n - 2);
				Mat D(n - 1, 1);
				for (i = 0; i < n - 2; i++)
				{
					A[0][i] = miu[0][i + 1];
				}
				for (i = 0; i < n - 1; i++)
				{
					B[0][i] = 2;
				}
				for (i = 0; i < n - 2; i++)
				{
					C[0][i] = lamda[0][i];
				}
				D[0][0] = d[1][0] - miu[0][0] * M[0][0];
				for (i = 1; i < n - 2; i++)
				{
					D[i][0] = d[i + 1][0];
				}
				D[n - 2][0] = d[n - 1][0] - lamda[0][n - 2] * M[n][0];
				Mat MM = Tridiagonal(A, B, C, D);
				for (i = 1; i < n; i++)
				{
					M[i][0] = MM[i - 1][0];
				}
			}
			if (bond == 2)
			{
				d[0][0] = 6 / h[0][0] * ((y[0][1] - y[0][0]) / h[0][0] - left_bond);
				d[n][0] = 6 / h[0][n - 1] * (right_bond - (y[0][n] - y[0][n - 1]) / h[0][n - 1]);
				Mat A(1, n);
				Mat B(1, n + 1);
				Mat C(1, n);
				for (i = 0; i < n - 1; i++)
				{
					A[0][i] = miu[0][i];
				}
				A[0][n - 1] = 1;
				for (i = 0; i < n + 1; i++)
				{
					B[0][i] = 2;
				}
				C[0][0] = 1;
				for (i = 1; i < n; i++)
				{
					C[0][i] = lamda[0][i - 1];
				}
				M = Tridiagonal(A, B, C, d);
			}
			if (bond == 3)
			{
				miu[0][n - 1] = h[0][n - 1] / (h[0][n - 1] + h[0][0]);
				lamda[0][n - 1] = 1 - miu[0][n - 1];
				d[n][0] = 6 / (h[0][n - 1] + h[0][0]) * ((y[0][1] - y[0][n]) / h[0][0] - (y[0][n] - y[0][n - 1]) / h[0][n - 1]);
				Mat B(1, n);
				Mat D(n, 1);
				Mat MM(n, 1);
				for (i = 0; i < n; i++)
				{
					B[0][i] = 2;
				}
				for (i = 0; i < n; i++)
				{
					D[i][0] = d[i + 1][0];
				}
				MM = Tridiagonal_Circle(miu, B, lamda, D);
				M[0][0] = MM[n - 1][0];
				for (i = 1; i < n + 1; i++)
				{
					M[i][0] = MM[i - 1][0];
				}
			}
			for (i = 0; i < n; i++)
			{
				S[i][0] = (M[i + 1][0] - M[i][0]) / (6 * h[0][i]);
				S[i][1] = (M[i][0] * x[0][i + 1] - M[i + 1][0] * x[0][i]) / (2 * h[0][i]);
				S[i][2] = (M[i + 1][0] * x[0][i] * x[0][i] - M[i][0] * x[0][i + 1] * x[0][i + 1]) / 2;
				S[i][2] += ((y[0][i + 1] - h[0][i] * h[0][i] * M[i + 1][0] / 6) - (y[0][i] - h[0][i] * h[0][i] * M[i][0] / 6));
				S[i][2] /= h[0][i];
				S[i][3] = (M[i][0] * x[0][i + 1] * x[0][i + 1] * x[0][i + 1] - M[i + 1][0] * x[0][i] * x[0][i] * x[0][i]) / 6;
				S[i][3] += ((y[0][i] - h[0][i] * h[0][i] * M[i][0] / 6) * x[0][i + 1] - (y[0][i + 1] - h[0][i] * h[0][i] * M[i + 1][0] / 6) * x[0][i]);
				S[i][3] /= h[0][i];
			}
			return S;
		}
	}
}