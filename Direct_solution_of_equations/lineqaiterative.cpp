// This file contains some function to solve the linear equations by iterative method.
//Author: Chen Jinlong
//Time: 21:39 2020/10/16

#include <iostream>
#include "lineqaiterative.h"

using namespace matrix;

namespace lineqa
{
	void Jacobi(const Mat& A, const Mat& b, const double& error, Mat& x)
	{
		if (A.row() == 0 || A.row() != A.column() || b.column() != 1 || x.column() != 1 || A.row() != b.row() || A.row() != x.row())
		{
			throw "Fail to solve equations: Parameter error.";
			exit(1);
		}
		else
		{
			Mat B = A, g = b;
			int i, j;
			int n = A.row();
			for (i = 0; i < n; i++)
			{
				for (j = 0; j < n; j++)
					B[i][j] = -B[i][j] / A[i][i];
				B[i][i] = 0;
				g[i][0] = g[i][0] / A[i][i];
			}
			Mat xx;
			int k = 0;
			while (true)
			{
				xx = x;
				x = B * x + g;
				if ((x - xx).norm_inf() < error)
					break;
				k++;
				if (k == INT_MAX)
				{
					std::cout << "Defocusing iteration!" << std::endl;
					break;
				}
			}
		}
	}

	void GaussSeidel(const Mat& A, const Mat& b, const double& error, Mat& x)
	{
		if (A.row() == 0 || A.row() != A.column() || b.column() != 1 || x.column() != 1 || A.row() != b.row() || A.row() != x.row())
		{
			throw "Fail to solve equations: Parameter error.";
			exit(1);
		}
		else
		{
			int i, j;
			int n = A.row();
			Mat xx;
			int k = 0;
			while (true)
			{
				xx = x;
				for (i = 0; i < n; i++)
				{
					double sum = 0;
					for (j = 0; j < n; j++)
					{
						if (i != j)
							sum += (A[i][j] * x[j][0]);
					}
					x[i][0] = (b[i][0] - sum) / A[i][i];
				}
				if ((x - xx).norm_inf() < error)
					break;
				k++;
				if (k == INT_MAX)
				{
					std::cout << "Defocusing iteration!" << std::endl;
					break;
				}
			}
		}
	}

	void SOR(const Mat& A, const Mat& b, const double& omega, const double& error, Mat& x)
	{
		if (A.row() == 0 || A.row() != A.column() || b.column() != 1 || x.column() != 1 || A.row() != b.row() || A.row() != x.row())
		{
			throw "Fail to solve equations: Parameter error.";
			exit(1);
		}
		else
		{
			int i, j;
			int n = A.row();
			Mat xx;
			int k = 0;
			while (true)
			{
				xx = x;
				for (i = 0; i < n; i++)
				{
					double sum = 0;
					for (j = 0; j < n; j++)
					{
							sum += (A[i][j] * x[j][0]);
					}
					x[i][0] += omega * (b[i][0] - sum) / A[i][i];
				}
				if ((x - xx).norm_inf() < error)
					break;
				k++;
				if (k == INT_MAX)
				{
					std::cout << "Defocusing iteration!" << std::endl;
					break;
				}
			}
		}
	}
}