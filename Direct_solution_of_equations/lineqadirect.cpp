//This file contains some function to solve the linear equations directly.
//Author: Chen Jinlong
//Time: 13:47 2020/9/26
#include <iostream>
#include <vector>
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


	/************************************************************
	结果采用LU分解的紧凑格式存储
	The results is stored in compact format with LU decomposition
	*************************************************************/
	Mat Doolittle(const Mat& A)
	{
		if (A.row() == 0 || A.row() != A.column())
		{
			throw "Fail to solve equations: Parameter error.";
			exit(1);
		}
		else
		{
			int n = A.row();
			Mat G(n, n);
			int i,j, k;
			for (i = 0; i < n; i++)
			{
				for (j = i; j < n; j++)
				{
					double sum = 0;
					for (k = 0; k < i; k++)
					{
						sum += (G[i][k] * G[k][j]);
					}
					G[i][j] = A[i][j] - sum;
				}
				for (j = i + 1; j < n; j++)
				{
					double sum = 0;
					for (k = 0; k < i; k++)
					{
						sum += (G[j][k] * G[k][i]);
					}
					G[j][i] = (A[j][i] - sum) / G[i][i];
				}
			}
			return G;
		}
	}

	Mat Doolittle_Solve(const Mat& A, const Mat& b)
	{
		int n = A.row();
		int eqa_num = b.column();
		Mat G = Doolittle(A);
		if (A.row() != b.row())
		{
			throw "Fail to solve equations: Parameter error.";
			exit(1);
		}
		else
		{
			Mat x(n, eqa_num), y(n, eqa_num);
			int i, j, k;
			for (k = 0; k < eqa_num; k++)
			{
				for (i = 0; i < n; i++)
				{
					y[i][k] = b[i][k];
					for (j = 0; j < i; j++)
					{
						y[i][k] -= (G[i][j] * y[j][k]);
					}
				}
				for (i = n; i > 0; i--)
				{
					double s = y[i - 1][k];
					for (j = i; j < n; j++)
					{
						s = s - G[i - 1][j] * x[j][k];
					}
					x[i - 1][k] = s / G[i - 1][i - 1];
				}
			}		
			return x;
		}
	}

	Mat Cholesky(const Mat& A)
	{
		if (A.row() == 0 || A.row() != A.column() || (!A.symmetric()))
		{
			throw "Fail to solve equations: Parameter error.";
			exit(1);
		}
		else
		{
			int n = A.row();
			Mat G(n, n);
			int i, j, k;
			for (j = 0; j < n; j++)
			{
				double sum = 0;
				for (k = 0; k < j; k++)
					sum += (G[j][k] * G[j][k]);
				G[j][j] = sqrt(A[j][j] - sum);
				for (i = j + 1; i < n; i++)
				{
					double sum = 0;
					for (k = 0; k < j; k++)
						sum += (G[i][k] * G[j][k]);
					G[i][j] = (A[i][j] - sum) / G[j][j];
				}
			}
			return G;
		}
	}

	Mat Cholesky_Solve(const Mat& A, const Mat& b)
	{
		int n = A.row();
		int eqa_num = b.column();
		Mat G = Cholesky(A), GT;
		GT = G.T();
		if (A.row() != b.row())
		{
			throw "Fail to solve equations: Parameter error.";
			exit(1);
		}
		else
		{
			Mat x(n, eqa_num), y(n, eqa_num);
			int i, j, k;
			for (k = 0; k < eqa_num; k++)
			{
				for (i = 0; i < n; i++)
				{
					double sum = 0;
					for (j = 0; j < i; j++)
					{
						sum += (G[i][j] * y[j][k]);
					}
					y[i][k] = (b[i][k] - sum) / G[i][i];
				}
				for (i = n; i > 0; i--)
				{
					double s = y[i - 1][k];
					for (j = i; j < n; j++)
					{
						s = s - GT[i - 1][j] * x[j][k];
					}
					x[i - 1][k] = s / GT[i - 1][i - 1];
				}
			}
			return x;
		}
	}

	Mat LDL(const Mat& A)
	{
		if (A.row() == 0 || A.row() != A.column() || (!A.symmetric()))
		{
			throw "Fail to solve equations: Parameter error.";
			exit(1);
		}
		else
		{
			int n = A.row();
			Mat G(n, n);
			int i, j, k;
			for (j = 0; j < n; j++)
			{
				double sum = 0;
				for (k = 0; k < j; k++)
					sum += (G[j][k] * G[j][k] * G[k][k]);
				G[j][j] = A[j][j] - sum;
				for (i = j + 1; i < n; i++)
				{
					double sum = 0;
					for (k = 0; k < j; k++)
						sum += (G[i][k] * G[j][k] * G[k][k]);
					G[i][j] = (A[i][j] - sum) / G[j][j];
				}
			}
			return G;
		}
	}

	Mat LDL_Solve(const Mat& A, const Mat& b)
	{
		int n = A.row();
		int eqa_num = b.column();
		Mat G = LDL(A);
		if (A.row() != b.row())
		{
			throw "Fail to solve equations: Parameter error.";
			exit(1);
		}
		else
		{
			Mat x(n, eqa_num), y(n, eqa_num);
			int i, j, k;
			for (k = 0; k < eqa_num; k++)
			{
				for (i = 0; i < n; i++)
				{
					double sum = 0;
					for (j = 0; j < i; j++)
					{
						sum += (G[i][j] * y[j][k]);
					}
					y[i][k] = b[i][k] - sum;
				}
				for (i = n; i > 0; i--)
				{
					//double s = y[i - 1][k] / G[i - 1][i - 1];
					double sum = 0;
					for (j = i; j < n; j++)
					{
						sum += (G[j][i - 1] * x[j][k]);
					}
					x[i - 1][k] = y[i - 1][k] / G[i - 1][i - 1] - sum;
				}
			}
			return x;
		}
	}

	/************************************
	G上三角部分存储R矩阵的上三角部分
	G的主对角线及下三角部分存储Householder变换矩阵中间值w
	d存储R矩阵的主对角线部分
	alpha存储Householder变换的每一步的数量sigma(sigma-akk)
	************************************/
	void QR_Householder(const Mat& A, Mat& G, Mat& d, Mat& alpha)
	{
		if (A.row() == 0 || A.column() == 0 || A.row() < A.column())
		{
			throw "Fail to solve equations: Parameter error.";
			exit(1);
		}
		else
		{
			int m = A.row();
			int n = A.column();
			G = A;
			d = Mat(n, 1);
			alpha = Mat(n, 1);
			int i, j, k;
			double sigma;
			for (k = 0; k < n - 1; k++)
			{
				double sum = 0;
				for (i = k; i < m; i++)
					sum += (G[i][k] * G[i][k]);
				sigma = sqrt(sum);
				if (G[k][k] >= 0)
					sigma = -sigma;
				d[k][0] = sigma;
				alpha[k][0] = sigma * (sigma - G[k][k]);
				G[k][k] = G[k][k] - sigma;			
				for (j = k + 1; j < n; j++)
				{
					double beta;
					double sum = 0;
					for (i = k + 1; i < m; i++)
						sum += G[i][k] * G[i][j];
					beta = (G[k][k] * G[k][j] + sum) / alpha[k][0];
					G[k][j] -= (beta * G[k][k]);
					for (i = k + 1; i < m; i++)
						G[i][j] -= (beta * G[i][k]);
				}
			}
			if (m == n)
				d[n - 1][0] = G[n - 1][n - 1];
			else
			{
				double sum = 0;
				for (i = n - 1; i < m; i++)
					sum += (G[i][n - 1] * G[i][n - 1]);
				sigma = sqrt(sum);
				if (G[n - 1][n - 1] >= 0)
					sigma = -sigma;
				d[n - 1][0] = sigma;
				alpha[n - 1][0] = sigma * (sigma - G[n - 1][n - 1]);
				G[n - 1][n - 1] = G[n - 1][n - 1] - sigma;
			}
		}		
	}

	void QR_Householder_Show(const Mat& A)
	{
		using std::cout;
		using std::endl;
		Mat G, d, alpha;
		QR_Householder(A, G, d, alpha);
		int m = G.row();
		int n = G.column();
		Mat H(m, m);
		Mat R(m, n);
		int i, j, k;
		Mat w(m, 1);
		for (j = 0; j < m; j++)
			w[j][0] = G[j][0];
		H = I(m) - 1 / alpha[0][0] * w * w.T();
		for (i = 1; i < n - 1; i++)
		{
			Mat w(m - i, 1);
			Mat Hi(m - i, m - i);
			for (j = i; j < m; j++)
				w[j - i][0] = G[j][i];
			Hi = I(m - i) - 1 / alpha[i][0] * w * w.T();
			Mat Ht = I(m);
			for (j = i; j < m; j++)
				for (k = i; k < n; k++)
					Ht[j][k] = Hi[j - i][k - i];
			H = Ht * H;
		}
		if (m > n)
		{
			Mat w(m - n + 1, 1);
			Mat Hn(m - n + 1, m - n + 1);
			for (j = n - 1; j < m; j++)
				w[j - n + 1][0] = G[j][n - 1];
			Hn = I(m - n + 1) - 1 / alpha[n - 1][0] * w * w.T();
			Mat Ht = I(m);
			for (j = n - 1; j < m; j++)
				for (k = n - 1; k < n; k++)
					Ht[j][k] = Hn[j - n + 1][k - n + 1];
			H = Ht * H;
		}
		Mat Q = H.T();
		for (i = 0; i < n; i++)
			R[i][i] = d[i][0];
		for (i = 0; i < n; i++)
			for (j = i + 1; j < m; j++)
				R[i][j] = G[i][j];
		cout << endl;
		cout << "Q=" << endl;
		cout << Q;
		cout << "R=" << endl;
		cout << R;
	}

	void QR_Householder_Show(const Mat& G, const Mat& d, const Mat& alpha)
	{
		using std::cout;
		using std::endl;
		if (G.row() < G.column() || (G.column() != d.row()) || (G.column() != alpha.row()))
		{
			throw "Parameter error.";
			exit(1);
		}
		else
		{
			int m = G.row();
			int n = G.column();
			Mat H(m, m);
			Mat R(m, n);
			int i, j, k;
			Mat w(m, 1);
			for (j = 0; j < m; j++)
				w[j][0] = G[j][0];
			H = I(m) - 1 / alpha[0][0] * w * w.T();
			for (i = 1; i < n - 1; i++)
			{
				Mat w(m - i, 1);
				Mat Hi(m - i, m - i);
				for (j = i; j < m; j++)
					w[j - i][0] = G[j][i];
				Hi = I(m - i) - 1 / alpha[i][0] * w * w.T();
				Mat Ht = I(m);
				for (j = i; j < m; j++)
					for (k = i; k < n; k++)
						Ht[j][k] = Hi[j - i][k - i];
				H = Ht * H;
			}
			if (m > n)
			{
				Mat w(m - n + 1, 1);
				Mat Hn(m - n + 1, m - n + 1);
				for (j = n - 1; j < m; j++)
					w[j - n + 1][0] = G[j][n - 1];
				Hn = I(m - n + 1) - 1 / alpha[n - 1][0] * w * w.T();
				Mat Ht = I(m);
				for (j = n - 1; j < m; j++)
					for (k = n - 1; k < n; k++)
						Ht[j][k] = Hn[j - n + 1][k - n + 1];
				H = Ht * H;
			}
			Mat Q = H.T();
			for (i = 0; i < n; i++)
				R[i][i] = d[i][0];
			for (i = 0; i < n; i++)
				for (j = i + 1; j < m; j++)
					R[i][j] = G[i][j];
			cout << endl;
			cout << "Q=" << endl;
			cout << Q;
			cout << "R=" << endl;
			cout << R;
		}
	}

	Mat QR_Householder_Solve(const Mat& A, const Mat& b)
	{
		if (A.row() != b.row() || A.row() == 0 || A.column() == 0 || A.row() < A.column() || b.column() != 1)
		{
			throw "Fail to solve equations: Parameter error.";
			exit(1);
		}
		else
		{
			Mat G, d, alpha;
			int m = A.row();
			int n = A.column();
			G = Mat(A, b);
			d = Mat(n, 1);
			alpha = Mat(n, 1);
			int i, j, k;
			double sigma;
			for (k = 0; k < n - 1; k++)
			{
				double sum = 0;
				for (i = k; i < m; i++)
					sum += (G[i][k] * G[i][k]);
				sigma = sqrt(sum);
				if (G[k][k] >= 0)
					sigma = -sigma;
				d[k][0] = sigma;
				alpha[k][0] = sigma * (sigma - G[k][k]);
				G[k][k] = G[k][k] - sigma;
				for (j = k + 1; j < n + 1; j++)
				{
					double beta;
					double sum = 0;
					for (i = k + 1; i < m; i++)
						sum += G[i][k] * G[i][j];
					beta = (G[k][k] * G[k][j] + sum) / alpha[k][0];
					G[k][j] -= (beta * G[k][k]);
					for (i = k + 1; i < m; i++)
						G[i][j] -= (beta * G[i][k]);
				}
			}
			if (m == n)
				d[n - 1][0] = G[n - 1][n - 1];
			else
			{
				double sum = 0;
				for (i = n - 1; i < m; i++)
					sum += (G[i][n - 1] * G[i][n - 1]);
				sigma = sqrt(sum);
				if (G[n - 1][n - 1] >= 0)
					sigma = -sigma;
				d[n - 1][0] = sigma;
				alpha[n - 1][0] = sigma * (sigma - G[n - 1][n - 1]);
				G[n - 1][n - 1] = G[n - 1][n - 1] - sigma;
				double beta;
				sum = 0;
				for (i = n; i < m; i++)
					sum += G[i][n - 1] * G[i][n];
				beta = (G[n - 1][n - 1] * G[n - 1][n] + sum) / alpha[n - 1][0];
				G[n - 1][n] -= (beta * G[n - 1][n - 1]);
				for (i = n; i < m; i++)
					G[i][n] -= (beta * G[i][n - 1]);
			}
			Mat x(n, 1);
			for (k = n; k > 0; k--)
			{
				double sum = 0;
				for (j = k; j < n; j++)
				{
					sum += (G[k - 1][j] * x[j][0]);
				}
				x[k - 1][0] = (G[k - 1][n] - sum) / d[k - 1][0];
			}
			return x;
		}
	}
}