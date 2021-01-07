//This file contains some function to solve function intergral.
//Author: Chen Jinlong
//Time: 19:24 2021/1/7

#include "funcintegral.h"
#include <cmath>

namespace function
{
	double Romberg(double (*func)(double), double a, double b, double error)
	{
		double T[2], S[2], C[2], R[2];
		//k=0
		T[0] = (b - a) / 2 * (func(a) + func(b));
		//k=1
		T[1] = T[0] / 2 + (b - a) / 2 * func(a + (b - a) / 2);
		S[0] = T[1] + (T[1] - T[0]) / 3;
		//k=2
		T[0] = T[1];
		T[1] = T[0] / 2 + (b - a) / 4 * (func(a + (b - a) / 4) + func(a + 3 * (b - a) / 4));
		S[1] = T[1] + (T[1] - T[0]) / 3;
		C[0] = S[1] + (S[1] - S[0]) / 15;
		//k=3
		T[0] = T[1];
		S[0] = S[1];
		T[1] = T[0] / 2 + (b - a) / 8 * (func(a + (b - a) / 8) + func(a + 3 * (b - a) / 8) + func(a + 5 * (b - a) / 8) + func(a + 7 * (b - a) / 8));
		S[1] = T[1] + (T[1] - T[0]) / 3;
		C[1] = S[1] + (S[1] - S[0]) / 15;
		R[0] = C[1] + (C[1] - C[0]) / 63;
		int pow2n = 16;
		//k=4及以后（程序中为0起）
		int k;
		for (k = 0; k < INT_MAX; k++)
		{
			T[0] = T[1];
			S[0] = S[1];
			C[0] = C[1];
			double sum = 0;
			int i;
			for (i = 1; i <= (pow2n / 2); i++)
				sum += func(a + (2 * i - 1) * (b - a) / pow2n);
			T[1] = T[0] / 2 + (b - a) / pow2n * sum;
			S[1] = T[1] + (T[1] - T[0]) / 3;
			C[1] = S[1] + (S[1] - S[0]) / 15;
			R[1] = C[1] + (C[1] - C[0]) / 63;
			if (fabs(R[1] - R[0]) < error)
				return R[1];
			pow2n *= 2;
			R[0] = R[1];
		}
		return R[1];
	}
}