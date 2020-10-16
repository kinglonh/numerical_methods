#include <iostream>
#include "matrix.h"
#include "lineqaiterative.h"

int main()
{
	using namespace matrix;
	using namespace lineqa;
	using namespace std;
	Mat A(3, 3), b(3,1), x(3, 1);
	double num1[9] = { 10, -1, -2, -1, 10, -2, -1, -1, 5 };
	double num2[3] = { 72, 83, 42 };
	double num3[3] = { 0, 0, 0 };
	A = num1;
	b = num2;
	x = num3;
	Jacobi(A, b, 0.0001, x);
	cout << x << endl;
	return 0;
}