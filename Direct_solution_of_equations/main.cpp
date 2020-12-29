#include <iostream>
#include "matrix.h"
#include "lineqaiterative.h"

int main()
{
	using namespace matrix;
	using namespace lineqa;
	using namespace std;
	Mat A(3, 3), b(3,1), x(3, 1);
	double num1[9] = { 2,0,1,0,1,0,1,0,2 };
	double num2[3] = { 3,1,3 };
	double num3[3] = { 0, 0, 0 };
	A = num1;
	b = num2;
	x = num3;
	Conjugate_Gradient(A, b, 0.000001, x);
	cout << x << endl;
	return 0;
}