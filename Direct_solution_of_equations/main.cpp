#include <iostream>
#include "matrix.h"
#include "lineqadirect.h"

int main()
{
	using namespace matrix;
	using namespace lineqa;
	using namespace std;
	Mat A(2, 2), b(2,1);
	double num1[4] = { 0.0003, 3, 1, 1 };
	double num2[2] = { 2.0001, 1 };
	A = num1;
	b = num2;
	Mat x;
	x = Gauss(A, b);
	cout << x << endl;
	return 0;
}