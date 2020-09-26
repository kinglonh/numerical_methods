#include <iostream>
#include "matrix.h"
#include "lineqadirect.h"

int main()
{
	using namespace matrix;
	using namespace lineqa;
	using namespace std;
	Mat A(3, 3), b(3,1);
	double num1[9] = { 1, 2, 1, 2, 2, 3, -1, -3 ,0 };
	double num2[3] = { 0, 3, 2 };
	A = num1;
	b = num2;
	Mat x;
	x = Gauss_MaxColumn(A, b);
	cout << x << endl;
	return 0;
}