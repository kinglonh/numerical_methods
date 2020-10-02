#include <iostream>
#include "matrix.h"
#include "lineqadirect.h"

int main()
{
	using namespace matrix;
	using namespace lineqa;
	using namespace std;
	Mat A(4, 4), b(4,2);
	double num1[16] = { 9, 18, 9, -27, 18, 45, 0, -45, 9, 0, 126, 9, -27, -45, 9, 135 };
	double num2[8] = { 1,1, 2,2, 16,3, 8,4 };
	A = num1;
	b = num2;
	Mat x;
	x = LDL_Solve(A, b);
	cout << x << endl;
	return 0;
}