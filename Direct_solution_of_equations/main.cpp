#include <iostream>
#include "matrix.h"
#include "lineqadirect.h"

int main()
{
	using namespace matrix;
	using namespace lineqa;
	using namespace std;
	Mat A(1, 5), B(1, 5), C(1, 5),b(5,1), x(5, 1),a(5,5);
	double numA[5] = { 1,2,-3,4,-5 };
	double numB[5] = { 1,3,4,7,6 };
	double numC[5] = { 1,2,1,2,1 };
	double numb[5] = { 5,9,2,19,-4 };
	double numa[25] = { 1,1,0,0,1,2,3,2,0,0,0,-3,4,1,0,0,0,4,7,2,1,0,0,-5,6 };
	//double numx[3] = { 0, 0, 0 };
	A = numA;
	B = numB;
	C = numC;
	b = numb;
	a = numa;
	//x = numx;
	x = Tridiagonal_Circle(A, B, C, b);
	cout << x << endl;
	cout << Doolittle(a) << endl;
	cout << Doolittle_Solve(a,b) << endl;
	return 0;
}