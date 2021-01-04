#include <iostream>
#include "matrix.h"
#include "funcinterpolation.h"

int main()
{
	using namespace matrix;
	using namespace function;
	using namespace std;
	Mat x(1,5),y(1,5);
	double numx[5] = { 0.25,0.30,0.39,0.45,0.50 };
	double numy[5] = { 0.5000,0.5477,0.6245,0.6708,0.5000 };
	//double numC[5] = { 1,2,1,2,1 };
	//double numb[5] = { 5,9,2,19,-4 };
	//double numa[25] = { 1,1,0,0,1,2,3,2,0,0,0,-3,4,1,0,0,0,4,7,2,1,0,0,-5,6 };
	//double numx[3] = { 0, 0, 0 };
	x = numx;
	y = numy;
	cout << Cubic_Spline(x,y,3) << endl;
	return 0;
}