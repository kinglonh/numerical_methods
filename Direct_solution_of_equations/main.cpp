#include <iostream>
#include "matrix.h"
#include "lineqadirect.h"

int main()
{
	using namespace matrix;
	using namespace lineqa;
	using namespace std;
	Mat A(20, 20);
	//double num1[6] = { 1, 4, 2, 5, 3, 6 };
	//double num2[3] = { 7, 8, 9 };
	//A = num1;
	//b = num2;
	int i, j;
	for(i = 0; i < 20; i++)
		for (j = 0; j < 20; j++)
		{
			if (i == j)
				A[i][j] = i + 1;
			else
				A[i][j] = min(i, j) + 1;
		}
	//Mat x;
	//x = QR_Householder_Solve(A, b);
	//x = Doolittle(A);
	//QR_Householder_Show(A);
	cout << LDL(A) << endl;
	return 0;
}