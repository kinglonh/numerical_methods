//This file contains some function to solve function intergral.
//Author: Chen Jinlong
//Time: 19:02 2021/1/7

#ifndef FUNCINTEGRAL
#define FUNCINTEGRAL

namespace function
{
	//Romberg积分法,func为被积函数,a为积分下界,b为积分上界
	double Romberg(double (*func)(double), double a, double b, double error);
}
#endif // !FUNCINTEGRAL

