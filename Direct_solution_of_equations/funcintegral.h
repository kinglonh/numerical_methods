//This file contains some function to solve function intergral.
//Author: Chen Jinlong
//Time: 19:02 2021/1/7

#ifndef FUNCINTEGRAL
#define FUNCINTEGRAL

namespace function
{
	//Romberg���ַ�,funcΪ��������,aΪ�����½�,bΪ�����Ͻ�
	double Romberg(double (*func)(double), double a, double b, double error);
}
#endif // !FUNCINTEGRAL

