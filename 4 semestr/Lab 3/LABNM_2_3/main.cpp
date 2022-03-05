#include <iostream>
#include <iomanip>
#include <math.h>
using namespace std;
#define A -2.9
#define B 0.4

double f(double x)
{
	return pow(x,5)-6.2*pow(x,3)+3.5*pow(x,2)-7*x-2.1;
}
double CavalieriSimpsonFormula(double a, double b, double eps)
{
	double I_prev = eps + 1, I_cur = 0;//I-предыдущее вычисленное значение интеграла, I1-новое, с большим N.
	for (int N = 2; (N <= 4) || (fabs(I_cur - I_prev) > eps); N *= 2)
	{
		double h, sum2 = 0, sum4 = 0, sum = 0;
		h = (b - a) / (2 * N);//Шаг интегрирования.
		for (int i = 1; i <= 2 * N - 1; i += 2)
		{
			sum4 += f(a + h * i);//Значения с нечётными индексами, которые нужно умножить на 4.
			sum2 += f(a + h * (i + 1));//Значения с чётными индексами, которые нужно умножить на 2.
		}
		sum = f(a) + 4 * sum4 + 2 * sum2 - f(b);//Отнимаем значение f(b) так как ранее прибавили его дважды. 
		I_prev = I_cur;
		I_cur = (h / 3) * sum;
	}
	return I_cur;
}
int main()
{
	for (double eps = 0.1; eps >= pow(10, -8); eps /= 10)
	{
		double Integral = CavalieriSimpsonFormula(A, B, eps);
		if (eps != 0.1)
		{
			cout << "\n\nError: " << eps;
			cout << "\nIntegral: " << setprecision(17) << Integral;
		}
		else
		{
			cout << "Error: " << eps;
			cout << "\nIntegral: " << setprecision(17) << Integral;
		}
	}
	
	return 0;
}