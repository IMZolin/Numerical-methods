#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <math.h>
using namespace std;
#define A -2.9
#define B 0.4
#define N 4
typedef struct struct_integral{
	vector<double> integral;
	vector<double> err;
	vector<int> number_nodes;
}struct_integral_t;
double f(double x);
double fIntegrated(double x);
double CavalieriSimpsonFormula(double a, double b, double eps);
void PrintVectorToFile(string filename, vector<double> vec);
struct_integral* FullSimpsonFormula(double a, double b, double eps);
double SimpsonFormula(double a, double b, double eps, int* n);
int main()
{
	vector<double> integral;
	vector<int> number_nodes;
	vector<double> epsilon;
	vector <double> error;
	vector<double> iterations;
	vector<double> f_integr;
	double t = 60.925837499999993;
	
	for (double eps = 0.1; eps >= pow(10, -8); eps /= 10)
	{
		int n = 0;
		integral.push_back(SimpsonFormula(A, B, eps, &n));
		epsilon.push_back(eps);
		f_integr.push_back(fIntegrated(B) - fIntegrated(A));
		error.push_back(fabs(integral.back() - f_integr.back()));
		iterations.push_back((int)log(n));
		/*integral.push_back(CavalieriSimpsonFormula(A, B, eps));
		epsilon.push_back(eps);
		f_integr.push_back(fIntegrated(B) - fIntegrated(A));
		error.push_back(fabs(integral.back() - f_integr.back()));*/

		if (eps != 0.1)
		{
			cout << "\n\nEpsilon: " << eps;
			cout << "\nIntegral: " << setprecision(17) << integral.back();
			cout << "\nExact integral: " << setprecision(17) << f_integr.back();
			cout << "\nError: " << setprecision(17) << error.back();
			cout << "\nIterations: " << iterations.back();
		}
		else
		{
			cout << "Epsilon: " << eps;
			cout << "\nIntegral: " << setprecision(17) << integral.back();
			cout << "\nExact integral: " << setprecision(17) << f_integr.back();
			cout << "\nError: " << setprecision(17) << error.back();
			cout << "\nIterations: " << iterations.back();
		}
	}

	string filename = "f_integral.txt";
	PrintVectorToFile(filename, integral);
	filename = "f_epsilon.txt";
	PrintVectorToFile(filename, epsilon);
	filename = "f_error.txt";
	PrintVectorToFile(filename, error);
	filename = "f_iterrations.txt";
	PrintVectorToFile(filename, iterations);
	return 0;
}

double f(double x)
{
	return pow(x,5)-6.2*pow(x,3)+3.5*pow(x,2)-7*x-2.1;
}

double fIntegrated(double x) {
	return pow(x, 6) / 6 - 31 * pow(x, 4) / 20 + 7 * pow(x, 3) / 6 - 3.5 * pow(x, 2) - 2.1 * x;
}
double SimpsonFormula(double a, double b,double eps, int* n)
{
	double I_prev = eps + 1, I_cur = 0;//I_prev-предыдущее вычисленное значение интеграла, I_cur-новое, с большим N.
	for (*n = 1; (*n <= N) || (fabs(I_cur - I_prev) > eps); *n *= 2)
	{
		double h, sum2 = 0, sum4 = 0, sum = 0;
		h = (b - a) / (2 * *n);//Шаг интегрирования.
		for (int i = 1; i <= 2 * *n - 1; i += 2)
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
double CavalieriSimpsonFormula(double a, double b, double eps)
{
	double I_prev = eps + 1, I_cur = 0;//I_prev-предыдущее вычисленное значение интеграла, I_cur-новое, с большим N.
	for (int n = 1; (n <= N) || (fabs(I_cur - I_prev) > eps); n *= 2)
	{
		double h, sum2 = 0, sum4 = 0, sum = 0;
		h = (b - a) / (2 * n);//Шаг интегрирования.
		for (int i = 1; i <= 2 * n - 1; i += 2)
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

struct_integral* FullSimpsonFormula(double a, double b, double eps)
{
	vector<double> err;
	vector<int>number_nodes;
	vector<double> integral;
	struct_integral_t* integr = new struct_integral_t[1];
	integr->err = err;
	integr->number_nodes = number_nodes;
	integr->integral = integral;
	int counter = 0;
	double I = CavalieriSimpsonFormula(a, b, eps);
	double I_prev = CavalieriSimpsonFormula(a, b, eps);
	counter++;
	integral.push_back(I);
	err.push_back((abs(I - I_prev) / (2 * N - 1)));
	while ((abs(I - I_prev) / (2 *N - 1)) > eps)
	{
		err.push_back((abs(I - I_prev) / (2 * N - 1)));
		I_prev = I;
		counter++;
		I = CavalieriSimpsonFormula(a, b, eps);
		integral.push_back(I);
		number_nodes.push_back(counter);
	}
	//return number_nodes;
	return integr;
}

void PrintVectorToFile(string filename, vector<double> vec)
{
	ofstream file;
	file.open(filename);
	if (!file.is_open())
	{
		cout << "Uh oh.File wasn't opened. Check if the name of file " << filename << " correct.\n" << endl;
		exit(1);
	}
	for (auto& v : vec)
	{
		file << setprecision(17)<<v;
		file << ' ';
	}
	file.close();
}

