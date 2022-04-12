#include <iostream>
#include <math.h>
#include <cmath>
#include <vector>
#include <iomanip>
#include <fstream>
using namespace std;
# define M_PI 3.14159265358979323846
#define A M_PI / 2
#define B 2 * M_PI
#define F(A) M_PI / 2
#define N 30
#define P 3
double f(double x, double y);
double Answer(double x);
vector<double> createSteadyGrid(double a, double b, int n);
double Runge_Kutta_method(int m, double h, double x, double y);
void PrintVector(vector<double> vec);
void PrintVectorToFile(string filename, vector<double> vec);
double Euler(int m, double h, double x, double y);
void Cycle(double y, vector<double> x, double eps, string filename_error, string filename_iter, int wtf);
void TestRunge_Kutta_method(double h, string filename_runge, string filename_stepLength, string filename_x);
int main(void)
{
	vector<double> x = createSteadyGrid(A, B, N);
	vector<double> answ;
	for (int i = 0; i < N; i++)
	{
		answ.push_back(Answer(x[i]));
	}
	string filename_x = "files/f_x.txt";
	string filename_y = "files/f_y.txt";
	PrintVectorToFile(filename_x, x);
	PrintVectorToFile(filename_y, answ);
	string filename_eps = "files/f_eps.txt";
	vector<double> epsilon;
	string filename_iter = "files/f_iter.txt";
	string filename_error = "files/f_error.txt";
	for (double eps = 0.001; eps > 0.00000001; eps /= 10) {
		epsilon.push_back(eps);
		printf("\neps=%e\n", eps);
		double y = F(A);
		//Cycle(y, x, eps, filename_error, filename_iter, 0);
	}
	PrintVectorToFile(filename_eps,epsilon);
	return 0;
}

void Cycle(double y, vector<double> x, double eps, string filename_error, string filename_iter, int wtf) {
	vector<double> iter;
	vector<double> error;
	for (int i = 0; i < N - 1; i++) {
		double h = (B - A) / (N - 1);
		int m = 1;
		double prev;
		double next = Runge_Kutta_method(m, h, x[i], y);
		do {
			h /= 2;
			m *= 2;
			prev = next;
			next = Runge_Kutta_method(m, h, x[i], y);
		} while (fabs((next - prev)) / (2 ^ P - 1) >= eps);

		y = next;
		if (i == 0 && wtf == 0) {
			iter.push_back((int)log2(m));
			
		}
		if (wtf == 0) {
			printf("\nx=%.15lf\ny=%.15lf\nnext=%.15lf\n", x[i + 1], Answer(x[i + 1]), next);
		}
		printf("error=%.15lf\n", fabs(next - Answer(x[i + 1])));
		error.push_back(fabs(next - Answer(x[i + 1])));
	}
	PrintVectorToFile(filename_iter, iter);
	PrintVectorToFile(filename_error, error);
}


void PrintVector(vector<double> vec)
{
	for (auto& v : vec)
	{
		cout << setprecision(17) << v;
		cout << ' ';
	}
}

double Runge_Kutta_method(int m, double h, double x, double y)
{
	double y_prev = y;
	double k1, k2, k3;
	//double x = A;
	//double y = F(A)
	for (int i = 0; i < m; i++)
	{
		//x = a + i * h;
		k1 = f(x, y);
		k2 = f(x + h / 3, y + h * k1 / 3);
		k3 = f(x + 2 * h / 3, y + 2 * h * k2 / 3);
		y_prev = y;
		y = y_prev + (h / 4)*(k1+3*k3);
		x = x + h;
	}
	return y;
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
		file << setprecision(17) << v;
		file << ' ';
	}
	file.close();
}

double Answer(double x)
{
	return x * sin(x);
}

double f(double x, double y)
{
	return y / x + x * cos(x);
}

vector<double> createSteadyGrid(double a, double b, int n) {
	vector<double> x;

	for (int i = 0; i < n; i++) {
		x.push_back(a + (b - a) * i / (n-1));
	}
	return x;
}

double Euler(int m, double h, double x, double y) {
	double halfY = y;
	for (int i = 0; i < m; i++) {
		halfY = y + h / 2 * f(x, y);
		y = y + h * f(x + h / 2, halfY);
		x += h;
	}

	return y;
}

void TestRunge_Kutta_method(double h, string filename_runge, string filename_stepLength, string filename_x)
{
	double y = F(A);
	int i = 0;
	vector<double> runge;
	vector<double> stepLength;
	vector<double>x_h;
	for (double x = A; x <= B; x += h, i++) {
		y = Runge_Kutta_method(1, h, x, y);
		x_h.push_back(x + h);
		runge.push_back(y);
	}
	PrintVectorToFile(filename_x, x_h);
	PrintVectorToFile(filename_runge, runge);
	stepLength.push_back(i);
	PrintVectorToFile(filename_stepLength, stepLength);
}