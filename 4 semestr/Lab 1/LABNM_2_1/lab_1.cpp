#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <math.h>
using namespace std;
# define M_PI 3.14159265358979323846
#define A 0.0
#define B 1.0
#define N 30
#define MAX_NODES 10
double f(double x);
double f5Der(double x);
vector<double> ChebyshevGrid(int n, double a, double b);
vector<double> UniformGrid(int k, double a, double b);
vector<double> GridFunc(vector<double> x_h, int k);
double LagrangePolynom(double x, vector<double> x_h, vector<double> y_h, int k);
void PrintVectorToFile(string filename, vector<double> vec);
vector<double> calcTheoreticError(vector<double> x, vector<double> x_h, int n);
vector<double> calcActualticError(vector<double> polynom, vector<double>func);
vector<double> GridForMaxError(vector<double> x_h, int k);
double MaxElementVector(vector<double> v);
void PrintVector(vector<double> vec);

int main()
{
	double a = A;
	double b = B;
	vector<double> polynom;
	vector<double> polynom_max;
	double pol = 0;
	double max_theoretic = 0;
	vector<double> x = ChebyshevGrid(N, a, b);
	vector<double> y = GridFunc(x, N);
	vector<double>teoretic_error;
	vector<double>actual_error;
	vector<double> max_error;
	string filename;
	filename = "f_x.txt";
	PrintVectorToFile(filename, x);
	filename = "f_y.txt";
	PrintVectorToFile(filename, y);

	for (int i = 6; i <= MAX_NODES; i+=2)
	{
		vector<double> x_h = ChebyshevGrid(i, a, b);
		//vector<double> x_h_max = GridForMaxError(x_h, i);
		filename = "f_x_h.txt";
		PrintVectorToFile(filename, x_h);
		vector<double> y_h = GridFunc(x_h, i);
		//vector<double> y_h_max = GridFunc(x_h_max, i);
		filename = "f_y_h.txt";
		
		PrintVectorToFile(filename, y_h);
		for (int j = 0; j <= N; j++)
		{
			pol = LagrangePolynom(x[j], x_h, y_h, i);
			polynom.push_back(pol);
			//pol = LagrangePolynom(x[j], x_h_max, y_h_max, i);
			//polynom_max.push_back(pol);
		}
		cout <<"\nNumb of nodes: " <<i;
		cout << "\nApproximation from polynom\n";
		PrintVector(polynom);
		cout << "\nActual error\n";
		actual_error = calcActualticError(polynom, y);
		PrintVector(actual_error);
		cout << "\nTheoretic error\n";
		teoretic_error = calcTheoreticError(x, x_h, i);
		max_theoretic = MaxElementVector(teoretic_error);
		PrintVector(teoretic_error);
		cout << "\nMAX"<<max_theoretic;

		//cout << "\nMax error\n";
		//max_error = calcMaxError(polynom_max, y_h_max);
		//PrintVector(max_error);
		cout << "\n";
		
		filename = "f_polynom.txt";
		PrintVectorToFile(filename, polynom);
		filename = "f_actual_error.txt";
		PrintVectorToFile(filename, actual_error);
		filename = "f_theoritic_error.txt";
		PrintVectorToFile(filename, teoretic_error);
		//filename = "f_max_error.txt";
		//PrintVectorToFile(filename, max_error);
	}
	return 0;
}

double f(double x) {
	return log10(x) + 7.0 / (2 * x + 6);
}

double f5Der(double x)
{
	//return 24 / log(10) * pow(x, 5) - 26880 / pow(2 * x + 6, 6);
	return 7.214;
}

vector<double> GridForMaxError(vector<double> x_h, int k) {
	vector<double> x_h_max;
	/*if (y_h.empty())
		cout << "Error allocation memory for vector";
		exit(1);
		f(x_h[i])*/
	for (int i = 0; i < k-1; i++) {
		x_h_max.push_back((x_h[i]+x_h[i+1])/2);
	}
	return x_h_max;
}

double MaxElementVector(vector<double> v)
{
	double max = v[0];
	for (int i = 1; i < N; i++)
	{
		if (max < v[i])
		{
			max = v[i];
		}
	}
	return max;
}

vector<double> calcActualticError(vector<double> polynom, vector<double>func)
{
	vector<double> error;
	for (int i = 0; i < N; i++)
	{
		error.push_back(polynom[i] - func[i]);
	}
	return error;
}

vector<double> calcTheoreticError(vector<double> x, vector<double> x_h, int n) {
	vector<double> res;
	/*if (res.empty()) {
		printf("Erorr allocating memory to theory\n");
		exit(1);
	}*/
	for (int i = 0; i < N; i++) {
		double root = 1;
		for (int j = 0; j < n; j++) {
			root *= x[i] - x_h[j];
		}
		double fact = 1;
		for (int j = 1; j < n + 2; j++) {
			fact *= j;
		}
		res.push_back(fabs(f5Der(x[i])) * fabs(root) / fact);
	}
	return res;
}
vector<double> ChebyshevGrid(int n, double a, double b) {
	vector<double> x_h;
/*	if (x_h.empty())
	{
		cout << "Error allocation memory for vector";
		exit(1);
	}*/	
	double tmp = 0;
	for (int i = 0; i < n+1; i++) {
		tmp = cos((M_PI * (2 * i + 1)) / (2 * (n + 1)));
		x_h.push_back(((b - a) * tmp / 2) + ((a + b) / 2));
	}
	return x_h;
}

vector<double> UniformGrid(int k, double a, double b) {
	int i = 0;
	double h = (b - a) / (k - 1); 
	double x0 = a;
	vector<double> x_h;
	/*if (x_h.empty())
	{
		cout << "Error allocation memory for vector";
		exit(1);
	}*/
	for (i = 0; i < k; i++)
	{
		x_h.push_back(x0 + i * h);
	}
	return x_h;
}

vector<double> GridFunc(vector<double> x_h, int k) {
	vector<double> y_h;
	/*if (y_h.empty())
		cout << "Error allocation memory for vector";
		exit(1);
		f(x_h[i])*/
	for (int i = 0; i < k; i++) {
		y_h.push_back(f(x_h[i]));
	}
	return y_h;
}

double LagrangePolynom(double x, vector<double> x_h, vector<double> y_h, int k)
{
	double lagrange_pol = 0;
	double basics_pol;
	vector<double> polynom;
	//for (int s = 0; s < k; s++)
	//{
	//	
		for (int i = 0; i < k; i++)
		{
			basics_pol = 1;
			for (int j = 0; j < k; j++)
			{
				if (j != i)
					basics_pol *= (x - x_h[j]) / (x_h[i] - x_h[j]);
			}
			lagrange_pol += basics_pol * y_h[i];
		}	
	//	polynom.push_back(lagrange_pol);
	//}
	//
	//return polynom;
	return lagrange_pol;
}

void PrintVector(vector<double> vec)
{
	for (auto& v : vec)
	{
		cout << setprecision(17) << v;
		cout << ' ';
	}
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

