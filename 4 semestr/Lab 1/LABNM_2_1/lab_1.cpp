#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <math.h>
using namespace std;
# define M_PI 3.14159265358979323846
#define A 0.0
#define B 1.0
#define N 100
#define MAX_NODES 5
#define MAX_ERROR_NODES 30
double f(double x);
double f5Der(double x);
vector<double> ChebyshevGrid(int n, double a, double b);
vector<double> UniformGrid(int k, double a, double b);
vector<double> GridFunc(vector<double> x_h, int k);
double LagrangePolynom(double x, vector<double> x_h, vector<double> y_h, int k);
void PrintVectorToFile(string filename, vector<double> vec);
vector<double> calcTheoreticError(vector<double> x, vector<double> x_h, int n);
vector<double> calcActualticError(vector<double> polynom, vector<double>func, int n);
vector<double> GridForMaxError(vector<double> x_h, int k);
void PrintVector(vector<double> vec);
double calcMaxError(vector<double> y_h, vector<double> polynom);
vector<double> GridMaxFunc(vector<double> x_h, int k);
double LagrangePolynomMax(double x, vector<double> x_h, vector<double> y_h, int k);
int main()
{
	double a = A;
	double b = B;

	double pol = 0;
	double max = 0;
	vector<double> x = ChebyshevGrid(N, a, b);
	vector<double> y = GridFunc(x, N);
	vector<double> nodes;
	vector<double> nodes_max;
	string filename;
	vector<double> max_error;
	filename = "f_x.txt";
	PrintVectorToFile(filename, x);
	filename = "f_y.txt";
	PrintVectorToFile(filename, y);
	ofstream file;
	for (int i =3; i <= MAX_ERROR_NODES; i+=2)
	{
		vector<double>teoretic_error;
		vector<double>actual_error;
		vector<double>actual_error_nodes;
		vector<double> polynom;
		vector<double> x_h_max;
		vector<double> y_h_max;
		vector<double> polynom_max;
		vector<double> x_h = ChebyshevGrid(i, a, b);
		
		if (i == 3)
		{
			filename = "f_x_h1.txt";
		}
		else if (i == 4)
		{
			filename = "f_x_h2.txt";
		}
		else if (i == 5)
		{
			filename = "f_x_h3.txt";
		}
		
		PrintVectorToFile(filename, x_h);
		vector<double> y_h = GridFunc(x_h, i);
		if (i == 3)
		{
			filename = "f_y_h1.txt";
		}
		else if (i == 4)
		{
			filename = "f_y_h2.txt";
		}
		else if (i == 5)
		{
			filename = "f_y_h3.txt";
		}
		
		PrintVectorToFile(filename, y_h);
		for (int j = 0; j <= N; j++)
		{
			pol = LagrangePolynom(x[j], x_h, y_h, i);
			polynom.push_back(pol);
		}
		cout <<"\nNumb of nodes: " <<i;
		//cout << "\nApproximation from polynom\n";
		//PrintVector(polynom);
		//cout << "\nActual error\n";
		actual_error = calcActualticError(polynom, y, N);
		//PrintVector(actual_error);
		//cout << "\nTheoretic error\n";
		teoretic_error = calcTheoreticError(x, x_h, i);

		x_h_max = GridForMaxError(x_h, i);
		y_h_max = GridMaxFunc(x_h_max, i);
		cout << "\nY\n";
		PrintVector(y_h_max);
		for (int j = 0; j < i; j+=1)
		{
			pol = LagrangePolynomMax(x_h_max[j], x_h, y_h, j);//N - длинна x или y
			polynom_max.push_back(pol);
		}
		cout << "\nPOLYNOM\n";
		PrintVector(polynom_max);
		max = calcMaxError(y_h_max, polynom_max);
		max_error.push_back(max);		
		nodes.push_back(i);
		//PrintVector(teoretic_error);
		cout << "\nMAX\n";
		PrintVector(max_error);

		cout << "\n";
		if (i == 3)
		{
			filename = "f_polynom1.txt";
			
		}
		else if (i == 4)
		{
			filename = "f_polynom2.txt";
		}
		else if (i == 5)
		{
			filename = "f_polynom3.txt";
		}
		PrintVectorToFile(filename, polynom);

		if (i == 3)
		{
			filename = "f_actual_error1.txt";
		}
		else if (i == 4)
		{
			filename = "f_actual_error2.txt";
		}
		else if (i == 5)
		{
			filename = "f_actual_error3.txt";
		}
		PrintVectorToFile(filename, actual_error);

		if (i == 3)
		{
			filename = "f_theoritic_error1.txt";
		}
		else if (i == 4)
		{
			filename = "f_theoritic_error2.txt";
		}
		else if (i == 5)
		{
			filename = "f_theoritic_error3.txt";
		}
		PrintVectorToFile(filename, teoretic_error);

		filename = "f_max_error.txt";
		PrintVectorToFile(filename, max_error);

		filename = "f_nodes.txt";
		PrintVectorToFile(filename, nodes);
		cout << "Nodes\n";
		PrintVector(nodes);
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
	for (int i = 0; i < k; i++) {
		x_h_max.push_back((x_h[i]+x_h[i+1])/2);
	}
	return x_h_max;
}

vector<double> GridMaxFunc(vector<double> x_h, int k) {
	vector<double> y_h;
	for (int i = 0; i < k; i++) {
		y_h.push_back(f(x_h[i]));
	}
	return y_h;
}

double calcMaxError(vector<double> y_h, vector<double> polynom)
{
	double max = fabs(polynom[0] - y_h[0]);
	for (int i = 1; i < y_h.size(); i++)
	{
		if (max < fabs(polynom[i] - y_h[i]))
		{
			max = fabs(polynom[i] - y_h[i]);
		}
	}
	return max;
}

double LagrangePolynomMax(double x, vector<double> x_h, vector<double> y_h, int k)
{
	double lagrange_pol = 0;
	double basics_pol;
	vector<double> polynom;
	for (int i = 0; i < k+1; i++)
	{
		basics_pol = 1;
		for (int j = 0; j < k+1; j++)
		{
			if (j != i)
				basics_pol *= (x - x_h[j]) / (x_h[i] - x_h[j]);
		}
		lagrange_pol += basics_pol * y_h[i];
	}
	return lagrange_pol;
}
vector<double> calcActualticError(vector<double> polynom, vector<double>func, int n)
{
	vector<double> error;
	for (int i = 0; i < n+1; i++)
	{
		error.push_back(polynom[i] - func[i]);
	}
	return error;
}

vector<double> calcTheoreticError(vector<double> x, vector<double> x_h, int n) {
	vector<double> res;
	for (int i = 0; i < N+1; i++) {
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
	for (i = 0; i < k; i++)
	{
		x_h.push_back(x0 + i * h);
	}
	return x_h;
}

vector<double> GridFunc(vector<double> x_h, int k) {
	vector<double> y_h;
	for (int i = 0; i < k+1; i++) {
		y_h.push_back(f(x_h[i]));
	}
	return y_h;
}

double LagrangePolynom(double x, vector<double> x_h, vector<double> y_h, int k)
{
	double lagrange_pol = 0;
	double basics_pol;
	vector<double> polynom;
	for (int i = 0; i < k+1; i++)
	{
		basics_pol = 1;
		for (int j = 0; j < k+1; j++)
		{
			if (j != i)
				basics_pol *= (x - x_h[j]) / (x_h[i] - x_h[j]);
		}
		lagrange_pol += basics_pol * y_h[i];
	}	
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

