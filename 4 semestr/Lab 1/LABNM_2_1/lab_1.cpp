#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <math.h>
using namespace std;
# define M_PI 3.14159265358979323846
#define A 1
#define B 5
#define N 100
#define DEGREE_POLYNOMIAL 2
#define MAX_NODES DEGREE_POLYNOMIAL+1
#define MAX_ERROR_NODES 100
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
	string filename_x = "files/f_x.txt";
	string filename_y = "files/f_y.txt";
	vector<double> max_error;
	PrintVectorToFile(filename_x, x);
	PrintVectorToFile(filename_y, y);
	ofstream file;
	for (int i = 2; i < MAX_NODES; i += 2)
	{
		vector<double>teoretic_error;
		vector<double>actual_error;
		vector<double> polynom;
		vector<double> x_midpoints;
		vector<double> y_midpoints;
		vector<double> polynom_max;
		vector<double> x_h_max;
		vector<double> y_h_max;
		vector<double> check;

		vector<double> x_h = ChebyshevGrid(i, a, b);
		vector<double> y_h = GridFunc(x_h, i);
		for (int j = 0; j < N + 1; j++)
		{
			pol = LagrangePolynom(x[j], x_h, y_h, i);//уменьшил количество узлов на 1
			actual_error.push_back((y[j] - pol));
			polynom.push_back(pol);
		}
		teoretic_error = calcTheoreticError(x, x_h, i);

		x_midpoints = GridForMaxError(x_h, i);
		y_midpoints = GridMaxFunc(x_midpoints, i);
		for (int j = 0; j < i; j++)
		{
			pol = LagrangePolynom(x_midpoints[j], x_h, y_h, i);//уменьшил количество узлов на 1
			//actual_error.push_back((y[j] - pol));
			polynom_max.push_back(pol);
		}
		max_error.push_back(calcMaxError(y_midpoints, polynom_max));
		int k = i + 1;
		nodes.push_back(k);

		string filename_x_h;
		string filename_y_h;
		string filename_polynom;
		string filename_actual_error;
		string filename_th_error;
		string filename_nodes;
		string filename_max_error;
		string filename_x;
		string filename_y;
		if (i == 2 || i == 4 || i == 6 || i >=20)
		{
			if (i == 2)
			{
				filename_x_h = "files/f_x_h1.txt";
				filename_y_h = "files/f_y_h1.txt";
				filename_polynom = "files/f_polynom1.txt";
				filename_actual_error = "files/f_actual_error1.txt";
				filename_th_error = "files/f_th_error1.txt";
				filename_nodes = "files/f_nodes1.txt";
				filename_max_error = "files/f_max_error1.txt";
			}
			else if (i == 4)
			{
				filename_x_h = "files/f_x_h2.txt";
				filename_y_h = "files/f_y_h2.txt";
				filename_polynom = "files/f_polynom2.txt";
				filename_actual_error = "files/f_actual_error2.txt";
				filename_th_error = "files/f_th_error2.txt";
				filename_nodes = "files/f_nodes2.txt";
				filename_max_error = "files/f_max_error2.txt";
			}
			else if (i == 6)
			{
				filename_x_h = "files/f_x_h3.txt";
				filename_y_h = "files/f_y_h3.txt";
				filename_polynom = "files/f_polynom3.txt";
				filename_actual_error = "files/f_actual_error3.txt";
				filename_th_error = "files/f_th_error3.txt";
				filename_nodes = "files/f_nodes3.txt";
				filename_max_error = "files/f_max_error3.txt";

			}
			else if (i >= 20)
			{
				filename_x_h = "files/f_x_h30.txt";
				filename_y_h = "files/f_y_h30.txt";
				filename_polynom = "files/f_polynom30.txt";
				filename_actual_error = "files/f_actual_error30.txt";
				filename_th_error = "files/f_th_error30.txt";
				filename_nodes = "files/f_nodes30.txt";
				filename_max_error = "files/f_max_error30.txt";
			}
			PrintVectorToFile(filename_x_h, x_h);
			PrintVectorToFile(filename_y_h, y_h);
			PrintVectorToFile(filename_polynom, polynom);
			PrintVectorToFile(filename_actual_error, actual_error);
			PrintVectorToFile(filename_th_error, teoretic_error);
			PrintVectorToFile(filename_max_error, max_error);
			PrintVectorToFile(filename_nodes, nodes);
		}
		filename_x = "files/f_m_nodes.txt";
		filename_y = "files/f_m_pol_nodes.txt";
		PrintVectorToFile(filename_x, x_midpoints);
		PrintVectorToFile(filename_y, polynom_max);
		//filename_nodes = "files/f_nodes3.txt";
		//filename_max_error = "files/f_max_error.txt";
		//PrintVectorToFile(filename_max_error, max_error);

		
	}
	return 0;
}

double f(double x) {
	//return log10(x) + 7.0 / (2 * x + 6);
	return pow((2.0 * x + log(x)),0.5);
}

double f5Der(double x)
{
	return 24 / log(10) * pow(x, 5) - 26880 / pow(2 * x + 6, 6);
	//return 7.214;
}

vector<double> GridForMaxError(vector<double> x_h, int k) {
	vector<double> x_h_max;
	for (int i = 0; i < k; i++) {
		x_h_max.push_back((x_h[i] + x_h[i + 1]) / 2);
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
	for (int i = 0; i < k + 1; i++)//k?
	{
		basics_pol = 1;
		for (int j = 0; j < k + 1; j++)//k?
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
	for (int i = 0; i < n + 1; i++)
	{
		error.push_back(polynom[i] - func[i]);
	}
	return error;
}

vector<double> calcTheoreticError(vector<double> x, vector<double> x_h, int n) {
	vector<double> res;
	for (int i = 0; i < N + 1; i++) {
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
	for (int i = 0; i < n + 1; i++) {//n?
		tmp = cos((M_PI * (2 * i + 1)) / (2 * (n + 1)));//n?
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
	for (int i = 0; i < k + 1; i++) {//x_h.size()
		y_h.push_back(f(x_h[i]));
	}
	return y_h;
}

double LagrangePolynom(double x, vector<double> x_h, vector<double> y_h, int k)
{
	double lagrange_pol = 0;
	double basics_pol;
	vector<double> polynom;
	for (int i = 0; i < x_h.size(); i++)//k?
	{
		basics_pol = 1;
		for (int j = 0; j < x_h.size(); j++)//k?
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

