#include <stdio.h>
#include <conio.h>
#include <math.h>
#include <stdlib.h>
#pragma warning(disable : 4996)
#define M_PI 3.14159265358979323846
#define START 0
#define FINISH M_PI / 2
double p_x(double x);
double q_x(double x);
double f_x(double x);
double y_exact(double x);

void solveMatrix(int n, double* a, double* c, double* b, double* f, double* x);
void MKR(double a, double b, int n1, double* p, double* q, double* f, double alfa0, double beta0, double* solution);
void ThomasAlgorithm(int n, double* under, double* diag, double* up, double* f, double* sol);
void Test();
int main(void)
{
	double a = START;
	double b = FINISH;
	int n = 12; //  dots amount
	double* p = (double*)malloc((n - 2)*sizeof(double));
	double* q = (double*)malloc((n - 2) * sizeof(double));
	double* f = (double*)malloc(n* sizeof(double));
	double alfa = 0;
	double beta = 1;
	double* solut = (double*)malloc(n * sizeof(double));
	FILE* f_y1 = fopen("f_y1.txt", "w");
	FILE* f_sol1 = fopen("f_sol1.txt", "w");
	FILE* f_x1 = fopen("f_x1.txt", "w");
	FILE* f_error1 = fopen("f_error1.txt", "w");
	FILE* f_y2 = fopen("f_y2.txt", "w");
	FILE* f_sol2 = fopen("f_sol2.txt", "w");
	FILE* f_x2 = fopen("f_x2.txt", "w");
	FILE* f_error2 = fopen("f_error2.txt", "w");
	FILE* f_perturb_error = fopen("f_perturb_error.txt", "w");
	FILE* f_perturbation = fopen("f_perturbation.txt", "w");
	FILE* f_actual_error = fopen("f_actual_error.txt", "w");
	FILE* f_steps = fopen("f_steps.txt", "w");
	MKR(a, b, n, p, q, f, alfa, beta, solut);
	printf("\n\n RESULTS \n\n");
	for (int i = 0; i < n; ++i)
	{
		printf("%lf\n", solut[i]);
	}

	double* y_example = (double*)malloc(4 * n * sizeof(double));
	double* x1 = (double*)malloc(2 * n * sizeof(double));
	double* sol1 = (double*)malloc(2 * n * sizeof(double));
	double* x2 = (double*)malloc(n * sizeof(double));
	double* sol2 = (double*)malloc(n * sizeof(double));

	double* p2 = (double*)malloc((n - 2) * sizeof(double));
	double* q2 = (double*)malloc((n - 2) * sizeof(double));
	double* f2 = (double*)malloc(n * sizeof(double));

	double* p1 = (double*)malloc((2 * n - 2) * sizeof(double));
	double* q1 = (double*)malloc((2*n - 2) * sizeof(double));
	double* f1 = (double*)malloc(2*n * sizeof(double));

	MKR(a, b, 2 * n - 1, p1, q1, f1, alfa, beta, sol1);
	MKR(a, b, n, p2, q2, f2, alfa, beta, sol2);

	printf("\n Example: \n");
	for (int i = 0; i < 4 * n; ++i)
	{
		y_example[i] = y_exact(a + (b - a) / (4 * n - 1) * i);
		printf("%lf ", y_example[i]);
		//fprintf(f_y, "%.15lf ", y_example[i]);
	}


	double* error1 = (double*)malloc((2 * n - 1) * sizeof(double));
	double* error2 = (double*)malloc(n * sizeof(double));


	printf("\n Solution 1: \n");
	for (int i = 0; i < 2 * n - 1; ++i)
	{
		error1[i] = fabs(sol1[i] - y_exact(a + (b - a) / (2 * n - 2) * i));
		printf("%lf ", error1[i]);
		fprintf(f_y1, "%.15lf ", y_exact(a + (b - a) / (2 * n - 2) * i));
		fprintf(f_sol1, "%.15lf ", sol1[i]);
		fprintf(f_error1, "%.15lf ", error1[i]);

		x1[i] = a + (b - a) / (2 * n - 2) * i;
		fprintf(f_x1, "%.15lf ", x1[i]);
	}


	printf("\n Solution 2: \n");
	for (int i = 0; i < n; ++i)
	{
		error2[i] = fabs(sol2[i] - y_exact(a + (b - a) / (n - 1) * i));
		printf("%lf ", error2[i]);
		fprintf(f_y2, "%.15lf ", y_exact(a + (b - a) / (n - 1) * i));
		fprintf(f_sol2, "%.15lf ", sol2[i]);
		fprintf(f_error2, "%.15lf ", error2[i]);

		x2[i] = a + (b - a) / (n - 1) * i;
		fprintf(f_x2, "%.15lf ", x2[i]);
	}



	//////-----second part------//


	int k = 100000;
	double* sol3 = (double*)malloc(k * sizeof(double));
	double* p3 = (double*)malloc((2 * k - 2) * sizeof(double));
	double* q3 = (double*)malloc((2 * k - 2) * sizeof(double));
	double* f3 = (double*)malloc((2 * k) * sizeof(double));
	double* actual_error = (double*)malloc(6 * sizeof(double));
	int step[] = { 10, 100, 1000, 10000, 20000, 100000 };
	double h_[] = { M_PI / 30, M_PI / 300, M_PI / 3000, M_PI / 30000, M_PI / 60000, M_PI / 300000 };
	double Const[6];
	for (int i = 0; i < 6; i += 1)
	{
		MKR(a, b, step[i], p3, q3, f3, alfa, beta, sol3);
		actual_error[i] = 0;
		for (int j = 0; j < step[i]; ++j)
		{
			if (actual_error[i] < fabs(y_exact(a + (M_PI / (2*(step[i] - 1))) * j) - sol3[j]))
				actual_error[i] = fabs(y_exact(a + (M_PI / (2 * (step[i] - 1))) * j) - sol3[j]);
			
			if (h_[i] == (M_PI / (2 * (step[i] - 1))))
			{
				printf("\ntrue\n");
			}
		}
		Const[i] = actual_error[i] / (h_[i]*h_[i]);
		printf("Const %lf\n", Const[i]);
		printf("%le ", actual_error[i]);
		fprintf(f_actual_error, "%.15lf ", actual_error[i]);
		fprintf(f_steps, "%.15lf ", (M_PI / (2 * (step[i] - 1))));
	}


	//////-------third part-------//

	printf("\n Perturb: \n");
	int l = 10000;
	double* sol4 = (double*)malloc(l * sizeof(double));
	double* p4 = (double*)malloc((2 * l - 2) * sizeof(double));
	double* q4 = (double*)malloc((2 * l - 2) * sizeof(double));
	double* f4 = (double*)malloc((2 * l) * sizeof(double));
	double* perturb_error = malloc(10 * sizeof(double));
	double* perturbation = malloc(10 * sizeof(double));
	for (int i = 0; i < 10; i++)
	{
		MKR(a + pow(10, -9 + i), b, 1000, p4, q4, f4, alfa, beta, sol4);
		perturbation[i] = pow(10, -9 + i);
		perturb_error[i] = 0;
		for (int j = 0; j < 1000; ++j)
		{
			if (perturb_error[i] < fabs(y_exact(a + (b - a) / (1000 - 1) * j) - sol4[j]))
				perturb_error[i] = fabs(y_exact(a + (b - a) / (1000 - 1) * j) - sol4[j]);
		}
		printf("%le ", perturb_error[i]);
		fprintf(f_perturbation, "%.15lf ", perturbation[i]);
		fprintf(f_perturb_error, "%.15lf ", perturb_error[i]);
	}
	Test();
	return 0;
}
double p_x(double x)
{
	return -tan(x);
}
double q_x(double x)
{
	return 3;
}

double f_x(double x)
{
	return sin(x);
}

double y_exact(double x)
{
	return sin(x);
}

void solveMatrix(int n, double* a, double* c, double* b, double* f, double* x)
{
	double m;
	for (int i = 1; i < n; i++)
	{
		m = a[i] / c[i - 1];
		c[i] = c[i] - m * b[i - 1];
		f[i] = f[i] - m * f[i - 1];
	}

	x[n - 1] = f[n - 1] / c[n - 1];

	for (int i = n - 2; i >= 0; i--)
	{
		x[i] = (f[i] - b[i] * x[i + 1]) / c[i];
	}
}

void ThomasAlgorithm(int n, double* under, double* diag, double* up, double* f, double* sol)
{
	double* delta = malloc((n) * sizeof(double));
	double* lambda = malloc((n) * sizeof(double));
	delta[0] = -up[0] / diag[0];
	//printf("\n");
	//printf("delta: %lf\n", delta[0]);
	lambda[0] = f[0] / diag[0];
	printf("lambda: %lf\n", lambda[0]);
	for (int i = 1; i < n - 1; i++)
	{
		//printf("\ni: %d\n", i);
		delta[i] = -up[i] / (under[i] * delta[i - 1] + diag[i]);
		//printf("delta: %lf\n", delta[i]);
		lambda[i] = (f[i] - under[i] * lambda[i - 1]) / (under[i] * delta[i - 1] + diag[i]);
		//printf("lambda: %lf\n", lambda[i]);
	}
	delta[n - 1] = 0;
	lambda[n - 1] = (f[n - 1] - under[n - 1] * lambda[n - 2] )/ (under[n - 1] * delta[n - 2] + diag[n-1]);
	//sol[n - 1] = lambda[n - 1];
	for (int i = n - 1; i >= 0; i--)
	{
		sol[i] = delta[i]*sol[i + 1] + lambda[i];
	}
}

void MKR(double a, double b, int n1, double* p, double* q, double* f, double alfa0, double beta0, double* solution)
{
	double h = (b - a) / (n1 - 1);
	for (int i = 1; i < n1 - 1; ++i)  //  except 0 and n
	{
		p[i] = p_x(a + h * i);
		//printf("pi: %lf\n", p[i]);
		q[i] = q_x(a + h * i);
		//printf("qi: %lf\n", q[i]);
		f[i] = f_x(a + h * i)*pow(h,2);
		//printf("fi: %lf\n", f[i]);
	}

	double* up = malloc((n1 - 1) * sizeof(double));
	double* diag = malloc((n1) * sizeof(double));
	double* under = malloc((n1) * sizeof(double));

	diag[0] = 1;
	f[0] = alfa0;
	up[0] = 0;

	diag[n1 - 1] = 1;
	f[n1 - 1] = beta0;
	under[n1 - 1] = 0;

	printf("\n");

	for (int i = 1; i < n1 - 1; ++i)
	{
		under[i] = (1 - h / 2 * p[i]);
		//printf("under: %lf\n", under[i]);
		diag[i] = (-2 + h * h * q[i]);
		//printf("diag: %lf\n", diag[i]);
		up[i] = (1 + h / 2 * p[i]);
		//printf("up: %lf\n", up[i]);
	}

	ThomasAlgorithm(n1, under, diag, up, f, solution);


}

void Test()
{
	int n = 5;
	double h = (FINISH - START) / (n - 1);
	double x = START;
	double alfa = 0;
	double beta = 1;
	double* p = (double*)malloc((n - 1) * sizeof(double));
	double* q = (double*)malloc((n - 1) * sizeof(double));
	double* f = (double*)malloc(n * sizeof(double));
	double* sol = (double*)malloc(n * sizeof(double));

	printf("\n TEST: \n");
	MKR(START, FINISH, n, p, q, f, alfa, beta, sol);

	
	printf("\nSolution:\n");
	for (int i = 0; i < n; ++i)
	{
		printf("%lf\n", sol[i]);
	}
	printf("\nFunction\n");
	for (int i = 0; i < n; ++i)
	{
		printf("%lf\n", y_exact(x));
		x += h;
	}
	printf("\nError\n");
	x = START;
	for (int i = 0; i < n; ++i)
	{
		printf("%lf\n", sol[i] - y_exact(x));
		x += h;
	}
}