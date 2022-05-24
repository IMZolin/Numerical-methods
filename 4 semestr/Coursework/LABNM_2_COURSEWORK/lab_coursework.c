#include <stdio.h>
#include <conio.h>
#include <math.h>
#include <stdlib.h>
#pragma warning(disable : 4996)
#define M_PI 3.14159265358979323846
#define START 0
#define FINISH M_PI / 2
#define N 10
#define Y_O 0
#define Y_DEF_O 1
#define A 0
#define B 1
//The Cauchy problem (lab 5)
double function(double x, double y, double z);
double Runge_Kutta(double a, double b, int n, double y0, double z0, double* p);
void AccuracyRunge(double* y1, double* y2, FILE* f_runge1, FILE* f_runge2, FILE* f_error_runge1, FILE* f_error_runge2);
//Boundary value problem (lab 7)
double p_x(double x);
double q_x(double x);
double f_x(double x);
double y_exact(double x);
void MKR(double a, double b, int n1, double* p, double* q, double* f, double alfa0, double beta0, double* solution);
void ThomasAlgorithm(int n, double* under, double* diag, double* up, double* f, double* sol);
void AccuracyMKR(double* y1, double* y2, FILE* f_mkr1, FILE* f_mkr2, FILE* f_error_mkr1, FILE* f_error_mkr2);

int main(void)
{
	FILE* f_x = fopen("f_x.txt", "w");
	FILE* f_y = fopen("f_y.txt", "w");
	FILE* f_runge1 = fopen("f_runge1.txt", "w");
	FILE* f_runge2 = fopen("f_runge2.txt", "w");
	FILE* f_mkr1 = fopen("f_mkr1.txt", "w");
	FILE* f_mkr2 = fopen("f_mkr2.txt", "w");
	FILE* f_x1 = fopen("f_x1.txt", "w");
	FILE* f_x2 = fopen("f_x2.txt", "w");
	FILE* f_error_runge1 = fopen("f_error_runge1.txt", "w");
	FILE* f_error_runge2 = fopen("f_error_runge2.txt", "w");
	FILE* f_error_mkr1 = fopen("f_error_mkr1.txt", "w");
	FILE* f_error_mkr2 = fopen("f_error_mkr2.txt", "w");
	double* x = (double*)malloc(10 * N + 1 * sizeof(double));
	double* y = (double*)malloc(10 * N + 1 * sizeof(double));
	double* y1 = (double*)malloc((N + 1) * sizeof(double));
	double* x1 = (double*)malloc((N + 1) * sizeof(double));
	double* y2 = (double*)malloc((2 * N + 1) * sizeof(double));
	double* x2 = (double*)malloc((2 * N + 1) * sizeof(double));
	for (int i = 0; i < N; i++)
	{
		double h = (FINISH - START) / N;
		x1[i] = START + i * h;
		y1[i] = y_exact(x1[i]);
	}
	for (int i = 0; i < 2 * N + 1; i++)
	{
		double h = (FINISH - START) / (2 * N + 1);
		x2[i] = START + i * h;
		y2[i] = y_exact(x2[i]);
	}
	for (int i = 0; i < 10 * N + 1; i++)
	{
		double h = (FINISH - START) / (2 * N + 1);
		x[i] = START + i * h;
		fprintf(f_x, "%lf ", x[i]);
		y[i] = y_exact(x[i]);
		fprintf(f_y,"%lf ", y[i]);
	}
	//AccuracyRunge(y1, y2, f_runge1, f_runge2, f_error_runge1, f_error_runge2);
	AccuracyMKR(y1, y2, f_mkr1, f_mkr2, f_error_mkr1, f_error_mkr2);
	fclose(f_y);
	fclose(f_x1);
	fclose(f_x2);
	fclose(f_runge1);
	fclose(f_runge2);
	fclose(f_mkr1);
	fclose(f_mkr2);
	fclose(f_error_runge1);
	fclose(f_error_runge2);
	fclose(f_error_mkr1);
	fclose(f_error_mkr2);
	return 0;
}

double function(double x, double y, double z)//z'
{
	return -p_x(x) * z - q_x(x) * y + f_x(x);
}

double y_exact(double x)
{
	return sin(x);
}

double Runge_Kutta(double a, double b, int n, double y0, double z0, double* p)  //  n - amount of sections
{
    double h = (b - a) / n;
	double* x = (double*)malloc(n * sizeof(double));
	double k1, k2, k3, q1, q2, q3;
    double sigma[4] = { 0 };
    double gamma[4] = { 0 };
    double t = y0;

    for (int i = 0; i < n + 1; ++i)
    {
        x[i] = a + h * i;
        *(p + i) = y0;
        t = y0;

        //printf("%lf %lf\n", y0, z0);
		k1 = function(x[i], y0, z0);
		q1 = z0;
		k2 = function(x[i] + h / 3, y0 + h * k1 / 3, z0 + h * k1 / 3);
		q2 = z0 + h * k1 / 3;
		k3 = function(x[i] + 2*h / 3, y0 + 2*h * k2 / 3, z0 + h * k2 / 3);;
		q3 = z0 + 2*h * k2 / 3;

		z0 = z0 + h * (k1 + 3 * k3) / 4;
        y0 = y0 + h * (k1 + 3 * k3) / 4;
    }

    return t;
}

void AccuracyRunge(double* y1, double* y2, FILE* f_runge1, FILE* f_runge2, FILE* f_error_runge1, FILE* f_error_runge2)
{
	double* error1 = (double*)malloc((N + 1) * sizeof(double));
	double* error2 = (double*)malloc((2*N + 1) * sizeof(double));
	double* runge1 = (double*)malloc((N + 1) * sizeof(double));
	double* runge2 = (double*)malloc((2 * N + 1) * sizeof(double));
	Runge_Kutta(START, FINISH, N, Y_O, Y_DEF_O, runge1);
	Runge_Kutta(START, FINISH, 2*N, Y_O, Y_DEF_O, runge2);
	for (int i = 0; i < N + 1; i++)
	{
		printf("%le ", error1[i]);
		error1[i] = fabs(runge1[i] - y1[i]);
		fprintf(f_runge1, "%lf ", runge1[i]);
		fprintf(f_error_runge1, "%le ", error1[i]);
	}
	for (int i = 0; i < 2*N + 1; i++)
	{
		printf("%le ", error2[i]);
		error2[i] = fabs(runge2[i] - y2[i]);
		fprintf(f_runge1, "%lf ", runge2[i]);
		fprintf(f_error_runge1, "%le ", error2[i]);
	}
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
	lambda[n - 1] = (f[n - 1] - under[n - 1] * lambda[n - 2]) / (under[n - 1] * delta[n - 2] + diag[n - 1]);
	for (int i = n - 1; i >= 0; i--)
	{
		sol[i] = delta[i] * sol[i + 1] + lambda[i];
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
		f[i] = f_x(a + h * i) * pow(h, 2);
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

void AccuracyMKR(double* y1, double* y2, FILE* f_mkr1, FILE* f_mkr2, FILE* f_error_mkr1, FILE* f_error_mkr2)
{
	double* error1 = (double*)malloc((N) * sizeof(double));
	double* error2 = (double*)malloc((2 * N - 1) * sizeof(double));
	double* mkr1 = (double*)malloc((N + 1) * sizeof(double));
	double* mkr2 = (double*)malloc((2 * N + 1) * sizeof(double));
	double* p1 = malloc((N - 2) * sizeof(double));
	double* q1 = malloc((N - 2) * sizeof(double));
	double* f1 = malloc((N) * sizeof(double));
	double* p2 = malloc((2*N - 2) * sizeof(double));
	double* q2 = malloc((2*N - 2) * sizeof(double));
	double* f2 = malloc((2*N) * sizeof(double));
	MKR(START, FINISH, N, p1, q1, f1, A, B, mkr1);
	MKR(START, FINISH, 2 * N - 1, p2, q2, f2, A, B, mkr2);
	for (int i = 0; i < N; i++)
	{
		error1[i] = fabs(mkr1[i] - y1[i]);
		printf("%le ", error1[i]);
		fprintf(f_mkr1, "%lf ", mkr1[i]);
		fprintf(f_error_mkr1, "%le ", error1[i]);
	}
	for (int i = 0; i < 2 * N - 1; i++)
	{
		error2[i] = fabs(mkr2[i] - y2[i]);
		printf("%le ", error2[i]);
		fprintf(f_mkr2, "%lf ", mkr2[i]);
		fprintf(f_error_mkr2, "%le ", error2[i]);
	}
}