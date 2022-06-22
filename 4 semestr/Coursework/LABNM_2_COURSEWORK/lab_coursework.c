#include <stdio.h>
#include <conio.h>
#include <math.h>
#include <stdlib.h>
#pragma warning(disable : 4996)
#define M_PI 3.14159265358979323846
#define START 0
#define FINISH M_PI / 2
#define N 12
#define Y_O 0
#define Y_DEF_O 1
#define A 0
#define B 1

//The Cauchy problem (lab 5)
double function(double x, double y, double z);
void Runge_Kutta(double a, double b, int n, double y0, double z0, double* y);
void accuracyRunge(double* y1, double* y2, FILE* f_runge1, FILE* f_runge2, FILE* f_error_runge1, FILE* f_error_runge2);
void OutrageRunge(double* y3, FILE* f_outrage_runge, FILE* f_outrage_runge_error);
void StepsErrorRunge(double n, double h, FILE* f_steps_runge_error);
//Boundary value problem (lab 7)
double p_x(double x);
double q_x(double x);
double f_x(double x);
double y_exact(double x);
void MKR(double a, double b, int n1, double* p, double* q, double* f, double alfa0, double beta0, double* solution);
void thomasAlgorithm(int n, double* under, double* diag, double* up, double* f, double* sol);
void accuracyMKR(double* y1, double* y2, FILE* f_mkr1, FILE* f_mkr2, FILE* f_error_mkr1, FILE* f_error_mkr2);
void OutrageMKR(double* y3, FILE* f_outrage_mkr, FILE* f_outrage_mkr_error);
void StepsErrorMKR(double n, double h, FILE* f_steps_mkr_error);
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
	FILE* f_outrage_mkr = fopen("f_outrage_mkr.txt","w");
	FILE* f_outrage_runge = fopen("f_outrage_runge.txt", "w");
	FILE* f_outrage_mkr_error = fopen("f_outrage_mkr_error.txt", "w");
	FILE* f_outrage_runge_error = fopen("f_outrage_runge_error.txt", "w");
	FILE* f_steps = fopen("f_steps.txt", "w");
	FILE* f_y_h1 = fopen("f_y_h1.txt", "w");
	FILE* f_y_h2 = fopen("f_y_h2.txt", "w");
	FILE* f_steps_mkr = fopen("f_steps_mkr.txt", "w");
	FILE* f_steps_runge = fopen("f_steps_runge.txt", "w");
	double* x = (double*)malloc((10 * N + 1) * sizeof(double));
	double* y = (double*)malloc((10 * N + 1) * sizeof(double));
	double* x1 = (double*)malloc((N) * sizeof(double));
	double* y1 = (double*)malloc((N) * sizeof(double));
	double* z1 = (double*)malloc((N) * sizeof(double));
	double* x2 = (double*)malloc((2 * N) * sizeof(double));
	double* y2 = (double*)malloc((2 * N) * sizeof(double));
	double* x3 = (double*)malloc((1000) * sizeof(double));
	double* y3 = (double*)malloc((1000) * sizeof(double));
	//2-nd part
	double* x4 = (double*)malloc((1000000) * sizeof(double));
	double* y4 = (double*)malloc((1000000) * sizeof(double));
	//3-rd part
	double n[10] = { 10, 50, 100, 500, 1000, 5000, 10000, 20000,50000, 100000 };
	double* h = (double*)malloc(10 * sizeof(double));
	double* y_h1 = (double*)malloc(10 * sizeof(double));
	double* y_h2 = (double*)malloc(10 * sizeof(double));
	//1-st part
	for (int i = 0; i < N; i++)
	{
		double h = (FINISH - START) / (N-1);
		x1[i] = START + i * h;
		y1[i] = y_exact(x1[i]);
		fprintf(f_x1, "%lf ", x1[i]);
	}
	for (int i = 0; i < 2 * N; i++)
	{
		double h = (FINISH - START) / (2 * N - 1);
		x2[i] = START + i * h;
		y2[i] = y_exact(x2[i]);
		fprintf(f_x2, "%lf ", x2[i]);
	}
	for (int i = 0; i < 10 * N + 1; i++)
	{
		double h = (FINISH - START) / (10 * N + 1);
		x[i] = START + i * h;
		fprintf(f_x, "%lf ", x[i]);
		y[i] = y_exact(x[i]);
		fprintf(f_y,"%lf ", y[i]);
	}
	//2-nd part
	for (int i = 0; i < 1000; i++)
	{
		double h = (FINISH - START) / (1000 - 1);
		x3[i] = START + i * h;
		y3[i] = y_exact(x3[i]);
	}
	for (int i = 0; i < 1000000; i++)
	{
		double h = (FINISH - START) / (1000000 - 1);
		x4[i] = START + i * h;
		y4[i] = y_exact(x4[i]);
	}
	//1-st part
	accuracyRunge(y1, y2, f_runge1, f_runge2, f_error_runge1, f_error_runge2);
	accuracyMKR(y1, y2, f_mkr1, f_mkr2, f_error_mkr1, f_error_mkr2);
	//2-nd part
	OutrageMKR(y3, f_outrage_mkr, f_outrage_mkr_error);
	OutrageRunge(y4, f_outrage_runge, f_outrage_runge_error);
	//3-rd extra part
	for (int i = 0; i < 10; i++)
	{
		h[i] = (FINISH - START) / (n[i] - 1);
		y_h1[i] = h[i] * h[i];
		y_h2[i] = h[i] * h[i] * h[i];
		fprintf(f_steps, "%le ", h[i]);
		fprintf(f_y_h1, "%le ", y_h1[i]);
		fprintf(f_y_h2, "%le ", y_h2[i]);
		StepsErrorRunge(n[i], h[i], f_steps_runge);
		StepsErrorMKR(n[i], h[i], f_steps_mkr);
	}
	free(x);
	free(y);
	free(x1);
	free(x2);
	free(y1);
	free(y2);
	free(x3);
	free(y3);
	free(x4);
	free(y4);
	free(h);
	free(y_h1);
	free(y_h2);
	fclose(f_x);
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
	fclose(f_outrage_mkr);
	fclose(f_outrage_runge);
	fclose(f_outrage_mkr_error);
	fclose(f_outrage_runge_error);
	fclose(f_steps_mkr);
	fclose(f_steps_runge);
	fclose(f_steps);
	fclose(f_y_h1);
	fclose(f_y_h2);
	return 0;
}

double function(double x, double y, double z)//z'
{
	return -p_x(x) * z - q_x(x) * y + f_x(x);
}

double gfunction(double x, double y, double z)
{
	return z;
}

double y_exact(double x)
{
	return sin(x);
}

void Runge_Kutta(double a, double b, int n, double y0, double z0, double* y)  //  n - amount of sections
{
    double h = (b - a) / (n-1);
	double* x = (double*)malloc(n * sizeof(double));
	double k1, k2, k3, q1, q2, q3, k4, q4;
    double t = y0;
	double s = z0;
    for (int i = 0; i < n; ++i)
    {
        x[i] = a + h * i;
        *(y + i) = y0;
        t = y0;
		s = z0;
		//printf("%lf\n", y0);
		/*k1 = function(x[i], y0, z0);
		q1 = gfunction(x[i], y0, z0);
		k2 = function(x[i] + h / 2, y0 + q1 / 2, z0 + k1 / 2);
		q2 = gfunction(x[i] + h / 2, y0 + q1 / 2, z0 + k1 / 2);
		k3 = function(x[i] + h / 2, y0 - q1 + 2 * q2, z0 - k1 + 2 * k2);
		q3 = gfunction(x[i] + h / 2, y0 - q1 + 2 * q2, z0 - k1 + 2 * k2);

		z0 = z0 + h * (k1 + 4 * k2 + k3) / 6;
        y0 = y0 + h * (q1 + 4 * q2 + q3) / 6;*/


		k1 = function(x[i], y0, z0);
		q1 = z0;

		k2 = function(x[i] + h / 2, y0 + h * k1 / 2, z0 + h * q1 / 2);
		q2 = z0 + h * q1 / 2;

		k3 = function(x[i] + h / 2, y0 + h * k2 / 2, z0 + h * q2 / 2);
		q3 = z0 + h * q2 / 2;

		k4 = function(x[i] + h, y0 + h * k3, z0 + h * q3);
		q4 = z0 + h * q3;

		z0 = z0 + h * (k1 + 2 * k2 + 2 * k3 + k4) / 6;
		y0 = y0 + h * (q1 + 2 * q2 + 2 * q3 + q4) / 6;
    }
	/*if (param == 0)
	{
		return t;
	}
	else {
		return s;
	}*/
}

void accuracyRunge(double* y1, double* y2, FILE* f_runge1, FILE* f_runge2, FILE* f_error_runge1, FILE* f_error_runge2)
{
	double* error1 = (double*)malloc((N) * sizeof(double));
	double* error2 = (double*)malloc((2*N) * sizeof(double));
	double* runge1 = (double*)malloc((N) * sizeof(double));
	double* runge2 = (double*)malloc((2 * N) * sizeof(double));
	Runge_Kutta(START, FINISH, N, Y_O, Y_DEF_O, runge1);
	Runge_Kutta(START, FINISH, 2*N, Y_O, Y_DEF_O, runge2);
	printf("RUNGE\n");
	printf("Error 1\n");
	for (int i = 0; i < N; i++)
	{
		//printf("%le\n", error1[i]);
		error1[i] = fabs(runge1[i] - y1[i]);
		fprintf(f_runge1, "%lf ", runge1[i]);
		fprintf(f_error_runge1, "%le ", error1[i]);
	}
	printf("Error 2\n");
	for (int i = 0; i < 2*N; i++)
	{
		//printf("%le\n", error2[i]);
		error2[i] = fabs(runge2[i] - y2[i]);
		fprintf(f_runge2, "%lf ", runge2[i]);
		fprintf(f_error_runge2, "%le ", error2[i]);
	}
}
void OutrageRunge(double* y3,FILE* f_outrage_runge, FILE* f_outrage_runge_error)
{
	int l = 10000000;
	double* sol = (double*)malloc(l * sizeof(double));
	double* perturb_error = malloc(10 * sizeof(double));
	double* perturbation = malloc(10 * sizeof(double));
	printf("\nOutrage Runge:\n");
	for (int i = 0; i < 10; i++)
	{
		Runge_Kutta(START + pow(10, -9 + i), FINISH, 1000000, Y_O, Y_DEF_O, sol);
		perturbation[i] = pow(10, -9 + i);
		perturb_error[i] = 0;
		for (int j = 0; j < 10000; ++j)
		{
			if (perturb_error[i] < fabs(y3[j] - sol[j]))
				perturb_error[i] = fabs(y3[j] - sol[j]);
			/*if(i == 0)
				printf("%lf %lf\n", y3[j], sol[j]);*/
		}
		printf("%le\n", perturb_error[i]);
		fprintf(f_outrage_runge, "%le ", perturbation[i]);
		fprintf(f_outrage_runge_error, "%le ", perturb_error[i]);
	}
}

void StepsErrorRunge(double n, double h, FILE* f_steps_runge_error)
{
	double* sol = (double*)malloc(n * sizeof(double));
	double* y_ex = (double*)malloc(n * sizeof(double));
	double* x = (double*)malloc(n * sizeof(double));
	double* error = (double*)malloc(n * sizeof(double));
	double max_error = 0;
	Runge_Kutta(START, FINISH, n,Y_O, Y_DEF_O, sol);
	for (int i = 0; i < n; i++)
	{
		x[i] = (FINISH - START) / (n - 1);
		y_ex[i] = y_exact(x[i]);
		error[i] = fabs(y_ex[i] - sol[i]);
		if (max_error < error[i])
			max_error = error[i];
	}
	fprintf(f_steps_runge_error, "%le ", max_error);
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

void thomasAlgorithm(int n, double* under, double* diag, double* up, double* f, double* sol)
{
	double* delta = malloc((n) * sizeof(double));
	double* lambda = malloc((n) * sizeof(double));
	delta[0] = -up[0] / diag[0];
	//printf("\n");
	//printf("delta: %lf\n", delta[0]);
	lambda[0] = f[0] / diag[0];
	//printf("lambda: %lf\n", lambda[0]);
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
	thomasAlgorithm(n1, under, diag, up, f, solution);
}

void accuracyMKR(double* y1, double* y2, FILE* f_mkr1, FILE* f_mkr2, FILE* f_error_mkr1, FILE* f_error_mkr2)
{
	double* error1 = (double*)malloc((N) * sizeof(double));
	double* error2 = (double*)malloc((2 * N) * sizeof(double));
	double* mkr1 = (double*)malloc((N) * sizeof(double));
	double* mkr2 = (double*)malloc((2 * N) * sizeof(double));
	double* p1 = malloc((N - 2) * sizeof(double));
	double* q1 = malloc((N - 2) * sizeof(double));
	double* f1 = malloc((N) * sizeof(double));
	double* p2 = malloc((2*N - 2) * sizeof(double));
	double* q2 = malloc((2*N - 2) * sizeof(double));
	double* f2 = malloc((2*N) * sizeof(double));
	MKR(START, FINISH, N, p1, q1, f1, A, B, mkr1);
	MKR(START, FINISH, 2 * N, p2, q2, f2, A, B, mkr2);
	printf("MKR\n");
	printf("Error 1\n");
	for (int i = 0; i < N; i++)
	{
		error1[i] = fabs(mkr1[i] - y1[i]);
		//printf("%le\n", error1[i]);
		fprintf(f_mkr1, "%lf ", mkr1[i]);
		fprintf(f_error_mkr1, "%le ", error1[i]);
	}
	printf("Error 2\n");
	error2[0] = 0.0;
	for (int i = 0; i < 2 * N; i++)
	{
		error2[i] = fabs(mkr2[i] - y2[i]);
		//printf("%le\n", error2[i]);
		fprintf(f_mkr2, "%lf ", mkr2[i]);
		fprintf(f_error_mkr2, "%le ", error2[i]);
	}
}

void OutrageMKR(double* y3, FILE* f_outrage_mkr, FILE* f_outrage_mkr_error)
{
	int l = 10000;
	double* sol = (double*)malloc(l * sizeof(double));
	double* p = (double*)malloc((2 * l - 2) * sizeof(double));
	double* q = (double*)malloc((2 * l - 2) * sizeof(double));
	double* f = (double*)malloc((2 * l) * sizeof(double));
	double* perturb_error = malloc(10 * sizeof(double));
	double* perturbation = malloc(10 * sizeof(double));
	printf("Outrage MKR:\n");
	for (int i = 0; i < 10; i++)
	{
		MKR(START + pow(10, -9 + i), FINISH, 1000, p, q, f, A, B, sol);
		perturbation[i] = pow(10, -9 + i);
		perturb_error[i] = 0;
		for (int j = 0; j < 1000; ++j)
		{
			if (perturb_error[i] < fabs(y3[j] - sol[j]))
				perturb_error[i] = fabs(y3[j] - sol[j]);
		}
		printf("%le ", perturb_error[i]);
		fprintf(f_outrage_mkr, "%le ", perturbation[i]);
		fprintf(f_outrage_mkr_error, "%le ", perturb_error[i]);
	}
}

void StepsErrorMKR(double n, double h, FILE* f_steps_mkr_error)
{
	double* sol = (double*)malloc(n * sizeof(double));
	double* y_ex = (double*)malloc(n * sizeof(double));
	double* x = (double*)malloc(n * sizeof(double));
	double* p = (double*)malloc((2 * n - 2) * sizeof(double));
	double* q = (double*)malloc((2 * n - 2) * sizeof(double));
	double* f = (double*)malloc((2 * n) * sizeof(double));
	double* error = (double*)malloc(n * sizeof(double));
	double max_error = 0;
	MKR(START, FINISH, n, p, q, f, A, B, sol);
	for (int i = 0; i < n; i++)
	{
		x[i] = (FINISH - START) / (n - 1);
		y_ex[i] = y_exact(x[i]);
		error[i] = fabs(y_ex[i] - sol[i]);
		if (max_error < error[i])
			max_error = error[i];
	}
	fprintf(f_steps_mkr_error, "%le ", max_error);
}