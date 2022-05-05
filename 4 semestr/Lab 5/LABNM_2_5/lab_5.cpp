#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#pragma warning(disable : 4996)
# define M_PI 3.14159265358979323846
#define A M_PI / 2
#define B 2 * M_PI
#define N 30
#define F(A) M_PI / 2
//#define F(A) 1.0
#define P 2
#define E 0.000000001

double fDiff(double x);
double f(double x, double y);
double Runge_Kutta_method(int m, double h, double x, double y);
double* createSteadyGrid(double a, double b, int n);
void cycle(double y, double* x, double eps, FILE* f_error, FILE* f_iter, int wtf);
void RungePoints(double h, FILE* f_euler, FILE* f_stepLength, FILE* f_x, FILE* f_error);
void RungePoints2(double h, FILE* f_euler, FILE* f_stepLength, FILE* f_x, FILE* f_error);
int main(int argc, char* argv[]) {
	FILE* f_y = fopen("f_y.txt", "w");
	FILE* f_x = fopen("f_x.txt", "w");
	FILE* f_euler = fopen("f_euler.txt", "w");
	FILE* f_error = fopen("f_error.txt", "w");
	FILE* f_eps = fopen("f_eps.txt", "w");
	FILE* f_iter = fopen("f_iter.txt", "w");
	FILE* f_wtfError = fopen("f_wtfError.txt", "w");
	FILE* f_wtfDer = fopen("f_wtfDer.txt", "w");
	FILE* f_maxError = fopen("f_maxError.txt", "w");
	FILE* f_stepLength = fopen("f_stepLength.txt", "w");
	if (!f_y || !f_x || !f_euler || !f_error || !f_eps || !f_iter || !f_wtfDer || !f_maxError) {
		printf("Error opening file\n");
		return 1;
	}

	double* x = createSteadyGrid(A, B, N);
	for (int i = 0; i < N; i++) {
		fprintf(f_x, "%.15lf ", x[i]);
		fprintf(f_y, "%.15lf ", fDiff(x[i]));
	}

	RungePoints((B - A) / 23, f_euler, f_stepLength, f_x, f_error);
	RungePoints2((B - A) / 46, f_euler, f_stepLength, f_x, f_error);
	/*double eps = 0.000000001;
	double y = F(A);
	for (double c = 0.1; c >= 0.00000000000001; c /= 10) {
		printf("Hui");
		fprintf(f_wtfDer, "%.15lf ", c);
		cycle(y - c, x, eps, f_wtfError, f_iter, 1);
	}*/

	for (double eps = 0.001; eps > 0.000000000000001; eps /= 10) {
		fprintf(f_eps, "%.15lf ", eps);
		printf("\n\neps=%e\n", eps);
		double y = F(A);
		cycle(y, x, eps, f_maxError, f_iter, 0);
		if (eps == 0.00001) {
			for (double c = 0.1; c >= 0.00000000000001; c /= 10) {
				fprintf(f_wtfDer, "%.15lf ", c);
				cycle(y - c, x, eps, f_wtfError, f_iter, 1);
			}
		}
	}

	fclose(f_eps);
	fclose(f_error);
	fclose(f_x);
	fclose(f_y);
	fclose(f_euler);
	fclose(f_iter);
	fclose(f_wtfError);
	fclose(f_wtfDer);
	fclose(f_stepLength);
	fclose(f_maxError);
	free(x);
	return 0;
}

void cycle(double y, double* x, double eps, FILE* f_error, FILE* f_iter, int wtf) {
	for (int i = 0; i < N - 1; i++) {
		double h = (B - A) / (N - 1.0);
		int m = 1;
		double prev;
		double next = Runge_Kutta_method(m, h, x[i], y);
		do {
			h /= 2;
			m *= 2;
			prev = next;
			next = Runge_Kutta_method(m, h, x[i], y);

		} while (fabs((next - prev)) / (pow(2, 3) - 1) >= eps);

		y = next;
		if (i == 0 && wtf == 0) {
			fprintf(f_iter, "%d ", (int)log2(m));
		}
		if (wtf == 0) {
			printf("\nx=%.15lf\ny=%.15lf\nnext=%.15lf\n", x[i + 1], fDiff(x[i + 1]), next);
		}
		printf("error=%.15lf\n", fabs(next - fDiff(x[i + 1])));
		fprintf(f_error, "%.15lf ", fabs(next - fDiff(x[i + 1])));
	}
}

double* createSteadyGrid(double a, double b, int n) {
	double* x = (double*)malloc(sizeof(double) * n);
	if (!x) {
		printf("Error allocating memory to x\n");
		exit(1);
	}

	for (int i = 0; i < n; i++) {
		x[i] = a + (b - a) * i / (n - 1.0);
	}

	return x;
}


double fDiff(double x) {
	return  x * sin(x);
}

/*double f(double x, double y) {
	//return (4 * x + 2 * y) / (2 * x + 1);
	return 2 * y / (x * log(x)) + 1 / x;
}*/

double f(double x, double y) {

	return y / x + x * cos(x);
}

/*double f(double x, double y) {
	//return (4 * x + 2 * y) / (2 * x + 1);
	double res;
	res = (2 * (log(x) - 1) + 1) / x;
	return res;
}*/

/*double euler(int m, double h, double x, double y) {
	double halfY = y;
	for (int i = 0; i < m; i++) {
		halfY = y + h / 2 * f(x, y);
		y = y + h * f(x + h / 2, halfY);
		x += h;
	}

	return y;
}*/

double Runge_Kutta_method(int m, double h, double x, double y) {
	double k1, k2, k3 = y;
	double x2, x3, y2, y3;
	for (int i = 0; i < m; i++) {
		k1 = f(x, y);
		x2 = x + h / 2;
		y2 = y + h * k1 / 2;
		k2 = f(x2, y2);
		y3 = y - h * k1 + 2 * h * k2;
		x3 = x + h;
		k3 = f(x3, y3);
		y = y + h * (k1 + 4 * k2 + k3) / 6;
		x += h;
	}

	return y;
}

void RungePoints(double h, FILE* f_euler, FILE* f_stepLength, FILE* f_x, FILE* f_error) {
	double y = F(A);
	int i = 0;

	//y = Runge_Kutta_method(1, h, A, y);
	fprintf(f_x, "%.15lf ", A);
	fprintf(f_euler, "%.15lf ", y);
	fprintf(f_error, "%.15lf ", fabs(0));
	for (double x = A; x < B; x += h, i++) {
		y = Runge_Kutta_method(1, h, x, y);
		fprintf(f_x, "%.15lf ", x+h);
		fprintf(f_euler, "%.15lf ", y);
		fprintf(f_error, "%.15lf ", fabs(y - fDiff(x+h)));
	}

	fprintf(f_stepLength, "%d ", i);
}

void RungePoints2(double h, FILE* f_euler, FILE* f_stepLength, FILE* f_x, FILE* f_error) {
	double y = F(A);
	int i = 0;

	//y = Runge_Kutta_method(1, h, A, y);
	fprintf(f_x, "%.15lf ", A);
	fprintf(f_euler, "%.15lf ", y);
	fprintf(f_error, "%.15lf ", fabs(0));
	for (double x = A; x < B-h; x += h, i++) {
		y = Runge_Kutta_method(1, h, x, y);
		fprintf(f_x, "%.15lf ", x + h);
		fprintf(f_euler, "%.15lf ", y);
		fprintf(f_error, "%.15lf ", fabs(y - fDiff(x + h)));
	}

	fprintf(f_stepLength, "%d ", i);
}


