#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#pragma warning(disable : 4996)
# define M_PI 3.14159265358979323846
#define A M_PI / 2
#define B 2 * M_PI
#define N 30
#define F(A) M_PI / 2
#define P 2

double fDiff(double x);
double f(double x, double y);
double Runge_Kutta_method(int m, double h, double x, double y);
double* createSteadyGrid(double a, double b, int n);
void cycle(double y, double* x, double eps, FILE* f_error, FILE* f_iter, int wtf);
void RungePoints(double h, FILE* f_euler, FILE* f_stepLength, FILE* f_x, FILE* f_error);
void Outrage(FILE* f_eps2, FILE* f_outrage, FILE* f_outrage_error, FILE* f_iter2);
void Accuracy(FILE* f_eps, FILE* f_maxError, FILE* f_iter);
int main(int argc, char* argv[]) {
	FILE* f_y = fopen("f_y.txt", "w");
	FILE* f_x = fopen("f_x.txt", "w");
	FILE* f_x1 = fopen("f_x1.txt", "w");
	FILE* f_x2 = fopen("f_x2.txt", "w");
	FILE* f_runge = fopen("f_runge1.txt", "w");
	FILE* f_runge2 = fopen("f_runge2.txt", "w");
	FILE* f_error = fopen("f_error1.txt", "w");
	FILE* f_error2 = fopen("f_error2.txt", "w");
	FILE* f_eps = fopen("f_eps.txt", "w");
	FILE* f_eps2 = fopen("f_eps2.txt", "w");
	FILE* f_iter = fopen("f_iter.txt", "w");
	FILE* f_iter2 = fopen("f_iter2.txt", "w");
	FILE* f_wtfError = fopen("f_wtfError.txt", "w");
	FILE* f_wtfDer = fopen("f_wtfDer.txt", "w");
	FILE* f_maxError = fopen("f_maxError.txt", "w");
	FILE* f_stepLength = fopen("f_stepLength.txt", "w");
	FILE* f_steps = fopen("f_steps.txt", "w");
	FILE* f_actual_error = fopen("f_actual_error.txt", "w");
	FILE* f_outrage = fopen("f_outrage.txt", "w");
	FILE* f_outrage_error = fopen("f_outrage_error.txt", "w");
	if (!f_y || !f_x || !f_runge || !f_error || !f_eps || !f_iter || !f_wtfDer || !f_maxError) {
		printf("Error opening file\n");
		return 1;
	}

	double* x = createSteadyGrid(A, B, N);
	for (int i = 0; i < N; i++) {
		fprintf(f_x, "%.15lf ", x[i]);
		fprintf(f_y, "%.15lf ", fDiff(x[i]));
	}
	//1-st part
	RungePoints((B - A) / 23, f_runge, f_stepLength, f_x1, f_error);
	RungePoints((B - A) / 46, f_runge2, f_stepLength, f_x2, f_error2);
	//2-nd part
	Accuracy(f_eps, f_maxError, f_iter);

	//3-rd part
	Outrage(f_eps2, f_outrage, f_outrage_error, f_iter2);

	//last part
	double* actual_error = (double*)malloc(8 * sizeof(double));
	double length = B - A;
	//double h[] = { length / 10, length / 50, length / 100, length / 500, length/1000, length / 5000, length/10000, length / 50000 };
	double h[] = { length / 10, length / 50, length / 100, length / 500, length / 1000, length / 5000, length / 10000, length / 50000 };
	double runge = F(A);
	printf("\n");
	for (int i=0; i<8; i++)
	{
		fprintf(f_steps, "%.15lf ", h[i]);
		actual_error[i] = 0.0;
		runge = F(A);
		for (double x = A; x < B; x += h[i])
		{
			runge = Runge_Kutta_method(1, h[i], x, runge);
			if(actual_error[i] < fabs(fDiff(x + h[i]) - runge))
				actual_error[i] = fabs(fDiff(x + h[i]) - runge);
			
			if (i == 0)
			{
				//printf("Runge %le\n", runge);
				//printf("Func %le\n", fDiff(x + h[i]));
				//printf("Error %le\n", actual_error[i]);
			}
		}
		//printf("FINISH\n");
		//printf("h %le\n", h[i]);
		//printf("Error %le\n", actual_error[i]);
		fprintf(f_actual_error, "%.15lf ", actual_error[i]);
	}

	fclose(f_eps);
	fclose(f_error);
	fclose(f_x);
	fclose(f_y);
	fclose(f_runge);
	fclose(f_iter);
	fclose(f_wtfError);
	fclose(f_wtfDer);
	fclose(f_stepLength);
	fclose(f_maxError);
	fclose(f_steps);
	fclose(f_actual_error);
	fclose(f_outrage_error);
	fclose(f_outrage);
	free(x);
	return 0;
}

void Accuracy(FILE* f_eps, FILE* f_maxError, FILE* f_iter)
{
	double* x = createSteadyGrid(A, B, N);
	for (double eps = 0.001; eps > 0.00000000000001; eps /= 10) {
		fprintf(f_eps, "%.15lf ", eps);
		printf("\n\neps=%e\n", eps);
		double y = F(A);
		cycle(y, x, eps, f_maxError, f_iter, 0);
	}
}

void Outrage(FILE* f_eps2, FILE* f_outrage, FILE* f_outrage_error, FILE* f_iter2)
{
	double* x = createSteadyGrid(A, B, N);
	for (double eps = 0.001; eps > 0.00000000000001; eps /= 10) {
		fprintf(f_eps2, "%.15lf ", eps);
		printf("\n\neps=%e\n", eps);
		double y = F(A);
		if (eps == 0.00001) {
			for (double c = 0.1; c >= 0.00000000000001; c /= 10) {
				fprintf(f_outrage, "%.15lf ", c);
				cycle(y - c, x, eps, f_outrage_error, f_iter2, 1);
			}
		}
	}
}

void cycle(double y, double* x, double eps, FILE* f_error, FILE* f_iter, int wtf) {
	double error = 0.0;
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
		//printf("m %d\n", (int)log2(m));
		//printf("h %le\n", h);
		y = next;
		if (i == 0 && wtf == 0) {
			fprintf(f_iter, "%d ", (int)log2(m));
		}
		if (wtf == 0) {
			printf("\nx=%.15lf\ny=%.15lf\nnext=%.15lf\n", x[i + 1], fDiff(x[i + 1]), next);
		}
		printf("error=%le\n", fabs(next - fDiff(x[i + 1])));
		if (error < fabs(next - fDiff(x[i + 1])))
			error = fabs(next - fDiff(x[i + 1]));
	}
	fprintf(f_error, "%.15lf ", error);
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

double f(double x, double y) {

	return y / x + x * cos(x);
}

double Runge_Kutta_method(int m, double h, double x, double y) {
	double k1, k2, k3 = y;
	double x2, x3, y2, y3;
	for (int i = 0; i < m; i++) {
		k1 = f(x, y);
		x2 = x + h / 3;
		y2 = y + h * k1 / 3;
		k2 = f(x2, y2);
		y3 = y + 2 * h * k2/3;
		x3 = x + 2 * h/ 3;
		k3 = f(x3, y3);
		y = y + h * (k1 + 3 * k3) / 4;
		x += h;
	}
	return y;
}

void RungePoints(double h, FILE* f_runge, FILE* f_stepLength, FILE* f_x, FILE* f_error) {
	double y = F(A);
	int i = 0;

	fprintf(f_x, "%.15lf ", A);
	fprintf(f_runge, "%.15lf ", y);
	fprintf(f_error, "%.15lf ", fabs(0));
	for (double x = A; x < B; x += h, i++) {
		y = Runge_Kutta_method(1, h, x, y);
		fprintf(f_x, "%.15lf ", x+h);
		fprintf(f_runge, "%.15lf ", y);
		fprintf(f_error, "%.15lf ", fabs(y - fDiff(x+h)));
	}

	fprintf(f_stepLength, "%d ", i);
}



