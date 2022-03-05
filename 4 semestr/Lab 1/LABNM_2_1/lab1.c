#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#pragma warning(disable : 4996)
#define N 30
#define A 0.0
#define B 1.0
#define PI 3.14
typedef struct {
	double* gridX;
	double* gridY;
	double* gridDer;
	int n;
} grid_t;

void printVectorToFile(FILE* f, double* vector, int n);
double f(double x);
double fDer(double x);
grid_t createGrid(int gridN, double a, double b);
void printVector(double* vector, int n);
double* createErmitPolynomValues(grid_t grid, double* x, int n);
double* createX(double a, double b, int n);
double* createFunctionValues(double* x, int n);
double* calcTheoreticError(double* lagrangeValues,double* x, grid_t grid, int n);
double f5Der(double x);
double* Lagrangepolynomial(double* x, grid_t grid, int n);

int main(int argc, char* argv[]) {
	double* x = createX(A, B, N);
	FILE* f_gridN = fopen("f_gridN.txt", "w");
	FILE* f_gridX = fopen("f_gridX.txt", "w");
	FILE* f_gridY = fopen("f_gridY.txt", "w");
	FILE* f_p = fopen("f_p.txt", "w");
	FILE* f_y = fopen("f_y.txt", "w");
	FILE* f_x = fopen("f_x.txt", "w");
	FILE* f_err = fopen("f_err.txt", "w");
	FILE* f_theory = fopen("f_theory.txt", "w");
	if (!f_gridN || !f_gridX || !f_gridY || !f_x || !f_y || !f_err || !f_p || !f_theory) {
		printf("Error opening file\n");
		exit(1);
	}
	fprintf(f_x, "%d ", N);
	printVectorToFile(f_x, x, N);
	double* functionValues = createFunctionValues(x, N);
	printVectorToFile(f_y, functionValues, N);
	for (int gridN = 4; gridN < 31; gridN += 2) {
		printf("\n\nn = %d\n", gridN);
		fprintf(f_gridN, "%d ", gridN);

		grid_t grid = createGrid(gridN, A, B);
		double* lagrangeValues = Lagrangepolynomial(x, grid, N);
		if (gridN == 4) {
			double* theoryError = calcTheoreticError(lagrangeValues, x, grid, N);
			printf("\ntheoryError = ");
			printVector(theoryError, N);
			printVectorToFile(f_theory, theoryError, N);
		}
		printf("\ngridX = ");
		printVector(grid.gridX, grid.n);
		printVectorToFile(f_gridX, grid.gridX, grid.n);

		printf("\ngridY = ");
		printVector(grid.gridY, grid.n);
		printVectorToFile(f_gridY, grid.gridY, grid.n);

		printf("\ngriDy' = ");
		printVector(grid.gridDer, grid.n);

		printf("\nx = ");
		printVector(x, N);

		//double* lagrangeValues = Lagrangepolynomial(x, grid, N);
		printf("\nP(x) = ");
		printVector(lagrangeValues, N);
		printVectorToFile(f_p, lagrangeValues, N);

		printf("\nf(x) = ");
		printVector(functionValues, N);

		double* error = (double*)malloc(sizeof(double) * N);
		for (int i = 0; i < N; i++) {
			error[i] = fabs(functionValues[i] - lagrangeValues[i]);
		}
		printf("\nError = ");
		printVector(error, N);
		printVectorToFile(f_err, error, N);
		free(lagrangeValues);
		free(error);
	}
	fclose(f_gridN);
	fclose(f_gridX);
	fclose(f_gridY);
	fclose(f_p);
	fclose(f_err);
	fclose(f_y);
	fclose(f_x);
	fclose(f_theory);
	free(x);
	free(functionValues);
	return 0;
}

double f5Der(double x) {
	return -exp(x);
}

double* calcTheoreticError(double* lagrangeValues, double* x, grid_t grid, int n) {
	double* res = (double*)malloc(sizeof(double) * n);
	if (!res) {
		printf("Erorr allocating memory to theory\n");
		exit(1);
	}
	for (int i = 0; i < n; i++) {
		double root = 1;
		for (int j = 0; j < grid.n; j++) {
			root *= x[i] - grid.gridX[j];
		}
		double fact = 1;
		for (int j = 1; j < grid.n + 2; j++) {
			fact *= j;
		}
		res[i] = fabs(f5Der(4)) * fabs(root) / fact;
	}
	return res;
}

double* createFunctionValues(double* x, int n) {
	double* res = (double*)malloc(sizeof(double) * n);
	if (!res) {
		printf("Error allocating memory to x\n");
		exit(1);
	}
	for (int i = 0; i < n; i++) {
		res[i] = f(x[i]);
	}
	return res;
}

double* createX(double a, double b, int n) {
	double* x = (double*)malloc(sizeof(double) * n);
	if (!x) {
		printf("Error allocating memory to x\n");
		exit(1);
	}
	double tmp = 0;
	for (int i = 0; i < n; i++) {
		tmp = cos((PI * (2 * i + 1)) / (2 * (n + 1)));
		x[i] = ((b - a) * tmp / 2) + ((a + b) / 2);
	}
	return x;
}

double f(double x) {
	return log10(x) + 7.0 / (2 * x + 6);
}

double fDer(double x) {
	return 3 * x * x - exp(x);
}

grid_t createGrid(int n, double a, double b) {
	grid_t grid;
	grid.n = n;
	grid.gridX = (double*)malloc(sizeof(double) * n);
	grid.gridX = createX(a,b,n);
	grid.gridY = (double*)malloc(sizeof(double) * n);
	grid.gridY = createFunctionValues(grid.gridX,n);
	grid.gridDer = (double*)malloc(sizeof(double) * n);
	if (!grid.gridY || !grid.gridX || !grid.gridDer) {
		printf("Error allocating memory to grid\n");
		exit(1);
	}
}

void printVectorToFile(FILE* f, double* vector, int n) {
	for (int i = 0; i < n; i++) {
		fprintf(f, "%.15lf ", vector[i]);
	}
}

void printVector(double* vector, int n) {
	for (int i = 0; i < n; i++) {
		printf("%.15lf ", vector[i]);
	}

	printf("\n");
}

double* createErmitPolynomValues(grid_t grid, double* x, int n) {
	double* polynom = (double*)malloc(n * sizeof(double));
	if (!polynom) {
		printf("Error allocating memory to polynom\n");
		exit(1);
	}

	for (int s = 0; s < n; s++) {
		double sum = 0;
		for (int j = 0; j < grid.n; j++) {
			double innerSum = 0;
			for (int k = 0; k < grid.n; k++) {
				if (j != k) {
					innerSum += (x[s] - grid.gridX[j]) / (grid.gridX[j] - grid.gridX[k]);
				}
			}

			double product = 1;
			for (int i = 0; i < grid.n; i++) {
				if (i != j) {
					product *= (x[s] - grid.gridX[i]) / (grid.gridX[j] - grid.gridX[i]) * (x[s] - grid.gridX[i]) / (grid.gridX[j] - grid.gridX[i]);
				}
			}

			sum += ((x[s] - grid.gridX[j]) * grid.gridDer[j] + (1 - 2 * innerSum) * grid.gridY[j]) * product;
		}

		polynom[s] = sum;
	}

	return polynom;
}

double* Lagrangepolynomial(double* x, grid_t grid, int n)
{
	double* polynom = (double*)malloc(n * sizeof(double));
	if (!polynom) {
		printf("Error allocating memory to polynom\n");
		exit(1);
	}
	double basics_pol = 1;
	double sum = 0;
	for (int s = 0; s < n; s++)
	{
		sum = 0;
		for (int i = 0; i < grid.n; i++)
		{
			basics_pol = 1;
			for (int j = 0; j < grid.n; j++)
			{
				if (j != i)
				{
					basics_pol *= (x[j] - grid.gridX[j]) / (grid.gridX[i] - grid.gridX[j]);
					//printf("%lf ", basics_pol);
				}
			}
			sum += grid.gridY[i] * basics_pol;
		}
		polynom[s] = sum;
	}
	return polynom;
}