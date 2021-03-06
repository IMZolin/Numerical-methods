#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <locale.h>
#define MATRIX_DIMENSION 10
#pragma warning(disable : 4996)

typedef struct {
	double** elements;
	int m;
	int n;
} matrix_t;

matrix_t CreateMatrix(int m, int n);
void DestroyMatrix(matrix_t* matrix);
void PrintMatrix(matrix_t matrix);
matrix_t GetIdentityMatrix(int size);
void FillDiagonalMatrix(matrix_t* matrix);
matrix_t SubtractOrAddMatrix(matrix_t a_matrix, matrix_t b_matrix, int operation);
matrix_t TransparentMatrix(matrix_t matrix);
matrix_t MultiplyMatrixes(matrix_t a_matrix, matrix_t b_matrix);
double SquaredNorm(matrix_t matrix);
double MaxNorm(matrix_t matrix);
matrix_t MultiplyMatrixByNumber(matrix_t matrix, double number);
matrix_t GetOrtoMatrix();
matrix_t GetMatrix(int size);
matrix_t GetInvMatrix(matrix_t matrix);
double InfNorm(matrix_t matrix);
double Get_m(matrix_t matrix);
matrix_t SeidelMethod(matrix_t A, matrix_t b, double epsilon, int* counter);

int main(void)
{
	setlocale(LC_CTYPE, "Russian");
	srand((unsigned int)time(NULL));

	FILE* f_iterations = fopen("iterations.txt", "wt");
	FILE* f_actual_error = fopen("actual_error.txt", "wt");
	FILE* f_discrepancy = fopen("discrepancy.txt", "wt");
	FILE* f_epsilon = fopen("epsilon.txt", "wt");

	double C_norm, alpha;
	matrix_t A;
	int counter = 0;
	double m;
	do
	{
		if (counter > 0)
			DestroyMatrix(&A);
		A = GetMatrix(MATRIX_DIMENSION);
		alpha = 1.0 / MaxNorm(A);

		matrix_t E = GetIdentityMatrix(MATRIX_DIMENSION);
		matrix_t C = SubtractOrAddMatrix(E, MultiplyMatrixByNumber(A, alpha), -1);
		C_norm = sqrt(SquaredNorm(C));
		m = Get_m(C);
		DestroyMatrix(&E);
		DestroyMatrix(&C);
		counter++;
	} while (C_norm >= 1);
	printf("%lf %lf\n\n", C_norm, alpha);
	PrintMatrix(A);
	for (double epsilon = 0.001; epsilon > 0.000000000001; epsilon *= 0.1) //?? 10^-3 ?? 10^-12 ? ????? 10^-1
	{

		matrix_t x_exact = CreateMatrix(MATRIX_DIMENSION, 1);

		for (int i = 0; i < x_exact.m; i++)
		{
			x_exact.elements[i][0] = i + 1.0;
		}
		matrix_t b = MultiplyMatrixes(A, x_exact);

		matrix_t x = SeidelMethod(A, b, epsilon, &counter);
		PrintMatrix(x);

		double actual_error = InfNorm(SubtractOrAddMatrix(x, x_exact, -1));   // ????? ??????????? ??????
		double discrepancy = InfNorm(SubtractOrAddMatrix(MultiplyMatrixes(A, x), b, -1));  // ????? ???????

		fprintf(f_epsilon, "%e ", epsilon);
		fprintf(f_iterations, "%d ", counter);
		fprintf(f_actual_error, "%e ", actual_error);
		fprintf(f_discrepancy, "%e ", discrepancy);
		printf("\n\n");
	}

	fclose(f_iterations);
	fclose(f_actual_error);
	fclose(f_discrepancy);
	fclose(f_epsilon);

	return 0;
}

matrix_t CreateMatrix(int m, int n) {
	matrix_t matrix;
	matrix.m = m;
	matrix.n = n;
	matrix.elements = (double**)calloc(m, sizeof(double*));
	if (matrix.elements == NULL) {
		printf("Erorr allocating memory\n");
		exit(1);
	}
	for (int i = 0; i < m; i++) {
		matrix.elements[i] = (double*)calloc(n, sizeof(double));
		if (matrix.elements[i] == NULL) {
			printf("Error allocating memory\n");
			exit(1);
		}
	}
	return matrix;
}

void DestroyMatrix(matrix_t* matrix) {
	for (int i = 0; i < matrix->m; i++) {
		free(matrix->elements[i]);
	}
	free(matrix->elements);
}

void PrintMatrix(matrix_t matrix)
{
	for (int i = 0; i < matrix.m; i++)
	{
		for (int j = 0; j < matrix.n; j++)
		{
			if (matrix.elements[i][j] > 0)
				printf(" ");
			printf(" %.15lf ", matrix.elements[i][j]);
		}
		printf("\n");
	}
	printf("\n");
}

matrix_t GetIdentityMatrix(int size) {
	matrix_t matrix = CreateMatrix(size, size);
	for (int i = 0; i < size; i++) {
		matrix.elements[i][i] = 1;
	}
	return matrix;
}

void FillDiagonalMatrix(matrix_t* matrix) {
	for (int i = 0; i < matrix->m; i++) {
		if (i == 0) {
			matrix->elements[i][i] = 1;
		}
		else if (i == matrix->m - 1) {
			matrix->elements[i][i] = 10;
		}
		else {
			double elem = (double)rand();
			while (elem > 1) {
				elem /= 10;
			}
			matrix->elements[i][i] = elem * 10;
		}
	}
}

matrix_t SubtractOrAddMatrix(matrix_t a_matrix, matrix_t b_matrix, int operation) {
	matrix_t res = CreateMatrix(a_matrix.m, a_matrix.n);
	for (int i = 0; i < res.m; i++) {
		for (int j = 0; j < res.n; j++) {
			res.elements[i][j] = a_matrix.elements[i][j] + operation * b_matrix.elements[i][j];
		}
	}
	return res;
}

matrix_t TransparentMatrix(matrix_t matrix) {
	matrix_t res = CreateMatrix(matrix.n, matrix.m);
	for (int i = 0; i < matrix.n; i++) {
		for (int j = 0; j < matrix.m; j++) {
			res.elements[i][j] = matrix.elements[j][i];
		}
	}
	return res;
}

matrix_t MultiplyMatrixes(matrix_t a_matrix, matrix_t b_matrix) {
	matrix_t res = CreateMatrix(a_matrix.m, b_matrix.n);
	for (int i = 0; i < a_matrix.m; i++) {
		for (int j = 0; j < b_matrix.n; j++) {
			for (int k = 0; k < a_matrix.n; k++) {
				res.elements[i][j] += a_matrix.elements[i][k] * b_matrix.elements[k][j];
			}
		}
	}
	return res;
}

matrix_t MultiplyMatrixByNumber(matrix_t matrix, double number) {
	matrix_t res = CreateMatrix(matrix.m, matrix.n);
	for (int i = 0; i < res.m; i++) {
		for (int j = 0; j < res.n; j++) {
			res.elements[i][j] = number * matrix.elements[i][j];
		}
	}
	return res;
}

double SquaredNorm(matrix_t matrix) {
	double res = 0;
	for (int i = 0; i < matrix.m; i++) {
		res += matrix.elements[i][0] * matrix.elements[i][0];
	}
	return res;
}

double MaxNorm(matrix_t matrix) {
	double res = matrix.elements[0][0];
	for (int i = 1; i < matrix.m; i++) {
		if (matrix.elements[i][0] > res) {
			res = matrix.elements[i][0];
		}
	}
	return res;
}
matrix_t GetOrtoMatrix() {
	matrix_t matrix = CreateMatrix(MATRIX_DIMENSION, 1);
	for (int i = 0; i < MATRIX_DIMENSION; i++) {
		matrix.elements[i][0] = (double)rand();
	}
	return SubtractOrAddMatrix(GetIdentityMatrix(MATRIX_DIMENSION), MultiplyMatrixByNumber(MultiplyMatrixes(matrix, TransparentMatrix(matrix)), 2 / SquaredNorm(matrix)), -1);
}

matrix_t GetMatrix(int size) {
	matrix_t matrix = CreateMatrix(size, size);
	FillDiagonalMatrix(&matrix);
	matrix_t orto = GetOrtoMatrix();
	matrix = MultiplyMatrixes(orto, matrix);
	return MultiplyMatrixes(matrix, TransparentMatrix(orto));
}




matrix_t filMatrix_X_k_1(matrix_t Matrix, matrix_t Matrix2, int N) {
	FILE* f;
	for (int i = 0; i < N; i++) {
		Matrix.elements[i][0] = Matrix2.elements[i][0];
	}
	return Matrix;
}

int converge(matrix_t xk, matrix_t xkp, int n, double eps, double t)
{
	double norma = 0;
	for (int i = 0; i < n; i++)
		norma += (xk.elements[i][0] - xkp.elements[i][0]) * (xk.elements[i][0] - xkp.elements[i][0]);
	return (sqrt(norma) < eps * (1 - t) / t);
}
int convergeNorm1(matrix_t xk, matrix_t xkp, int n, double eps, double t)
{
	double norma = 0;
	for (int i = 0; i < n; i++)
		norma += fabs(xk.elements[i][0] - xkp.elements[i][0]);
	return (norma < eps* (1 - t) / t);
}

double InfNorm(matrix_t matrix)
{
	double norm = 0;

	for (int i = 0; i < matrix.m; i++)
	{
		double sum = 0;

		for (int j = 0; j < matrix.n; j++)
		{
			sum += fabs(matrix.elements[i][j]);
		}

		if (i == 0)
			norm = sum;
		else if (sum > norm)
			norm = sum;
	}

	return norm;
}

double Get_m(matrix_t matrix)
{
	double alpha;
	double beta;
	double m = 0;

	for (int i = 0; i < matrix.m; i++)
	{
		alpha = 1.0 / MaxNorm(matrix);

		beta = 0;
		for (int j = i; j < matrix.m; j++)
		{
			beta += fabs(matrix.elements[i][j]);
		}

		if (i == 0)
			m = beta / (1 - alpha);
		else if (beta / (1 - alpha) > m)
			m = beta / (1 - alpha);

	}

	return m;
}

matrix_t SeidelMethod(matrix_t A, matrix_t b, double epsilon, int* counter)
{
	*counter = 0;
	double alpha = 1.0 / InfNorm(A);

	matrix_t g = MultiplyMatrixByNumber(b, alpha);
	matrix_t E = GetIdentityMatrix(MATRIX_DIMENSION);
	matrix_t C = SubtractOrAddMatrix(E, MultiplyMatrixByNumber(A, alpha), -1);
	matrix_t x = CreateMatrix(MATRIX_DIMENSION, 1);
	matrix_t xk = CreateMatrix(MATRIX_DIMENSION, 1);
	for (int i = 0; i < xk.m; i++)
		xk.elements[i][0] = 100.0;
	double C_norm = MaxNorm(C);
	double m = Get_m(C);
	//printf("Norm of C:%lf\n", C_norm);
	//printf("Alfa %lf\n", alpha);
	//printf("Mu is %lf\n", m);
	do
	{
		for (int i = 0; i < xk.m; i++)
		{
			x.elements[i][0] = xk.elements[i][0];
		}

		for (int i = 0; i < x.m; i++)
		{
			xk.elements[i][0] = g.elements[i][0];
			int t = i;
			for (int j = 0; j < C.n; j++)
			{
				if (t > 0)
				{
					xk.elements[i][0] += C.elements[i][j] * xk.elements[j][0];  
					//?????? i - t ????? ?????? j
					t--;
				}
				else
				{
					xk.elements[i][0] += C.elements[i][j] * x.elements[j][0];
				}

			}

		}
		(*counter)++;
		//printf("%e\n", epsilon * (1.0 - C_norm) / C_norm);

	} while (InfNorm(SubtractOrAddMatrix(xk, x, -1)) > epsilon * (1.0 - C_norm) / C_norm);

	return xk;
}


