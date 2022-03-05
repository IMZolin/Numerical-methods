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
void FillDiagonalMatrix(matrix_t* matrix, double condition_number);
matrix_t SubtractOrAddMatrix(matrix_t a_matrix, matrix_t b_matrix, int operation);
matrix_t TransparentMatrix(matrix_t matrix);
matrix_t MultiplyMatrixes(matrix_t a_matrix, matrix_t b_matrix);
double Norm(matrix_t matrix);
double MaxNorm(matrix_t matrix);
matrix_t MultiplyMatrixByNumber(matrix_t matrix, double number);
matrix_t GetOrtoMatrix();
matrix_t GetMatrix(double condition_number, int size);
matrix_t RotationMethod(matrix_t A, matrix_t b);

int main(void)
{
	setlocale(LC_CTYPE, "Russian");
	srand((unsigned int)time(NULL));

	FILE* f_condition = fopen("condition.txt", "wt");
	FILE* f_actual_error = fopen("actual_error.txt", "wt");
	FILE* f_discrepancy = fopen("discrepancy.txt", "wt");
	FILE* f_delta_b = fopen("delete_b.txt", "wt");
	FILE* f_y = fopen("y.txt", "wt");

	for (double cond_number = 10; cond_number < pow(10, 9); cond_number *= 10)
	{
		cond_number = (double)(int)cond_number;
		matrix_t A = GetMatrix(cond_number, MATRIX_DIMENSION);
		//printf("-----\n");
		//PrintMatrix(A);

		matrix_t x_known = CreateMatrix(MATRIX_DIMENSION, 1);
		for (int i = 0; i < x_known.m; i++)
		{
			x_known.elements[i][0] = i + 1.0;
		}
		matrix_t b = MultiplyMatrixes(A, x_known);
		//PrintMatrix(TransparentMatrix(b));
		matrix_t mul = MultiplyMatrixes(A, x_known);
		//PrintMatrix(TransparentMatrix(x_known));
		
		matrix_t x = RotationMethod(A, b);
	
		//PrintMatrix(TransparentMatrix(mul));

		double actual_error = MaxNorm(SubtractOrAddMatrix(x, x_known, -1));
		double nevyazka = MaxNorm(SubtractOrAddMatrix(MultiplyMatrixes(A, x), b, -1));
		matrix_t delta_b = CreateMatrix(b.m, b.n);
	for (int i = 0; i < delta_b.m; i++) {
		delta_b.elements[i][0] = b.elements[i][0] + (rand() % 100 + 1.0) / 100.0;
	}
	matrix_t y = RotationMethod(A, b);
		printf("\n");
		

		fprintf(f_condition, "%d ", (int)cond_number);
		fprintf(f_actual_error, "%e ", actual_error);
		fprintf(f_discrepancy, "%e ", nevyazka);
		for (int i = 0; i < delta_b.m; i++)
		{
			fprintf(f_delta_b, "%lf", delta_b.elements[i]);
		}
		
		/*printf("Condition number %.0lf\n\n x* = ", cond_number);
		PrintMatrix(TransparentMatrix(x_known));
		printf("x = ");
		PrintMatrix(TransparentMatrix(x));
		printf("y = ");
		PrintMatrix(TransparentMatrix(y));
		printf("%.15lf <= %.15lf\n\n",Norm(SubtractOrAddMatrix(x, y, -1)) / Norm(x), cond_number * Norm(delta_b) / Norm(b));*/
		//printf("x = ");
		//fprintf("%d",f_delta_b);
		//PrintMatrix(TransparentMatrix(x));
		//printf("y = ");
		//PrintMatrix(TransparentMatrix(y));
		DestroyMatrix(&A);
		DestroyMatrix(&x);
		DestroyMatrix(&x_known);
		DestroyMatrix(&b);
	}

	fclose(f_condition);
	fclose(f_actual_error);
	fclose(f_discrepancy);
	fclose(f_delta_b);
	fclose(f_y);
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
			printf(" %.1lf ", matrix.elements[i][j]);
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

void FillDiagonalMatrix(matrix_t* matrix, double condition_number) {
	for (int i = 0; i < matrix->m; i++) {
		if (i == 0) {
			matrix->elements[i][i] = 1;
		}
		else if (i == matrix->m - 1) {
			matrix->elements[i][i] = condition_number;
		}
		else {
			double elem = (double)rand();
			while (elem > 1) {
				elem /= 10;
			}
			matrix->elements[i][i] = elem * condition_number;
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

double Norm(matrix_t matrix) {
	double res = 0;
	for (int i = 0; i < matrix.m; i++) {
		res += matrix.elements[i][0] * matrix.elements[i][0];
	}
	return res;
}

double MaxNorm(matrix_t matrix) {
	double res = fabs(matrix.elements[0][0]);
	for (int i = 1; i < matrix.m; i++) {
		if (fabs(matrix.elements[i][0]) > res) {
			res = fabs(matrix.elements[0][0]);
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

matrix_t GetOrtoMatrix() {
	matrix_t matrix = CreateMatrix(MATRIX_DIMENSION, 1);
	for (int i = 0; i < MATRIX_DIMENSION; i++) {
		matrix.elements[i][0] = (double)rand();
	}
	return SubtractOrAddMatrix(GetIdentityMatrix(MATRIX_DIMENSION), MultiplyMatrixByNumber(MultiplyMatrixes(matrix, TransparentMatrix(matrix)), 2 / Norm(matrix)), -1);
}

matrix_t GetMatrix(double condition_number, int size) {
	matrix_t matrix = CreateMatrix(size, size);
	FillDiagonalMatrix(&matrix, condition_number);
	matrix_t orto = GetOrtoMatrix();
	matrix = MultiplyMatrixes(orto, matrix);
	return MultiplyMatrixes(matrix, TransparentMatrix(orto));
}

matrix_t RotationMethod(matrix_t A, matrix_t B)
{
	matrix_t X = CreateMatrix(B.m, B.n);
	double C, S = 0;
	double tmp = 0;
	for (int i = 0; i < MATRIX_DIMENSION; i++)
	{
		for (int j = i + 1; j < MATRIX_DIMENSION; j++)
		{
			C = A.elements[i][i] / sqrt(A.elements[i][i] * A.elements[i][i] + A.elements[j][i] * A.elements[j][i]);
			S = A.elements[j][i] / sqrt(A.elements[i][i] * A.elements[i][i] + A.elements[j][i] * A.elements[j][i]);
			A.elements[i][i] = C * A.elements[i][i] + S * A.elements[j][i];
			for (int k = i + 1; k < MATRIX_DIMENSION; k++)
			{
				tmp = A.elements[i][k];
				A.elements[i][k] = C * A.elements[i][k] + S * A.elements[j][k];
				A.elements[j][k] = -S * tmp + C * A.elements[j][k];
			}
			tmp = B.elements[i][0];
			B.elements[i][0] = C * B.elements[i][0] + S * B.elements[j][0];
			B.elements[j][0] = -S * tmp + C * B.elements[j][0];
		}
	}
	//PrintMatrix(A);
	for (int i = MATRIX_DIMENSION - 1; i >= 0; i--)
	{
		double right_sum = 0.0;
		for (int j = i + 1; j <= MATRIX_DIMENSION - 1; j++)
		{
			right_sum += A.elements[i][j] * X.elements[j][0];
		}
		X.elements[i][0] = (B.elements[i][0] - right_sum) / A.elements[i][i];
	}
	return X;
}
