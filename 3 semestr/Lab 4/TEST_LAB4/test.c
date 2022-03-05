#pragma once
#pragma warning(disable: 4996)
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#define MATRIX_DIMENSION 3

typedef struct {
	double** elements;
	int m;
	int n;
	double condition_number;
} matrix_t;

matrix_t createMatrix(int m, int n);
void freeMatrix(matrix_t* matrix);
void printMatrix(matrix_t matrix);
void fillDiagonalMatrix(matrix_t* matrix);
matrix_t subtractOrAddMatrix(matrix_t a_matrix, matrix_t b_matrix, int operation);
matrix_t transparentMatrix(matrix_t matrix);
matrix_t multiplyMatrixes(matrix_t a_matrix, matrix_t b_matrix);
double findNorm(matrix_t matrix);
matrix_t getIdentityMatrix(int size);
matrix_t multiplyMatrixByNumber(matrix_t matrix, double number);
double MaxNorm(matrix_t matrix);
double MaxMatrixNorm(matrix_t matrix);
matrix_t getOrtoMatrix();
matrix_t getFinalMatrix(int size);
matrix_t LRMethod(matrix_t matrix, double e, int* iterations);
matrix_t* inverseIteration(matrix_t matrix, double* eigenValues, double eigenValuesKnown, double e, int* iterations);
matrix_t solveSystem(matrix_t A, matrix_t b);
double average(matrix_t x, matrix_t y);
void LUmethod(matrix_t* a, matrix_t* u, matrix_t* l);
void readMatrixes(matrix_t* matrix, matrix_t* eigenValues, matrix_t* eigenVectors);
double ScalarMul(matrix_t v1, matrix_t v2);
double ScalarProductMethod(matrix_t A, double e, int* iterations);
void readVectors(matrix_t* vector1, matrix_t* vector2);
void WriteVector(matrix_t W, matrix_t X);
int main(int argc, char* argv[]) {
	srand((unsigned int)time(NULL));
	matrix_t eigenValuesKnown = createMatrix(MATRIX_DIMENSION, 1);
	double eigenValues;
	double eigenValues2;
	matrix_t afterMethod;
	matrix_t* eigenVectors = NULL;
	matrix_t* eigenVectors2 = NULL;

	FILE* f_eigenValuesNorm = fopen("eigenValuesNorm.txt", "w");
	FILE* f_eigenVectorsNorm = fopen("eigenVectorsNorm.txt", "w");
	
	FILE* f_nevyazka = fopen("nevyazka.txt", "w");
	
	FILE* f_iterations = fopen("iterations.txt", "w");
	FILE* f_epsilon = fopen("epsilon.txt", "w");
	FILE* eigenVector1 = fopen("eigenVector1.txt", "r");
	FILE* eigenVector2 = fopen("eigenVector1.txt", "r");
	if (!f_eigenValuesNorm || !f_eigenVectorsNorm || !f_nevyazka || !f_iterations || !f_epsilon) {
		printf("Error opening files\n");
		exit(1);
	}
	matrix_t matrixMatlab = createMatrix(MATRIX_DIMENSION, MATRIX_DIMENSION);
	matrix_t eigenValuesMatlab = createMatrix(MATRIX_DIMENSION, 1);
	matrix_t* eigenVectorsMatlab = (matrix_t*)malloc(MATRIX_DIMENSION * sizeof(matrix_t));
	if (!eigenVectorsMatlab) {
		printf("Error allocating memory\n");
		exit(1);
	}
	readMatrixes(&matrixMatlab, &eigenValuesMatlab, eigenVectorsMatlab);
	printf("Start matrix\n");
	printMatrix(matrixMatlab);
	int iterations = 0;
	//for (double e = 0.1; e >= pow(10, -14); e /= 10) {
	
	{
		int inverseIterations = 0;
		eigenValues = ScalarProductMethod(matrixMatlab, 0.1, &iterations);



		
		printf("eigenValues %.4lf\n", eigenValues);
		double eigenValuesNorm = fabs(eigenValuesMatlab.elements[1][0] - eigenValues);
		printf("Error %lf", eigenValuesNorm);
		fprintf(f_eigenValuesNorm, "%e ", eigenValuesNorm);
		fprintf(f_iterations, "%d ", iterations);
		
	}
		//printf("EPSILON %e\n", e);
		
		
		//fprintf(f_epsilon, "%e ", e);
	//}
	//

	fclose(f_epsilon);
	fclose(f_iterations);
	fclose(f_eigenValuesNorm);
	fclose(f_eigenVectorsNorm);
	fclose(f_nevyazka);
	fclose(eigenVector1);
	fclose(eigenVector2);
	free(eigenVectorsMatlab);
	free(eigenVectors);
	freeMatrix(&matrixMatlab);
	freeMatrix(&eigenValuesMatlab);
	//freeMatrix(&eigenValues);
	return 0;
}

void readMatrixes(matrix_t* matrix, matrix_t* eigenValues, matrix_t* eigenVectors) {
	FILE* f_matrix = fopen("matrix.txt", "r");
	FILE* f_eigenVectors = fopen("eigenVectorsKnown.txt", "r");
	FILE* f_eigenValues = fopen("eigenValuesKnown.txt", "r");
	if (!f_matrix || !f_eigenValues || !f_eigenVectors) {
		printf("Error opening files\n");
		exit(1);
	}
	printf("Eigenvectorsknown\n");
	for (int i = 0; i < matrix->m; i++) {
		eigenVectors[i] = createMatrix(MATRIX_DIMENSION, 1);
		for (int j = 0; j < matrix->n; j++) {
			if (j == matrix->n - 1) {
				fscanf(f_matrix, "%lf\n", &matrix->elements[i][j]);
				fscanf(f_eigenVectors, "%lf\n", &eigenVectors[i].elements[j][0]);
			}
			else {
				fscanf(f_matrix, "%lf ", &matrix->elements[i][j]);
				fscanf(f_eigenVectors, "%lf ", &eigenVectors[i].elements[j][0]);
			}
		}
		printMatrix(transparentMatrix(eigenVectors[i]));
		fscanf(f_eigenValues, "%lf ", &eigenValues->elements[i][0]);
	}
	fclose(f_matrix);
	fclose(f_eigenVectors);
	fclose(f_eigenValues);
}

matrix_t createMatrix(int m, int n) {
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

void freeMatrix(matrix_t* matrix) {
	for (int i = 0; i < matrix->m; i++) {
		free(matrix->elements[i]);
	}
	free(matrix->elements);
}

void printMatrix(matrix_t matrix) {
	for (int i = 0; i < matrix.m; i++) {
		for (int j = 0; j < matrix.n; j++) {
			printf(" %.4lf ", matrix.elements[i][j]);
		}
		printf("\n");
	}
	printf("\n");
}

matrix_t getIdentityMatrix(int size) {
	matrix_t matrix = createMatrix(size, size);
	for (int i = 0; i < size; i++) {
		matrix.elements[i][i] = 1;
	}
	return matrix;
}

void fillDiagonalMatrix(matrix_t* matrix) {
	for (int i = 0; i < matrix->m; i++) {
		for (int j = i; j < matrix->n; j++) {
			if (i == j) {
				matrix->elements[i][i] = i + 1.0;
			}
			else {
				matrix->elements[i][j] = rand() % 100;
			}
		}
	}
}

matrix_t subtractOrAddMatrix(matrix_t a_matrix, matrix_t b_matrix, int operation) {
	matrix_t res = createMatrix(a_matrix.m, a_matrix.n);
	for (int i = 0; i < res.m; i++) {
		for (int j = 0; j < res.n; j++) {
			res.elements[i][j] = a_matrix.elements[i][j] + operation * b_matrix.elements[i][j];
		}
	}
	return res;
}

matrix_t transparentMatrix(matrix_t matrix) {
	matrix_t res = createMatrix(matrix.n, matrix.m);
	for (int i = 0; i < matrix.n; i++) {
		for (int j = 0; j < matrix.m; j++) {
			res.elements[i][j] = matrix.elements[j][i];
		}
	}
	return res;
}

matrix_t multiplyMatrixes(matrix_t a_matrix, matrix_t b_matrix) {
	matrix_t res = createMatrix(a_matrix.m, b_matrix.n);
	for (int i = 0; i < a_matrix.m; i++) {
		for (int j = 0; j < b_matrix.n; j++) {
			for (int k = 0; k < a_matrix.n; k++) {
				res.elements[i][j] += a_matrix.elements[i][k] * b_matrix.elements[k][j];
			}
		}
	}
	return res;
}

double findNorm(matrix_t matrix) {
	double res = 0;
	for (int i = 0; i < matrix.m; i++) {
		res += matrix.elements[i][0] * matrix.elements[i][0];
	}
	return res;
}

matrix_t multiplyMatrixByNumber(matrix_t matrix, double number) {
	matrix_t res = createMatrix(matrix.m, matrix.n);
	for (int i = 0; i < res.m; i++) {
		for (int j = 0; j < res.n; j++) {
			res.elements[i][j] = number * matrix.elements[i][j];
		}
	}
	return res;
}

double MaxNorm(matrix_t matrix) {
	double res = matrix.elements[0][0];
	for (int i = 1; i < matrix.m; i++) {
		if (fabs(matrix.elements[i][0]) > res) {
			res = fabs(matrix.elements[i][0]);
		}
	}
	return res;
}

double MaxMatrixNorm(matrix_t matrix) {
	double res = 0;
	for (int j = 0; j < matrix.m; j++) {
		for (int i = j + 1; i < matrix.n; i++) {
			if (fabs(matrix.elements[i][j]) > res) {
				res = fabs(matrix.elements[i][j]);
			}
		}
	}
	return res;
}

matrix_t getOrtoMatrix() {
	matrix_t matrix = createMatrix(MATRIX_DIMENSION, 1);
	for (int i = 0; i < MATRIX_DIMENSION; i++) {
		matrix.elements[i][0] = (double)rand();
	}
	return subtractOrAddMatrix(getIdentityMatrix(MATRIX_DIMENSION), multiplyMatrixByNumber(multiplyMatrixes(matrix, transparentMatrix(matrix)), 2 / findNorm(matrix)), -1);
}

matrix_t getFinalMatrix(int size) {
	matrix_t matrix = createMatrix(size, size);
	fillDiagonalMatrix(&matrix);
	matrix_t orto = getOrtoMatrix();
	matrix = multiplyMatrixes(transparentMatrix(orto), matrix);
	return multiplyMatrixes(matrix, orto);
}

matrix_t LRMethod(matrix_t matrix, double e, int* iterations) {
	matrix_t u = createMatrix(matrix.m, matrix.n);
	matrix_t l;
	matrix_t res = matrix;
	do {
		l = getIdentityMatrix(matrix.m);
		LUmethod(&res, &u, &l);
		res = multiplyMatrixes(u, l);
		(*iterations)++;
	} while (MaxMatrixNorm(res) > e);
	/*matrix_t eigenValues = createMatrix(MATRIX_DIMENSION, 1);
	for (int i = 0; i < MATRIX_DIMENSION; i++) {
		eigenValues.elements[i][0] = res.elements[i][i];
	}*/
	freeMatrix(&u);
	freeMatrix(&l);
	return res;
}

void SaveVectors(matrix_t W, matrix_t X, double lambda)
{
	matrix_t W1 = createMatrix(MATRIX_DIMENSION, 1);
	for (int i = 0; i < MATRIX_DIMENSION; i++) {
		W1.elements[i][0] = W.elements[i][0] + lambda*X.elements[i][0];
	}

	double normaW1 = sqrt(findNorm(W1));


	matrix_t W2 = createMatrix(MATRIX_DIMENSION, 1);
	for (int i = 0; i < MATRIX_DIMENSION; i++) {
		W2.elements[i][0] = W.elements[i][0] - lambda*X.elements[i][0];
	}
	double normaW2 = sqrt(findNorm(W2));

	//FILE* eigenVector1 = fopen("eigenVector1.txt", "w+");
	printf("Vector1\n");
	for (int i = 0; i < MATRIX_DIMENSION; i++)
	{
		printf("%.4lf\n", W1.elements[i][0]/normaW1);
	}
	//fclose(eigenVector1);
	printf("Vector2\n");
	//FILE* eigenVector2 = fopen("eigenVector2.txt", "w+");
	for (int i = 0; i < MATRIX_DIMENSION; i++)
	{
		printf("%.4lf\n", W2.elements[i][0]/normaW2);
	}
	//fclose(eigenVector2);

	freeMatrix(&W1);
	freeMatrix(&W2);
}

void WriteVector(matrix_t W, matrix_t X)
{
	matrix_t W1 = createMatrix(MATRIX_DIMENSION, 1);
	for (int i = 0; i < MATRIX_DIMENSION; i++) {
		W1.elements[i][0] = W.elements[i][0] + X.elements[i][0];
	}

	double normaW1 = sqrt(findNorm(W1));


	matrix_t W2 = createMatrix(MATRIX_DIMENSION, 1);
	for (int i = 0; i < MATRIX_DIMENSION; i++) {
		W2.elements[i][0] = W.elements[i][0] - X.elements[i][0];
	}
	double normaW2 = sqrt(findNorm(W2));

	FILE* eigenVector1 = fopen("eigenVector1.txt", "w+");
	for (int i = 0; i < MATRIX_DIMENSION; i++)
	{
		fprintf(eigenVector1, "%.15lf\n", W1.elements[i][0] / normaW1);
	}
	fclose(eigenVector1);

	FILE* eigenVector2 = fopen("eigenVector2.txt", "w+");
	for (int i = 0; i < MATRIX_DIMENSION; i++)
	{
		fprintf(eigenVector2, "%.15lf\n", W2.elements[i][0] / normaW2);
	}
	fclose(eigenVector2);

	freeMatrix(&W1);
}

void WriteVector2(matrix_t W, matrix_t X, double lambda)
{
	matrix_t W1 = createMatrix(MATRIX_DIMENSION, 1);
	for (int i = 0; i < MATRIX_DIMENSION; i++) {
		W1.elements[i][0] = W.elements[i][0] + lambda*X.elements[i][0];
	}

	double normaW1 = sqrt(findNorm(W1));


	matrix_t W2 = createMatrix(MATRIX_DIMENSION, 1);
	for (int i = 0; i < MATRIX_DIMENSION; i++) {
		W2.elements[i][0] = W.elements[i][0] - lambda*X.elements[i][0];
	}
	double normaW2 = sqrt(findNorm(W2));

	FILE* eigenVector1 = fopen("eigenVector1.txt", "w+");
	for (int i = 0; i < MATRIX_DIMENSION; i++)
	{
		fprintf(eigenVector1, "%.15lf\n", W1.elements[i][0] / normaW1);
	}
	fclose(eigenVector1);

	FILE* eigenVector2 = fopen("eigenVector2.txt", "w+");
	for (int i = 0; i < MATRIX_DIMENSION; i++)
	{
		fprintf(eigenVector2, "%.15lf\n", W2.elements[i][0] / normaW2);
	}
	fclose(eigenVector2);

	freeMatrix(&W1);
}


double coverge(matrix_t A, matrix_t W, double lamda, int size, matrix_t X) {
	matrix_t AW1_res;
	matrix_t W1 = createMatrix(size, 1);
	for (int i = 0; i < size; i++) {
		W1.elements[i][0] = W.elements[i][0] + X.elements[i][0];
	}

	AW1_res = multiplyMatrixes(A, W1);
	double normaW1 = sqrt(findNorm(W1));

	//меняем W1
	for (int i = 0; i < size; i++) {
		W1.elements[i][0] = W1.elements[i][0] * lamda;
	}
	matrix_t AW1_lamdaW1 = subtractOrAddMatrix(AW1_res, W1, -1);
	double normaAW1_lamdaW1 = sqrt(findNorm(AW1_lamdaW1));
	freeMatrix(&AW1_res);
	freeMatrix(&W1);
	freeMatrix(&AW1_lamdaW1);
	/*printf("chislitel %lf\n", normaAW1_lamdaW1);
	printf("znam %lf\n", normaW1);*/

	return normaAW1_lamdaW1 / normaW1;

}


double ScalarMul(matrix_t v1, matrix_t v2)
{
	double result = 0;
	for (int i = 0; i < v1.m; i++)
	{
		result += v1.elements[i][0] * v2.elements[i][0];
	}
	return result;
}

double ScalarProductMethod(matrix_t A, double e, int* iterations)
{
	double prev_lambda = 0;
	double cur_lambda = 0;
	matrix_t w_prev = createMatrix(MATRIX_DIMENSION, 1);
	matrix_t w_prevprev = createMatrix(MATRIX_DIMENSION, 1);
	double scal_prod;
	matrix_t w_cur = createMatrix(MATRIX_DIMENSION, 1);
	for (int i = 0; i < w_prev.m; i++) {
		w_cur.elements[i][0] = 0.5774;
		w_prev.elements[i][0] = 0.5774;
		w_prevprev.elements[i][0] = 0.5774;
	}
	matrix_t tmp;
	matrix_t tmp_next;
	matrix_t tmp_next_next;
	double norm = sqrt(findNorm(w_cur));
	double isCoverge = 100;
	while ((*iterations) != 3)
	{
		printf("Iterations %d\n", *iterations);
		
		prev_lambda = cur_lambda;
		cur_lambda = 0;
		//tmp = w_cur;
		tmp = multiplyMatrixes(A, w_cur);
		printf("After Ax\n");
		printMatrix(tmp);
		tmp_next = multiplyMatrixes(A, tmp);
		printf("After Ay\n");
		printMatrix(tmp_next);
		tmp_next_next = multiplyMatrixes(A, tmp_next);
		printf("Norm y %.4lf\n", sqrt(findNorm(tmp)));
		//freeMatrix(&w_cur);
		cur_lambda = sqrt(ScalarMul(w_cur, tmp_next) / ScalarMul(w_cur, w_cur));
		for (int i = 0; i < MATRIX_DIMENSION; i++)
		{
			w_cur.elements[i][0] = tmp.elements[i][0] / sqrt(findNorm(tmp));
		}
		printf("x_cur\n");
		printMatrix(w_cur);


		//cur_lambda = ScalarMul(tmp, tmp_next_next);

		//for (int i = 0; i < MATRIX_DIMENSION; i++)
		//{
		//	w_prevprev.elements[i][0] = w_prev.elements[i][0];
		//	w1.elements[i][0] = tmp_next.elements[i][0] + cur_lambda * tmp.elements[i][0];
		//	w2.elements[i][0] = tmp_next.elements[i][0] - cur_lambda * tmp.elements[i][0];
		//	/*w1.elements[i][0] = tmp_next_next.elements[i][0] + cur_lambda * tmp_next.elements[i][0];
		//	w2.elements[i][0] = tmp_next_next.elements[i][0] - cur_lambda * tmp_next.elements[i][0];*/
		//}
		/*printf("After norm\n");
		printMatrix(w_cur);
		for (int i = 0; i < MATRIX_DIMENSION; i++)
		{
			w_prevprev.elements[i][0] = w_prev.elements[i][0];
			w_prev.elements[i][0] = w_cur.elements[i][0];
		}
		matrix_t W1 = createMatrix(MATRIX_DIMENSION, 1);
		for (int i = 0; i < MATRIX_DIMENSION; i++) {
			W1.elements[i][0] = w_prev.elements[i][0] + cur_lambda*w_prevprev.elements[i][0];
		}*/
		/*printf("w1 \n");
		printMatrix(W1);

		printf("w2 \n");
		for (int i = 0; i < MATRIX_DIMENSION; i++) {
			W1.elements[i][0] = w_prev.elements[i][0] - cur_lambda * w_prevprev.elements[i][0];
		}
		printMatrix(W1);*/
		printf("Lambda %.4lf\n", cur_lambda);
		SaveVectors(tmp_next, tmp, cur_lambda);
		printf("Error %.4lf\n", fabs(cur_lambda - prev_lambda));
		(*iterations)++;
	}
		
	
	//WriteVector(w_prev, w_prevprev);
	freeMatrix(&w_prev);
	freeMatrix(&w_prevprev);
	freeMatrix(&w_cur);

	printf("Iter: %d\n", *iterations);
	return cur_lambda;
}

void LUmethod(matrix_t* a, matrix_t* u, matrix_t* l) {
	for (int i = 0; i < a->m; i++) {
		for (int j = 0; j < a->n; j++) {
			if (i <= j) {
				double sum = 0;
				for (int k = 0; k < i; k++) {
					sum += l->elements[i][k] * u->elements[k][j];
				}
				u->elements[i][j] = a->elements[i][j] - sum;
			}
			else {
				double sum = 0;
				for (int k = 0; k < j; k++) {
					sum += l->elements[i][k] * u->elements[k][j];
				}
				l->elements[i][j] = (a->elements[i][j] - sum) / u->elements[j][j];
			}
		}
	}
}

matrix_t* inverseIteration(matrix_t matrix, double* eigenValues, double eigenValuesKnown, double e, int* iterations) {
	matrix_t* eigenVectors = (matrix_t*)malloc(1 * sizeof(matrix_t));
	if (!eigenVectors) {
		printf("Error allocating memory\n");
		exit(1);
	}
	matrix_t A;
	matrix_t y;
	matrix_t x_0 = createMatrix(matrix.m, 1);
	x_0.elements[0][0] = 1;

	eigenVectors[0] = createMatrix(matrix.m, 1);
	int j = 0;
	do {
		A = subtractOrAddMatrix(matrix, multiplyMatrixByNumber(getIdentityMatrix(matrix.m), *eigenValues), -1);
		if (j == 0) {
			y = solveSystem(A, x_0);
			*eigenValues += average(x_0, y);
		}
		else {
			y = solveSystem(A, eigenVectors[0]);
			*eigenValues += average(eigenVectors[0], y);
		}
		eigenVectors[0] = multiplyMatrixByNumber(y, 1 / sqrt(findNorm(y)));
		j++;

		(*iterations)++;

	} while (fabs(*eigenValues - eigenValuesKnown) > e);

	freeMatrix(&y);
	return eigenVectors;
}

matrix_t solveSystem(matrix_t A, matrix_t b) {
	matrix_t y = createMatrix(b.m, b.n);
	matrix_t x = createMatrix(b.m, b.n);
	matrix_t u = createMatrix(A.m, A.n);
	matrix_t l = getIdentityMatrix(A.m);
	LUmethod(&A, &u, &l);
	for (int i = 0; i < y.m; i++) {
		double sum = 0;
		for (int k = 0; k < i; k++) {
			sum += l.elements[i][k] * y.elements[k][0];
		}
		y.elements[i][0] = b.elements[i][0] - sum;
	}
	for (int i = x.m - 1; i >= 0; i--) {
		double sum = 0;
		for (int k = i + 1; k < x.m; k++) {
			sum += u.elements[i][k] * x.elements[k][0];
		}
		x.elements[i][0] = (y.elements[i][0] - sum) / u.elements[i][i];
	}
	freeMatrix(&u);
	freeMatrix(&l);
	freeMatrix(&y);
	return x;
}

double average(matrix_t x, matrix_t y) {
	double res = 0;
	for (int i = 0; i < x.m; i++) {
		res += x.elements[i][0] / y.elements[i][0];
	}
	return res / x.m;
}