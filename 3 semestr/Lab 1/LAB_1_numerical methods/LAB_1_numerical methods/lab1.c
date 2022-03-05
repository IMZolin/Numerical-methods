#include<stdio.h>
#include<math.h>
#include<conio.h>
//Основные функции 
double Polynomial(double x);
double Transcendental(double x);
void BisectionMethod(double (*func)(double x), double a, double b, double e);
void FixedPointIterations(double (*func)(double x), double a, double b, double q, double e);
double PolynomialFi(double x);
double TranscendentalFi(double x);
void Run();
//Для вывода в Матлаб
void IterationsFixedPointIterations(double (*fi)(double x), double a, double b, double q, double e);
void IterationsBisectionMethod(double (*func)(double x), double a, double b, double e);
void RunIterationsBisectionMethod(double (*func)(double x), double a, double b);
void RunIterationsFixedPointIterations(double (*func)(double x), double a, double b, double q);
void RootsBisectionMethod(double (*func)(double x), double a, double b, double e);
void RootsFixedPointIterations(double (*fi)(double x), double a, double b, double q, double e);
void RunRootsBisectionMethod(double(*func)(double x), double a, double b);
void RunRootsFixedPointIterations(double (*fi)(double x), double a, double b, double q);
void BisectionMethodFunction(double(*func)(double x), double a, double b);
void RunForMatlab();

int main(int argc, char* argv[])
{
	//Run();
	RunForMatlab();
	return 0;
}

/*
*@brief Алгебраическая функция
*@param x число элементов в массиве
*/
double Polynomial(double x)
{
	return 2 * pow(x, 4) - pow(x, 2) - 10;
}

/*
*@brief Трансцедентная функция
*@param x число элементов в массиве
*/
double Transcendental(double x)
{
	return x + log10(1 + x) - 1.5;
}

/*
*@brief Функция Фи алгебраической функции
*@param x
*/
double PolynomialFi(double x)
{
	// M1 = max(f'(x)) = 60
	//alfa = 2/M1 = 1/30
	//func_fi = x - 1/30*func(x) (func_fi = x - alfa*f(x))
	//q = 0.8
	return -2 * pow(x, 4) / 30 + pow(x, 2) / 30 + x + 1.0 / 3;
}

/*
*@brief Функция Фи трансцедентной функции
*@param x
*/
double TranscendentalFi(double x)
{
	//M1 = max(f'(a), f'(b)) = f'(a) т.к. f'(x) убывает
	//M1 = 1.255
	//alfa = 2/M1 = 1.593
	//func_fi = x - 1.643*func(x) (func_fi = x - alfa*func(x))
	//q = -0.593 -1.593/(ln(10)*(x+1) (q = 1 - alfa*f'(x))
	//return -0.643 * x - 1.643 * log10(1 + x) + 2.465;
	return -0.59303 * x - 1.59303 * log10(1 + x) + 2.3895;
}

double Testfunc(double x)
{
	return pow(x, 2) + 10 * x;
}

double TestfuncFi(double x)
{
	return -1 * pow(x, 2) / 5 - x / 3;
}
 
/*
*@brief Метод половинного деления
*@param f(x) функция
*@param a левая граница интервала (а, b)
*@param b чправая граница интервала (а, b)
*@param e машинное эпсилон
*/
void BisectionMethod(double (*func)(double x), double a, double b, double e)
{
	double c;
	int i = 0;
	while (fabs(b - a) > 2 * e)
	{
		c = (a + b) / 2;
		if (func(a) * func(c) < 0)
		{
			b = c;
		}
		else
		{
			a = c;
		}
		i++;
	}
	printf("Root : %.15lf; Iterations: %d\n", (a + b) / 2, i);
	//printf("\nIterations: %d\n", i);
}

/*
*@brief Метод простых итераций
*@param f(x) функция
*@param a левая граница интервала (а, b)
*@param b чправая граница интервала (а, b)
*@param e машинное эпсилон
*/
void FixedPointIterations(double (*fi)(double x), double a, double b, double q, double e)
{
	int i = 0;
	double x_k_1;
	double x_k = a;
	do
	{
		x_k_1 = x_k;
		x_k = fi(x_k_1);
		i++;
	} while (fabs(x_k - x_k_1) >= (1 - q) / q * e);
	printf("Root : %.15lf; Iterations: %d\n", x_k, i);
}

//Обычный вывод результатов
void Run()
{
	double e = 0.1;
	while (e > 0.000001)
	{
		printf("\nEpsilon: %lf", e);
		printf("\nTHE BISECTION METHOD:\n");
		printf("Algebraic equation:\n");
		BisectionMethod(&Polynomial, 1, 2, e);
		printf("Transcendental equation:\n");
		BisectionMethod(&Transcendental, 0.7, 1.7, e);
		printf("\nFIXED-POINT ITERATIONS METHOD\n");
		printf("Algebraic equation:\n");
		FixedPointIterations(&PolynomialFi, 1, 2, 0.8, e);
		printf("Transcendental equation:\n");
		FixedPointIterations(&TranscendentalFi, 0.7, 1.7, 0.85, e);
		e /= 10;
	}
}

//Для вывода в Матлаб

//Итерации метода половинного деления
void IterationsBisectionMethod(double (*func)(double x), double a, double b, double e)
{
	double c;
	int i = 0;
	while (fabs(b - a) > 2 * e)
	{
		c = (a + b) / 2;
		if (func(a) * func(c) < 0)
		{
			b = c;
		}
		else
		{
			a = c;
		}
		i++;
	}
	printf(" %d", i);
}

//Вывод итераций метода простых итераций
void RunIterationsBisectionMethod(double (*func)(double x), double a, double b)
{
	double e = 0.1;
	printf("\n[");

	while (e > 0.000001)
	{
		IterationsBisectionMethod(func, a, b, e);
		e /= 10;
	}
	printf("];");
}

//Итерации метода простых итераций
void IterationsFixedPointIterations(double (*fi)(double x), double a, double b, double q, double e)
{
	int i = 0;
	double x_k_1;
	double x_k = a;
	do
	{
		x_k_1 = x_k;
		x_k = fi(x_k_1);
		i++;
	} while (fabs(x_k - x_k_1) >= (1 - q) / q * e);
	printf(" %d", i);
}

//Вывод итераций метода половинного деления
void RunIterationsFixedPointIterations(double (*fi)(double x), double a, double b, double q)
{
	double e = 0.1;
	printf("\n[");

	while (e > 0.000001)
	{
		IterationsFixedPointIterations(fi, a, b, q, e);
		e /= 10;
	}
	printf("];");
}

//Корни метода половинного деления
void RootsBisectionMethod(double (*func)(double x), double a, double b, double e)
{
	double c;
	int i = 0;
	while (fabs(b - a) > 2 * e)
	{
		c = (a + b) / 2;
		if (func(a) * func(c) < 0)
		{
			b = c;
		}
		else
		{
			a = c;
		}
		i++;
	}
	printf(" %.15lf", (a + b) / 2);
}

//Корни метода простых итераций
void RootsFixedPointIterations(double (*fi)(double x), double a, double b, double q, double e)
{
	int i = 0;
	double x_k_1;
	double x_k = a;
	do
	{
		x_k_1 = x_k;
		x_k = fi(x_k_1);
		i++;
	} while (fabs(x_k - x_k_1) >= (1 - q) / q * e);
	printf(" %.15lf", x_k);
}

//Вывод корни метода половинного деления
void RunRootsBisectionMethod(double(*func)(double x), double a, double b)
{
	double e = 0.1;
	printf("\n[");

	while (e > 0.000001)
	{
		RootsBisectionMethod(func, a, b, e);
		e /= 10;
	}
	printf("];");
}
//Вывод корни метода половинного деления
void RunRootsFixedPointIterations(double (*fi)(double x), double a, double b, double q)
{
	double e = 0.1;
	printf("\n[");

	while (e > 0.000001)
	{
		RootsFixedPointIterations(fi, a, b, q, e);
		e /= 10;
	}
	printf("];");
}

void AproximationsBisectionMethod(double(*func)(double x), double a, double b)
{
	double e = 0.1;
	for (int i = 0; i < 5; i++)
	{
		double c;
		int i = 0;
		while (fabs(b - a) > 2 * e)
		{
			c = (a + b) / 2;
			//printf("(%.15lf)", (b - a) / 2.0);
			printf(" %.15lf ", (a + b) / 2);
			//printf(" %.lf ", func(a) * func(c));
			if (func(a) * func(c) < 0)
			{
				b = c;
				printf(" %.15lf ", (a + b) / 2);
				
			}
			else
			{
				a = c;
				printf(" %.15lf ", (a + b) / 2);
				
			}
			//printf(" %.15lf ", (a + b) / 2);
			i++;
		}
	}
}

void BisectionMethodFunction(double(*func)(double x), double a, double b)
{
	double e = 0.1;
	for (int i = 0; i < 5; i++)
	{
		double c;
		int i = 0;
		printf(" %.lf ", func(a) * func(b));
		while (fabs(b - a) > 2 * e)
		{
			c = (a + b) / 2;
			//printf(" %.lf ", func(a) * func(c));
			if (func(a) * func(c) < 0)
			{
				b = c;
				//printf(" %.lf ", func(a) * func(c));
				//printf(" %.15lf ", (a + b) / 2);
			}
			else
			{
				a = c;
				//printf(" %.lf ", func(a) * func(c));
				//printf(" %.15lf ", (a + b) / 2);
			}
			printf(" %.lf ", func(a) * func(c));
			i++;
		}
	}
}


void AproximationsFixedPointIterations(double(*fi)(double x), double a, double b, double q)
{
	int i = 0;
	double e = 0.1;
	double x_k_1;
	double x_k = a;
	do
	{
		x_k_1 = x_k;
		x_k = fi(x_k_1);
		i++;
		printf(" %.15lf", x_k);
	} while (fabs(x_k - x_k_1) >= (1 - q) / q * e);
	//printf("Root : %.15lf; Iterations: %d\n", x_k, i);
}

//Вывод результатов для матлаба
void RunForMatlab()
{
	printf("Iterations:");
	printf("\nTHE BISECTION METHOD:");
	RunIterationsBisectionMethod(&Polynomial,1,2);
	RunIterationsBisectionMethod(&Transcendental,0.7,1.7);
	printf("\nFIXED-POINT ITERATIONS METHOD");
	RunIterationsFixedPointIterations(&PolynomialFi, 1, 2, 0.8);
	RunIterationsFixedPointIterations(&TranscendentalFi, 0.7, 1.7, 0.85);
	printf("\nTHE BISECTION METHOD:");
	RunRootsBisectionMethod(&Polynomial, 1, 2);
	RunRootsBisectionMethod(&Transcendental, 0.7, 1.7);
	printf("\nFIXED-POINT ITERATIONS METHOD");
	RunRootsFixedPointIterations(&PolynomialFi, 1, 2, 0.8);
	RunRootsFixedPointIterations(&TranscendentalFi, 0.7, 1.7, 0.85);
	printf("\nAproximations:\n");
	AproximationsBisectionMethod(&Polynomial, 1, 2);
	printf("\n");
	AproximationsBisectionMethod(&Transcendental, 0.7, 1.7);
	printf("\n");
	printf("\nTEST:");
	RunIterationsBisectionMethod(&Testfunc, -1, 2);
	RunRootsBisectionMethod(&Testfunc, -1, 2);
	RunIterationsFixedPointIterations(&TestfuncFi, -1, 2, 1.0 / 3);
	RunRootsFixedPointIterations(&TestfuncFi, -1, 2, 1.0 / 3);
	printf("\n");
	AproximationsBisectionMethod(&Testfunc, -1, 2);
	printf("\n");
	//BisectionMethodFunction(&Testfunc, -2, 1);
	
}


