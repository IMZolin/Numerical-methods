#include <stdio.h>
#include <conio.h>
#include <math.h>
#include <stdlib.h>
#define PI 3.1415926535897932384
#define N 10


void File_Print_MKR(double A[], double B[], double C[], double D[], double E[], double ST[], double ST_C[])
{
    FILE* fout;
    fout = fopen("C:\\Users\\Le\\Desktop\\AppMath\\Numerical Methods\\CourseWork\\CW_MKR.txt", "w");


    for(int i = 0; i < 4*N; ++i)
    {
        fprintf(fout, "%le ", A[i]);
    }

    fprintf(fout, "\n");

    for(int i = 0; i < 2*N; ++i)
    {
        fprintf(fout, "%le ", B[i]);
    }

    fprintf(fout, "\n");

    for(int i = 0; i < N; ++i)
    {
        fprintf(fout, "%le ", C[i]);
    }

    fprintf(fout, "\n");

    for(int i = 0; i < 2*N; ++i)
    {
        fprintf(fout, "%le ", D[i]);
    }

    fprintf(fout, "\n");

    for(int i = 0; i < N; ++i)
    {
        fprintf(fout, "%le ", E[i]);
    }

    fprintf(fout, "\n");

    for(int i = 0; i < 14; ++i)
    {
        fprintf(fout, "%le ", ST[i]);
    }

    fprintf(fout, "\n");

    for(int i = 0; i < 14; ++i)
    {
        fprintf(fout, "%le ", ST_C[i]);
    }


    fprintf(fout, "\n");
    fclose(fout);

}


double ExactFun(double x)
{
    return 1 + 1/x;
}

double p_x(double x)
{
    return -1/(x*x*(x+1));
    //return 0.4*x*x;
}

double q_x(double x)
{
    return -2/(x*x*(x+1));
    //return 5*x;
}

double f_x(double x)
{
    return 1/(pow(x,4)*(x+1));
    //return 10;
}

/**
	 * n - число уравнений (строк матрицы)
	 * b - диагональ, лежащая над главной (нумеруется: [0;n-2])
	 * c - главная диагональ матрицы A (нумеруется: [0;n-1])
	 * a - диагональ, лежащая под главной (нумеруется: [1;n-1])
	 * f - правая часть (столбец)
	 * x - решение, массив x будет содержать ответ
*/
void solveMatrix (int n, double *a, double *c, double *b, double *f, double *x)
{
    double m;
    for (int i = 1; i < n; i++)
    {
        m = a[i]/c[i-1];
        c[i] = c[i] - m*b[i-1];
        f[i] = f[i] - m*f[i-1];
    }

    x[n-1] = f[n-1]/c[n-1];

    for (int i = n - 2; i >= 0; i--)
    {
        x[i]=(f[i]-b[i]*x[i+1])/c[i];
    }
}

/**
     * a left board
     * b right board
     * n dots amount
     * p
 */
void MKR(double a, double b, int n1, double* p, double* q, double* f, double alfa0, double beta0, double* solution)
{
    double h = (b-a)/(n1-1);
    for(int i = 1; i < n1-1; ++i)  //  except 0 and n
    {
        p[i] = p_x(a + h*i);
        //printf("pi: %lf\n", p[i]);
        q[i] = q_x(a + h*i);
        //printf("qi: %lf\n", q[i]);
        f[i] = f_x(a + h*i);
        //printf("fi: %lf\n", f[i]);
    }

    double* up = malloc((n1 - 1)*sizeof(double));
    double* diag = malloc((n1)*sizeof(double));
    double* under = malloc((n1)*sizeof(double));

    diag[0] = 1;
    f[0] = alfa0;
    up[0] = 0;

    diag[n1-1] = 1;
    f[n1-1] = beta0;
    under[n1-1] = 0;



    for(int i = 1; i < n1-1; ++i)
    {
        under[i] = (1 - h/2*p[i])/(h*h);
        //printf("under: %lf\n", under[i]);
        diag[i] = (-2 + h*h*q[i])/(h*h);
        //printf("diag: %lf\n", diag[i]);
        up[i] = (1 + h/2*p[i])/(h*h);
        //printf("up: %lf\n", up[i]);
    }

    solveMatrix(n1, under, diag, up, f, solution);

}

double solveMatrix_STEP(int n, double *a, double *c, double *b, double *f, double l, double r, double* absc)
{
    double m;
    double* x = malloc(n*sizeof(double));
    double tmp = 0;
    for (int i = 1; i < n; i++)
    {
        m = a[i]/c[i-1];
        c[i] = c[i] - m*b[i-1];
        f[i] = f[i] - m*f[i-1];
    }

    x[n-1] = f[n-1]/c[n-1];
    tmp = fabs(x[n-1] - ExactFun(r));

    for (int i = n - 2; i >= 0; i--)
    {
        x[i]=(f[i]-b[i]*x[i+1])/c[i];
        if(tmp < fabs(x[i] - ExactFun(l + (r-l)/(n-1)*i)))
        {
            tmp = fabs(x[i] - ExactFun(l + (r-l)/(n-1)*i));
            *absc = l + (r-l)/(n-1)*i;
        }
    }
    return tmp;
}

double MKR_STEP(double a, double b, int n1, double* p, double* q, double* f, double alfa0, double beta0, double* absc)
{
    double h = (b-a)/(n1-1);
    for(int i = 1; i < n1-1; ++i)  //  except 0 and n
    {
        p[i] = p_x(a + h*i);
        //printf("pi: %lf\n", p[i]);
        q[i] = q_x(a + h*i);
        //printf("qi: %lf\n", q[i]);
        f[i] = f_x(a + h*i);
        //printf("fi: %lf\n", f[i]);
    }

    double* up = malloc((n1 - 1)*sizeof(double));
    double* diag = malloc((n1)*sizeof(double));
    double* under = malloc((n1)*sizeof(double));

    diag[0] = 1;
    f[0] = alfa0;
    up[0] = 0;

    diag[n1-1] = 1;
    f[n1-1] = beta0;
    under[n1-1] = 0;



    for(int i = 1; i < n1-1; ++i)
    {
        under[i] = (1 - h/2*p[i])/(h*h);
        //printf("under: %lf\n", under[i]);
        diag[i] = (-2 + h*h*q[i])/(h*h);
        //printf("diag: %lf\n", diag[i]);
        up[i] = (1 + h/2*p[i])/(h*h);
        //printf("up: %lf\n", up[i]);
    }

    return solveMatrix_STEP(n1, under, diag, up, f, a, b, absc);
}


int main()
{
    double a = 0.2;
    double b = 1;
    int n = N; //  dots amount
    double* p = malloc((n-2)*sizeof(double));
    double* q = malloc((n-2)*sizeof(double));
    double* f = malloc((n)*sizeof(double));
    double alfa = 6;
    double beta = 2;
    double* solut = malloc(n*sizeof(double));
    double* absc = malloc(sizeof(double));

    printf("MAX MISTAKE AND ITS COORDINATE WITH %d DOTS:\n", n);
    printf("%lf ",MKR_STEP(a,b,n,p,q,f,alfa,beta, absc));
    printf("%lf", *absc);

    // putting into the loop

    printf("\n");
    int v = 15;
    double* STEP = malloc(20*sizeof(double));
    double* STEP_COORDINATE = malloc(20*sizeof(double));
    double* pSTEP = malloc((pow(2,v))*sizeof(double));
    double* qSTEP = malloc((pow(2,v))*sizeof(double));
    double* fSTEP = malloc((pow(2,v))*sizeof(double));
    for(int i = 1; i < v; ++i)
    {
        STEP[i-1] = MKR_STEP(a,b,pow(2,i),pSTEP,qSTEP,fSTEP,alfa,beta, &(STEP_COORDINATE[i-1]));
        printf("%le %lf\n", STEP[i-1], STEP_COORDINATE[i-1]);
    }






    double* A = malloc(4*n*sizeof(double));
    double* B = malloc(2*n*sizeof(double));
    double* C = malloc(n*sizeof(double));

    double* pC = malloc((n-2)*sizeof(double));
    double* qC = malloc((n-2)*sizeof(double));
    double* fC = malloc(n*sizeof(double));

    double* pB = malloc((2*n-2)*sizeof(double));
    double* qB = malloc((2*n-2)*sizeof(double));
    double* fB = malloc((2*n)*sizeof(double));

    MKR(a,b,2*n-1,pB,qB,fB,alfa,beta,B);
    MKR(a,b,n,pC,qC,fC,alfa,beta,C);


    for (int i = 0; i < 4*n; ++i)
    {
        A[i] = ExactFun(a + (b-a)/(4*n - 1)*i);
    }


    double* B_Mistake = malloc((2*n-1)*sizeof(double));
    double* C_Mistake = malloc(n*sizeof(double));

    for(int i = 0; i < 2*n-1; ++i)
    {
        B_Mistake[i] = fabs(B[i] - ExactFun(a + (b-a)/(2*n - 2)*i));
    }

    for(int i = 0; i < n; ++i)
    {
        C_Mistake[i] = fabs(C[i] - ExactFun(a + (b-a)/(n - 1)*i));
    }





    File_Print_MKR(A, B, C, B_Mistake, C_Mistake, STEP, STEP_COORDINATE);

    return 0;
}

