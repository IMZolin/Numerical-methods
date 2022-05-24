#include <stdio.h>
#include <conio.h>
#include <math.h>
#include <stdlib.h>
#define PI 3.1415926535897932384
#define N 9

void File_Print_RK(double A[], double B[], double C[], double D[], double E[], double ST[], double ST_C[])
{
    FILE* fout;
    fout = fopen("C:\\Users\\Le\\Desktop\\AppMath\\Numerical Methods\\CourseWork\\CW_RK.txt", "w");


    for(int i = 0; i < 4*N + 1; ++i)
    {
        fprintf(fout, "%le ", A[i]);
    }

    fprintf(fout, "\n");

    for(int i = 0; i < N + 1; ++i)
    {
        fprintf(fout, "%le ", B[i]);
    }

    fprintf(fout, "\n");

    for(int i = 0; i < 2*N + 1; ++i)
    {
        fprintf(fout, "%le ", C[i]);
    }

    fprintf(fout, "\n");

    for(int i = 0; i < N + 1; ++i)
    {
        fprintf(fout, "%le ", D[i]);
    }

    fprintf(fout, "\n");

    for(int i = 0; i < 2*N + 1; ++i)
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

double Function(double x, double y, double z)
{
    return -p_x(x)*z - q_x(x)*y + f_x(x);
    //return 10;
}

double Runge_Kutta(double a, double b, int n, double y0, double z0, double* p)  //  n - amount of sections
{
    double h = (b - a)/n;
    double x[n+1];
    double sigma[4] = {0};
    double gamma[4] = {0};
    double t = y0;

    for(int i = 0; i < n + 1; ++i)
    {
        x[i] = a + h*i;
        *(p + i) = y0;
        t = y0;

        //printf("%lf %lf\n", y0, z0);

        sigma[0] = Function(x[i], y0, z0);
        gamma[0] = z0;

        sigma[1] = Function(x[i] + h/2, y0 + h*gamma[0]/2, z0 + h*sigma[0]/2);
        gamma[1] = z0 + h*sigma[0]/2;

        sigma[2] = Function(x[i] + h/2, y0 + h*gamma[1]/2, z0 + h*sigma[1]/2);
        gamma[2] = z0 + h*sigma[1]/2;

        sigma[3] = Function(x[i] + h, y0 + h*gamma[2], z0 + h*sigma[2]);
        gamma[3] = z0 + h*sigma[2];

        z0 = z0 + h*(sigma[0] + 2*sigma[1] + 2*sigma[2] + sigma[3])/6;
        y0 = y0 + h*(gamma[0] + 2*gamma[1] + 2*gamma[2] + gamma[3])/6;


    }

    return t;
}

double Runge_Kutta_2(double a, double b, int n, double y0, double z0)  //  n - amount of sections
{
    double h = (b - a)/n;
    double x[n+1];
    double sigma[4] = {0};
    double gamma[4] = {0};
    double t = y0;

    for(int i = 0; i < n + 1; ++i)
    {
        x[i] = a + h*i;
        t = y0;

        //printf("%lf %lf\n", y0, z0);

        sigma[0] = Function(x[i], y0, z0);
        gamma[0] = z0;

        sigma[1] = Function(x[i] + h/2, y0 + h*gamma[0]/2, z0 + h*sigma[0]/2);
        gamma[1] = z0 + h*sigma[0]/2;

        sigma[2] = Function(x[i] + h/2, y0 + h*gamma[1]/2, z0 + h*sigma[1]/2);
        gamma[2] = z0 + h*sigma[1]/2;

        sigma[3] = Function(x[i] + h, y0 + h*gamma[2], z0 + h*sigma[2]);
        gamma[3] = z0 + h*sigma[2];

        z0 = z0 + h*(sigma[0] + 2*sigma[1] + 2*sigma[2] + sigma[3])/6;
        y0 = y0 + h*(gamma[0] + 2*gamma[1] + 2*gamma[2] + gamma[3])/6;
    }

    return t;
}

double Runge_Kutta_3(double a, double b, int n, double y0, double z0)  //  n - amount of sections
{
    double h = (b - a)/n;
    double x[n+1];
    double sigma[4] = {0};
    double gamma[4] = {0};
    double t = y0;

    for(int i = 0; i < n + 1; ++i)
    {
        x[i] = a + h*i;
        t = z0;

        //printf("%lf %lf\n", y0, z0);

        sigma[0] = Function(x[i], y0, z0);
        gamma[0] = z0;

        sigma[1] = Function(x[i] + h/2, y0 + h*gamma[0]/2, z0 + h*sigma[0]/2);
        gamma[1] = z0 + h*sigma[0]/2;

        sigma[2] = Function(x[i] + h/2, y0 + h*gamma[1]/2, z0 + h*sigma[1]/2);
        gamma[2] = z0 + h*sigma[1]/2;

        sigma[3] = Function(x[i] + h, y0 + h*gamma[2], z0 + h*sigma[2]);
        gamma[3] = z0 + h*sigma[2];

        z0 = z0 + h*(sigma[0] + 2*sigma[1] + 2*sigma[2] + sigma[3])/6;
        y0 = y0 + h*(gamma[0] + 2*gamma[1] + 2*gamma[2] + gamma[3])/6;
    }

    return t;
}

int Accuracy(double a, double b, double e, int n, double* p1, double* p2, double y0, double z0)
{
    double l = a;
    double r = b;
    int veshka = 1;
    int veshka_tmp = 1;
    int i = 0;
    for(int k = 0; k < n; ++k)
    {
        r = a + (b-a)/n*(k+1);
        double x = Runge_Kutta_2(l, r, pow(2, veshka), y0, z0);
        double y = Runge_Kutta_2(l, r, 2*pow(2, veshka), y0, z0);
        while(fabs((x - y)/15) > e)
        {
            veshka++;
            x = y;
            y = Runge_Kutta_2(l, r, 2*pow(2, veshka), y0, z0);
        }
        z0 = Runge_Kutta_3(l, r, pow(2, veshka), y0, z0);
        l = r;
        y0 = x;
        *(p1 + i + 1) = r;
        *(p2 + i + 1) = y;
        //printf("%lf %lf\n", r, y);
        i++;
        if(veshka > veshka_tmp)
            veshka_tmp = veshka;
        veshka = 1;
    }
    return veshka_tmp;
}

int main()
{
    int n = N;
    double a = 0.2;
    double b = 1;
    double y0 = 6;
    double z0 = -25;
    double* A = malloc((4*n + 1)*sizeof(double));
    if(NULL == A)
        printf("MEMORY IS NOT GIVEN");
    double* B = malloc((n + 1)*sizeof(double));
    if(NULL == B)
        printf("MEMORY IS NOT GIVEN");
    double* C = malloc((2*n + 1)*sizeof(double));
    if(NULL == C)
        printf("MEMORY IS NOT GIVEN");

    double* B_mistake = malloc((n + 1)*sizeof(double));
    if(NULL == B_mistake)
        printf("MEMORY IS NOT GIVEN");
    double* C_mistake = malloc((2*n + 1)*sizeof(double));
    if(NULL == C_mistake)
        printf("MEMORY IS NOT GIVEN");


    Runge_Kutta(a, b, n, y0, z0, B);
    Runge_Kutta(a, b, 2*n, y0, z0, C);

    for(int i = 0; i < 4*n + 1; ++i)
    {
        double h = (b - a)/(4*n);
        A[i] = ExactFun(a + h*i);
    }

    for(int i = 0; i < n+1; ++i)
    {
        B_mistake[i] = fabs(A[4*i] - B[i]);
    }

    for(int i = 0; i < 2*n+1; ++i)
    {
        C_mistake[i] = fabs(A[2*i] - C[i]);
    }


    double P1[1000];
    double P2[1000];
    //Accuracy(a, b, 0.00001, 10, P1, P2, y0, z0);

    //printf("%lf", Runge_Kutta_3(a, b, 20, y0, z0));

    // Second part:
    double* STEP = malloc(20*sizeof(double));
    if(NULL == STEP)
        printf("MEMORY IS NOT GIVEN");
    double* STEP_COORDINATE = malloc(20*sizeof(double));
    if(NULL == STEP_COORDINATE)
        printf("MEMORY IS NOT GIVEN");
    int v = 15;
    double tmp = 0;
    double* X = malloc((pow(2,v) + 1)*sizeof(double));
    if(NULL == X)
        printf("MEMORY IS NOT GIVEN");
    for(int i = 1; i < v; ++i)
    {
        tmp = 0;
        Runge_Kutta(a, b, pow(2,i), y0, z0, X);
        for(int k = 0; k < pow(2,i) + 1; ++k)
        {
            if(tmp < fabs(X[k] - ExactFun(a + (b-a)/(pow(2,i))*k)))
            {
                tmp = fabs(X[k] - ExactFun(a + (b-a)/(pow(2,i))*k));
                STEP_COORDINATE[i-1] = a + (b-a)/(pow(2,i))*k;
            }
        }
        STEP[i-1] = tmp;
        printf("%le %lf\n", STEP[i-1], STEP_COORDINATE[i-1]);
    }



    File_Print_RK(A, B, C, B_mistake, C_mistake, STEP, STEP_COORDINATE);


    return 0;
}

