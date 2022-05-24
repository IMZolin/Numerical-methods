#include <stdio.h>
#include <conio.h>
#include <math.h>
#include <stdlib.h>
#define PI 3.1415926535897932384

void File_Print(double A[], double B[], double C[], double D[], double E[], double F[], double G[], int H[], double J[], double K[], double L[], double M[])
{
    FILE* fout;
    fout = fopen("C:\\Users\\Le\\Desktop\\AppMath\\Numerical Methods\\Lab_5\\Attempt_1\\L_5_text.txt", "w");

    int N = 1000;

    for(int i = 0; i < N; ++i)
    {
        fprintf(fout, "%le ", A[i]);
    }

    fprintf(fout, "\n");

    for(int i = 0; i < N; ++i)
    {
        fprintf(fout, "%le ", B[i]);
    }

    fprintf(fout, "\n");

    for(int i = 0; i < N; ++i)
    {
        fprintf(fout, "%le ", C[i]);
    }

    fprintf(fout, "\n");

    for(int i = 0; i < N; ++i)
    {
        fprintf(fout, "%le ", D[i]);
    }

    fprintf(fout, "\n");

    for(int i = 0; i < N; ++i)
    {
        fprintf(fout, "%le ", E[i]);
    }

    fprintf(fout, "\n");

    for(int i = 0; i < N; ++i)
    {
        fprintf(fout, "%lf ", F[i]);
    }

    fprintf(fout, "\n");

    for(int i = 0; i < N; ++i)
    {
        fprintf(fout, "%le ", G[i]);
    }

    fprintf(fout, "\n");

    for(int i = 0; i < N; ++i)
    {
        fprintf(fout, "%d ", H[i]);
    }

    fprintf(fout, "\n");

    for(int i = 0; i < N; ++i)
    {
        fprintf(fout, "%le ", J[i]);
    }

    fprintf(fout, "\n");

    for(int i = 0; i < N; ++i)
    {
        fprintf(fout, "%le ", K[i]);
    }

    fprintf(fout, "\n");

    for(int i = 0; i < N; ++i)
    {
        fprintf(fout, "%le ", L[i]);
    }

    fprintf(fout, "\n");

    for(int i = 0; i < N; ++i)
    {
        fprintf(fout, "%le ", M[i]);
    }




    fprintf(fout, "\n");
    fclose(fout);

}

double Function(double x, double y)
{
    if (x == 1)
        return -1;
    else
        return (2*y + log(x))/(x*log(x));
}

double ExactFun(double x)
{
    return log(x)*log(x) - log(x);
}

double Runge_Kutta(double a, double b, int n, double y, double* p2)  //  n - amount of sections
{
    double h = (b - a)/n;
    double x[n+1];  //  x - absc for functions
    //double y = -0.08622614;
    double sigma[4] = {0};

    for(int i = 0; i < n + 1; ++i)
    {
        x[i] = a + h*i;
        *(p2 + i) = y;

        sigma[0] = Function(x[i], y);
        sigma[1] = Function(x[i] + h/3, y + h*sigma[0]/3);
        sigma[2] = Function(x[i] + 2*h/3, y - h*sigma[0]/3 + h*sigma[1]);
        sigma[3] = Function(x[i] + h, y + h*sigma[0] - h*sigma[1] + h*sigma[2]);

        y = y + h*(sigma[0] + 3*sigma[1] + 3*sigma[2] + sigma[3])/8;
    }
    return *(p2 + n);
}

int Accuracy(double a, double b, double e, double* px, double* py, int n)
{
    double tmp1[2];
    double tmp2[3];
    int start = 1;
    int start_tmp = 0;
    double l = a;
    double r = b;
    int i = 0;
    double y0 = -0.08622614;

    for(int k = 0; k < n; ++k)
    {
        r = a + (b-a)/n*(k+1);
        while (l != r)
        {
            double x = Runge_Kutta(l, r, 2, y0, tmp2);
            double y = Runge_Kutta(l, r, 1, y0, tmp1);
            while (fabs(x-y)/15 > e)
            {
                start++;
                r = l + (r - l) / 2;
                x = Runge_Kutta(l, r, 2, y0, tmp2);
                y = Runge_Kutta(l, r, 1, y0, tmp1);

            }
            *(px + i + 1) = r;
            *(py + i + 1) = y;
            l = r;
            r = a + (b-a)/n*(k+1);
            i++;
            y0 = y;
        }
        if(start > start_tmp)
            start_tmp = start;
        start = 1;
    }
    return start_tmp;
}

double Runge_Kutta_2(double a, double b, int n, double y)  //  n - amount of sections
{
    double h = (b - a)/n;
    double x[n+1];  //  x - absc for functions
    //double y = -0.08622614;
    double sigma[4] = {0};
    double t = y;

    for(int i = 0; i < n + 1; ++i)
    {
        x[i] = a + h*i;
        t = y;
        sigma[0] = Function(x[i], y);
        sigma[1] = Function(x[i] + h/3, y + h*sigma[0]/3);
        sigma[2] = Function(x[i] + 2*h/3, y - h*sigma[0]/3 + h*sigma[1]);
        sigma[3] = Function(x[i] + h, y + h*sigma[0] - h*sigma[1] + h*sigma[2]);

        y = y + h*(sigma[0] + 3*sigma[1] + 3*sigma[2] + sigma[3])/8;
    }
    //return fabs((sigma[1] - sigma[2])/(sigma[0] - sigma[1]));
    return t;
}

int Accuracy_2(double a, double b, double e, int n, double* p1, double* p2, double y)
{
    double l = a;
    double r = b;
    double y0 = y;
    int veshka = 1;
    int veshka_tmp = 1;
    int i = 0;
    for(int k = 0; k < n; ++k)
    {
        r = a + (b-a)/n*(k+1);
        double x = Runge_Kutta_2(l, r, pow(2, veshka), y0);
        double y = Runge_Kutta_2(l, r, 2*pow(2, veshka), y0);
        double z = (x - y)/15;
        while(fabs((x - y)/15) > e)
        {
            veshka++;
            x = y;
            y = Runge_Kutta_2(l, r, 2*pow(2, veshka), y0);
            z = (x - y)/15;
        }
        l = r;
        y0 = x;
        *(p1 + i + 1) = r;
        *(p2 + i + 1) = y;
        i++;
        if(veshka > veshka_tmp)
            veshka_tmp = veshka;
        veshka = 1;
    }
    return veshka_tmp;
}

int main()
{

    int n = 200;
    double a = 1.001; //  change to 1
    int b = 3;
    double y0 = ExactFun(a);

    double A[4*n+1];  //  Exact fun (4 - s zapasom, tak kak budem brat n and 2n)
    double B[n+1];  //  n dots
    double C[2*n+1];  //  n dots

    double B_mistake[n + 1];
    double C_mistake[2*n + 1];

    for(int i = 0; i < 4*n + 1; ++i)
    {
        double h = (b - a)/(4*n);
        A[i] = ExactFun(a + h*i);
    }

    Runge_Kutta(a, b, n, y0, B);
    Runge_Kutta(a, b, 2*n, y0,C);

    double Perturb[2*n + 1];
    Runge_Kutta(a, b, 2*n, y0, Perturb);

    double Perturb2[2*n + 1];
    Runge_Kutta(a - 0.001, b, 2*n, y0, Perturb2);


    //---------------//

    double P1[1000] = {0};
    double P2[1000] = {0};
    double Iter[15] = {0};
    double Accur[15] = {0};
    P1[0] = a;
    P2[0] = y0;
    int iter;
    for(int l = 0; l < 14; ++l)  //  change 5 to 14
    {
        Iter[l] = Accuracy_2(a, b, pow(10, -l-1), n, P1, P2, y0);
        iter = 0;
        while(P1[iter] != b)
        {
            if(Accur[l] < fabs(P2[iter] - ExactFun(P1[iter])))
                Accur[l] = fabs(P2[iter] - ExactFun(P1[iter]));
            iter++;
        }
        if(Accur[l] < fabs(P2[iter] - ExactFun(P1[iter])))
            Accur[l] = fabs(P2[iter] - ExactFun(P1[iter]));
        printf("%le ", Accur[l]);
    }

    //---------------//

    double PP1[1000] = {0};
    double PP2[1000] = {0};
    PP1[0] = 2;
    PP2[0] = y0;
    int w = 0;
    double Perturm_FM[14] = {0};

    for(int i = 0; i < 14; ++i)
    {
        w = 0;
        Accuracy_2(a + pow(10, i - 15),b, pow(10,-10), n, PP1, PP2, y0); //  change -2 to -6
        while(PP1[w] != b)
        {
            w++;
            if(Perturm_FM[i] < fabs(PP2[w] - ExactFun(PP1[w])))
                Perturm_FM[i] = fabs(PP2[w] - ExactFun(PP1[w]));
        }
        if(Perturm_FM[i] < fabs(PP2[w] - ExactFun(PP1[w])))
            Perturm_FM[i] = fabs(PP2[w] - ExactFun(PP1[w]));
        printf("\n OK: %le\n", Perturm_FM[i]);
    }


    //--------------//




    for(int i = 0; i < n+1; ++i)
    {
        B_mistake[i] = fabs(A[4*i] - B[i]);
    }

    for(int i = 0; i < 2*n+1; ++i)
    {
        C_mistake[i] = fabs(A[2*i] - C[i]);
    }


    double px[1000] = {0};
    double py[1000] = {0};
    double FM[13];
    px[0] = a;
    py[0] = -0.08622614;
    /*
    //printf("%d\n",Accuracy(a, b, pow(10, -7), px, py, 1));
    for(int k = 0; k < 13; ++k)
    {
        Accuracy_2(a, b, pow(10, -k-1), 1, px, py);
        int i = 0;
        while (px[i] != 3) {
            //printf("%lf ", px[i]);
            //printf("%le\n", py[i] - ExactFun(px[i]));
            i++;
        }
        FM[k] = fabs(py[i] - ExactFun(px[i]));
        printf("%le\n", FM[k]);
        //printf("%lf ", px[i]);
        //printf("%le\n", py[i] - ExactFun(px[i]));
    }
     */

    double px_1[1000] = {0};  //  с запасом
    double py_1[1000] = {0};
    int Steps[14];
    double Fact_Mist[14];
    /*
    for(int i = 0; i < 14; ++i)
    {
       Steps[i] = Accuracy(a, b, pow(10, -i-1), px_1, py_1, 1);
       printf("%d ", Steps[i]);
    }
     */


    for(int i = 0; i < 14; ++i)
        printf("%lf, ", Iter[i]);


    File_Print(A,B,C, B_mistake, C_mistake, px, py, Iter, Accur, Perturb, Perturb2, Perturm_FM);

    return 0;
}

