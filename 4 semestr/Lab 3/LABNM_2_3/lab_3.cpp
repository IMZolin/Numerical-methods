#include <iostream>
#include <math.h>
#include <cmath>
#include <vector>
#include <iomanip>
#include <fstream>
using namespace std;

double f(double x) {
    return pow(x, 5) - 6.2 * pow(x, 3) + 3.5 * pow(x, 2) - 7 * x - 2.1;
}

double f4(double x) {
    return pow(x, 6) / 6 - 31 * pow(x, 4) / 20 + 7 * pow(x, 3) / 6 - 3.5 * pow(x, 2) - 2.1 * x;
}

double simpson_integral(double a, double b, int m) {
    const double h = (b - a) / (m * 2);
    double s1 = 0, s2 = 0;

    for (int i = 2; i <= (2 * m - 2); i += 2) { // четные
        s2 += f(a + i * h);
    }
    for (int i = 1; i <= (2 * m - 1); i += 2) { // нечетные
        s1 += f(a + i * h);
    }
    return h / 3 * (f(a) + f(b) + 4 * s1 + 2 * s2);
}

void PrintVectorToFile(string filename, vector<double> vec)
{
    ofstream file;
    file.open(filename);
    if (!file.is_open())
    {
        cout << "Uh oh.File wasn't opened. Check if the name of file " << filename << " correct.\n" << endl;
        exit(1);
    }
    for (auto& v : vec)
    {
        file << setprecision(17) << v;
        file << ' ';
    }
    file.close();
}

int main() {
    double a = -2.9; // Левая граница интегрирования
    double b = 0.4; // Правая граница интегрирования
    double I, I_prev;
    double integr13 = 0.1393832154470942;
    double integr11 = 1.493648265624854;
    double t = 60.925837499999993;
    t = f4(b) - f4(a);
    vector<double> epsilon;
    vector<double> iter;
    vector<double> error;
    vector<double> h;
    string filename;
    for (double eps = 0.1; eps >= 10e-8; eps *= 0.1) {
        int m = 1; //начальное внутренние число точек
        I_prev = 0;
        I = 0;
        I = simpson_integral(a, b, m); //первое приближение для интеграла
        do {
            I_prev = I;     //второе приближение
            m = 2 * m;  //увеличение числа шагов в два раза,
                        //т.е. уменьшение значения шага в два раза
            I = simpson_integral(a, b, m);
        } while (fabs(I - I_prev) / 15 > eps);  //сравнение приближений с заданной точностью
        cout << eps << "  " << I << "  " <<log2(2 * m) << endl;
        iter.push_back((int)log2(2*m));
        h.push_back((b - a) / (2*m));
        // Подсчет практической ошибки
        double error_pract = fabs(t - I);
        cout << " E = " << error_pract << " m = " << log2(2 * m) << endl;
        filename = "f_epsilon.txt";
        epsilon.push_back(eps);
        PrintVectorToFile(filename, epsilon);
        
        filename = "f_iterations.txt";
        PrintVectorToFile(filename, iter);
        error.push_back(fabs(I - I_prev) / 15);
        filename = "f_error.txt";
        PrintVectorToFile(filename, error);
        filename = "f_h.txt";
        PrintVectorToFile(filename, h);
    }


    return 0;
}