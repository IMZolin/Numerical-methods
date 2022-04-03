#include <iostream>
#include <math.h>
#include <cmath>
#include <vector>
#include <iomanip>
#include <fstream>
using namespace std;

double f(double x) {
    //return (x);
    return pow(x, 5) - 6.2 * pow(x, 3) + 3.5 * pow(x, 2) - 7 * x - 2.1;
}

double f4(double x) {
    //return (x);
    return (16 * pow(x, 4) * exp(-pow(x, 2)) - 48 * pow(x, 2) * exp(-pow(x, 2)) + 12 * exp(-pow(x, 2)));
}

double simpson_integral(double a, double b, int m) {
    const double h = (b - a) / (m * 2);
    double k1 = 0, k2 = 0;

    for (int i = 2; i <= (2 * m - 2); i += 2) { // четные
        k2 += f(a + i * h);
    }
    for (int i = 1; i <= (2 * m - 1); i += 2) { // нечетные
        k1 += f(a + i * h);
    }


    return h / 3 * (f(a) + f(b) + 4 * k1 + 2 * k2);
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
    double s1, s;
    double integr13 = 0.1393832154470942;
    double integr11 = 1.493648265624854;
    double t = 60.925837499999993;
    vector<double> epsilon;
    vector<double> iter;
    vector<double> error;
    string filename;
    for (double eps = 0.1; eps >= 10e-8; eps *= 0.1) {

        int m = 1; //начальное внутренние число точек
        s = 0;
        s1 = 0;
        s1 = simpson_integral(a, b, m); //первое приближение для интеграла
        int k = 0;
        do {
            k++;
            s = s1;     //второе приближение
            m = 2 * m;  //увеличение числа шагов в два раза,
                        //т.е. уменьшение значения шага в два раза
            s1 = simpson_integral(a, b, m);
        } while (fabs(s1 - s) / 15 > eps);  //сравнение приближений с заданной точностью
        cout << eps << "  " << s1 << "  " << 2 * m << endl;

        // Подсчет практической ошибки
        double error_pract = fabs(t - s1);
        cout << " E = " << error_pract << " m = " << 2 * m << endl;
        filename = "f_epsilon.txt";
        epsilon.push_back(eps);
        PrintVectorToFile(filename, epsilon);
        iter.push_back((int)log(m));
        filename = "f_iterations.txt";
        PrintVectorToFile(filename, iter);
        error.push_back(error_pract);
        filename = "f_error.txt";
        PrintVectorToFile(filename, error);
    }


    return 0;
}