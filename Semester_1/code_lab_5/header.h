#pragma once

#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <cmath>

using namespace std;

/* ЗАДАЕМ УСЛОВИЯ  */

/* Исходные функции */
double f(double x, int key) {
    switch (key) {
        case 1:
            // x in [-1, 1]
            // Корни: 0.1, 0.22, 0,55, 0,7, 0.75
            return (x - 0.1) * (x - 0.22) * (x - 0.55) * (x - 0.7) * (x - 0.75);
        case 2:
            // x in [-1, 10]
            // Корень: 0
            // Приближение x0 = 8 - метод ньютона ломается
            return sqrt(x + 1) - 1;
        case 3:
            // x in [0, 1]
            // Корень: 0.2
            // С аналитической производной метод ньютона зацикливается при х0 = 0.
            return 35 * x * x * x - 67 * x * x - 3 * x + 3;
    };
}

/* Аналитические производные исходных функций */
double df(double x, int key) {
    switch (key) {
        case 1:
            return 0.121495 + 5 * x * (-0.30238 + x * (1.1907 + (-1.856 + x) * x));
        case 2:
            return 1 / (2 * sqrt(1 + x));
        case 3:
            return -3 + x * (-134 + 105 * x);
    };
}

/* Исходные системы функций */
double g(double x, double y, int key) {
    switch (key)
    {
        // x in [-10, 10], y in [-10, 10]
        // Корни: [4, -1], [-4, 1]
        case 1:
            return x * x - y * y - 15;
        case 2:
            return x * y + 4;

        // x in [-10, 10], y in [-10, 10]
        // Корни: [-3, 1], [1, -3], [1, 2], [2, 1]
        case 3:
            return x * x + y * y + x + y - 8;
        case 4:
            return x * x + y * y + x * y - 7;
    };
}

vector<double> dg(double x, double y, int key) {
    switch (key)
    {
        case 1:
            return{ 2 * x, -2 * y }; //df1/dx df1/dy
        case 2:
            return { y, x }; //df2/dx df2/dy
        case 3:
            return{ 2 * x + 1, 2 * y + 1 };
        case 4:
            return { 2 * x + y, 2 * y + x };
    };
}

/* Аналитические производные исходных систем функций */
//df1/dx df1/dy
//df2/dx df2/dy

/* Функция вывода матрицы на экран */
template <typename T>
void print(const vector<vector<T>>& matrix);

/* Функция вывода вектора на экран */
template <typename T>
void print(const vector<T>& vec);

/* Функция численного вычисления производной функции в точке */
double df(double x, int key, double eps);

/* Функция локализации корней */
vector<vector<double>> locale_roots(int key, double l, double r, int div);

/* Функция решения уравнения методом Бисекции */
double Method_Bisection(int key, vector<double> line, double eps);

/* Функция для начального приближения решения методом Хорд */
double Method_Hord(int key, double xprev, double eps);

/* Функция решения уравнения методом Ньютона */

double Method_Newton(int key, vector<double> line, double eps, double x0, bool flag, bool flag2);

/* Функция решения уравнения методом Ньютона c модификацией против зацикливания */
double Method_Newton2(int key, vector<double> line, double eps, double alpha, double x0, bool flag, bool flag2);

/* Функция решения уравнения модифицированным методом Ньютона*/
double Method_Newton_mod(int key, vector<double> line, double x0, double eps, double alpha, bool flag, bool flag2);

/* Функция решения системы уравнений методом Ньютона */

vector<double> Method_Newton_sys(int key, vector<double> x0, double eps, bool flag);

vector<double> Method_Newton_sys2(int key1, int key2,double L1, double L2, double eps, bool flag, bool flag2);