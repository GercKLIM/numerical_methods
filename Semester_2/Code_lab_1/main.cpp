/*  Лабораторная работа №1 "Методы численного решения ОДУ"
 *  Разработчик программы: Климов О.Д. ФН2-51Б 2024г.
 *
 * */

#include <iostream>
#include <iostream>
#include <vector>
#include <iomanip>
#include <fstream>

#include "algebra.cpp"
#include "methods.cpp"

using namespace std;


int main() {

    // Равномерная сетка
    double h = 0.01;            // Шаг
    vector<double> Diapazon = {0, 10}; // Отрезок

    // Начальное условие
    double t = 0;
    vector<double> u0 = {0, 1};

    // Метод Эйлера явный
    Method_Euler_implicit(*ODU_0, u0, Diapazon, h);

    // Метод Эйлера неявный
    Method_Euler_explicit(*ODU_0, u0, Diapazon, h);

    // Метод симметричной схемы
    Method_Symmetric_scheme(*ODU_0, u0, Diapazon, h);

    // Метод Рунге-Кутты 2-x-стадийный
    Method_Runge_Kutta_2ord(*ODU_0, u0, Diapazon, h);

    // Метод Рунге-Кутты 4-x-стадийный
    Method_Runge_Kutta_4ord(*ODU_0, u0, Diapazon, h);

    // Метод Адамса_Башфорта 4-го порядка
    Method_Adamsa_bashforma(*ODU_0, u0, Diapazon, h);

    // Метод Прогноз-Коррекция
    Method_Predictor_corrector(*ODU_0, u0, Diapazon, h);

    cout << "Complete!" << endl;
    return 0;
}
