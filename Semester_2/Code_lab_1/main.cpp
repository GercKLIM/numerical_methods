/*  Лабораторная работа №1 "Методы численного решения ОДУ"
 *  Разработчик программы: Климов О.Д. ФН2-51Б 2024г.
 *
 * */

#include <iostream>
#include <vector>
#include <chrono>

#include "algebra.cpp"
#include "methods.cpp"

using namespace std;


int main() {

    // Равномерная сетка
    double h = 1e-1;            // Шаг
    vector<double> Diapazon = {0, 1}; // Отрезок

    // Начальное условие
    double t = 0;
    vector<double> u0 = {0, 1};



    /* ### Методы ### */

    // Метод Эйлера явный
    //auto start_time = chrono::high_resolution_clock::now();
    Method_Euler_explicit(*ODU_0, u0, Diapazon, h);
    //auto end_time = chrono::high_resolution_clock::now();
    //chrono::duration<double> execution_time = end_time - start_time;

    //cout << "Time: " << execution_time.count() << " seconds" << endl;

    // Метод Эйлера неявный
    Method_Euler_implicit(*ODU_0, u0, Diapazon, h);

    // Метод симметричной схемы
    Method_Symmetric_scheme(*ODU_0, u0, Diapazon, h);

    // Метод Рунге-Кутты 2-x-стадийный
    Method_Runge_Kutta_2ord(*ODU_0, u0, Diapazon, h);

    // Метод Рунге-Кутты 2-x-стадийный c автоматическим шагом
    Method_Runge_Kutta_2ord_auto(*ODU_0, u0, Diapazon, h);

    // Метод Рунге-Кутты 4-x-стадийный
    Method_Runge_Kutta_4ord(*ODU_0, u0, Diapazon, h);

    // Метод Рунге-Кутты 4-x-стадийный c автоматическим шагом
    Method_Runge_Kutta_4ord_auto(*ODU_0, u0, Diapazon, h);

    // Метод Адамса_Башфорта 4-го порядка
    Method_Adamsa_bashforma(*ODU_0, u0, Diapazon, h);

    // Метод Прогноз-Коррекция
    Method_Predictor_corrector(*ODU_0, u0, Diapazon, h);

    cout << "Complete!" << endl;
    return 0;
}
