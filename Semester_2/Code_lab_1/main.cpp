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

    Method_Euler(t, u0, Diapazon, h);


    //cout << ODU(0, u0) << endl;
    cout << "Complete!" << endl;
    return 0;
}
