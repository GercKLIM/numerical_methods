#pragma once

#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <cmath>

using namespace std;



/* Метод Эйлера явный */
void Method_Euler_explicit(vector<double> (*ODU)(const vector<double>&), const vector<double> u0, const vector<double>& diapazon, double h);

/* Метод Эйлера неявный */
void Method_Euler_implicit(vector<double> (*ODU)(const vector<double>&), const vector<double> u0, const vector<double>& diapazon, double h);

/* Метод Симетричной схемы 2-х шаговый */
void Method_Symmetric_scheme(vector<double> (*ODU)(const vector<double>&), const vector<double> u0, const vector<double>& diapazon, double h);

/* Метод Рунге-Кутты 2-го порядка точности */
void Method_Runge_Kutta_2ord(vector<double> (*ODU)(const vector<double>&), const vector<double> u0, const vector<double>& diapazon, double h);

/* Метод Рунге-Кутты 4-го порядка точности */
void Method_Runge_Kutta_3ord(vector<double> (*ODU)(const vector<double>&), const vector<double> u0, const vector<double>& diapazon, double h);

/* Метод Адамса-Башформа 4-го порядка точности */
void Method_Adamsa_bashforma(vector<double> (*ODU)(const vector<double>&), const vector<double> u0, const vector<double>& diapazon, double h);

/* Метод "Прогноз-Коррекция" 4-го порядка точности */
void Method_Predictor_corrector(vector<double> (*ODU)(const vector<double>&), const vector<double> u0, const vector<double>& diapazon, double h);