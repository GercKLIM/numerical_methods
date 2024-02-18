#pragma once

#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <cmath>

using namespace std;



/* Метод Эйлера явный */
void Method_Euler(vector<double> (*ODU)(vector<double>), const vector<double> u0, const vector<double>& diapazon, double h);

/* Метод Эйлера неявный */
void Method_Euler_implicit(vector<double> (*ODU)(const vector<double>&), const vector<double> u0, const vector<double>& diapazon, double h);

/* Метод Рунге-Кутты 4-го порядка точности */

/* Метод Адамса-Башформа 4-го порядка точности */

/* Метод "Прогноз-Коррекция" 4-го порядка точности */