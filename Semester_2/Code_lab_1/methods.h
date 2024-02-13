#pragma once

#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <cmath>

using namespace std;

/* Переопределение вывода для vector */
template<typename T>
ostream& operator<<(ostream& os, const vector<T>& vec);

/* Метод Рунге-Кутты 4-го порядка точности */

/* Метод Адамса-Башформа 4-го порядка точности */

/* Метод "Прогноз-Коррекция" 4-го порядка точности */