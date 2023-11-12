#pragma once

/* methods.h и methods.cpp - файлы,
 * в которых соответственно объявляются и реализовываются
 * функции методов необходимых в лабораторной работе.
 * */


#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <cmath>


using namespace std;


/* ### Функции лабы 1 ### */

/* Функция для LU-разложения с частичным выбором */
template <typename T>
void lu_decomposition(vector<vector<T>>& A, vector<vector<T>>& L, vector<vector<T>>& U);


/* Функция для вычисления нормы вектора невязки */
template <typename T>
T norm_vector_nevazki(vector<vector<T>> A, vector<T> b, vector<T> x, const int n);


/* Функция для решения СЛАУ прямым методом Гаусса */
template <typename T>
vector<T> method_Gaussa(vector<vector<T>>& matrix, vector<T>& vec);


/* Функция для QR-разложения матрицы */
template <typename T>
void QR_decomposition(const vector<vector<T>>& A, vector<vector<T>>& Q, vector<vector<T>>& R);


/* Функция для решения СЛАУ методом QR-разложения */
template <typename T>
vector<T> method_QR(vector<vector<T>>& A, vector<T>& b);


/* Функция для оценки изменения числа обуcловленности от возмущения вектора правой части */
template <typename T>
void min_change_cond(vector<vector<T>> matrix, vector<T> vec, vector<T> mod);


/* ### Функций лабы 2 ### */

/* Функция решения СЛАУ методом Простой Итерации */
template <typename T>
vector<T> My_Solve_Simple_Iter(vector<vector<T>> A, vector<T> b, T eps, T tau, T MaxIterations);

/* Функция решения СЛАУ методом Якоби */

/* Функция решения СЛАУ методом Зейделя */

/* Функция решения СЛАУ методом Релаксации */

/* Функция решения трехдиагональной СЛАУ большой размерности методом Зейделя */

/* Функция решения трехдиагональной СЛАУ большой размерности методом Релаксации */

/* Функция представления матрицы С в виде: C = C_l + C_d + D_u */

/* Функция исследования итерационного параметра */

/* Функция исследования сходимости при различных начальных приближениях */

