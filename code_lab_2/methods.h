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


/* Функция преобразования матрицы в сумму из Нижнетреугольной, Диагональной, Верхнетреугольной */
template<typename T>
void LDU_decomposotion(vector<vector<T>> A, vector<vector<T>>& L, vector<vector<T>>& D, vector<vector<T>>& U);

/* Функция решения СЛАУ методом Простой Итерации */
template<typename T>
vector<T> method_SimpleIteration(vector<vector<T>> A, vector<T> b, vector<T> x0, T tau, T eps, int MaxIter);


/* Функция решения СЛАУ методом Якоби */
template<typename T>
vector<T> method_Yacobi(vector<vector<T>> A, vector<T> b, vector<T> x0, T eps, int MaxIter);

/* Функция решения СЛАУ методом Зейделя */
template<typename T>
vector<T> method_Zeidel(vector<vector<T>> A, vector<T> b, vector<T> x0, T eps, int MaxIter);

/* Функция решения СЛАУ методом Релаксации */
template<typename T>
vector<T> method_Relax(vector<vector<T>> A, vector<T> b, vector<T> x0, T w, T eps, int MaxIter);

/* Функция решения трехдиагональной СЛАУ большой размерности методом Зейделя */
template <typename T>
vector<T> method_Zeidel_diag(vector<T> A, vector<T> B, vector<T> C, vector<T> b, vector<T> x0, T eps, int maxIterations);

/* Функция решения трехдиагональной СЛАУ большой размерности методом Релаксации */
template <typename T>
vector<T> method_Relax_diag(vector<T> A, vector<T> B, vector<T> C, vector<T> b, vector<T> x0, T w, T eps, int MaxIter);
/* Функция представления матрицы С в виде: C = C_l + C_d + D_u */

/* Функция исследования итерационного параметра */

/* Функция исследования сходимости при различных начальных приближениях */
