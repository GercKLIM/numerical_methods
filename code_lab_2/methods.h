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
void lu_decomposition(const vector<vector<T>>& A, vector<vector<T>>& L, vector<vector<T>>& U);


/* Функция для вычисления нормы вектора невязки */
template <typename T>
T norm_vector_nevazki(const vector<vector<T>>& A, const vector<T>& b, const vector<T>& x, const int& n);


/* Функция для решения СЛАУ прямым методом Гаусса */
template <typename T>
vector<T> method_Gaussa(const vector<vector<T>>& matrix, const vector<T>& vec);


/* Функция для QR-разложения матрицы */
template <typename T>
void QR_decomposition(const vector<vector<T>>& A, vector<vector<T>>& Q, vector<vector<T>>& R);


/* Функция для решения СЛАУ методом QR-разложения */
template <typename T>
vector<T> method_QR(const vector<vector<T>>& A, const vector<T>& b);


/* Функция для оценки изменения числа обуcловленности от возмущения вектора правой части */
template <typename T>
void min_change_cond(const vector<vector<T>>& matrix, const vector<T>& vec, const vector<T>& mod);


/* ### Функций лабы 2 ### */

/* Структура, с помощью которой будет выводится результат метода */
template<typename T>
struct Result {
    vector<T> solve;
    int iterations;
    vector<vector<T>> C;
    T norm_C;
};

/* Функция преобразования матрицы в сумму из Нижнетреугольной, Диагональной, Верхнетреугольной */
template<typename T>
void LDU_decomposotion(const vector<vector<T>>& A, vector<vector<T>>& L, vector<vector<T>>& D, vector<vector<T>>& U);


/* Функция решения СЛАУ методом Простой Итерации */
template<typename T>
Result<T> method_SimpleIteration(const vector<vector<T>>& A, const vector<T>& b, const vector<T>& x0, const T& tau, const T& eps, const int& p, const int& MaxIter);


/* Функция решения СЛАУ методом Якоби */
template<typename T>
Result<T> method_Yacobi(const vector<vector<T>>& A, const vector<T>& b, const vector<T>& x0, const T& eps, const int& p, const int& MaxIter);


/* Функция решения СЛАУ методом Релаксации */
template<typename T>
Result<T> method_Relax(const vector<vector<T>>& A, const vector<T>& b, const vector<T>& x0, const T& w, const T& eps, const int& MaxIter);


/* Функция решения СЛАУ методом Зейделя */
template<typename T>
Result<T> method_Zeidel(const vector<vector<T>>& A, const vector<T>& b, const vector<T>& x0, const T& eps, const int& MaxIter);


/* Функция решения трехдиагональной СЛАУ большой размерности методом Зейделя */
template <typename T>
vector<T> method_Zeidel_diag(const vector<T>& A, const vector<T>& B, const vector<T>& C, const vector<T>& b, const vector<T>& x0, const T& eps, const T& maxIterations);


/* Функция решения трехдиагональной СЛАУ большой размерности методом Релаксации */
template <typename T>
vector<T> method_Relax_diag(const vector<T>& A, const vector<T>& B, const vector<T>& C, const vector<T>& b, const vector<T>& x0, const T& w, const T& eps, const T& MaxIter);


/* Функция для вычисления нормы вектора невязки трехдиагональной СЛАУ */
template <typename T>
T norm_vector_nevazki(const vector<T>& A, const vector<T>& B, const vector<T>& C, const vector<T>& b, const vector<T>& solution, const int& n);

