#pragma once
#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <cmath>
#include "source.cpp"
#include "method_gaussa.cpp"
#include "method_QR.cpp"


using namespace std;


/* Функция импорта матрицы из текстового файла*/
template <typename T>
vector<vector<T>> importSLAU(const string& filename);


/* Функция вывода матрицы на экран */

template <typename T>
void print(vector<vector<T>> matrix);

template <typename T>
void print(vector<T> vec);


/* Функция для получения матрицы из СЛАУ */
template <typename T>
vector<vector<T>> SLAU_to_matrix(vector<vector<T>> SLAU);


/* Функция для получения векторая из СЛАУ */
template <typename T>
vector<T> SLAU_to_vec(vector<vector<T>> SLAU);


/* Функция для транспонирования матрицы */
template <typename T>
vector<vector<T>> transpon(vector<vector<T>> matrix);

/* Функция для создания единичной матрицы размера n x n */
template <typename T>
vector<vector<T>> create_identity_matrix(int n);

/* Функция для LU-разложения с частичным выбором */
template <typename T>
void lu_decomposition(vector<vector<T>>& A, vector<vector<T>>& L, vector<vector<T>>& U);


/* Функция для вычисления обратной матрицы методом Гаусса-Жордана */
template <typename T>
vector<vector<T>> inverseMatrix(vector<vector<double>> matrix);


/* Функция для решения СЛАУ методом QR-разложения
template <typename T>
vector<T> method_QR(vector<vector<T>> A, vector<T> b);
*/

/* Функция для умножения матриц */
template <typename T>
vector<vector<T>> MatrixMultiply(vector<vector<T>> A, vector<vector<T>> B);

/* Функция округления чисел в матрицах */
template <typename T>
vector<vector<T>> Matrix_round(vector<vector<T>> A, double eps);

/* Функция для вычисления 1-нормы матрицы */
template <typename T>
T norm_1(vector<vector<T>> A);

/* Функция для вычисления 2-нормы матрицы */
template <typename T>
T norm_2(vector<vector<T>> A);

/* Функция для вычисления оо-нормы матрицы */
template <typename T>
T norm_oo(vector<vector<T>> A);


/* Функция для вычисления числа обусловленности матрицы c нормой 1*/
template <typename T>
T cond_1(vector<vector<T>> matrix);

/* Функция для вычисления числа обусловленности матрицы c нормой 2*/
template <typename T>
T cond_2(vector<vector<T>> matrix);

/* Функция для вычисления числа обусловленности матрицы c нормой oo*/
template <typename T>
T cond_oo(vector<vector<T>> matrix);


/* Функция для вычисления нормы-1 вектора невязки */
template <typename T>
T norm_vector_nevazki(vector<vector<T>> A, vector<T> b, vector<T> x, const int n);


/* Функция для скалярного умножения векторов */
template <typename T>
T dot_vec(vector<T> a, vector<T> b);

/* Функция для нормы-1 вектора */
template <typename T>
T norm_1(vector<T> vec);

/* Функция для нормы-2 вектора */
template <typename T>
T norm_2(vector<T> vec);

/* Функция для нормы-oo вектора */
template <typename T>
T norm_oo(vector<T> vec);


/* Функция для cложения векторов */
template <typename T>
vector<T> vec_sum(vector<T> vec1, vector<T> vec2);
