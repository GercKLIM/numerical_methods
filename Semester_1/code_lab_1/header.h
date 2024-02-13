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
void print(const vector<vector<T>>& matrix);

/* Функция вывода вектора на экран */
template <typename T>
void print(const vector<T>& vec);


/* Функция для получения матрицы из СЛАУ */
template <typename T>
vector<vector<T>> SLAU_to_matrix(const vector<vector<T>>& SLAU);


/* Функция для получения векторая из СЛАУ */
template <typename T>
vector<T> SLAU_to_vec(const vector<vector<T>>& SLAU);


/* Функция для транспонирования матрицы */
template <typename T>
vector<vector<T>> transpon(const vector<vector<T>>& matrix);


/* Функция для создания единичной матрицы размера n x n */
template <typename T>
vector<vector<T>> create_identity_matrix(int n);


/* Функция для LU-разложения с частичным выбором */
template <typename T>
void lu_decomposition(const vector<vector<T>>& A, vector<vector<T>>& L, vector<vector<T>>& U);

/* Функция поворота матрицы вправо */
template <typename T>
vector<vector<T>> MatrixRotateRight(const vector<vector<T>>& A);

/* Функция поворота матрицы влево */
template <typename T>
vector<vector<T>> MatrixRotateLeft(const vector<vector<T>>& A);

/* Функция для вычисления обратной матрицы методом Гаусса-Жордана */
template <typename T>
vector<vector<T>> inverseMatrix(const vector<vector<T>>& matrix);


/* Операции с матрицами и векторами */

/* Функция для скалярного умножения векторов */
template <typename T>
T dot_vec(const vector<T>& a, const vector<T>& b);

/* Функция для умножения матриц */
template <typename T>
vector<vector<T>> MatrixMultiply(const vector<vector<T>>& A, const vector<vector<T>>& B);


/* Функция округления чисел в матрицах */
template <typename T>
vector<vector<T>> Matrix_round(const vector<vector<T>>& A, const double& eps);


/* Функция для вычисления 1-нормы матрицы */
template <typename T>
T norm_1(const vector<vector<T>>& A);

/* Функция для вычисления 2-нормы матрицы */
template <typename T>
T norm_2(const vector<vector<T>>& A);

/* Функция для вычисления оо-нормы матрицы */
template <typename T>
T norm_oo(const vector<vector<T>>& A);


/* Функция для вычисления числа обусловленности матрицы c нормой 1*/
template <typename T>
T cond_1(const vector<vector<T>>& matrix);

/* Функция для вычисления числа обусловленности матрицы c нормой 2*/
template <typename T>
T cond_2(const vector<vector<T>>& matrix);

/* Функция для вычисления числа обусловленности матрицы c нормой oo*/
template <typename T>
T cond_oo(const vector<vector<T>>& matrix);


/* Функция для вычисления нормы-1 вектора невязки */
template <typename T>
T norm_vector_nevazki(const vector<vector<T>>& A, const vector<T>& b, const vector<T>& x, const int& n);


/* Функция для нормы-1 вектора */
template <typename T>
T norm_1(const vector<T>& vec);
/* Функция для нормы-2 вектора */
template <typename T>
T norm_2(const vector<T>& vec);

/* Функция для нормы-oo вектора */
template <typename T>
T norm_oo(const vector<T>& vec);


/* Функция для cложения векторов */
template <typename T>
vector<T> vec_sum(const vector<T>& vec1, const vector<T>& vec2);
