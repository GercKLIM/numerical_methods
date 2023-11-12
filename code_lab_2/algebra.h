#pragma once

/* algebra.h и algebra.cpp - основательные файлы,
 * в которых соответственно объявляются и реализовываются функции
 * для работы с матрицами, векторами и т.д.
 * */


#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <cmath>

using namespace std;

/* *** Начальные функции для испорта/экспорта данных *** */

/* Функция импорта матрицы из текстового файла*/
template <typename T>
vector<vector<T>> importSLAU(const string& filename);

/* Функция вывода матрицы на экран */
template <typename T>
void print(vector<vector<T>> matrix);


/* Функция вывода вектора на экран */
template <typename T>
void print(vector<T> vec);


/* Функция для получения матрицы из СЛАУ */
template <typename T>
vector<vector<T>> SLAU_to_matrix(vector<vector<T>> SLAU);


/* Функция для получения векторая из СЛАУ */
template <typename T>
vector<T> SLAU_to_vec(vector<vector<T>> SLAU);



/* *** Функции математики векторов *** */

/* Операция cложения векторов */
template <typename T>
vector<T> operator+(vector<T> vec1, vector<T> vec2);

/* Операция вычитания векторов */
template <typename T>
vector<T> operator-(vector<T> vec1, vector<T> vec2);

/* Операция почленного умножения векторов */
template <typename T>
vector<T> operator*(vector<T> vec1, vector<T> vec2);

/* Операция умножения вектора на число */
template <typename T>
vector<T> operator*(T c, vector<T> vec2);

template <typename T>
vector<T> operator*(vector<T> vec2, T c);



/* Операция почленного деления векторов */
template <typename T>
vector<T> operator/(vector<T> vec1, vector<T> vec2);

/* Функция для скалярного умножения векторов */
template <typename T>
T dot(vector<T> vec1, vector<T> vec2);


/* Функция для нормы-1 вектора */
template <typename T>
T norm_1(vector<T> vec);


/* Функция для нормы-2 вектора */
template <typename T>
T norm_2(vector<T> vec);


/* Функция для нормы-oo вектора */
template <typename T>
T norm_oo(vector<T> vec);

/* *** Функции математики матриц *** */

/* Матричное умножение */
template <typename T>
vector<vector<T>> operator*(vector<vector<T>> A, vector<vector<T>> B);

/* Операция для умножения матрицы на число */
template <typename T>
vector<vector<T>> operator*(vector<vector<T>> A,  T scalar);

template <typename T>
vector<vector<T>> operator*(T scalar, vector<vector<T>>);

template <typename T>
vector<vector<T>> operator*(vector<vector<T>> A,  const T scalar);

/* Функция поэлементного сложения матриц */
template <typename T>
vector<vector<T>> operator+(vector<vector<T>> A, vector<vector<T>> B);


/* Функция поэлементного вычитания матриц */
template <typename T>
vector<vector<T>> operator-(vector<vector<T>> A, vector<vector<T>> B);


/* Функция для умножения матрицы на вектор */
template <typename T>
vector<T> operator*(vector<vector<T>> matrix, vector<T> vec);

/* Функция для транспонирования матрицы */
template <typename T>
vector<vector<T>> transpose(const vector<vector<T>>& A);


/* Функция для создания единичной матрицы размера n x n */
template <typename T>
vector<vector<T>> create_identity_matrix(int n);


/* Функция для поэлементного умножения матриц */
template <typename T>
vector<vector<T>> Multyply(vector<vector<T>> A, vector<vector<T>> B);


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


/* Функция поворота матрицы вправо */
template <typename T>
vector<vector<T>> RotateRight(vector<vector<T>> A);


/* Функция поворота матрицы влево */
template <typename T>
vector<vector<T>> RotateLeft(vector<vector<T>> A);


// Функция для обратной матрицы с проверкой на вырожденность
template <typename T>
vector<vector<T>> inverseMatrix(vector<vector<T>> A);

