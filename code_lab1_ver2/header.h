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


/* Функция вывода матрицы на экран
 * (работает вне зависимости от размера матрицы)
 * */
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


/* Функция для вычисления обратной матрицы методом Гаусса-Жордана */
template <typename T>
vector<vector<T>> inverseMatrix(vector<vector<double>> matrix);


/* Функция для решения СЛАУ прямым методом Гаусса */
template <typename T>
vector<double> method_Gaussa(const vector<vector<T>>& A, const vector<T>& b);

/* Функция для вычисления скалярного произведения векторов */
template <typename T>
T dotProduct(const vector<T>& a, const vector<T>& b);

/* Функция для получения длины вектора */
template <typename T>
double vectorLength(vector<T> v);

/* Функция для ортогонализации векторов методом Грама-Шмидта */
template <typename T>
vector<T> gramSchmidt(const vector<T>& v, const vector<vector<T>>& basis);


/* Функция для QR-разложения матрицы */
template <typename T>
void qrDecomposition(vector<vector<T>>& A, vector<vector<T>>& Q, vector<vector<T>>& R);


/* Функция для решения СЛАУ методом QR-разложения */
template <typename T>
vector<T> method_QR(vector<vector<T>> A, vector<T> b);


/* Функция для вывода погрешности решения СЛАУ для определенного варианта
template <typename T>
T eval_eps(vector<T> true_solve, vector<T> numerical_solve);
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

/* Функция для вычисления оо-нормы матрицы */
template <typename T>
T norm_oo(vector<vector<T>> A);


/* Функция для вычисления числа обусловленности матрицы */
template <typename T>
T cond(vector<vector<T>> A);

/* Функция для вычисления нормы вектора невязки */
template <typename T>
T norm_vector_nevazki(vector<T> true_solve, vector<T> numerical_solve);

/* Функция для cложения векторов */
template <typename T>
vector<T> vec_sum(vector<T> vec1, vector<T> vec2);

/* Функция для оценки изменения числа обуcловленности от возмущения вектора правой части */

template <typename T>
T evaluate_change_cond(vector<vector<T>> matrix, vector<T> vec, vector<T> mod);
