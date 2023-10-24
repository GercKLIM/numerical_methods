#pragma once
#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include "source.cpp"
#include "method_gaussa.cpp"
#include "method_QR.cpp"
#include <math.h>

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
void qrDecomposition(const vector<vector<T>>& A, vector<vector<T>>& Q, vector<vector<T>>& R);


/* Функция для решения СЛАУ методом QR-разложения */
template <typename T>
vector<T> method_QR(const vector<vector<T>>& A, const vector<T>& b);
