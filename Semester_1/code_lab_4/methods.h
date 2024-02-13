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
vector<T> method_Gaussa(const vector<vector<T>>& matrix, const vector<T>& vec, const T& eps);
template <typename T>
vector<T> method_Gaussa2(const vector<vector<T>>& matrix, const vector<T>& vec, const T& eps);

/* Функция QR-разложения матрицы методом вращений */
template<typename T>
void QR_decomposition(const vector<vector<T>>& matrix, vector<vector<T>>& Q, vector<vector<T>>& R, const T& eps);


/* Функция для решения СЛАУ методом QR-разложения */
template <typename T>
vector<T> method_QR(const vector<vector<T>>& A, const vector<T>& b, const T& eps);


/* Функция для оценки изменения числа обуcловленности от возмущения вектора правой части */
template <typename T>
void min_change_cond(const vector<vector<T>>& matrix, const vector<T>& vec, const vector<T>& mod);


/* ### Функций лабы 2 ### */

/* Структура, с помощью которой будет выводится результат метода лаб 3 */
template<typename T>
struct MyResult2 {
    vector<T> solve;
    int iterations;
    vector<vector<T>> C;
    vector<T> y;
    T batch;

};


/* Функция преобразования матрицы в сумму из Нижнетреугольной, Диагональной, Верхнетреугольной */
template<typename T>
void LDU_decomposotion(const vector<vector<T>>& A, vector<vector<T>>& L, vector<vector<T>>& D, vector<vector<T>>& U);


/* Функция решения СЛАУ методом Простой Итерации */
template<typename T>
MyResult2<T> method_SimpleIteration(const vector<vector<T>>& A, const vector<T>& b, const vector<T>& x0, const T& tau, const T& eps, const int& p, const int& MaxIter);


/* Функция решения СЛАУ методом Якоби */
template<typename T>
MyResult2<T> method_Yacobi(const vector<vector<T>>& A, const vector<T>& b, const vector<T>& x0, const T& eps, const int& p, const int& MaxIter);


/* Функция решения СЛАУ методом Релаксации */
template<typename T>
MyResult2<T> method_Relax(const vector<vector<T>>& A, const vector<T>& b, const vector<T>& x0, const T& w, const T& eps, const int& p, const int& MaxIter);


/* Функция решения СЛАУ методом Зейделя */
template<typename T>
MyResult2<T> method_Zeidel(const vector<vector<T>>& A, const vector<T>& b, const vector<T>& x0, const T& eps, const int& p, const int& MaxIter);


/* Функция решения трехдиагональной СЛАУ большой размерности методом Зейделя */
template <typename T>
vector<T> method_Zeidel_diag(const vector<T>& A, const vector<T>& B, const vector<T>& C, const vector<T>& b, const vector<T>& x0, const T& eps, const T& maxIterations);


/* Функция решения трехдиагональной СЛАУ большой размерности методом Релаксации */
template <typename T>
vector<T> method_Relax_diag(const vector<T>& A, const vector<T>& B, const vector<T>& C, const vector<T>& b, const vector<T>& x0, const T& w, const T& eps, const T& MaxIter);


/* Функция для вычисления нормы вектора невязки трехдиагональной СЛАУ */
template <typename T>
T norm_vector_nevazki(const vector<T>& A, const vector<T>& B, const vector<T>& C, const vector<T>& b, const vector<T>& solution, const int& n);


/* Функция исследования итерационного параметра tau для метода простых итераций (Метод Золотого сечения)*/
template<typename T>
T SimpleIterations_method_matrix_norm_C(const vector<vector<T>>& A, const T& tau, const int& p);


// Метод золотого сечения для поиска минимума функции на заданном интервале [a, b]
template<typename T>
T golden_section_search_tau(const vector<vector<T>>& A, T a, T b, const int& p, const T& epsilon);


/* Функция исследования итерационного параметра W для метода Релаксации (Метод Золотого сечения)*/
template<typename T>
T golden_section_search_W(const vector<vector<T>>& A, T a, T b, const int& p, const T& eps);


/* Функция от которой ищется минимум в золотом сечении для релаксации */
template <typename T>
T C_matrix_for_relax(const vector<vector<T>>& A, const T& w, const int& p);


/* Функция исследования итерационного параметра W для метода Релаксации для трехдиагональной матрицы (Метод Золотого сечения)*/
template<typename T>
T golden_section_search_W(vector<T> A, vector<T> B, vector<T> C, vector<T> vec, vector<T> x, T EPS, int MaxIteration, T a, T b);


/* Функция априорной оценки */
template <typename T>
void aprior_eps(const vector<vector<T>>& C, const vector<T>& y, const vector<T>& x0, const int& p);


/* Функция апостериорной оценки */
template <typename T>
void aposter_eps(const vector<vector<T>>& C, T norm_delta, const int& p);



/* ### Функций лабы 4 ### */

/* Структура для вывода данных 4 лабы*/
template <typename T>
struct MyResult4{
    vector<T>  eigen;
    int iterations = 0;
    vector<vector<T>> eigens_vec;
    vector<vector<T>> R;
    vector<vector<vector<T>>> A_iter;
};

/* Функция приведения матрицы к форме Хессенберга методом вращений */
template <typename T>
vector<vector<T>> Hessenberg_decomposition(const vector<vector<T>>& matrix);


/* Функция нахождения собственных значений матрицы методом QR-разложения за одну итерацию */
template <typename T>
MyResult4<T> Eigen_method_QR(const vector<vector<T>>& A, const T& eps, const int& maxIterations);

/* Функция нахождения собственных значений матрицы методом QR-разложения */
template <typename T>
MyResult4<T> Eigen_method_QR(const vector<vector<T>>& matrix, const T& sigma, const T& eps, const int& maxIterations);

template <typename T>
MyResult4<T> Eigen_method_QR2(const vector<vector<T>>& matrix, const T& eps, const int& maxIterations);

template <typename T>
MyResult4<T> Eigen_method_QR3(const vector<vector<T>>& matrix, const T& eps, const int& maxIterations);
template <typename T>
MyResult4<T> Eigen_method_QR3(const vector<vector<T>>& matrix, const T&sigma,const T& eps, const int& maxIterations);

/* Функция нахождения собственных векторов матрицы методом Обратных Итераций */
template<typename T>
MyResult4<T> reverse_iteration(const vector<vector<T>>& matrix, const vector<T>& lambda, const T& eps, const int& maxIteration);

/* Функция нахождения собственных значений и собственных векторов методом Обратных Итераций
 * с использованием отношения Рэлея (Модификация метода Обратных Итераций) */

//template <typename T>
//MyResult4<T> reverse_iterator_with_reley(const vector<vector<T>>& matrix, const vector<vector<T>>& X0, const T eps, const int& maxIteration);

template <typename T>
MyResult4<T> reverse_iterator_with_reley(const vector<vector<T>>& matrix, const vector<vector<T>>& X0, const T eps, const int& maxIteration);