#pragma once

//const string File = "C:\Users\al-ba\Desktop\���� ����\3 ����\Lab_3\Lab_3";

template <typename T>
const T pi = 3.141592653589793;

// Заданные функции для интерполяции
template <typename T>
T func1(const T x) {
	return pow(x, 2);
}

template <typename T>
T func2(const T x) {
	return 1 / (1 + 25 * pow(x, 2));
}

template <typename T>
T func3(const T x) {
	return 1 / (1 + atan(1 + 10 * pow(x, 2)));
}

template <typename T>
T func4(const T x) {
	T mnogochlen = (4 * pow(x, 3) + 2 * pow(x, 2) - 4 * x + 2);
	T R1 = pow(mnogochlen, sqrt(2));
	T R2 = asin(1 / (5 + x - x * x));
	return R1 + R2 - 5;
}

template <typename T>
T func5(const T x) {
	T R1 = tan(pow(x, 2) * 5 + 3 * x - 2);
	T R2 = exp((pow(x, 3) + pow(x, 2) * 6 + 12 * x + 8) / (pow(x, 2) * 2 + 8 * x + 7));
	return R1 + R2 - 2.0;
}

template <typename T>
T func6(const T x) {
	T R1 = tan((-3 * pow(x, 2) - pow(5, 1 / 3) * x + log(2)) / sqrt(19));
	T R2 = log((-sqrt(2) * pow(x, 4) - sqrt(6) * pow(x, 3) + 4 * x + 4 * sqrt(3) - 1) / (x + sqrt(3)));

	return R1 + R2 - 2.2;
}

template <typename T>
T func7(const T x) {
	return 1;
}

template <typename T>
T func8(const T x) {
	if (x < 0) return -pow(x, 3);
	if (x >= 0) return pow(x, 3);
}

template <typename T>
T func9(const T x) {
	if ((x <= pi) && (x > 0)) return sin(x);
	if ((x >= -2 * pi) && (x <= 0)) return sin(2 * x);
}

template <typename T>
T func10(const T x) {
	return exp(x);
};


// Указатель на функцию, который передается методам

typedef double(*func) (double x);

// Массив функций
func funcD[10] = { func1,func2,func3,func4,func5,func6,func7,func8,func9,func10 };


// Функия возвращает
template <typename T>
T Lagrange(int n, T x, vector<pair<T, T>>& vec);


// Функция,
template <typename T>
void Progonka(int n, vector<T>& L, vector<T>& D, vector<T>& U, vector<T>& R, vector<T>& c);

template <typename T>
void SplineCoeff(int n, vector<pair<T, T>>& vec, vector<T>& a, vector<T>& b, vector<T>& c, vector<T>& d);

template <typename T>
void LeftTridiagRun(const int n, vector<T>& x, const vector<T>& a, const vector<T>& b, const vector<T>& c, const vector<T>& d);

template <typename T>
void CubicPolyCoeff(const int n, const vector<pair<T, T>>& vec, vector<T>& a, vector<T>& b, vector<T>& c, vector<T>& d);

