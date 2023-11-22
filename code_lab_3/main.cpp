#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

#define T double

using namespace std;

#include "Header.h"
void test(){
	int n = 20; //число узлов
	int N = 200; //число точек на графике
	T A = -2 * pi, B = pi;//отрезок интерполирования
	//T A = -1, B = 1;//отрезок интерполирования

	T x;
	vector<T> x_i, x_i2;
	vector<pair<T, T>> Cravnomer, Cravnomer2, CChebyshov;

	func f = funcD[8];

	ofstream file1, file2, file3;
	file1.open("CoordRavnomer.txt");
	file2.open("CoordChebyshev.txt");

	//формируем сетки
	for (int j = 0; j <= n; ++j)
	{
		//равномерная сетка
		x = A + (B - A) * j / n;
		Cravnomer.push_back(make_pair(x, f(x)));
		file1 << x << " " << f(x) << endl;

		//чебышевская сетка
		x = (A + B) / 2 + ((B - A) / 2) * cos(((2 * (n - j) + 1) * pi) / (2 * (n + 1)));
		CChebyshov.push_back(make_pair(x, f(x)));
		file2 << x << " " << f(x) << endl;
	}
	file1.close();
	file2.close();

	//точные значения функции. создание файла для рисования в Wolfram
	file1.open("Func.txt");
	T h = (B - A) / (N - 1);
	x = A;
	for (int j = 0; j < N; ++j)
	{
		file1 << x << " " << f(x) << endl;
		x = x + h;
	}

	file1.close();

	//интерполяция многочленом Лагранжа. создание файлов для рисования в Wolfram
	file1.open("LagrangeRavnomer.txt");
	file2.open("LagrangeChebyshev.txt");
	h = (B - A) / (N - 1);
	x = A;
	for (int j = 0; j < N; ++j)
	{
		file1 << x << " " << Lagrange(n, x, Cravnomer) << endl;
		file2 << x << " " << Lagrange(n, x, CChebyshov) << endl;
		x = x + h;

	}
	file1.close();
	file2.close();

	//интерполяция сплайном
	vector<T> a, b, c, d; //коэффициенты сплайна
	file1.open("SplineRavnomer.txt");
	file2.open("SplineChebyshev.txt");
	SplineCoeff(n, Cravnomer, a, b, c, d);
	h = (B - A) / (N - 1);
	x = A;
	for (int j = 0; j < N - 1; ++j)
	{
		file1 << x << " " << Spline(n, x, Cravnomer, a, b, c, d) << endl;
		x = x + h;
	}
	file1.close();

	//сплайн для чебышевской сетки
	x = A;
	SplineCoeff(n, CChebyshov, a, b, c, d);
	for (int j = 0; j < N - 1; ++j)
	{
		file2 << x << " " << Spline(n, x, CChebyshov, a, b, c, d) << endl;
		x = x + h;

	}
	file2.close();



	//сплайн для оценки(=4)
	//h = 1e-2;
	//n = int((B-A)/h);
	Cravnomer.clear(); Cravnomer2.clear();
	vector<T> a1, b1, c1, d1;
	n = 500;
	h = (B - A) / n;
	x = A;
	for (int j = 0; j <= n; ++j)
	{
		x = A + (B - A) * j / n;
		Cravnomer.push_back(make_pair(x, f(x)));
	}
	for (int j = 0; j <= 2 * n; ++j)
	{
		x = A + (B - A) * j / (2 * n);
		Cravnomer2.push_back(make_pair(x, f(x)));
	}
	file3.open("SplineRavn_h.txt"); //для сравнения (=4)
	file1.open("SplineRavn_0.5h.txt");
	SplineCoeff(n, Cravnomer, a, b, c, d);
	SplineCoeff(2 * n, Cravnomer2, a1, b1, c1, d1);
	x = A;
	for (int j = 0; j < n; ++j)
	{
		x_i.push_back(x + (h / 3));
		file3 << abs(f(x_i[j]) - Spline(n, x_i[j], Cravnomer, a, b, c, d)) << endl;
		file1 << abs(f(x_i[j]) - Spline(2 * n, x_i[j], Cravnomer2, a1, b1, c1, d1)) << endl;
		x = x + h;
	}
	file3.close();
	file1.close();
}

int main(){
    cout << ("Hello, World!");
}