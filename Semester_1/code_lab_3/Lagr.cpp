#include <iostream>
using namespace std;

template <typename T>
T lagrange1(double* x, double* y, short n, double _x)
{
	double result = 0.0;

	for (short i = 0; i < n; i++)
	{
		double P = 1.0;

		for (short j = 0; j < n; j++)
			if (j != i)
				P *= (_x - x[j]) / (x[i] - x[j]);

		result += P * y[i];
	}

	return result;
}

template <typename T>
T lagrange2(double* x, double* y, short n, double _x)
{
	double result = 0.0;

	double h = x[1] - x[0];

	for (short i = 0; i < n; i++)
	{
		double P = 1.0;

		for (short j = 0; j < n; j++)
			if (i != j)
				P *= (_x - x[0] - h * j) / h / (i - j);

		result += P * y[i];
	}

	return result;
}



