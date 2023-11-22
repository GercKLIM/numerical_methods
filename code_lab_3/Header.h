#pragma once

//const string File = "C:\Users\al-ba\Desktop\Лабы вычи\3 лаба\Lab_3\Lab_3";
const T pi = 3.141592653589793;

//интерполируемые функции
T func1(const T x) {
	return pow(x, 2);
};
T func2(const T x) {
	return 1 / (1 + 25 * pow(x, 2));
};
T func3(const T x) {
	return 1 / (1 + atan(1 + 10 * pow(x, 2)));
};
T func4(const T x) {
	T mnogochlen = (4 * pow(x, 3) + 2 * pow(x, 2) - 4 * x + 2);
	T R1 = pow(mnogochlen, sqrt(2));
	T R2 = asin(1 / (5 + x - x * x));
	return R1 + R2 - 5;
};
T func5(const T x) {
	T R1 = tan(pow(x, 2) * 5 + 3 * x - 2);
	T R2 = exp((pow(x, 3) + pow(x, 2) * 6 + 12 * x + 8) / (pow(x, 2) * 2 + 8 * x + 7));
	return R1 + R2 - 2.0;
};
T func6(const T x) {
	T R1 = tan((-3 * pow(x, 2) - pow(5, 1 / 3) * x + log(2)) / sqrt(19));
	T R2 = log((-sqrt(2) * pow(x, 4) - sqrt(6) * pow(x, 3) + 4 * x + 4 * sqrt(3) - 1) / (x + sqrt(3)));

	return R1 + R2 - 2.2;
};
T func7(const T x) {
	return 1;
};
T func8(const T x) {
	if (x < 0) return -pow(x, 3);
	if (x >= 0) return pow(x, 3);
};
T func9(const T x) {
	if ((x <= pi) && (x > 0)) return sin(x);
	if ((x >= -2 * pi) && (x <= 0)) return sin(2 * x);
};
T func10(const T x) {
	return exp(x);
};

typedef T(*func) (T x);

//для выбора функции
func funcD[10] = { func1,func2,func3,func4,func5,func6,func7,func8,func9,func10 };

//построение многочлена Лагранжа
T Lagrange(int n, T x, vector<pair<T, T>>& vec) {

	T w = 1.0, L = 0.0;
	for (int i = 0; i <= n; i++) {
		for (int j = 0; j <= n; j++) {
			if (j != i) {
				w *= (x - vec[j].first) / (vec[i].first - vec[j].first);
				//cout << w << endl;
			}
		}
		L = L + w * vec[i].second;
		w = 1.0;
	}
	return L;
}

void Progonka(int n, vector<T>& L, vector<T>& D, vector<T>& U, vector<T>& R, vector<T>& c)
{
	U.push_back(0);                 // U[n] = 0
	L.insert(L.begin(), 0);        // L[0] = 0

	vector<T> alpha, betta;
	T gamma;

	// Прямой ход
	gamma = D[0];
	alpha.push_back(-U[0] / gamma);
	betta.push_back(R[0] / gamma);
	for (int i = 1; i <= n - 1; ++i)
	{
		gamma = D[i] + L[i] * alpha[i - 1];
		alpha.push_back(-U[i] / gamma);
		betta.push_back((R[i] - L[i] * betta[i - 1]) / gamma);
	}
	gamma = D[n - 1] + L[n - 1] * alpha[n - 2];
	betta.push_back((R[n - 1] - L[n - 1] * betta[n - 2]) / gamma);

	// Обратный ход
	c.push_back(betta[n]);
	for (int i = n - 1; i >= 0; --i)
		c.insert(c.begin(), alpha[i] * c[0] + betta[i]);
}

void SplineCoeff(int n, vector<pair<T, T>>& vec, vector<T>& a, vector<T>& b, vector<T>& c, vector<T>& d) {

	a.clear();
	b.clear();
	c.clear();
	d.clear();

	vector<T> L, D, U, R;
	vector<T> h, g;

	D.push_back(1);
	R.push_back(0);   // c[0] = 0

	for (int i = 0; i < n; ++i)
	{
		a.push_back(vec[i].second);
		h.push_back(vec[i + 1].first - vec[i].first);                   // h[i] = x[i] - x[i-1]                i = 0, 1, ..., n-1
		g.push_back((vec[i + 1].second - vec[i].second) / h.back());      // g[i] = (y[i] - y[i-1]) / h[i]       i = 0, 1, ..., n-1
		if (i >= 1)
		{
			D.push_back(2 * (h[i - 1] + h[i]));            // D[i] = 2*(h[i - 1] + h[i])       i = 1, 2, ..., n-1
			R.push_back(3 * (g[i] - g[i - 1]));       // R[i] = 3*(g[i] - g[i - 1])  i = 1, 2, ..., n-1
		}
	}
	L = h;
	L[L.size() - 1] = 0;
	U = h;
	U[0] = 0;
	D.push_back(1);
	R.push_back(0);
	Progonka(n, L, D, U, R, c);

	for (int i = 0; i <= n - 1; ++i)
	{
		b.push_back(g[i] - (c[i + 1] + 2 * c[i]) * h[i] / 3);
		d.push_back((c[i + 1] - c[i]) / (3 * h[i]));
	}
	c.pop_back();
}

// Левая 3-х диагональная прогонка
void LeftTridiagRun(const int n, vector<T>& x, const vector<T>& a, const vector<T>& b, const vector<T>& c, const vector<T>& d) {

	vector<T> xi(n);
	vector<T> eta(n);
	for (auto i = 0; i < n - 1; ++i) {
		double denom = b[i] - c[i] * xi[i + 1];
		xi.push_back(a[i] / denom);
		eta.push_back((d[i] + c[i] * eta[i + 1]) / denom);
	}
	xi.push_back(a[n - 1] / b[n - 1]);
	eta.push_back(d[n - 1] / b[n - 1]);
	x.push_back(0); // (d[0] + c[0] * eta[1]) / (b[0] - c[0] * xi[1]);
	for (auto i = 1; i < n - 1; ++i) {
		x.push_back(xi[i] * x[i - 1] + eta[i]);
	}
	x.push_back(0);
}

void CubicPolyCoeff(const int n, const vector<pair<T, T>>& vec, vector<T>& a, vector<T>& b, vector<T>& c, vector<T>& d) {
	vector<T> h(n), g(n);
	vector<T> A(n), B(n), C(n), D(n);

	for (auto i = 0; i < n; ++i) {
		a.push_back(vec[i].second);
		h.push_back(vec[i + 1].first - vec[i].first);
		g.push_back((vec[i + 1].second - vec[i].second) / h[i]);
	}
	// Коэффициенты СЛАУ в прогонке

	for (auto i = 1; i < n; ++i) {
		A.push_back(h[i - 1]);
		B.push_back(-2 * (h[i - 1] + h[i]));
		C.push_back(h[i]);
		D.push_back(-3 * (g[i] - g[i - 1]));
	}
	LeftTridiagRun(n, c, A, B, C, D);
	c.push_back(0);
	for (auto i = 0; i < n; ++i) {
		b.push_back(g[i] - (c[i + 1] + 2 * c[i]) * h[i] / 3);
		d.push_back((c[i + 1] - c[i]) / (3 * h[i]));
	}
}


T Spline(int n, T x, vector<pair<T, T>>& vec, vector<T>& a, vector<T>& b, vector<T>& c, vector<T>& d) {

	T Si = 0, xpred = vec[0].first;
	for (int i = 0; i < n; ++i)
	{
		if ((x <= vec[i + 1].first) && (x > xpred))
		{
			Si = a[i] + b[i] * (x - xpred) + c[i] * pow((x - xpred), 2) + d[i] * pow((x - xpred), 3);
			break;
		}
		xpred = vec[i + 1].first;
	}
	return Si;
}

