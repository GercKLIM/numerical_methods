#include <vector>
using namespace std;
#include <math.h>


/* Функция для вычисления скалярного произведения векторов */
template <typename T>
T dotProduct(const vector<T>& a, const vector<T>& b) {
    T result = 0;
    for (int i = 0; i < a.size(); i++) {
        result += a[i] * b[i];
    }
    return result;
}


/* Функция для получения длины вектора */
template <typename T>
double vectorLength(const vector<T> v) {
    return sqrt(dotProduct(v, v));
}


/* Функция для ортогонализации векторов методом Грама-Шмидта */
template <typename T>
vector<T> gramSchmidt(const vector<T>& v, const vector<vector<T>>& basis) {
    vector<T> result = v;
    for (const vector<T>& b : basis) {
        double projection = dotProduct(v, b) / dotProduct(b, b);
        for (int i = 0; i < v.size(); i++) {
            result[i] -= projection * b[i];
        }
    }
    return result;
}


/* Функция для QR-разложения матрицы */
template <typename T>
void qrDecomposition(const vector<vector<T>>& A, vector<vector<T>>& Q, vector<vector<T>>& R) {
    int n = A.size();
    Q = vector<std::vector<T>>(n, std::vector<T>(n, 0));
    R = vector<std::vector<T>>(n, std::vector<T>(n, 0));

    vector<vector<T>> basis;

    for (int j = 0; j < n; j++) {
        vector<T> v;
        for (int i = 0; i < n; i++) {
            v.push_back(A[i][j]);
        }

        for (int i = 0; i < j; i++) {
            vector<T> b;
            for (int k = 0; k < n; k++) {
                b.push_back(Q[k][i]);
            }
            v = gramSchmidt(v, basis);
        }

        basis.push_back(v);
        double length = vectorLength(v);

        for (int i = 0; i < n; i++) {
            Q[i][j] = v[i] / length;
        }

        for (int i = j; i < n; i++) {
            R[j][i] = dotProduct(v, basis[i - j]);
        }
    }
}

template <typename T>
vector<vector<T>> qrDecomposition_Q(const vector<vector<T>> A){
    vector<vector<T>> Q, R;
    qrDecomposition(A, Q, R);
    return Q;
}

template <typename T>
vector<vector<T>> qrDecomposition_R(const vector<vector<T>> A){
    vector<vector<T>> Q, R;
    qrDecomposition(A, Q, R);
    return R;
}


/* Функция для решения СЛАУ методом QR-разложения */
template <typename T>
vector<T> method_QR(const vector<vector<T>>& A, const vector<T>& b) {
    int n = A.size();
    vector<vector<T>> Q, R;
    qrDecomposition(A, Q, R);

    // Решение системы уравнений Rx = Q^T * b
    vector<T> Qt_b(n, 0);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            Qt_b[i] += Q[j][i] * b[j];
        }
    }

    vector<T> x(n, 0);

    for (int i = n - 1; i >= 0; i--) {
        x[i] = Qt_b[i];
        for (int j = i + 1; j < n; j++) {
            x[i] -= R[i][j] * x[j];
        }
        x[i] /= R[i][i];
    }

    return x;
}