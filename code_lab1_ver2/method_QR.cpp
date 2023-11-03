#include <iostream>
#include <vector>
#include <cmath>


using namespace std;


template <typename T>
vector<vector<T>> transpose(const vector<vector<T>>& A) {
    int rows = A.size();
    int cols = A[0].size();
    vector<vector<T>> result(cols, vector<T>(rows));

    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            result[j][i] = A[i][j];
        }
    }

    return result;
}

template <typename T>
vector<vector<T>> multiply(const vector<vector<T>>& A, const vector<vector<T>>& B) {
    int rowsA = A.size();
    int colsA = A[0].size();
    int colsB = B[0].size();
    vector<vector<T>> result(rowsA, vector<T>(colsB, 0));

    for (int i = 0; i < rowsA; i++) {
        for (int j = 0; j < colsB; j++) {
            for (int k = 0; k < colsA; k++) {
                result[i][j] += A[i][k] * B[k][j];
            }
        }
    }

    return result;
}

template <typename T>
void QR_decomposition(const vector<vector<T>>& A, vector<vector<T>>& Q, vector<vector<T>>& R) {
    int n = A.size();
    Q = A;
    R = vector<vector<T>>(n, vector<T>(n, 0));

    for (int j = 0; j < n; j++) {
        for (int i = 0; i < j; i++) {
            T dotProduct = 0;
            for (int k = 0; k < n; k++) {
                dotProduct += Q[k][i] * A[k][j];
            }
            R[i][j] = dotProduct;
            for (int k = 0; k < n; k++) {
                Q[k][j] -= R[i][j] * Q[k][i];
            }
        }

        T norm = 0;
        for (int i = 0; i < n; i++) {
            norm += Q[i][j] * Q[i][j];
        }
        R[j][j] = sqrt(norm);
        for (int i = 0; i < n; i++) {
            Q[i][j] /= R[j][j];
        }
    }
}

template <typename T>
vector<T> method_QR(vector<vector<T>>& A, vector<T>& b) {
    int n = A.size();
    vector<vector<T>> Q, R;
    QR_decomposition(A, Q, R);

    // Решение системы Q^T * y = b
    vector<T> y(n, 0);
    for (int i = 0; i < n; i++) {
        T dotProduct = 0;
        for (int j = 0; j < n; j++) {
            dotProduct += Q[j][i] * b[j];
        }
        y[i] = dotProduct;
    }

    // Решение системы R * x = y методом обратной подстановки
    vector<T> x(n, 0);
    for (int i = n - 1; i >= 0; i--) {
        x[i] = y[i];
        for (int j = i + 1; j < n; j++) {
            x[i] -= R[i][j] * x[j];
        }
        x[i] /= R[i][i];
    }

    return x;
}