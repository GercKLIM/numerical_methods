#include <iostream>
#include <vector>
#include <cmath>


/* Функция для Q-разложения матрицы */
template <typename T>
vector<vector<T>> Q_decomposition(vector<vector<T>> A) {
    int n = A.size();
    vector<vector<T>> Q = A;

    for (int j = 0; j < n; j++) {
        // Нормализация j-го столбца
        T length = 0;
        for (int i = 0; i < n; i++) {
            length += Q[i][j] * Q[i][j];
        }
        length = std::sqrt(length);

        for (int i = 0; i < n; i++) {
            Q[i][j] /= length;
        }

        // Ортогонализация оставшихся столбцов
        for (int k = j + 1; k < n; k++) {
            T dotProduct = 0;
            for (int i = 0; i < A.size(); i++) {
                dotProduct += Q[i][j] * Q[i][k];
            }
            for (int i = 0; i < A.size(); i++) {
                Q[i][k] -= dotProduct * Q[i][j];
            }
        }
    }
    return Q;
}


/* Функция для R-разложения матрицы */
template <typename T>
vector<vector<T>> R_decomposition(vector<vector<T>> A) {
    int n = A.size();

    // Создаем копии матрицы и вектора
    vector<vector<T>> R = A;

    for (int i = 0; i < n; i++) {
        T a = R[i][i];
        if (a == 0) {
            printf("Error: Det(matrix) = 0 \n" );
            exit(1);
        }

        // Делаем текущий диагональный элемент равным 1
        for (int j = i; j < n; j++) {
            R[i][j] /= a;
        }

        // Обнуляем элементы под текущим диагональным элементом
        for (int k = i + 1; k < n; k++) {
            T factor = R[k][i];
            for (int j = i; j < n; j++) {
                R[k][j] -= factor * R[i][j];
            }
        }
    }
    return R;
}

template <typename T>
vector<T> method_QR(vector<vector<T>> A, vector<T> b) {
    vector<vector<T>> Q = Q_decomposition(A);
    vector<vector<T>> R = R_decomposition(A);

    int n = Q.size();
    vector<T> x(n, 0);
    for (int i = n - 1; i >= 0; i--) {
        T sum = 0;
        for (int j = i + 1; j < n; j++) {
            sum += R[i][j] * x[j];
        }
        x[i] = (b[i] - sum) / R[i][i];
    }

    return x;
}