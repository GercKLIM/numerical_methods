#include <iostream>
#include <vector>
#include <cmath>


using namespace std;


template <typename T>
vector<vector<T>> transpose(const vector<vector<T>>& A) {
    int rows = A.size();
    int cols = A[0].size();
    vector<vector<T>> result(cols, vector<T>(rows, 0.0));

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
/*0.863868 0.00421149 0.307373 0.399043
0.431934 -0.280064 -0.795667 -0.319235
0.259161 0.452735 0.301536 -0.798087
0 0.846509 -0.426042 0.319235*/


/* Функция QR-разложения матрицы методом вращений */
template<typename T>
void QR_decomposition(const vector<vector<T>>& matrix, vector<vector<T>>& Q, vector<vector<T>>& R, const T& eps){

    int n = matrix.size();
    R = matrix;                        // R - копия матрицы A
    Q = create_identity_matrix<T>(n);  // Q - единичная матрица
    T c, s;                            // Коэффициенты с и s

    for (int i = 0; i < n; ++i){
        int m = i;
        for (int k = i; k < n; ++k){
            if (fabs(R[k][i]) > fabs(R[m][i])){
                m = k;
            };
        }

        for (int k = 0; k < n; ++k){
            swap(R[m][k], R[i][k]);
            swap(Q[k][m], Q[k][i]);
        }

        if (fabs(R[i][i]) <= eps){
            cout << "Error in QR_decomposition" << endl;
            system("pause");
            exit(1);
        }

        for (int j = i + 1; j < n; ++j){
            c = (R[i][i]) / sqrt(R[i][i] * R[i][i] + R[j][i] * R[j][i]);
            s = (R[j][i]) / sqrt(R[i][i] * R[i][i] + R[j][i] * R[j][i]);

            for (int k = 0; k < n; ++k) {
                T aa = R[i][k];
                T ab = R[j][k];

                R[i][k] = c * aa + s * ab;
                R[j][k] = c * ab - s * aa;

                T qa = Q[k][i];
                T qb = Q[k][j];

                Q[k][i] = c * qa + s * qb;
                Q[k][j] = c * qb - s * qa;
            }
            R[j][i] = 0;
        }
    }
}

/* Функция QR-разложения матрицы методом Грамма-Шмидта */
//template <typename T>
//void QR_decomposition(const vector<vector<T>>& A, vector<vector<T>>& Q, vector<vector<T>>& R) {
//    int n = A.size();
//    Q = A;
//    R = vector<vector<T>>(n, vector<T>(n, 0));
//
//    for (int j = 0; j < n; j++)  {                    // Столбцы
//        for (int i = 0; i < j; i++) {               // Строки
//            T dotProduct = 0;
//            for (int k = 0; k < n; k++) {
//                dotProduct += Q[k][i] * A[k][j];
//            }
//            R[i][j] = dotProduct;
//            for (int k = 0; k < n; k++) {
//                Q[k][j] -= R[i][j] * Q[k][i];
//            }
//        }
//
//        T norm = 0;
//        for (int i = 0; i < n; i++) {
//            norm += Q[i][j] * Q[i][j];
//        }
//        R[j][j] = sqrt(norm);
//        for (int i = 0; i < n; i++) {
//            Q[i][j] /= R[j][j];
//        }
//    }
//}

template <typename T>
vector<T> method_QR(const vector<vector<T>>& A, const vector<T>& b, const T& eps) {
    int n = A.size();
    vector<vector<T>> Q, R;
    QR_decomposition(A, Q, R, eps);

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