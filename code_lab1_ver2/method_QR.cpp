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
void QR_decomposition(const vector<vector<T>>& matrix, vector<vector<T>>& Q, vector<vector<T>>& R){
    int rows = matrix.size();
    vector<vector<T>> A(matrix), G(Q);
    T accuracy = 10e-10;
    int cols;
    if (rows != 0) {
        cols = rows;

    } else {
        cout << "Error: Matrix A not quadratic" << endl;
        exit(1);
    }

    G.resize(rows);
    for (int i = 0; i < rows; i++){
        G[i].resize(cols);
        for (size_t j = 0; j < cols; j++){
            G[i][j] = 0.0;
        }
    }
    for (int i = 0; i < rows; i++){
        G[i][i] = 1.0;
    }
    for (int k = 0; k < rows; k++){
        for (int i = k + 1; i < rows; i++){
            if (abs(A[i][k]) >= accuracy){
                T c = A[k][k]/sqrt(A[k][k] * A[k][k] + A[i][k] * A[i][k]);
                T s = A[i][k]/sqrt(A[k][k] * A[k][k] + A[i][k] * A[i][k]);
                for (int j = 0; j < cols; j++){
                    T temp = G[k][j];
                    G[k][j] = c * G[k][j] + s * G[i][j];
                    G[i][j] = -s * temp + c * G[i][j];
                    if (abs(G[i][j]) < accuracy)
                        G[i][j] = 0.0;
                }
                for (int j = k; j < cols; j++){
                    T temp = A[k][j];
                    G[k][j] = c * G[k][j] + s * G[i][j];
                    A[i][j] = -s * temp + c * A[i][j];
                    if (abs(A[i][j]) < accuracy)
                        A[i][j] = 0.0;
                }
            }
        }
    }
    R = multiply(G, A);
    Q = transpose(G);

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
vector<T> method_QR(const vector<vector<T>>& A, const vector<T>& b) {
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