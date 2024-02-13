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

template <typename T>
vector<T> method_QR_withoutQR(const vector<vector<T>>& A, const vector<T>& b, const vector<vector<T>>& Q, const vector<vector<T>>& R, const T& eps) {

    int n = A.size();
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

// Функция для обратной матрицы с проверкой на вырожденность
template <typename T>
vector<vector<T>> inverseMatrix3(const vector<vector<T>>& A, const T& eps) {
    vector<vector<T>> E = create_identity_matrix<T>(A.size());
    vector<vector<T>> E_rotate = MatrixRotateLeft(E);
    vector<T> e(A.size());
    vector<vector<T>> X(A.size(), vector<T>(A.size(), 0));


    for (int i = 0; i < A.size(); i++){
        e = E_rotate[i];
        X[i] = method_QR(A, e, eps);

    }
    vector<vector<T>> A_inv = MatrixRotateLeft(X);
    return A_inv;
}

/* Функция для оценки изменения числа обуcловленности от возмущения вектора правой части */
template <typename T>
void min_change_cond2(const vector<vector<T>>& matrix, const vector<T>& vec, const vector<T>& mod, const T& eps) {
    /* Находим минимальное значение числа обусловленности */

    // Находим относительную погрешность vec
    T delta_b_1 = norm_1(mod) / norm_1(vec);
    T delta_b_2 = norm_1(mod) / norm_1(vec);
    T delta_b_oo = norm_1(mod) / norm_1(vec);

    // Находим относительную погрешность x
    T delta_x_1 = 0;
    T delta_x_2 = 0;
    T delta_x_oo = 0;

    vector<T> solve = method_Gaussa(matrix, vec, eps);

    vector<T> mod_vec = vec;
    vector<vector<T>> all_mod_vec = generateCombinations(mod);

    int n = matrix.size();
    vector<vector<T>> Q, R;
    QR_decomposition(matrix, Q, R, eps);

    for (int k = 0; k < all_mod_vec.size(); k++) {
        // Cоздаем модифицированный вектор правой части
        mod_vec = vec_sum(mod_vec, all_mod_vec[k]);
        // Ищем максимальное изменение нормы вектора изменения решения
        vector<T> mod_solve = method_QR_withoutQR(matrix, mod_vec, Q, R, eps);

        for (int i = 0; i < mod_solve.size(); i++) {
            mod_solve[i] = abs(mod_solve[i] - solve[i]);
        }
        delta_x_1 = (delta_x_1 <= norm_1(mod_solve)) ? norm_1(mod_solve) : delta_x_1;
        delta_x_2 = (delta_x_2 <= norm_2(mod_solve)) ? norm_2(mod_solve) : delta_x_2;
        delta_x_oo = (delta_x_oo <= norm_oo(mod_solve)) ? norm_oo(mod_solve) : delta_x_oo;
    }


    delta_x_1 /= norm_1(solve);
    delta_x_2 /= norm_1(solve);
    delta_x_oo /= norm_1(solve);

    T min_cond_1 = delta_x_1 / delta_b_1;
    T min_cond_2 = delta_x_2 / delta_b_2;
    T min_cond_oo = delta_x_oo / delta_b_oo;

    /* Находим максимальное значение числа обусловленности */

    vector<vector<T>> U(matrix);
    vector<vector<T>> L(matrix);

    lu_decomposition(matrix, L, U);

    L = transpon(L);
    T max_cond_1 = cond_1(L) * cond_1(U);
    T max_cond_2 = cond_2(L) * cond_2(U);
    T max_cond_oo = cond_oo(L) * cond_oo(U);


    cout << endl;
    cout << min_cond_1 << " <= cond_1(A) <= " << max_cond_1 << endl;
    cout << min_cond_2 << " <= cond_2(A) <= " << max_cond_2 << endl;
    cout << min_cond_oo << " <= cond_oo(A) <= " << max_cond_oo << endl;
    cout << endl;
}