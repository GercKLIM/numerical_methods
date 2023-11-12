#include "methods.h"
#include "algebra.h"
using namespace std;


// Функция для LU-разложения с частичным выбором
template <typename T>
void lu_decomposition(vector<vector<T>>& A, vector<vector<T>>& L, vector<vector<T>>& U) {
    int n = A.size();
    for (int i = 0; i < n; i++) {
        int pivot_row = i;
        T max_val = 0;
        for (int k = i; k < n; k++) {
            if (abs(A[k][i]) > max_val) {
                max_val = abs(A[k][i]);
                pivot_row = k;
            }
        }

        if (pivot_row != i) {
            swap(A[i], A[pivot_row]);
            swap(L[i], L[pivot_row]);
            swap(U[i], U[pivot_row]);
        }

        for (int j = i; j < n; j++) {
            U[i][j] = A[i][j];
            L[i][j] = i == j ? 1 : 0;
            for (int k = 0; k < i; k++) {
                U[i][j] -= L[i][k] * U[k][j];
            }
        }

        for (int j = i + 1; j < n; j++) {
            L[j][i] = A[j][i];
            for (int k = 0; k < i; k++) {
                L[j][i] -= L[j][k] * U[k][i];
            }
            L[j][i] /= U[i][i];
        }
    }
}

/* Функция для вычисления нормы вектора невязки */
template <typename T>
T norm_vector_nevazki(vector<vector<T>> A, vector<T> b, vector<T> x, const int n) {
    int s = A.size();
    vector<T> residual(s, 0);


    for (int i = 0; i < s; i++) {
        for (int j = 0; j < s; j++) {
            residual[i] += A[i][j] * x[j];
        }
        residual[i] = b[i] - residual[i];
    }
    T residual_norm;
    if (n == 1) {
        residual_norm = norm_1(residual);
        return residual_norm;
    }
    if (n == 2) {
        residual_norm = norm_2(residual);
        return residual_norm;
    }
    if (n == 0) {
        residual_norm = norm_oo(residual);
        return residual_norm;
    } else {
        cout << "Error: U stupid" << n << endl;
        exit(1);
    }

}

/* Функция для решения СЛАУ прямым методом Гаусса */

template <typename T>
vector<T> method_Gaussa(vector<vector<T>>& matrix, vector<T>& vec){
    int n = matrix.size();

    // Создаем копии матрицы и вектора
    vector<vector<T>> A(matrix);
    vector<T> b(vec);

    // Прямой ход
    for (int i = 0; i < n; i++) {
        // Поиск максимального элемента в текущем столбце и его индекса
        int maxRow = i;
        T maxVal = abs(A[i][i]);
        for (int k = i + 1; k < n; k++) {
            if (abs(A[k][i]) > maxVal) {
                maxVal = abs(A[k][i]);
                maxRow = k;
            }
        }

        if (maxVal < numeric_limits<T>::epsilon()) {
            printf("Error: Det(matrix) = 0 \n");
            exit(1);
        }

        // Обмен строк, если необходимо
        if (maxRow != i) {
            swap(A[i], A[maxRow]);
            swap(b[i], b[maxRow]);
        }

        // Делаем текущий диагональный элемент равным 1
        T a = A[i][i];
        for (int j = i; j < n; j++) {
            A[i][j] /= a;
        }
        b[i] /= a;

        // Обнуляем элементы под текущим диагональным элементом
        for (int k = i + 1; k < n; k++) {
            T factor = A[k][i];
            for (int j = i; j < n; j++) {
                A[k][j] -= factor * A[i][j];
            }
            b[k] -= factor * b[i];
        }
    }

    // Подстановка обратно в систему
    vector<T> x(n, 0);
    for (int i = n - 1; i >= 0; i--) {
        x[i] = b[i];
        for (int j = i + 1; j < n; j++) {
            x[i] -= A[i][j] * x[j];
        }
    }

    return x;

}



/* Функция для оценки изменения числа обуcловленности от возмущения вектора правой части */
template <typename T>
void min_change_cond(vector<vector<T>> matrix, vector<T> vec, vector<T> mod) {
    /* Находим минимальное значение числа обусловленности */

    // Находим относительную погрешность
    T delta_b_1 = norm_1(mod) / norm_1(vec);
    T delta_b_2 = norm_1(mod) / norm_1(vec);
    T delta_b_oo = norm_1(mod) / norm_1(vec);

    // Находим относительную погрешность x
    T delta_x_1 = 0;
    T delta_x_2 = 0;
    T delta_x_oo = 0;

    vector<T> solve = method_Gaussa(matrix, vec);
    vector<T> mod_vec;

    for (int epo = 0; epo < 50; epo++) {

        // создаем модифицированный вектор правой части
        mod_vec = vec;

        for (int i = 0; i < mod_vec.size(); i++) {
            mod_vec[i] += mod[i] * pow(-1, rand());
        }

        // Ищем максимальное изменение нормы вектора изменения решения
        vector<T> mod_solve = method_Gaussa(matrix, mod_vec);

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

    int n = matrix.size();

    vector<vector<T>> U(matrix);
    vector<vector<T>> L(matrix);

    lu_decomposition(matrix, L, U);

    L = transpose(L);
    T max_cond_1 = cond_1(L) * cond_1(U);
    T max_cond_2 = cond_2(L) * cond_2(U);
    T max_cond_oo = cond_oo(L) * cond_oo(U);


    cout << endl;
    cout << min_cond_1 << " <= cond_1(A) <= " << max_cond_1 << endl;
    cout << min_cond_2 << " <= cond_2(A) <= " << max_cond_2 << endl;
    cout << min_cond_oo << " <= cond_oo(A) <= " << max_cond_oo << endl;
    cout << endl;
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


/* ### Функций лабы 2 ### */




/* Функция решения СЛАУ методом Простой Итерации */
template <typename T>
vector<T> My_Solve_Simple_Iter(vector<vector<T>> A, vector<T> b, T eps, T tau, T MaxIterations) {

    vector<T> bufM(b.size());                                         // Вектор-обменник
    vector<vector<T>> E = create_identity_matrix<T>(A.size());     // Единичная матрица
    vector<vector<T>> BUF = MatrixMultiply(A, tau);                   // BUF = tau * A
    vector<vector<T>>  C = Matrix_plus(E, BUF);                 // C = E - BUF
    vector<T> bufX = b;                                               // bufX = b       (Вектор-обменник)
    vector<T> y = vecMultiply(bufX, tau);                    // y = tau * bufX (Вектор итерации)

    T normN = 0;                                                      // Норма вектора невязки
    T normC1 = norm_1(A);                                             // Норма матрицы С
    T p0 = 0;                                                         // Расстояние от x0 до x1
    T q = normC1;
    vector<T> x(b.size()); // x
    int counter = 0;                                                   // Счётчик итераций

    T CmpEps = (1 - normC1) * eps / normC1;                            // Сравнительная точность

    int Kest = 10;
    p0 = normN;
    Kest = round(log((1 - q) * eps / p0) / log(q));

    for (int i = 0; i < MaxIterations; i++){

        bufX = x;                                                // bufX = x(k)
        bufM = MatrixMultiply(C, x);                             // bufM = C * x(k)
        x = vec_sum(bufM, y);                         // x(k + 1) = bufM + y

        bufM = vec_minus(x, bufX);                         // bufM = x - bufX

        normN = norm_1(bufM);

        if (normN > eps){
            break;
        }
    }

    cout << "||C|| = " << normC1 << endl;

    cout << "K(est) = " << Kest << endl;

    cout << "Number of iterations: " << counter << endl;

    cout << "Error: " << norm_1(bufM) << endl;

    return x;

};

/* Функция решения СЛАУ методом Якоби */

template <typename T>
vector<T> method_Yacobi(vector<vector<T>> A, vector<T>& b, T tolerance, T maxIterations) {
    int n = A.size();
    vector<double> x(n, 0.0);  // начальное приближение

    for (int k = 0; k < maxIterations; ++k) {
        vector<double> newX(n, 0.0);

        for (int i = 0; i < n; ++i) {
            double sum = 0.0;
            for (int j = 0; j < n; ++j) {
                if (j != i) {
                    sum += A[i][j] * x[j];
                }
            }
            newX[i] = (b[i] - sum) / A[i][i];
        }

        // Проверка на сходимость
        double error = 0.0;
        for (int i = 0; i < n; ++i) {
            error += abs(newX[i] - x[i]);
        }

        if (error < tolerance) {
            return newX;  // Возвращаем результат, если достигнута сходимость
        }

        x = newX;
    }

    // Возвращаем последнее найденное решение (может быть неточным)
    return x;
}


/* Функция решения СЛАУ методом Зейделя */
template <typename T>
vector<T> method_Zeidela(vector<vector<T>> A, vector<T> b, T tolerance, T maxIterations = 1000) {
    int n = A.size();
    vector<T> x(n, 0.0);  // начальное приближение

    for (int k = 0; k < maxIterations; ++k) {
        vector<T> newX(n, 0.0);

        for (int i = 0; i < n; ++i) {
            T sum1 = 0.0;
            for (int j = 0; j < i; ++j) {
                sum1 += A[i][j] * newX[j];
            }

            T sum2 = 0.0;
            for (int j = i + 1; j < n; ++j) {
                sum2 += A[i][j] * x[j];
            }

            newX[i] = (b[i] - sum1 - sum2) / A[i][i];
        }

        // Проверка на сходимость
        T error = 0.0;
        for (int i = 0; i < n; ++i) {
            error += abs(newX[i] - x[i]);
        }

        if (error < tolerance) {
            return newX;  // Возвращаем результат, если достигнута сходимость
        }

        x = newX;
    }

    // Возвращаем последнее найденное решение (может быть неточным)
    return x;
}



/* Функция решения СЛАУ методом Релаксации */

/* Функция решения трехдиагональной СЛАУ большой размерности методом Зейделя */

/* Функция решения трехдиагональной СЛАУ большой размерности методом Релаксации */

/* Функция представления матрицы С в виде: C = C_l + C_d + D_u */

/* Функция исследования итерационного параметра */

/* Функция исследования сходимости при различных начальных приближениях */
