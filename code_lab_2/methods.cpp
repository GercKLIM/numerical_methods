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

/* Функция представления матрицы С в виде: C = C_l + C_d + D_u (Нижнетреугольной, Диагональной, Верхнетреугольной) */
template<typename T>
void LDU_decomposotion(vector<vector<T>> A, vector<vector<T>> &L, vector<vector<T>> &D, vector<vector<T>> &U){
    int n = A.size();

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (i < j) {
                U[i][j] = A[i][j];
            } else if (i > j) {
                L[i][j] = A[i][j];
            } else {
                D[i][j] = A[i][j];
            }
        }
    }
}

/* Функция исследования итерационного параметра tau для метода простых итераций (Метод Золотого сечения)*/


// Целевая функция, для которой мы ищем минимум

template<typename T>
T SimpleIterations_method_matrix_C(vector<vector<T>> A, T tau) {
    vector<vector<T>> E = create_identity_matrix<T>(A.size()); // Единичный вектор
    vector<vector<T>> C = -1.0 * (tau * A - E);                   // Матрица С
    return norm_1(C);
}

// Метод золотого сечения для поиска минимума функции на заданном интервале [a, b]
template<typename T>
T golden_section_search(vector<vector<T>> A, T a, T b, T epsilon) {

    const double golden_ratio = 1.618033988749895; // Золотое Сечение
    double x1, x2;

    // Начальные точки
    x1 = b - (b - a) / golden_ratio;
    x2 = a + (b - a) / golden_ratio;

    while (fabs(b - a) > epsilon) {
        if (SimpleIterations_method_matrix_C(A, x1) < SimpleIterations_method_matrix_C(A, x2)) {
            b = x2;
            x2 = x1;
            x1 = b - (b - a) / golden_ratio;
        } else {
            a = x1;
            x1 = x2;
            x2 = a + (b - a) / golden_ratio;
        }
    }
    return (a + b) / 2;
}

template<typename T>
vector<T> method_SimpleIteration(vector<vector<T>> A, vector<T> b, vector<T> x0, T eps, int MaxIter) {


    //T tau = golden_section_search(A, -1000.0, 1000.0, eps);
    T tau = 0.01;
    vector<vector<T>> E = create_identity_matrix<T>(A.size()); // Единичный вектор
    vector<vector<T>> C = -1.0 * (tau * A - E);                   // Матрица С
    vector<T> y = tau * b;                                        // Вектор y

    vector<T> xk = x0;
    vector<T> xk_new = xk;
    T C_norm = norm_2(C);
    cout << "Tau = " <<  tau << endl;
    cout << "norm(C) = " <<  C_norm << endl;

    for (int i = 0; i < MaxIter; ++i){

        xk_new = C * xk + y;

        // Критерий останова итерационного процесса
        vector<T> delta_stop = xk_new - xk;
        if (norm_2(delta_stop) <= ((1 - norm_2(C)) / norm_2(C)) * eps) {
            break;
        }
        xk = xk_new;
    }
    return xk;

}


/* Функция решения СЛАУ методом Якоби */
template<typename T>
vector<T> method_Yacobi(vector<vector<T>> A, vector<T> b, vector<T> x0, T eps, int MaxIter){

    vector<vector<T>> L(A.size(), vector<T>(A.size(), 0)),D(A.size(), vector<T>(A.size(), 0)), U(A.size(), vector<T>(A.size(), 0));
    LDU_decomposotion(A, L, D, U);

    vector<vector<T>> D_inv = inverseMatrix(D);

    vector<vector<T>> C = -1.0 * D_inv * (L + U);
    vector<T> y = D_inv * b;

    vector<T> xk = x0;
    vector<T> xk_new = xk;
    T C_norm = norm_2(C);
    cout << "norm(C) = " <<  C_norm << endl;

    for (int i = 0; i < MaxIter; ++i){

        xk_new = C * xk + y;

        // Критерий останова итерационного процесса
        vector<T> delta_stop = xk_new - xk;
        if (norm_2(delta_stop) <= ((1 - norm_2(C)) / norm_2(C)) * eps) {
            break;
        }
        xk = xk_new;
    }
    return xk;

}

/* Функция решения СЛАУ методом Релаксации */
template<typename T>
vector<T> method_Relax(vector<vector<T>> A, vector<T> b, vector<T> x0, T w, T eps, int MaxIter){
    int n = A.size();
    vector<T> x(n, 0);  // Начальное приближение
    cout << "W = " << w << endl;
    for (int k = 0; k < MaxIter; ++k) {
        std::vector<T> x_new(n, 0);

        for (int i = 0; i < n; ++i) {
            T sum1 = 0;
            T sum2 = 0;

            for (int j = 0; j < i; ++j) {
                sum1 += A[i][j] * x_new[j];
            }

            for (int j = i + 1; j < n; ++j) {
                sum2 += A[i][j] * x[j];
            }

            x_new[i] = (1 - w) * x[i] + (w / A[i][i]) * (b[i] - sum1 - sum2);
        }

        // Проверка на сходимость
        T max_error = 0;
        for (int i = 0; i < n; ++i) {
            T error = abs(x_new[i] - x[i]);
            if (error > max_error) {
                max_error = error;
            }
        }

        // Если достигнута необходимая точность, завершаем итерации
        if (max_error < eps) {
            return x_new;
        }

        x = x_new;
    }

    // Если не достигнута необходимая точность за максимальное число итераций
    cout << "Error: Relax don't converge :(" << endl;
    return x;
}

/* Функция решения СЛАУ методом Зейделя */
template<typename T>
vector<T> method_Zeidel(vector<vector<T>> A, vector<T> b, vector<T> x0, T eps, int MaxIter){
    vector<T> x = method_Relax(A, b, x0, 1.0, eps, MaxIter);
    return x;
}

/* Функция решения трехдиагональной СЛАУ большой размерности методом Зейделя */
template <typename T>
vector<T> method_Zeidel_diag(vector<T> A, vector<T> B, vector<T> C, vector<T> b, vector<T> x0, T eps, int maxIterations) {
    size_t n = A.size();
    vector<T> x = x0; // Начальное приближение

    for (int iter = 0; iter < maxIterations; ++iter) {
        T max_diff = 0;

        for (size_t i = 0; i < n; ++i) {
            T sum1 = 0.0;
            if (i > 0) {
                sum1 = A[i] * x[i - 1];
            }

            T sum2 = 0.0;
            if (i < n - 1) {
                sum2 = C[i] * x[i + 1];
            }

            T new_x_i = (1 - 1) * x[i] + (1 / B[i]) * (b[i] - sum1 - sum2);

            max_diff = max(max_diff, abs(new_x_i - x[i]));

            x[i] = new_x_i;
        }

        // Проверка на достижение необходимой точности
        if (max_diff < eps) {
            cout << "Converged after " << iter + 1 << " iterations\n";
            break;
        }
    }
    return x;
}

/* Функция решения трехдиагональной СЛАУ большой размерности методом Релаксации */
template <typename T>
vector<T> method_Relax_diag(vector<T> A, vector<T> B, vector<T> C, vector<T> b, vector<T> x0, T w, T eps, int MaxIter){
    size_t n = A.size();
    vector<T> x = x0; // Начальное приближение

    // Итерации метода релаксации
    for (int iter = 0; iter < MaxIter; ++iter) {
        T max_diff = 0;

        for (size_t i = 0; i < n; ++i) {
            T sum = 0;

            if (i > 0) {
                sum += A[i] * x[i - 1];
            }

            if (i < n - 1) {
                sum += C[i] * x[i + 1];
            }

            T new_x_i = (1 - w) * x[i] + (w / B[i]) * (b[i] - sum);

            max_diff = max(max_diff, abs(new_x_i - x[i]));

            x[i] = new_x_i;
        }

        // Проверка на достижение необходимой точности
        if (max_diff < eps) {
            cout << "Converged after " << iter + 1 << " iterations\n";
            break;
        }
    }
    return x;
}

/* Функция исследования сходимости при различных начальных приближениях */
