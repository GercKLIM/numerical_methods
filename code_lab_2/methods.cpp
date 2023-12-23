#include "methods.h"
#include "algebra.h"
using namespace std;


// Функция для LU-разложения с частичным выбором
template <typename T>
void lu_decomposition(const vector<vector<T>>& A, vector<vector<T>>& L, vector<vector<T>>& U) {
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
T norm_vector_nevazki(const vector<vector<T>>& A, const vector<T>& b, const vector<T>& x, const int& n) {
    int s = A.size();
    vector<T> residual(s, 0);


    for (int i = 0; i < s; i++) {
        for (int j = 0; j < s; j++) {
            residual[i] += A[i][j] * x[j];
        }
        residual[i] = b[i] - residual[i];
    }
    T residual_norm = norm(residual, n);
    return residual_norm;
}

/* Функция для решения СЛАУ прямым методом Гаусса */
template <typename T>
vector<T> method_Gaussa(const vector<vector<T>>& matrix, const vector<T>& vec, const T& eps){
    int n = matrix.size();

    // Создаем копии матрицы и вектора
    vector<vector<T>> A(matrix);
    vector<T> b(vec);

    // Прямой ход
    for (int i = 0; i < n; i++) {
        // Поиск максимального элемента в текущем столбце и его индекса
        int maxRow = i;
        T maxVal = fabs(A[i][i]);
        for (int k = i + 1; k < n; k++) {
            if (fabs(A[k][i]) > maxVal) {
                maxVal = fabs(A[k][i]);
                maxRow = k;
            }
        }

        if (maxVal < eps) {
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


/* Функция для оценки изменения числа обуcловленности от возмущения вектора правой части через метод Гаусса */
template <typename T>
void min_change_cond(const vector<vector<T>>& matrix,const  vector<T>& vec, const vector<T>& mod) {
    /* Находим минимальное значение числа обусловленности */

    // Находим относительную погрешность
    T delta_b_1 = norm(mod, 1) / norm(vec, 1);
    T delta_b_2 = norm(mod, 2) / norm(vec, 2);
    T delta_b_oo = norm(mod, 0) / norm(vec, 0);

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
        delta_x_1 = (delta_x_1 <= norm(mod_solve, 1)) ? norm(mod_solve, 1) : delta_x_1;
        delta_x_2 = (delta_x_2 <= norm(mod_solve, 2)) ? norm(mod_solve, 2) : delta_x_2;
        delta_x_oo = (delta_x_oo <= norm(mod_solve, 0)) ? norm(mod_solve, 0) : delta_x_oo;
    }

    delta_x_1 /= norm(solve, 1);
    delta_x_2 /= norm(solve, 2);
    delta_x_oo /= norm(solve, 0);

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


/* Функция решения СЛАУ методом QR */
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

/* Метод QR с известным разложением */
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

/* Функция для оценки изменения числа обуcловленности от возмущения вектора правой части через метод QR */
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

/* ### Функций лабы 2 ### */


/* Функция представления матрицы С в виде: C = C_l + C_d + D_u (Нижнетреугольной, Диагональной, Верхнетреугольной) */
template<typename T>
void LDU_decomposotion(const vector<vector<T>>& A, vector<vector<T>> &L, vector<vector<T>> &D, vector<vector<T>> &U){
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
template<typename T>
T SimpleIterations_method_matrix_norm_C(const vector<vector<T>>& A, const T& tau, const int& p) {
    vector<vector<T>> E = create_identity_matrix<T>(A.size()); // Единичный вектор
    vector<vector<T>> C = -(tau * A - E);                   // Матрица С
    return norm(C, p);
}

// Метод золотого сечения для поиска минимума функции на заданном интервале [a, b]
template<typename T>
T golden_section_search_tau(const vector<vector<T>>& A, T a, T b, const int& p, const T& epsilon) {

    const T golden_ratio = 1.618033988749895; // Золотое Сечение
    T x1, x2;

    // Начальные точки
    x1 = b - (b - a) / golden_ratio;
    x2 = a + (b - a) / golden_ratio;

    while (fabs(b - a) > epsilon) {
        if (SimpleIterations_method_matrix_norm_C(A, x1, p) < SimpleIterations_method_matrix_norm_C(A, x2, p)) {
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
Result<T> method_SimpleIteration(const vector<vector<T>>& A, const vector<T>& b, const vector<T>& x0, const T& tau, const T& eps, const int& p, const int& MaxIter) {

    Result<T> result;
    vector<vector<T>> E = create_identity_matrix<T>(A.size());   // Единичный вектор
    vector<vector<T>> C = -(tau * A - E);                           // Матрица С
    vector<T> y = tau * b;                                          // Вектор y
    vector<T> xk = x0;
    vector<T> xk_new = xk;

    result.C = C;
    result.y = y;

    for (int i = 0; i < MaxIter; ++i){

        xk_new = C * xk + y;

        // Критерий останова итерационного процесса
        vector<T> delta_stop = xk_new - xk;
        xk = xk_new;
        if ((norm_vector_nevazki(A, b, xk, p) <= eps) or (norm(delta_stop, p) <= (((1 - norm(C, p)) / norm(C, p)) * eps))) {
            result.solve = xk;
            result.iterations = i + 1;
            result.batch =  norm(delta_stop, p);

            return result;
        }
    }
    result.solve = xk;
    result.iterations = MaxIter;
    return result;
}


/* Функция решения СЛАУ методом Якоби */
template<typename T>
Result<T> method_Yacobi(const vector<vector<T>>& A, const vector<T>& b, const vector<T>& x0, const T& eps, const int& p, const int& MaxIter){

    Result<T> result;
    vector<vector<T>> L(A.size(), vector<T>(A.size(), 0)),D(A.size(), vector<T>(A.size(), 0)), U(A.size(), vector<T>(A.size(), 0));
    LDU_decomposotion(A, L, D, U);

    vector<vector<T>> D_inv = inverseMatrix(D);
    vector<vector<T>> C =  -D_inv * (L + U);
    vector<T> y = D_inv * b;

    result.C = C;
    result.y = y;

    vector<T> xk = x0;
    vector<T> xk_new = xk;
    T C_norm = norm(C, p);



    for (int i = 0; i < MaxIter; ++i){

        xk_new = C * xk + y;

        // Критерий останова итерационного процесса
        vector<T> delta_stop = xk_new - xk;
        xk = xk_new;
        if  ((norm_vector_nevazki(A, b, xk, 1) <= eps) or (norm(delta_stop, p) <= (((1 - norm(C, p)) / norm(C, p)) * eps))){
            result.solve = xk;
            result.iterations = i + 1;
            result.batch =  norm(delta_stop, p);
            return result;
        }
    }
    result.solve = xk;
    result.iterations = MaxIter;
    return result;

}
template <typename T>
T C_matrix_for_relax(const vector<vector<T>>& A, const T& w, const int& p){

    vector<vector<T>> L(A.size(), vector<T>(A.size(), 0)), D(A.size(), vector<T>(A.size(), 0)), U(A.size(), vector<T>(A.size(),0));
    LDU_decomposotion(A, L, D, U);
    vector<vector<T>> C = ((1 / w) * D - L);
    C = inverseMatrix(C);
    C = C * (((1 - w) / w) * D + U);
    return norm(C, p);
}

/* Функция исследования итерационного параметра W для метода Релаксации для трехдиагональной матрицы (Метод Золотого сечения)*/
template<typename T>
T golden_section_search_W(const vector<vector<T>>& A, T a, T b, const int& p, const T& eps) {

    const T golden_ratio = 1.618033988749895; // Золотое Сечение
    T x1, x2;

    // Начальные точки
    x1 = b - (b - a) / golden_ratio;
    x2 = a + (b - a) / golden_ratio;

    while (abs(b - a) > eps) {
        if (C_matrix_for_relax(A, x1, p) < C_matrix_for_relax(A, x2, p) ) {
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

/* Функция решения СЛАУ методом Релаксации через каноническую формулу */
template <typename T>
Result<T> method_Relax(const vector<vector<T>>& A, const vector<T>& b, const vector<T>& x0, const T& w, const T& eps, const int& p, const int& MaxIter) {
    Result<T> result;

    int n = A.size();
    vector<T> x = x0;  // Начальное приближение

    // Вычисление ненужной матрицы С для вывода
    vector<vector<T>> L(A.size(), vector<T>(A.size(), 0)), D(A.size(), vector<T>(A.size(), 0)), U(A.size(), vector<T>(A.size(),0));
    LDU_decomposotion(A, L, D, U);
    vector<vector<T>> C = ((1 / w) * D - L);
    C = inverseMatrix(C);
    vector<T> y = C * b;
    C = C * (((1 - w) / w) * D + U);

    result.C = C;
    result.y = y;

    for (int k = 0; k < MaxIter; ++k) {
        vector<T> x_new(n, 0);

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

        // Критерий останова итерационного процесса
        vector<T> delta_stop = x_new - x;
        x = x_new;
        if ((norm_vector_nevazki(A, b, x, p) <= eps) or
            (norm(delta_stop, p) <= (((1 - norm(C, p)) / norm(C, p)) * eps))) {
            result.solve = x;
            result.iterations = k + 1;
            result.batch =  norm(delta_stop, p);
            return result;
        }
    }
    // Если не достигнута необходимая точность за максимальное число итераций
    result.solve = x;
    result.iterations = MaxIter;
    return result;

}

/* Функция решения СЛАУ методом Релаксации через матричную формулу */
//template<typename T>
//Result<T> method_Relax(const vector<vector<T>>& A, const vector<T>& b, const vector<T>& x0, const T& w, const T& eps, const int& p, const int& MaxIter){
//    Result<T> result;
//
//    vector<vector<T>> L(A.size(), vector<T>(A.size(), 0)),D(A.size(), vector<T>(A.size(), 0)), U(A.size(), vector<T>(A.size(), 0));
//    LDU_decomposotion(A, L, D, U);
//
//    vector<vector<T>> C = ((1 / w) * D - L);
//    C = inverseMatrix(C);
//    vector<T> y = C * b;
//    C = C * (((1 - w) / w) * D + U);
//
//    vector<T> xk = x0;
//    vector<T> xk_new = xk;
//
//    result.C = C;
//
//    for (int i = 0; i < MaxIter; ++i){
//
//        xk_new = C * xk + y;
//
//        // Критерий останова итерационного процесса
//        vector<T> delta_stop = xk_new - xk;
//        xk = xk_new;
//        if ((norm_vector_nevazki(A, b, xk, p) <= eps) or (norm(delta_stop, p) <= (((1 - norm(C, p)) / norm(C, p)) * eps))) {
//            result.solve = xk;
//            result.iterations = i + 1;
//            return result;
//        }
//    }
//    result.solve = xk;
//    result.iterations = MaxIter;
//    return result;
//}


/* Функция решения СЛАУ методом Зейделя */
template <typename T>
Result<T> method_Zeidel(const vector<vector<T>>& A, const vector<T>& b, const vector<T>& x0, const T& eps, const int& p, const int& MaxIter){
    Result<T> result = method_Relax<T>(A, b, x0, 1, eps, p, MaxIter);
    return result;
}


/* Функция для вычисления нормы вектора невязки трехдиагональной СЛАУ */
template <typename T>
T norm_vector_nevazki(const vector<T>& A, const vector<T>& B, const vector<T>& C, const vector<T>& b, const vector<T>& solution, const int& n){
    // Вычисление невязки
    vector<T> residual(n);
    for (int i = 0; i < n; ++i) {
        residual[i] = b[i] - (A[i - 1] * solution[i - 1] + B[i] * solution[i] + C[i] * solution[i + 1]);
    }

    // Норма невязки
    T residual_norm = norm(residual, n);
    return residual_norm;
}


/* Функция решения трехдиагональной СЛАУ большой размерности методом Зейделя */
template <typename T>
vector<T> method_Zeidel_diag(const vector<T>& A, const vector<T>& B, const vector<T>& C, const vector<T>& b, const vector<T>& x0, const T& eps, const T& maxIterations) {
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
            cout << "Method_Zeidel_diag converged after " << iter + 1 << " iterations\n";
            return x;
        }
    }
    cout << "Method_Zeidel_diag DON'T converged after " << maxIterations << " iterations\n";
    return x;
}



/* Функция решения трехдиагональной СЛАУ большой размерности методом Релаксации */
template <typename T>
vector<T> method_Relax_diag(const vector<T>& A, const vector<T>& B, const vector<T>& C, const vector<T>& b, const vector<T>& x0, const T& w, const T& eps, const T& MaxIter){
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
            cout << "Method_Relax_diag converged after " << iter + 1 << " iterations\n";
            return x;
        }
    }
    cout << "Method_Relax_diag DON'T converged after " << MaxIter << " iterations\n";
    return x;
}


//* Функция исследования итерационного параметра W для метода Релаксации для трехдиагональной матрицы (Метод Золотого сечения)*/
template<typename T>
T golden_section_search_W(vector<T> A, vector<T> B, vector<T> C, vector<T> vec, vector<T> x, T EPS, int MaxIteration, T a, T b) {

    const T golden_ratio = 1.618033988749895; // Золотое Сечение
    T x1, x2;

    // Начальные точки
    x1 = b - (b - a) / golden_ratio;
    x2 = a + (b - a) / golden_ratio;

    while (abs(b - a) > EPS) {

        vector<T> sol1 = method_Relax_diag(A, B, C, vec, x, x1, EPS,MaxIteration);
        vector<T> sol2 = method_Relax_diag(A, B, C, vec, x, x2, EPS,MaxIteration);
        T norm_sol1 = norm_vector_nevazki(A, B, C, vec, sol1, 1);
        T norm_sol2 = norm_vector_nevazki(A, B, C, vec, sol2, 1);



        if (norm_sol1 < norm_sol2) {
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


/* Функция, которая делает диагональные элементы матрицы больше нуля */

template<typename T>
vector<vector<T>> make_plus_diaf_matrix(vector<vector<T>> matrix) {
    int rows = matrix.size();
    int cols = matrix[0].size();


    vector<vector<T>> result(rows, vector<T>(cols, 0));

    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            // Увеличиваем диагональные элементы
            if (i == j) {
                result[i][j] = matrix[i][j] + 1;
            } else {
                result[i][j] = matrix[i][j];
            }
        }
    }

    return result;
}



/* Функция априорной оценки */
template <typename T>
void aprior_eps(const vector<vector<T>>& C, const vector<T>& y, const vector<T>& x0, const int& p){
    T normC = norm(C, p);
    vector<T> delta_stop = C * x0 + y;
    T norm_delta_stop = norm(delta_stop, p);
    cout << "Aprior Eps <= " << abs((normC/ (1 - normC)) * norm_delta_stop) << endl;
}


/* Функция апостериорной оценки */
template <typename T>
void aposter_eps(const vector<vector<T>>& C, T norm_delta, const int& p){
    T normC = norm(C, p);
    cout << "Aposterior Eps <= " << abs((normC/ (1 - normC)) * norm_delta) << endl;
}


