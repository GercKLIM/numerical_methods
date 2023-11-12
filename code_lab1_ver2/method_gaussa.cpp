#include <iostream>
#include <vector>
using namespace std;

/* Функция для решения СЛАУ прямым методом Гаусса */

template <typename T>
vector<T> method_Gaussa(vector<vector<T>> matrix, vector<T> vec){
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

/* Функция поворота матрицы вправо */
template <typename T>
vector<vector<T>> MatrixRotateRight(vector<vector<T>> A){

    vector<vector<T>> A_rotate(A.size(), vector<T>(A.size(), 0));

    for (int i = 0; i < A.size(); ++i) {
        for (int j = 0; j < A.size(); ++j) {
            A_rotate[A.size() - 1 - j][i] = A[i][j];
        }
    }

    return A_rotate;

}

/* Функция поворота матрицы влево */
template <typename T>
vector<vector<T>> MatrixRotateLeft(vector<vector<T>> A){

    vector<vector<T>> A_rotate(A.size(), vector<T>(A.size(), 0));

    for (int i = 0; i < A.size(); ++i) {
        for (int j = 0; j < A.size(); ++j) {
            A_rotate[j][A.size() - 1 - i] = A[i][j];
        }
    }

    return A_rotate;
}

// Функция для создания единичной матрицы размера n x n
template <typename T>
vector<vector<T>> create_identity_matrix(int n) {
    vector<vector<T>> identity(n, vector<T>(n, 0));
    for (int i = 0; i < n; i++) {
        identity[i][i] = 1;
    }
    return identity;
}

// Функция для обратной матрицы с проверкой на вырожденность
template <typename T>
vector<vector<T>> inverseMatrix2(vector<vector<T>> A) {
    vector<vector<T>> E = create_identity_matrix<T>(A.size());
    vector<vector<T>> E_rotate = MatrixRotateLeft(E);
    vector<T> e(A.size());
    vector<vector<T>> X(A.size(), vector<T>(A.size(), 0));


    for (int i = 0; i < A.size(); i++){
        e = E_rotate[i];
        X[i] = method_Gaussa(A, e);

    }
    vector<vector<T>> A_inv = MatrixRotateLeft(X);
    return A_inv;
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