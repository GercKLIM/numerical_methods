#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>


using namespace std;



/* Функция импорта матрицы из текстового файла*/
template <typename T>
vector<vector<T>> importSLAU(const string& filename) {
    vector<vector<T>> matrix;
    vector<T> vec;
    ifstream file(filename);

    if (!file.is_open()) {
        printf("Error: not open file \n" );
        exit(1);
    }

    int size;
    file >> size;


    matrix.resize(size, vector<T>(size+1));

    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size+1; ++j) {
            T value;
            if (file >> value) {
                matrix[i][j] = value;
            }
        }
    }

    file.close();
    return matrix;
};



/* Функция вывода матрицы на экран */

template <typename T>
void print(vector<vector<T>> matrix) {
    for (vector<T> row : matrix) {
        for (T value: row) {
            cout << value << ' ';
        }
        cout << '\n';
    }
    cout << endl;
}

template <typename T>
void print(vector<T> vec) {
    for (T value : vec) {
        cout << value << ' ';
    }
    cout << endl;
}

/* Функция для получения матрицы из СЛАУ */
template <typename T>
vector<vector<T>> SLAU_to_matrix(vector<vector<T>> SLAU){
    vector<vector<T>> matrix;
    matrix.resize(SLAU.size(), vector<T>(SLAU.size()));

    for (int i = 0; i < SLAU.size(); i++) {
        for (int j = 0; j < SLAU.size(); j++) {
            matrix[i][j] = SLAU[i][j];
        }
    }
    return matrix;
}


/* Функция для получения вектора из СЛАУ */
template <typename T>
vector<T> SLAU_to_vec(vector<vector<T>> SLAU){
    int s = SLAU.size();
    vector<T> vec(s);

    for (int i = 0; i < SLAU.size(); i++) {
        vec[i] = SLAU[i][s];
    }
    return vec;
}


/* Функция для транспонирования матрицы */
template <typename T>
vector<vector<T>> transpon(vector<vector<T>> matrix){
    vector<vector<T>> TMatrix = matrix;
    int s = matrix.size();
    for (int i = 0; i < s; i++) {
        for (int j = 0; j < s; j++) {
            TMatrix[j][i] = matrix[i][j];
        }
    }
    return TMatrix;
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

// Функция для обратной матрицы с проверкой на вырожденность
template <typename T>
vector<vector<T>> inverseMatrix(vector<vector<T>> A) {
    int n = A.size();
    vector<vector<T>> L(n, vector<T>(n, 0));
    vector<vector<T>> U(n, vector<T>(n, 0));

    // Выполняем LU-разложение с частичным выбором
    lu_decomposition(A, L, U);

    // Создаем векторы для решения системы
    vector<vector<T>> identity = create_identity_matrix<T>(n);
    vector<vector<T>> A_inv(n, vector<T>(n, 0));

    // Решаем системы уравнений Ly = I и Ux = y для каждой строки I
    for (int i = 0; i < n; i++) {
        vector<T> y(n, 0);
        vector<T> x(n, 0);

        for (int j = 0; j < n; j++) {
            T sum = 0;
            for (int k = 0; k < j; k++) {
                sum += L[j][k] * y[k];
            }
            y[j] = (identity[i][j] - sum) / L[j][j];
        }

        for (int j = n - 1; j >= 0; j--) {
            T sum = 0;
            for (int k = j + 1; k < n; k++) {
                sum += U[j][k] * x[k];
            }
            x[j] = (y[j] - sum) / U[j][j];
        }

        A_inv[i] = x;
    }
    return A_inv;
}

/* Функция для умножения матриц */
template <typename T>
vector<vector<T>> MatrixMultiply(vector<vector<T>> A, vector<vector<T>> B){
    int m = A.size();    // Количество строк в матрице A
    int n = A[0].size(); // Количество столбцов в матрице A
    int p = B[0].size(); // Количество столбцов в матрице B

    if (n != B.size()) {
        printf("Error: impossible multiply matrix");
        exit(1);
    }

    vector<vector<T>> result(m, vector<T>(p, 0.0));

    for (int i = 0; i < m; i++) {
        for (int j = 0; j < p; j++) {
            for (int k = 0; k < n; k++) {
                result[i][j] += A[i][k] * B[k][j];
            }
        }
    }
    return result;
}


/* Функция округления чисел в матрицах */
template <typename T>
vector<vector<T>> Matrix_round(vector<vector<T>> A, double eps){
    vector<vector<T>> roundA = A;
    int size = A.size();

    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            roundA[i][j] = (round(A[i][j]) >= 0)? round(abs(A[i][j]) * (1 / eps)) / (1 / eps): -1 * round(abs(A[i][j]) * (1 / eps)) / (1 / eps);
        }
    }
    return roundA;
}


/* Функция для вычисления 1-нормы матрицы */
template <typename T>
T norm_1(vector<vector<T>> matrix){
    int rows = matrix.size();
    int cols = matrix[0].size();
    T norm = 0;

    for (int j = 0; j < cols; j++) {
        T columnSum = 0;
        for (int i = 0; i < rows; i++) {
            columnSum += std::abs(matrix[i][j]);
        }
        if (columnSum > norm) {
            norm = columnSum;
        }
    }

    return norm;
}

/* Функция для вычисления 2-нормы матрицы */
template <typename T>
T norm_2(vector<vector<T>> matrix){
    int n = matrix.size();
    T norm = 0;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            T abs_value = std::abs(matrix[i][j]);
            norm += abs_value * abs_value;
        }
    }
    norm = sqrt(norm);
    return norm;
}

/* Функция для вычисления оо-нормы матрицы */
template <typename T>
T norm_oo(vector<vector<T>> matrix){
    int rows = matrix.size();
    int cols = matrix[0].size();
    T norm = 0;

    for (int i = 0; i < rows; i++) {
        T rowSum = 0;
        for (int j = 0; j < cols; j++) {
            rowSum += std::abs(matrix[i][j]);
        }
        if (rowSum > norm) {
            norm = rowSum;
        }
    }

    return norm;
}


/* Функция для вычисления числа обусловленности матрицы c нормой 1*/
template <typename T>
T cond_1(vector<vector<T>> matrix){
    T n_1 = norm_1(matrix);
    if (n_1 == 0) {
        printf("Error: Det(A) = 0  =>  cond_1(A) = oo");
        return numeric_limits<T>::infinity();
    }
    vector<vector<T>> inverse_matrix = inverseMatrix(matrix);
    T n_2 = norm_1(inverse_matrix);
    T cond = n_1 * n_2;
    return cond;
}

/* Функция для вычисления числа обусловленности матрицы c нормой 2*/
template <typename T>
T cond_2(vector<vector<T>> matrix){
    T n_1 = norm_2(matrix);
    if (n_1 == 0) {
        printf("Error: Det(A) = 0  =>  cond_2(A) = oo");
        return numeric_limits<T>::infinity();
    }
    vector<vector<T>> inverse_matrix = inverseMatrix(matrix);
    T n_2 = norm_2(inverse_matrix);
    T cond = n_1 * n_2;
    return cond;
}

/* Функция для вычисления числа обусловленности матрицы с нормой oo*/
template <typename T>
T cond_oo(vector<vector<T>> matrix){
    T n_1 = norm_oo(matrix);
    if (n_1 == 0) {
        printf("Error: Det(A) = 0  =>  cond_oo(A) = oo");
        return numeric_limits<T>::infinity();
    }
    vector<vector<T>> inverse_matrix = inverseMatrix(matrix);
    T n_2 = norm_oo(inverse_matrix);
    T cond = n_1 * n_2;
    return cond;
}

/* Функция для сложения векторов */
template <typename T>
vector<T> vec_sum(vector<T> vec1, vector<T> vec2) {
    vector<T> pert_vec = vec1;
    for (int i = 0; i < vec1.size(); i++) {
        pert_vec[i] += vec2[i];
    }
    return pert_vec;
}


/* Функция для скалярного умножения векторов */
template <typename T>
T dot_vec(vector<T> a, vector<T> b){
    if (a.size != b.size()) {
        cout << "Error: different lenght of vectors" << endl;
        exit(1);
    }
    T sum = 0;
    for (int i = 0; i < a.size(); i++){
        sum += a[i] * b[i];
    }
    return sum;
}

/* Функция для нормы-1 вектора */
template <typename T>
T norm_1(vector<T> vec) {
    T norm = 0;
    for (const T& value : vec) {
        norm += abs(value);
    }
    return norm;
}

/* Функция для нормы-2 вектора */
template <typename T>
T norm_2(vector<T> a){
    T sum = 0;
    for (int i = 0; i < a.size(); i++){
        sum += a[i] * a[i];
    }
    return sum;
}

/* Функция для нормы-оо вектора */
template <typename T>
T norm_oo(vector<T> vec) {
    T norm = 0;
    for (const T& value : vec) {
        T abs_value = std::abs(value);
        if (abs_value > norm) {
            norm = abs_value;
        }
    }
    return norm;
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




