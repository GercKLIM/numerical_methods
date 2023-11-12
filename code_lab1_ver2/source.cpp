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
            y[j] = (i == j) ?  1 - sum : -sum;  // Изменение в этой строке
        }

        for (int j = n - 1; j >= 0; j--) {
            T sum = 0;
            for (int k = j + 1; k < n; k++) {
                sum += U[j][k] * x[k];
            }
            x[j] = (y[j] - sum) / U[j][j];
        }

        // Изменение в этой строке
        for (int j = 0; j < n; j++) {
            A_inv[j][i] = (x[j] == -0) ? -x[j] : x[j];
        }
    }

    /*vector<vector<T>> A_inv_rotated = A_inv;

    for (int i = 0; i < A_inv.size(); ++i) {
        for (int j = 0; j < A_inv.size(); ++j) {
            A_inv_rotated[A_inv.size() - 1 - j][i] = A_inv[i][j];
        }
    }

    return A_inv_rotated;
    */
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


