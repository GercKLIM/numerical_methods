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


/* Функция для вычисления обратной матрицы методом Гаусса-Жордана */
template <typename T>
vector<vector<T>> inverseMatrix(vector<vector<T>> matrix) {
    int n = matrix.size();

    // Создаем расширенную матрицу, объединяя исходную матрицу с единичной матрицей
    vector<vector<T>> augmentedMatrix(n, vector<T>(2 * n, 0));
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            augmentedMatrix[i][j] = matrix[i][j];
        }
        augmentedMatrix[i][i + n] = 1;
    }

    // Приведение к диагональному виду
    for (int i = 0; i < n; i++) {
        double pivot = augmentedMatrix[i][i];
        if (pivot == 0) {
            cout << "Определитель матрицы ноль" << std::endl;
            exit(1);
        }
        for (int j = 0; j < 2 * n; j++) {
            augmentedMatrix[i][j] /= pivot;
        }
        for (int k = 0; k < n; k++) {
            if (k != i) {
                double factor = augmentedMatrix[k][i];
                for (int j = 0; j < 2 * n; j++) {
                    augmentedMatrix[k][j] -= factor * augmentedMatrix[i][j];
                }
            }
        }
    }

    // Извлечение обратной матрицы из расширенной
    vector<vector<double>> inverse(n, vector<double>(n, 0));
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            inverse[i][j] = augmentedMatrix[i][j + n];
        }
    }

    return inverse;
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


/* Функция для вычисления числа обусловленности матрицы */
template <typename T>
T cond(vector<vector<T>> matrix){
    T n_1 = norm_1(matrix);
    if (n_1 == 0) {
        printf("Error: Det(A) = 0  =>  cond(A) = oo");
        return numeric_limits<T>::infinity();
    }
    vector<vector<T>> inverse_matrix = inverseMatrix(matrix);
    T n_2 = norm_1(inverse_matrix);
    T cond = n_1 * n_2;
    return cond;
}


/* Функция для вычисления нормы вектора невязки */
template <typename T>
T norm_vector_nevazki(vector<T> true_solve, vector<T> numerical_solve){
    vector<T> vec_nev = true_solve;
    for (int i = 0; i < true_solve.size(); ++i) {
        vec_nev[i] = abs(true_solve[i] - numerical_solve[i]);
    }
    T n_vn = 0;
    for (int i = 0; i < vec_nev.size(); i++) {
        n_vn += vec_nev[i] * vec_nev[i];
    }
    return n_vn;
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
