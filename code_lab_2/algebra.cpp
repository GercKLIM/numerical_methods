#include "algebra.h"
#include "methods.h"

using namespace std;

/* *** Начальные функции для испорта/экспорта данных *** */

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

/* Функция вывода вектора на экран */
template <typename T>
void print(vector<T> vec) {
    for (T value : vec) {
        cout << value << ' ';
    }
    cout << endl;
}

/* Функция вывода обрезанного вектора на экран */
template <typename T>
void print_short(vector<T> vec, const int n){

    for (int i = 0; i < n; ++i){
        cout << vec[i] << ' ';
    }
    cout << "..." << endl;
}

/* Функция вывода разделительной линии на экран */
void printline(const int& n){
    for (int i = 0; i < n; i ++){
        cout << "-";
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


/* *** Функции математики векторов *** */

/* Функция для сложения векторов */
template <typename T>
vector<T> operator+(vector<T> vec1, vector<T> vec2){
    vector<T> pert_vec = vec1;
    for (int i = 0; i < vec1.size(); i++) {
        pert_vec[i] += vec2[i];
    }
    return pert_vec;
}


/* Функция вычитания векторов */
template <typename T>
vector<T> operator-(vector<T> a, vector<T> b){
    // Проверка на возможность умножения
    if (a.size() != b.size()) {
        cout << "Error: size a != size b in substraction vectors." << endl;
        exit(1);
    }
    // Создание результирующего вектора
    vector<T> result(a.size(), 0);

    // Умножение матрицы на вектор
    for (int i = 0; i < a.size(); ++i) {
        result[i] += a[i] - b[i];
    }
    return result;

}

/* Операция почленного умножения векторов */
template <typename T>
vector<T> operator*(vector<T> vec1, vector<T> vec2){
    if (vec1.size() != vec2.size()) {
        cout << "Error: vector1 size != vector2 size in operator*." << endl;
        exit(1);
    }
    vector<T> result(vec1.size(), 0);
    for (int i = 0; i < vec1.size(); i++){
        result[i] = vec1[i] * vec2[i];
    }
    return result;
}

/* Операция умножения вектора на число */
template <typename T>
vector<T> operator*(T c, vector<T> vec){
    vector<T> result(vec.size(), 0);
    for (int i = 0; i < vec.size(); i++){
        result[i] = vec[i] * c;
    }
    return result;
}

template <typename T>
vector<T> operator*(vector<T> vec, T c){
    vector<T> result(vec.size(), 0);
    for (int i = 0; i < vec.size(); i++){
        result[i] = vec[i] * c;
    }
    return result;
}


/* Операция почленного деления векторов */
template <typename T>
vector<T> operator/(vector<T> vec1, vector<T> vec2){
    if (vec1.size() != vec2.size()) {
        cout << "Error: vector1 size != vector2 size in operator*." << endl;
        exit(1);
    }
    vector<T> result(vec1.size(), 0);
    for (int i = 0; i < vec1.size(); i++){
        result[i] = vec1[i] / vec2[i];
    }
    return result;
}

/* Функция для скалярного умножения векторов */
template <typename T>
T dot(vector<T> vec1, vector<T> vec2){
    if (vec1.size() != vec2.size()) {
        cout << "Error: vector1 size != vector2 size in operator*." << endl;
        exit(1);
    }
    T result;
    for (int i = 0; i < vec1.size(); i++){
        result[i] += vec1[i] * vec2[i];
    }
    return result;
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
T norm_2(vector<T> vec){
    T sum = 0;
    for (int i = 0; i < vec.size(); i++){
        sum += vec[i] * vec[i];
    }
    return sum;
}

/* Функция для нормы-оо вектора */
template <typename T>
T norm_oo(vector<T> vec) {
    T norm = 0;
    for (const T& value : vec) {
        T abs_value = abs(value);
        if (abs_value > norm) {
            norm = abs_value;
        }
    }
    return norm;
}



/* *** Функции математики матриц *** */

/* Операция для умножения матрицы на число */
template <typename T>
vector<vector<T>> operator*(vector<vector<T>> A,  T scalar){
    // Создание результирующей матрицы с теми же размерами
    vector<vector<T>> result(A.size(), vector<T>(A[0].size(), 0));

    // Умножение каждого элемента матрицы на число
    for (size_t i = 0; i < A.size(); ++i) {
        for (size_t j = 0; j < A[0].size(); ++j) {
            result[i][j] = A[i][j] * scalar;
        }
    }

    return result;
}

/* Операция для умножения  числа на матрицу */
template <typename T>
vector<vector<T>> operator*(T scalar, vector<vector<T>> A){
    // Создание результирующей матрицы с теми же размерами
    vector<vector<T>> result(A.size(), vector<T>(A[0].size(), 0));

    // Умножение каждого элемента матрицы на число
    for (size_t i = 0; i < A.size(); ++i) {
        for (size_t j = 0; j < A[0].size(); ++j) {
            result[i][j] = A[i][j] * scalar;
        }
    }

    return result;
}

/* Операция поэлементного сложения матриц */
template <typename T>
vector<vector<T>> operator+(vector<vector<T>> A, vector<vector<T>> B){
    // Проверка на совпадение размеров матриц
    if (A.size() != B.size() || A[0].size() != B[0].size()) {
        cout << "Error: size A != size B in addition matrix." << endl;
        exit(1);
    }

    // Создание результирующей матрицы с теми же размерами
    vector<vector<T>> result(A.size(), vector<T>(A[0].size(), 0));

    // Поэлементное сложение
    for (size_t i = 0; i < A.size(); ++i) {
        for (size_t j = 0; j < A[0].size(); ++j) {
            result[i][j] = A[i][j] + B[i][j];
        }
    }

    return result;
}

/* Операция поэлементного вычитания матриц */
template <typename T>
vector<vector<T>> operator-(vector<vector<T>> A, vector<vector<T>> B){
    // Проверка на совпадение размеров матриц
    if (A.size() != B.size() || A[0].size() != B[0].size()) {
        cout << "Error: size A != size B in substraction matrix." << endl;
        exit(1);
    }

    // Создание результирующей матрицы с теми же размерами
    vector<vector<T>> result(A.size(), vector<T>(A[0].size(), 0));

    // Поэлементное сложение
    for (size_t i = 0; i < A.size(); ++i) {
        for (size_t j = 0; j < A[0].size(); ++j) {
            result[i][j] = A[i][j] - B[i][j];
        }
    }
    return result;
}

/* Операция умножения матрицы на вектор */
template <typename T>
vector<T> operator*(vector<vector<T>> matrix, vector<T> vec) {
    // Проверка на возможность умножения
    if (matrix[0].size() != vec.size()) {
        cout << "Error: size A != size b in multiply Matrix By Vector." << endl;
        exit(1);
    }
    // Создание результирующего вектора
    vector<T> result(matrix.size(), 0);

    // Умножение матрицы на вектор
    for (int i = 0; i < matrix.size(); ++i) {
        for (size_t j = 0; j < matrix[0].size(); ++j) {
            result[i] += matrix[i][j] * vec[j];
        }
    }
    return result;
}


/* Матричное умножение */
template <typename T>
vector<vector<T>> operator*(vector<vector<T>> A, vector<vector<T>> B){
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

// Определение оператора отрицания для матрицы
template <typename T>
vector<vector<T>> operator-(vector<vector<T>> matrix) {
    int rows = matrix.size();
    int cols = matrix[0].size();
    vector<vector<T>> result(rows, vector<T>(cols, 0));
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            result[i][j] = -matrix[i][j];
        }
    }
    return result;
}

/* Функция для поэлементного умножения матриц */
template <typename T>
vector<vector<T>> Multyply(vector<vector<T>> A, vector<vector<T>> B){
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
            result[i][j] = A[i][j] * B[i][j];
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
            columnSum += abs(matrix[i][j]);
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


/* Функция поворота матрицы вправо */
template <typename T>
vector<vector<T>> RotateRight(vector<vector<T>> A){

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
vector<vector<T>> RotateLeft(vector<vector<T>> A){

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
vector<vector<T>> inverseMatrix(vector<vector<T>> A) {
    vector<vector<T>> E = create_identity_matrix<T>(A.size());
    vector<vector<T>> E_rotate = RotateLeft(E);
    vector<T> e(A.size());
    vector<vector<T>> X(A.size(), vector<T>(A.size(), 0));


    for (int i = 0; i < A.size(); i++){
        e = E_rotate[i];
        X[i] = method_Gaussa(A, e);

    }
    vector<vector<T>> A_inv = RotateLeft(X);
    return A_inv;
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
