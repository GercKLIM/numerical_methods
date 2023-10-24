#include <iostream>
#include <fstream>
#include <vector>
using namespace std;


/* Функция импорта матрицы из текстового файла*/
template <typename T>
vector<vector<T>> importSLAU(const string& filename) {
    vector<vector<T>> matrix;
    vector<T> vec;
    ifstream file(filename);

    if (!file.is_open()) {

        cout << "error with open file" << filename << endl;
        return matrix;
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
