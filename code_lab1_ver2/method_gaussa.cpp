#include <iostream>
#include <vector>
using namespace std;

/* Функция для решения СЛАУ прямым методом Гаусса */
template <typename T>
vector<T> method_Gaussa(vector<vector<T>>& matrix, vector<T>& vec){
    int n = matrix.size();

    // Создаем копии матрицы и вектора
    vector<vector<T>> A(matrix);
    vector<T> b(vec);

    // Прямой ход
    for (int i = 0; i < n; i++) {
        T a = A[i][i];
        if (a == 0) {
            printf("Error: Det(matrix) = 0 \n" );
            exit(1);
        }

        // Делаем текущий диагональный элемент равным 1
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
T evaluate_change_cond(vector<vector<T>> matrix, vector<T> vec, vector<T> mod){
    vector<double> solve = method_Gaussa(matrix, vec);
    vector<double> mod_vec = vec_sum(vec, mod);
    vector<double> mod_solve = method_Gaussa(matrix, mod_vec);

    T delta = 0; // Максимальное зименение в решении
    for (int i = 0; i < solve.size(); i++) {
        delta = (delta > abs(mod_solve[i] - solve[i]))? delta : abs(mod_solve[i] - solve[i]);
    }
    delta = delta / 0.01; // Если оценка возмущения близка к 1, это говорить о хорошей устойчивости системы

    T min_cond = delta * cond_1(matrix) ;
    return min_cond;
}
