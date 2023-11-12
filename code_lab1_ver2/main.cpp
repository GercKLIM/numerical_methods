#include <iostream>
#include <vector>
#include "header.h"

using namespace std;

/* Тест программы */
template <typename T>
void test_programm() {
    // Путь к файлу
    const string filename = "input_data/TEST/D5.txt";

    // Базовые функции
    vector<vector<T>> SLAU = importSLAU<T>(filename);         // Импорт СЛАУ из текстового файла
    vector<vector<T>> matrix = SLAU_to_matrix(SLAU);          // Получение матрицы из СЛАУ
    vector<T> vec = SLAU_to_vec(SLAU);                        // Получение вектора из СЛАУ
    vector<vector<T>> trans_matrix = transpon(matrix);        // Транспонирование матрицы
    vector<vector<T>> inverse_matrix = inverseMatrix2(matrix); // Обратная матрица


    cout << "Precision: FLOAT \n \n";
    printf("Input matrix: \nA = \n");
    print(matrix);
    printf("Input vec: \nb = ");
    print(vec);
    cout << endl;

    cout << "Norm-1(A) = " << norm_1(matrix) << endl;                    // Норма-1 матрицы
    cout << "Norm-oo(A) = " << norm_1(matrix) << endl;                   // Норма-oo матрицы
    cout << "Cond_1(A) = " << cond_1(matrix) << endl;                    // Число обусловленности через определение
    cout << "Cond_2(A) = " << cond_2(matrix) << endl;                    // Число обусловленности через определение
    cout << "Cond_oo(A) = " << cond_oo(matrix) << endl;                  // Число обусловленности через определение
    vector <T> mod = {0.01, 0.01, 0.01, 0.01};                           // число модификаций для:
    min_change_cond(matrix, vec, mod);                                   // Оценка числа обусловленности через изменение вектора правой части
    cout << "A^-1 = " << endl;
    print(inverse_matrix);

    // Решение СЛАУ методом Гаусса (прямым)
    vector<T> solve = method_Gaussa(matrix, vec);
    printf("Gauss Solve SLAU = ");
    print(solve);

    cout << endl;
    T n_nev1 = norm_vector_nevazki(matrix, vec, solve, 1); // Норма вектора незязки
    cout << "Norm_1(b - b1) = " << n_nev1 << endl;
    n_nev1 = norm_vector_nevazki(matrix, vec, solve, 2); // Норма вектора незязки
    cout << "Norm_2(b - b1) = " << n_nev1 << endl;
    n_nev1 = norm_vector_nevazki(matrix, vec, solve, 0); // Норма вектора незязки
    cout << "Norm_oo(b - b1) = " << n_nev1 << endl;
    cout << endl;
    // Решение СЛАУ методом QR-разложения
    vector<T> sol = method_QR(matrix, vec);
    cout << "QR solve SLAU= ";
    print(sol);

    vector<vector<T>> Q, R;
    QR_decomposition(matrix, Q, R);
    cout << "Q = \n";
    print(Q);
    cout << "R = \n";
    print(R);

    cout << endl;
    T n_nev2 = norm_vector_nevazki(matrix, vec, sol, 1); // Норма вектора незязки
    cout << "Norm_1(b - b1) = " << n_nev2 << endl;
    n_nev2 = norm_vector_nevazki(matrix, vec, sol, 2); // Норма вектора незязки
    cout << "Norm_2(b - b1) = " << n_nev2 << endl;
    n_nev2 = norm_vector_nevazki(matrix, vec, sol, 0); // Норма вектора незязки
    cout << "Norm_oo(b - b1) = " << n_nev2 << endl;
    cout << endl;

    print(matrix);
    print(inverse_matrix);
    vector<vector<T>> E = MatrixMultiply(matrix, inverse_matrix);       //
    vector<vector<T>> roundE = Matrix_round(E, 0.01);                 //  Проверка A * A^-1 = E
    cout << "E = " << endl;                                                   //
    print(roundE);                                                            //
}

int main() {
    test_programm<double>();
    return 0;
}
