#include <iostream>
#include <vector>
#include "header.h"

using namespace std;

int main() {
    const string filename = "input_data/TEST/D1.txt";                              // Путь к файлу

    // Базовые функции
    vector<vector<double>> SLAU = importSLAU<double>(filename);    // Импорт СЛАУ из текстового файла
    vector<vector<double>> matrix = SLAU_to_matrix(SLAU);          // Получение матрицы из СЛАУ
    vector<double> vec = SLAU_to_vec(SLAU);                        // Получение вектора из СЛАУ
    vector<vector<double>> trans_matrix = transpon(matrix);        // Транспонирование матрицы
    vector<vector<double>> inverse_matrix = inverseMatrix(matrix); // Обратная матрица

    printf("Input matrix: \n");
    print(matrix);
    printf("Input vec: \n");
    print(vec);

    // Решение СЛАУ методом Гаусса (прямым)
    vector<double> solve = method_Gaussa(matrix, vec);
    printf("Solve SLAU = \n");
    print(solve);

    //vector<vector<double>> E = MatrixMultiply(matrix, inverse_matrix); //
    //vector<vector<double>> roundE = Matrix_round(E, 0.01);             //  Проверка A * A^-1 = E
    // print(roundE);                                                    //

    //double n_1 = norm_1(matrix);                                       // 1-норма матрицы

    //double n_00 = norm_oo(matrix);                                     // 00-норма матрицы

    //double cond_matrix = cond(matrix);                                 // Число обусловленности матрицы


    // vector<double> true_solve_sys10 = {1, 2, 4, 20};                  // Образцовое решение СЛАУ SYS10
    // vector<double> true_solve_sys3 = {1, 2, 7, 0.000000111915};       // Образцовое решение СЛАУ SYS3
    // double n_nev = norm_vector_nevazki(solve, true_solve_sys3);       // Норма вектора незязки
    // cout << "Norm(b - b1) = " << n_nev << endl;

    vector <double> mod = {0.01, 0.01, 0.01, 0.01};                     // число модификаций
    double min_cond = evaluate_change_cond(matrix, vec, mod);           // Оценка числа обучловленности
    cout << "Cond(A) ~ " << min_cond << endl;                                           // через изменение вектора правой части



    /* Не работает, пока :( */
    // Решение СЛАУ методом QR-разложения

    //vector<vector<double>> Q = qrDecomposition_Q(matrix);
    //vector<vector<double>> R = qrDecomposition_R;
    //vector<double> solution = method_QR(matrix, vec);
    //print(Q);
    //print(R);
    //print(solution);

    return 0;
}
