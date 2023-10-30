#include <iostream>
#include <vector>
#include "header.h"

using namespace std;

/* Тест программы */
template <typename T>
void test_programm(){
    const string filename = "input_data/TEST/D1.txt";              // Путь к файлу

    // Базовые функции
    vector<vector<T>> SLAU = importSLAU<T>(filename);    // Импорт СЛАУ из текстового файла
    vector<vector<T>> matrix = SLAU_to_matrix(SLAU);          // Получение матрицы из СЛАУ
    vector<T> vec = SLAU_to_vec(SLAU);                        // Получение вектора из СЛАУ
    vector<vector<T>> trans_matrix = transpon(matrix);        // Транспонирование матрицы
    vector<vector<T>> inverse_matrix = inverseMatrix(matrix); // Обратная матрица

    printf("Input matrix: \nA = \n");
    print(matrix);
    printf("Input vec: \nb = ");
    print(vec);

    // Решение СЛАУ методом Гаусса (прямым)
    vector<T> solve = method_Gaussa(matrix, vec);
    printf("Gauss Solve SLAU = ");
    print(solve);


    // Решение СЛАУ методом QR-разложения
    vector<T> sol = method_QR(matrix, vec);
    vector<vector<T>> Q = Q_decomposition(matrix);
    cout << "Q = \n";
    print(Q);
    vector<vector<T>> R = R_decomposition(matrix);
    cout << "R = \n";
    print(Q);
    cout << "QR solve SLAU= ";
    print(sol);

    //vector<vector<double>> E = MatrixMultiply(matrix, inverse_matrix); //
    //vector<vector<double>> roundE = Matrix_round(E, 0.01);             //  Проверка A * A^-1 = E
    //print(roundE);                                                     //

    //double n_1 = norm_1(matrix);                                       // 1-норма матрицы

    //double n_00 = norm_oo(matrix);                                     // 00-норма матрицы

    //double cond_matrix = cond(matrix);                                 // Число обусловленности матрицы


    //vector<double> true_solve_sys10 = {1, 2, 4, 20};                   // Образцовое решение СЛАУ SYS10
    //vector<double> true_solve_sys3 = {1, 2, 7, 0.000000111915};        // Образцовое решение СЛАУ SYS3
    //double n_nev = norm_vector_nevazki(solve, true_solve_sys3);        // Норма вектора незязки
    //cout << "Norm(b - b1) = " << n_nev << endl;

    vector <T> mod = {0.01, 0.01, 0.01, 0.01};                     // число модификаций
    T min_cond = evaluate_change_cond(matrix, vec, mod);           // Оценка числа обусловленности
    cout << "Cond_1(A) ~ " << min_cond << endl;                           // через изменение вектора правой части

    cout << "Cond_1(A) = " << cond_1(matrix) << endl;                       // Число обусловленности через определение
    cout << "Cond_oo(A) = " << cond_oo(matrix) << endl;                       // Число обусловленности через определение
}

int main() {
    test_programm<double>();
    return 0;
}
