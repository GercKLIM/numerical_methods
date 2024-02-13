#include <iostream>
#include <vector>
#include "header.h"

using namespace std;

/* Тест программы */
template <typename T>
void test_programm() {
    cout << "Precision: DOUBLE" << endl;

    const string filename = "input_data/TEST/D6.txt";                     // Путь к файлу
    const T eps = 1e-10;                                                  // Погрешность
    //const T eps = numeric_limits<T>::epsilon();

    // Базовые функции
    vector<vector<T>> SLAU = importSLAU<T>(filename);                     // Импорт СЛАУ из текстового файла
    vector<vector<T>> matrix = SLAU_to_matrix<T>(SLAU);                   // Получение матрицы из СЛАУ
    vector<T> vec = SLAU_to_vec(SLAU);                                    // Получение вектора из СЛАУ
    vector<vector<T>> trans_matrix = transpon(matrix);                    // Транспонирование матрицы
    vector<vector<T>> inverse_matrix = inverseMatrix3(matrix, eps);    // Обратная матрица
    vector<vector<T>> R_dec = R_decomposition(matrix);                    // R - Верхнетреугольная матрица
    vector<T> solve1 = method_Gaussa(matrix, vec, eps);                   // x - Вектор решения СЛАУ методом Гаусса
    T n_nev1_1 = norm_vector_nevazki(matrix, vec, solve1, 1); // Норма-1 вектора незязки
    T n_nev2_1 = norm_vector_nevazki(matrix, vec, solve1, 2); // Норма-2 вектора незязки
    T n_nevoo_1 = norm_vector_nevazki(matrix, vec, solve1, 0);// Норма-oo вектора незязки

    vector<vector<T>> Q, R;
    QR_decomposition(matrix, Q, R, eps);
    vector<T> solve2 = method_QR(matrix, vec, eps);                 // x - Вектор решения СЛАУ методом Гаусса
    T n_nev1_2 = norm_vector_nevazki(matrix, vec, solve2, 1); // Норма-1 вектора незязки
    T n_nev2_2 = norm_vector_nevazki(matrix, vec, solve2, 2); // Норма-2 вектора незязки
    T n_nevoo_2 = norm_vector_nevazki(matrix, vec, solve2, 0);// Норма-oo вектора незязки

    vector<vector<T>> E = MatrixMultiply(matrix, inverse_matrix);   // Проверка E = A^-1 * A
    vector<vector<T>> roundE = Matrix_round(E, eps);                   // Округление E до Eps

    vector<T> true_D4 = {2, 1, -0.5, 0.5};
    //vector<T> true_D5 = {1, 1000, -2.0, 3.0};
    vector<T> true_D6 = {1, 1000, -20.0, 3.0};
    vector<T> eps_vec_1 = true_D6 - solve1;
    vector<T> eps_vec_2 = true_D6 - solve2;
    T eps_gauss_1 = norm_1(eps_vec_1);
    T eps_gauss_2 = norm_2(eps_vec_1);
    T eps_gauss_oo = norm_oo(eps_vec_1);
    T eps_QR_1 = norm_1(eps_vec_1);
    T eps_QR_2 = norm_2(eps_vec_1);
    T eps_QR_oo = norm_oo(eps_vec_1);


    printline(30);
    cout << "Method Gaussa\n" << endl;
    cout << "A = \n";
    print(matrix);
    cout << "b = ";
    print(vec);
    cout << endl;
    cout << "R = \n";
    print(R_dec);
    printf("Solve: \nx = ");
    print(solve1);
    cout << endl;
    cout << "Eps_1 = " <<  eps_gauss_1 << endl;
    cout << "Eps_2 = " <<  eps_gauss_2 << endl;
    cout << "Eps_oo = " <<  eps_gauss_oo << endl;
    cout << "Norm_1(b - b1) = " << n_nev1_1 << endl;
    cout << "Norm_2(b - b1) = " << n_nev2_1 << endl;
    cout << "Norm_oo(b - b1) = " << n_nevoo_1 << endl;

    cout << endl << endl;

    printline(30);
    cout << "Method QR\n" << endl;
    cout << "A = \n";
    print(matrix);
    cout << "b = ";
    print(vec);
    cout << endl;
    cout << "Q = \n";
    print(Q);
    cout << "R = \n";
    print(R);
    printf("Solve: \nx = ");
    print(solve2);
    cout << endl;
    cout << "Eps_1 = " <<  eps_QR_1 << endl;
    cout << "Eps_2 = " <<  eps_QR_2 << endl;
    cout << "Eps_oo = " <<  eps_QR_oo << endl;
    cout << "Norm_1(b - b1) = " << n_nev1_2 << endl;
    cout << "Norm_2(b - b1) = " << n_nev1_2 << endl;
    cout << "Norm_oo(b - b1) = " << n_nev1_2 << endl;

    cout << endl << endl;

    printline(30);
    cout << "Conditionality\n" << endl;
    cout << "Norm-1(A) = " << norm_1(matrix) << endl;                    // Норма-1 матрицы
    cout << "Norm-2(A) = " << norm_2(matrix) << endl;                    // Норма-1 матрицы
    cout << "Norm-oo(A) = " << norm_1(matrix) << endl << endl;           // Норма-oo матрицы
    cout << "Cond_1(A) = " << cond_1(matrix) << endl;                    // Число обусловленности через определение
    cout << "Cond_2(A) = " << cond_2(matrix) << endl;                    // Число обусловленности через определение
    cout << "Cond_oo(A) = " << cond_oo(matrix) << endl;                  // Число обусловленности через определение
    vector <T> mod = {0.01, 0.01, 0.01, 0.01};                           // число модификаций для:
    min_change_cond2(matrix, vec, mod, eps);                              // Оценка числа обусловленности через изменение вектора правой части
    cout << "E = " << endl;
    print(roundE);
    printline(30);

    vector<vector<T>> A_prim = matrix - Q * R;
    cout << "A - Q * R = " << endl;
    print(A_prim);

}
//template <typename T>
//void test_programm2() {
//    // Путь к файлу
//    const string filename = "input_data/TEST/D5.txt";
//
//    // Базовые функции
//    vector<vector<T>> SLAU = importSLAU<T>(filename);              // Импорт СЛАУ из текстового файла
//    vector<vector<T>> matrix = SLAU_to_matrix<T>(SLAU);            // Получение матрицы из СЛАУ
//    vector<T> vec = SLAU_to_vec(SLAU);                             // Получение вектора из СЛАУ
//    vector<vector<T>> trans_matrix = transpon(matrix);             // Транспонирование матрицы
//    vector<vector<T>> inverse_matrix = inverseMatrix2(matrix);  // Обратная матрица
//
//    //cout << "Precision: FLOAT \n \n";
//    printf("Input matrix: \nA = \n");
//    print(matrix);
//    printf("Input vec: \nb = ");
//    print(vec);
//    cout << endl;
//
//    cout << "Norm-1(A) = " << norm_1(matrix) << endl;                    // Норма-1 матрицы
//    cout << "Norm-oo(A) = " << norm_1(matrix) << endl;                   // Норма-oo матрицы
//    cout << "Cond_1(A) = " << cond_1(matrix) << endl;                    // Число обусловленности через определение
//    cout << "Cond_2(A) = " << cond_2(matrix) << endl;                    // Число обусловленности через определение
//    cout << "Cond_oo(A) = " << cond_oo(matrix) << endl;                  // Число обусловленности через определение
//    vector <T> mod = {0.01, 0.01, 0.01, 0.01};                           // число модификаций для:
//    min_change_cond(matrix, vec, mod);                                   // Оценка числа обусловленности через изменение вектора правой части
//    cout << "A^-1 = " << endl;
//    print(inverse_matrix);
//
//    // Решение СЛАУ методом Гаусса (прямым)
//    vector<T> solve = method_Gaussa(matrix, vec);
//    printf("Gauss Solve SLAU = ");
//    print(solve);
//
//    cout << endl;
//    T n_nev1 = norm_vector_nevazki(matrix, vec, solve, 1); // Норма вектора незязки
//    cout << "Norm_1(b - b1) = " << n_nev1 << endl;
//    n_nev1 = norm_vector_nevazki(matrix, vec, solve, 2);   // Норма вектора незязки
//    cout << "Norm_2(b - b1) = " << n_nev1 << endl;
//    n_nev1 = norm_vector_nevazki(matrix, vec, solve, 0);   // Норма вектора незязки
//    cout << "Norm_oo(b - b1) = " << n_nev1 << endl;
//    cout << endl;
//    // Решение СЛАУ методом QR-разложения
//    vector<T> sol = method_QR(matrix, vec);
//    cout << "QR solve SLAU= ";
//    print(sol);
//
//    vector<vector<T>> Q, R;
//    QR_decomposition(matrix, Q, R);
//    cout << "Q = \n";
//    print(Q);
//    cout << "R = \n";
//    print(R);
//
//    cout << endl;
//    T n_nev2 = norm_vector_nevazki(matrix, vec, sol, 1); // Норма вектора незязки
//    cout << "Norm_1(b - b1) = " << n_nev2 << endl;
//    n_nev2 = norm_vector_nevazki(matrix, vec, sol, 2); // Норма вектора незязки
//    cout << "Norm_2(b - b1) = " << n_nev2 << endl;
//    n_nev2 = norm_vector_nevazki(matrix, vec, sol, 0); // Норма вектора незязки
//    cout << "Norm_oo(b - b1) = " << n_nev2 << endl;
//    cout << endl;
//
//    print(matrix);
//    print(inverse_matrix);
//    vector<vector<T>> E = MatrixMultiply(matrix, inverse_matrix);       //
//    vector<vector<T>> roundE = Matrix_round(E, 1e-10);                //  Проверка A * A^-1 = E
//    cout << "E = " << endl;                                                   //
//    print(roundE);                                                            //
//}

/* мой тест программы */
//template <typename T>
//void test(){
//    const string filename = "input_data/TEST/D5.txt";
//    vector<vector<T>> SLAU = importSLAU<T>(filename);
//    vector<vector<T>> matrix = SLAU_to_matrix<T>(SLAU);
//    vector<T> vec = SLAU_to_vec(SLAU);
//    vector<vector<T>> trans_matrix = transpon(matrix);
//    vector<vector<T>> inverse_matrix = inverseMatrix2(matrix);
//    cout << "Cond_1(A) = " << cond_1(matrix) << endl;                    // Число обусловленности через определение
//    cout << "Cond_2(A) = " << cond_2(matrix) << endl;                    // Число обусловленности через определение
//    cout << "Cond_oo(A) = " << cond_oo(matrix) << endl;                  // Число обусловленности через определение
//    vector <T> mod = {0.01, 0.01, 0.01, 0.01};
//    min_change_cond(matrix, vec, mod);
//
//}

int main() {
    test_programm<float>();
    //test<double>();
    return 0;
}
