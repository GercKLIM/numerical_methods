#include <iostream>
#include <vector>

#include"methods.cpp"
#include"algebra.cpp"

using namespace std;

/* Тест программы */
template <typename T>
void test_programm() {
    // Путь к файлу
    const string filename = "input_data/TEST2/D2.txt";

    // Базовые функции
    vector<vector<T>> SLAU = importSLAU<T>(filename);            // Импорт СЛАУ из текстового файла
    vector<vector<T>> matrix = SLAU_to_matrix(SLAU);             // Получение матрицы из СЛАУ
    vector<T> vec = SLAU_to_vec(SLAU);                           // Получение вектора из СЛАУ
    vector<vector<T>> trans_matrix = transpose(matrix);       // Транспонирование матрицы
    vector<vector<T>> inverse_matrix = inverseMatrix(matrix); // Обратная матрица

    cout << "Precision: DOUBLE \n \n";
    printf("Input matrix: \nA = \n");
    print(matrix);
    printf("Input vec: \nb = ");
    print(vec);
    cout << endl;

    // Параметры методов
    vector<T> x0(vec.size(), 0);
    T EPS = 10e-3;
    long int MaxIteration = 100000;


    // Метод Простой Итерации
    cout << "Method Simple Iteration:" << endl;

    T tau = 0; //golden_section_search_tau<T>(matrix, -1000.0, 1000.0, EPS); // Поиск тау меотодом золотого сечения

    cout << "Tau = " << tau << endl;
    cout << "norm_1(C) = " << SimpleIterations_method_matrix_norm_C(matrix, tau) << endl;
    vector<T> sol1 = method_SimpleIteration(matrix, vec, x0, tau, EPS, MaxIteration);
    cout << "x = ";
    print(sol1);
    cout << "norm_1(b - b1) = " << norm_vector_nevazki(matrix, vec, sol1, 1) << endl;
    cout << endl;

    // Метод Якоби
    cout << "Method Yacobi:" << endl;
    vector<T> sol2 = method_Yacobi(matrix, vec, x0, EPS, MaxIteration);
    cout << "x = ";
    print(sol2);
    cout << "norm_1(b - b1) = " << norm_vector_nevazki(matrix, vec, sol2, 1) << endl;
    cout << endl;

    // Метод Зейделя
    cout << "Method Zeidela:" << endl;
    vector<T> sol3 = method_Zeidel(matrix, vec, x0, EPS, MaxIteration);
    cout << "x = ";
    print(sol3);
    cout << "norm_1(b - b1) = " << norm_vector_nevazki(matrix, vec, sol3, 1) << endl;
    cout << endl;

    // Метод Релаксации
    cout << "Method Relaxation:" << endl;
    T w = 1;
    vector<T> sol4 = method_Relax(matrix, vec, x0, w, EPS, MaxIteration);
    cout << "x = ";
    print(sol4);
    cout << "norm_1(b - b1) = " << norm_vector_nevazki(matrix, vec, sol4, 1) << endl;
    cout << endl;

    // Функция представления матрицы С в виде: C = C_l + C_d + D_u
    //vector<vector<T>> L(matrix.size(), vector<T>(matrix.size(), 0)), D(matrix.size(), vector<T>(matrix.size(), 0)), U(matrix.size(), vector<T>(matrix.size(), 0));
    //LDU_decomposotion(matrix, L, D, U);
    //cout << "A = " << endl;
    //print(L);
    //cout << " + " << endl;
    //print(D);
    //cout << " + " << endl;
    //print(U);
    //cout << " = " << endl;
    //print(matrix);


    // Функции для трехдиагональных матриц реализованы, но не протестированы

    // зададим 3-диагональную матрицу через векторы диагоналей
    int n = 210;
    vector<T> A(n, 1), B(n, 4), C(n, 1), b(n, 0);  // Диагонали и вектор правой части трехдиагональной матрицы
    vector<T> x_diag(n, 0); // true sol
    vector<T> x0_diag(n, 0); // Начальное приближение


    for (int i = 0; i < n; ++i) {
        b[i] = 10 - 2 * ((i+1) % 2);
        x_diag[i] = 2 - ((i+1) % 2);
    }
    b[0] = 6;
    b[n] = 9 - 3 * (n % 2);



    T EPS_diag = 10e-7;
    T MaxIteration_diag = 10e10;

    cout << "Method Zeidela 3-diag:" << endl;
    vector<T> sol5 = method_Zeidel_diag(A, B, C, b, x0_diag, EPS_diag, MaxIteration_diag);
    cout << "x = ";
    print_short(sol5, 10);
    cout << "norm_1(b - b1) = " << norm_vector_nevazki(A, B, C, b, sol5, 1) << endl;
    cout << endl;

    T W = 1;
    cout << "Method Relaxation 3-diag:" << endl;
    cout << "W = " << W << endl;
    vector<T> sol6 = method_Relax_diag(A, B, C, b, x0_diag, W, EPS_diag, MaxIteration_diag);
    cout << "x = ";
    print_short(sol6, 10);
    cout << "norm_1(b - b1) = " << norm_vector_nevazki(A, B, C, b, sol6, 1) << endl;
    cout << endl;

}


int main() {
    test_programm<double>();
    return 0;
}

