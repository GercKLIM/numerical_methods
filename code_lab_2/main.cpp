#include <iostream>
#include <vector>

#include"methods.cpp"
#include"algebra.cpp"

using namespace std;

/* Тест программы */
template <typename T>
void test_programm() {
    // Путь к файлу
    const string filename = "input_data/TEST2/D1.txt";

    // Базовые функции
    vector<vector<T>> SLAU = importSLAU<T>(filename);         // Импорт СЛАУ из текстового файла
    vector<vector<T>> matrix = SLAU_to_matrix(SLAU);          // Получение матрицы из СЛАУ
    vector<T> vec = SLAU_to_vec(SLAU);                        // Получение вектора из СЛАУ
    vector<vector<T>> trans_matrix = transpose(matrix);        // Транспонирование матрицы
    vector<vector<T>> inverse_matrix = inverseMatrix(matrix); // Обратная матрица


    cout << "Precision: DOUBLE \n \n";
    printf("Input matrix: \nA = \n");
    print(matrix);
    printf("Input vec: \nb = ");
    print(vec);
    cout << endl;

    // Параметры методов
    vector<T> x0(vec.size(), 0);
    T EPS = 10e-9;
    int MaxIteration = 10000;


    // Метод Простой Итерации
    cout << "Method Simple Iteration:" << endl;
    vector<T> sol1 = method_SimpleIteration(matrix, vec, x0, EPS, MaxIteration);
    cout << "x = ";
    print(sol1);
    cout << endl;

    // Метод Якоби
    cout << "Method Yacobi:" << endl;
    vector<T> sol2 = method_Yacobi(matrix, vec, x0, EPS, MaxIteration);
    cout << "x = ";
    print(sol2);
    cout << endl;

    // Метод Зейделя
    cout << "Method Zeidela:" << endl;
    vector<T> sol3 = method_Zeidel(matrix, vec, x0, EPS, MaxIteration);
    cout << "x = ";
    print(sol3);
    cout << endl;

    // Метод Релаксации
    cout << "Method Relaxation:" << endl;
    T w = 0.01;
    vector<T> sol4 = method_Relax(matrix, vec, x0, w, EPS, MaxIteration);
    cout << "x = ";
    print(sol4);
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
    vector<T> A(n, 1), B(n, 4), C(n, 1);
    vector<T> D(n, 0), x(n, 0), x0_diag(n, 0);

    T W = 10e-11;
    D[0] = 6;
    for (int i = 0; i <= n; ++i) {
        D[i] = 10 - 2 * (i / 2);
        x0_diag[i] = 2 - (i % 2);
    }

    vector<T> sol5 = method_Relax_diag(A, B, C, D, x0_diag, W, EPS, MaxIteration);
    cout << "Method Relaxation 3-diag:" << endl;
    cout << "x = ";
    print(sol5);

    vector<T> sol6 = method_Zeidel_diag(A, B, C, D, x0_diag, EPS, MaxIteration + 10e20);
    cout << "Method Zeidela 3-diag:" << endl;
    cout << "x = ";
    print(sol6);
    cout << endl;

    cout << "True x = ";
    print(x0_diag);
    cout << endl;
}


int main() {
    test_programm<double>();
    return 0;
}

