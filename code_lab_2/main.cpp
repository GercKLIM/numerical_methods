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

    /* Параметры методов */
    vector<T> x0 = {0, 0, 0, 0};
    T Tau = 0.01;
    T EPS = 10e-9;
    int MaxIteration = 10000;


    /* Метод Простой Итерации */
    cout << "Method Simple Iteration:" << endl;
    vector<T> sol1 = method_SimpleIteration(matrix, vec, x0, Tau, EPS, MaxIteration);
    cout << "x = ";
    print(sol1);
    cout << endl;


    /* Метод Якоби */
    cout << "Method Yacobi:" << endl;
    vector<T> sol2 = method_Yacobi(matrix, vec, x0, EPS, MaxIteration);
    cout << "x = ";
    print(sol2);
    cout << endl;

    /* Метод Зейделя */
    cout << "Method Zeidela:" << endl;
    vector<T> sol3= method_Zeidel(matrix, vec, x0, EPS, MaxIteration);
    cout << "x = ";
    print(sol3);
    cout << endl;

    /* Метод Релаксации */
    cout << "Method Relaxation:" << endl;
    T w = 0.01;
    vector<T> sol4= method_Relax(matrix, vec, x0, w, EPS, MaxIteration);
    cout << "x = ";
    print(sol4);
    cout << endl;

    /*
    vector<vector<T>> L(matrix.size(), vector<T>(matrix.size(), 0)), D(matrix.size(), vector<T>(matrix.size(), 0)), U(matrix.size(), vector<T>(matrix.size(), 0));
    LDU_decomposotion(matrix, L, D, U);
    cout << "A = " << endl;
    print(L);
    cout << " + " << endl;
    print(D);
    cout << " + " << endl;
    print(U);
    cout << " = " << endl;
    print(matrix);
    */

    /* Функции для трехдиагональных матриц реализованы, но не протестированы */
}

int main() {
    //test_programm<double>();
    test_programm<double>();
    return 0;
}

