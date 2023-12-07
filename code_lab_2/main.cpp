#include <iostream>
#include <vector>

#include"methods.cpp"
#include"algebra.cpp"

using namespace std;

/* Тест программы */
template <typename T>
void test_programm() {
    //cout << "Precision: DOUBLE \n \n";

    const string filename = "input_data/TEST2/D1.txt";           // Путь к файлу

    vector<vector<T>> SLAU = importSLAU<T>(filename);            // Импорт СЛАУ из текстового файла
    vector<vector<T>> matrix = SLAU_to_matrix(SLAU);             // Получение матрицы из СЛАУ
    vector<T> vec = SLAU_to_vec(SLAU);                           // Получение вектора из СЛАУ

    vector<T> x0(vec.size(), 0);                                 // Начальное приближение
    T EPS = 1e-6;                                               // Погрешность
    long int MaxIteration = 100;                                 // Максимальное количество итераций метода


    cout << "\nA = \n";
    print(matrix);
    cout << "b = ";
    print(vec);
    cout << endl;

    int p1 = 0, p2 = 0, p3 = 0, p4 = 0; // Тип нормы
    /* 1) Метод простых итераций */
    T tau = 0.05; //golden_section_search_tau<T>(matrix, 0.0, 10.0, p1,EPS); // Поиск тау методом золотого сечения
    Result<T> result1 = method_SimpleIteration(matrix, vec, x0, tau, EPS, 0, MaxIteration);

    /* 2) Метод Якоби */
    Result<T> result2 = method_Yacobi(matrix, vec, x0, EPS, p2, MaxIteration);

    /* 3) Метод Зейделя */
    Result<T> result3 = method_Zeidel(matrix, vec, x0, EPS, p3,MaxIteration);

    /* 4) Метод Релаксации */
    T W = 0.95; //golden_section_search_W(matrix, 0.0, 5.0, p4, EPS);
    Result<T> result4 = method_Relax(matrix, vec, x0, W, EPS, p4, MaxIteration);

    vector<T> true_D1 = {5.0, -7.0, 12.0, 4.0};
    vector<T> true_D2 = {10.0, -10.0, 12.0, 4.0};
    vector<T> eps_vec_1;
    vector<T> eps_vec_2;
    vector<T> eps_vec_3;
    vector<T> eps_vec_4;

    if (filename[18] == '1') {
        eps_vec_1 = true_D1 - result1.solve;
        eps_vec_2 = true_D1 - result2.solve;
        eps_vec_3 = true_D1 - result3.solve;
        eps_vec_4 = true_D1 - result4.solve;
    } else if (filename[18] == '2') {
        eps_vec_1 = true_D2 - result1.solve;
        eps_vec_2 = true_D2 - result2.solve;
        eps_vec_3 = true_D2 - result3.solve;
        eps_vec_4 = true_D2 - result4.solve;
    }

    /* Вывод результатов */

    // Метод Простой Итерации
    printline(30);
    cout << "Method Simple Iteration:\n" << endl;
    cout << "Tau = " << tau << endl;
    cout << "x = ";
    print(result1.solve);
    cout << "Converge on iterations = " << result1.iterations << endl;
    cout << "C = " << endl;
    print(result1.C);
    cout << "norm_" << p1 << "(C) = " << norm(result1.C, p1) << endl;
    cout << "norm_" << p1 << "(b - b1) = " << norm_vector_nevazki(matrix, vec, result1.solve, p1) << endl;
    cout << endl;
    aprior_eps(result1.C, result1.y, x0, p1);
    aposter_eps(result1.C, result1.batch, p1);
    cout << "EPS = " << norm(eps_vec_1, p1) << endl;
    printline(30);

    // Метод Якоби
    cout << "Method Yacobi:" << endl;
    cout << "x = ";
    print(result2.solve);
    cout << "Converge on iterations = " << result2.iterations << endl;
    cout << "C = " << endl;
    print(result2.C);
    cout << "norm_" << p2 << "(C) = " << norm(result2.C, p2) << endl;
    cout << "norm_" << p2 << "(b - b1) = " << norm_vector_nevazki(matrix, vec, result2.solve, p2) << endl;
    cout << endl;
    aprior_eps(result2.C, result2.y, x0, p2);
    aposter_eps(result2.C, result2.batch, p1);
    cout << "EPS = " << norm(eps_vec_2, p2) << endl;
    printline(30);

    // Надо править - сделать по феншую
    // Метод Зейделя
    cout << "Method Zeidela:" << endl;
    cout << "x = ";
    print(result3.solve);
    cout << "Converge on iterations = " << result3.iterations << endl;
    cout << "C = " << endl;
    print(result3.C);
    cout << "norm_" << p3 << "(C) = " << norm(result3.C, p3) << endl;
    cout << "norm_" << p3 << "(b - b1) = " << norm_vector_nevazki(matrix, vec, result3.solve, p3) << endl;
    cout << endl;
    aprior_eps(result3.C, result3.y, x0, p3);
    aposter_eps(result3.C, result3.batch, p3);
    cout << "EPS = " << norm(eps_vec_3, p3) << endl;
    printline(30);

    // Метод Релаксации
    cout << "Method Relaxation:" << endl;
    cout << "W = " << W << endl << endl;
    cout << "x = ";
    print(result4.solve);
    cout << "Converge on iterations = " << result4.iterations << endl;
    cout << "C = " << endl;
    print(result4.C);
    cout << "norm_" << p4 << "(C) = " << norm(result1.C, p4) << endl;
    cout << "norm_" << p4 << "(b - b1) = " << norm_vector_nevazki(matrix, vec, result4.solve, p4) << endl;
    cout << endl;
    aprior_eps(result3.C, result3.y, x0, p3);
    aposter_eps(result3.C, result3.batch, p3);
    cout << "EPS = " << norm(eps_vec_3, p3) << endl;
    printline(30);

    // Функция представления матрицы С в виде: C = C_l + C_d + D_u
//    vector<vector<T>> L(matrix.size(), vector<T>(matrix.size(), 0)), D(matrix.size(), vector<T>(matrix.size(), 0)), U(matrix.size(), vector<T>(matrix.size(), 0));
//    LDU_decomposotion(matrix, L, D, U);
//    cout << "A = " << endl;
//    print(L);
//    cout << " + " << endl;
//    print(D);
//    cout << " + " << endl;
//    print(U);
//    cout << " = " << endl;
//    print(matrix);
}


// Все работает :)
/* Тестирование 3-диагональных методов */
template <typename T>
void test_3diad(){

    /* Создание 3-диагональной матрицы */

    int n = 210;       // Размерность матрицы
    vector<T> A(n, 1); // Нижняя диагональ
    vector<T> B(n, 4); // Центральная диагональ
    vector<T> C(n, 1); // Верхняя диагональ
    vector<T> b(n, 0); // Вектор правой части матрицы

    vector<T> x_diag(n, 0);  // Истинное решение 3-диагональной СЛАУ
    vector<T> x0_diag(n, 0); // Начальное приближение

    for (int i = 0; i < n; ++i) {
        b[i] = 10 - 2 * ((i+1) % 2);
        x_diag[i] = 2 - ((i+1) % 2);
    }
    b[0] = 6;
    b[n] = 9 - 3 * (n % 2);


    /* Параметры методов */
    T EPS_diag = 10e-7;             // Точность
    T MaxIteration_diag = 1e10;    // Максимальное число итераций


    cout << "Method Zeidela 3-diag:" << endl;
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
    test_3diad<double>();
    cout << "Complete!";
    return 0;
}

