#include <iostream>
#include <vector>

#include"methods.cpp"
#include"algebra.cpp"

using namespace std;

/* Тест программы */
template <typename T>
void test_programm() {
    //cout << "Precision: DOUBLE \n \n";

    const string filename = "input_data/TEST4/D1_1.txt";           // Путь к файлу

    vector<vector<T>> SLAU = importSLAU<T>(filename);            // Импорт СЛАУ из текстового файла
    vector<vector<T>> A = SLAU_to_matrix(SLAU);                  // Получение матрицы из СЛАУ
    //vector<T> vec = SLAU_to_vec(SLAU);                           // Получение вектора из СЛАУ
    int n = A.size();
    T EPS = 1e-10;                                               // Погрешность
    vector<T> true_D1 = {1.0, 2.0, 3.0, 4.0};
    long int MaxIteration = 100;                                 // Максимальное количество итераций метода


    cout << "\nA = \n";
    print(A);

    vector<vector<T>> A_hess = Hessenberg_decomposition2(A, EPS);
    cout << "A_hess = " << endl;
    print(A_hess);

    printline(30);
    cout << "Method QR-decompositions: " << endl;
    // Нахождение собственных чисел методом QR-разложения
    MyResult4<T> result1 = Eigen_method_QR(A, EPS, MaxIteration);
    cout << "lambda = ";
    print(result1.eigen);

    cout << "Number of iterations = " << result1.iterations << endl;


    cout << "TEST: ";
    vector<vector<T>> E = create_identity_matrix<T>(n);
    vector<vector<T>> A_test1 = A - result1.eigen[3] * E;
    T test1 = det(A_test1);
    cout << " Det(A - lamda * E) = " << test1 << endl;


    printline(30);
    cout << "Method rewerse iteration: " << endl;

    vector<vector<T>> eigen_vec = reverse_iteration(A, result1.eigen, EPS, MaxIteration);
    for (int i = 0; i < n; i++){
        cout << "vec_l" << i + 1 << " = ";
        print_vec(eigen_vec[i]);
    }

    printline(30);
    cout << "Method rewerse iteration with Reley: " << endl;
    vector<vector<T>> X0 = eigen_vec - EPS*E;
    MyResult4<T> result2 = reverse_iterator_with_reley(A, X0, EPS, MaxIteration);

    cout << "lambda = ";
    print(result2.eigen);
    for (int i = 0; i < n; i++){
        cout << "vec_l" << i + 1 << " = ";
        print_vec(result2.eigens_vec[i]);
    }
}


int main() {
    test_programm<double>();
    cout << "Complete!";
    return 0;
}

