#include <iostream>
#include <vector>

#include"methods.cpp"
#include"algebra.cpp"

using namespace std;

/* Тест программы */
template <typename T>
void test_programm() {
    //cout << "Precision: DOUBLE \n \n";

    const string filename = "input_data/TEST4/D1.txt";          // Путь к файлу
    vector<vector<T>> A = importMatrix<T>(filename);            // Импорт СЛАУ из текстового файла
    int n = A.size();                                           // Размер матрицы
    vector<vector<T>> E = create_identity_matrix<T>(n);         // Единичная матрица
    T EPS = 1e-4;                                               // Погрешность
    vector<T> true_D1 = {0.997313, 2.00425, 2.98702, 4.01142};  // Правильное решение теста D1
    long int MaxIteration = 1000;                               // Максимальное количество итераций метода






    cout << "\nA = \n";
    print(A);
    printline(30);
    cout << "Method QR-decompositions: " << endl << endl;

    /* Матрица A в форме Хессенберга */
    if (false) {
        vector<vector<T>> A_hess = Hessenberg_decomposition2(A, EPS);
        A = A_hess;
        cout << "A_hess = " << endl;
        print(A_hess);
    }

    /* 1) Метод QR-разложения */
    T sigma = 0;
    MyResult4<T> result1 = Eigen_method_QR3(A, EPS, MaxIteration);
    T error_QR = testeps(result1.eigen, true_D1, 2);
    T eps_QR = test_eigen(A, result1.eigen);

    cout << "lambda = ";
    print(result1.eigen);
    cout << "Number of iterations = " << result1.iterations << endl;
    cout << "Error = " << error_QR << endl;
    cout << "det(A - lambda * E) = " << eps_QR << endl;
//    cout << "R = " << endl;
//    print(result1.R);

    if (false) {
        double len = (result1.A_iter).size();
        for (int i = 0; i < len; i++) {
            print(result1.A_iter[i]);
        }
    }

    printline(30);
    cout << "Method rewerse iteration: " << endl << endl;

    /* 2) Метод Обратных Итераций */
    //vector<T> eigen = {1, 2, 3, 4};
    MyResult4<T> result2 = reverse_iteration(A, result1.eigen, EPS, MaxIteration);
    T eps_RI = test_eigen_vec(A, result2.eigens_vec, result2.eigen);

    for (int i = 0; i < n; i++){
        cout << "vec_l" << i + 1 << " = ";
        print_vec(result2.eigens_vec[i]);
    }
    cout << "Number of iterations = " << result2.iterations << endl;
    cout << "|A * e - lambda * e| = " << eps_RI << endl;

    printline(30);
    cout << "Method rewerse iteration with Reley: " << endl << endl;

    /* 3) Метод Обратных Итераций с отношением Рэлея */
    vector<vector<T>> X0 = result2.eigens_vec;
    MyResult4<T> result3 = reverse_iterator_with_reley(A, X0, EPS, MaxIteration);
    T error_Reley = testeps(result3.eigen, true_D1, 2);
    T eps_Reley = test_eigen(A, result3.eigen);
    T eps_RIReley = test_eigen_vec(A, result3.eigens_vec, result3.eigen);

    cout << "lambda = ";
    print(result3.eigen);
    for (int i = 0; i < n; i++){
        cout << "vec_l" << i + 1 << " = ";
        print_vec(result3.eigens_vec[i]);
    }
    cout << endl << "Number of iterations = " << result3.iterations << endl;
    cout << "Error = " << error_Reley << endl;
    //cout << "det(A - lambda * E) = " << eps_Reley << endl;
    cout << "|A * e - lambda * e| = " << eps_RIReley << endl;
}

//template <typename T>
//void test(){
//    const string filename = "input_data/TEST4/D1.txt";          // Путь к файлу
//    vector<vector<T>> A = importMatrix<T>(filename);            // Импорт СЛАУ из текстового файла
//    int n = A.size();                                           // Размер матрицы
//    vector<vector<T>> E = create_identity_matrix<T>(n);         // Единичная матрица
//    T EPS = 1e-4;                                               // Погрешность
//    vector<T> true_D1 = {0.997313, 2.00425, 2.98702, 4.01142};  // Правильное решение теста D1
//    long int MaxIteration = 1000;
//
//    T sigma = 0;
//    MyResult4<T> result1 = Eigen_method_QR3(A, sigma, EPS, MaxIteration);
//    //T error_QR = testeps(result1.eigen, true_D1, 2);
//    //T eps_QR = test_eigen(A, result1.eigen);
//
//    cout << "lambda = ";
//    print(result1.eigen);
//    //cout << "Number of iterations = " << result1.iterations << endl;
//    //cout << "Error = " << error_QR << endl;
//    //cout << "det(A - lambda * E) = " << eps_QR << endl;
//
//    if (true) {
//        double len = (result1.A_iter).size();
//        for (int i = 0; i < len; i++) {
//            print(result1.A_iter[i]);
//        }
//    }
//
//}


int main() {
    test_programm<double>();
    cout << "Complete!";
    return 0;
}

