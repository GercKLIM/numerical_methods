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
}




int main() {
    test_programm<double>();
    cout << "Complete!";
    return 0;
}

