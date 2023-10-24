#include <iostream>
#include <vector>
#include "header.h"

using namespace std;

int main() {

    setlocale(LC_ALL, "Russian");
    const string filename = "input_data/SYS1/DATA10.txt";                              // Путь к файлу

    // Базовые функции
    vector<vector<double>> SLAU = importSLAU<double>(filename);    // Импорт СЛАУ из текстового файла
    vector<vector<double>> matrix = SLAU_to_matrix(SLAU);          // Получение матрицы из СЛАУ
    vector<double> vec = SLAU_to_vec(SLAU);                        // Получение вектора из СЛАУ
    vector<vector<double>> trans_matrix = transpon(matrix);        // Транспонирование матрицы
    vector<vector<double>> inverse_matrix = inverseMatrix(matrix); // Обратная матрица

    //print(matrix);
    //print(vec);


    // Решение СЛАУ методом Гаусса (прямым)
    vector<double> solve = method_Gaussa(matrix, vec);
    print(solve);



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
