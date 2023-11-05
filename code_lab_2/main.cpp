#include <iostream>
#include <vector>
#include "header.h"

using namespace std;

/* Тест программы */
template <typename T>

void test_programm() {
    // Путь к файлу
    const string filename = "input_data/TEST/D2.txt";

    // Базовые функции
    vector<vector<T>> SLAU = importSLAU<T>(filename);         // Импорт СЛАУ из текстового файла
    vector<vector<T>> matrix = SLAU_to_matrix(SLAU);          // Получение матрицы из СЛАУ
    vector<T> vec = SLAU_to_vec(SLAU);                        // Получение вектора из СЛАУ

}

int main() {
    test_programm<double>();
    return 0;
}

