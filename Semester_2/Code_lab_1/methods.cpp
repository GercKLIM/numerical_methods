#include "methods.h"

using namespace std;

/* Переопределение вывода для vector */
template<typename T>
ostream& operator<<(ostream& os, const vector<T>& vec) {
    os << "[";
    for (int i = 0; i < vec.size(); ++i) {
        os << vec[i];
        if (i != vec.size() - 1) {
            os << ", ";
        }
    }
    os << "]";
    return os;
}



/* Метод Эйлера явный */



/* Метод Рунге-Кутты 4-го порядка точности */

/* Метод Адамса-Башформа 4-го порядка точности */

/* Метод "Прогноз-Коррекция" 4-го порядка точности */