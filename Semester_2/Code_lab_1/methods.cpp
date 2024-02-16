#include "methods.h"
#include "tests.cpp"
using namespace std;


/* Метод Эйлера явный */

void Method_Euler(const double t, const vector<double> u0, const vector<double>& diapazon, double h){

    // Открытие файла для записи
    ofstream data("data/data1.txt");

    int RANG = u0.size();
    int n = (diapazon[1] - diapazon[0]) / h; // Количество разбиений отрезка
    vector<double> u_old = u0, u_new = u0;


    // Цикл по шагу
    for (int i = 0; i < n; i++){
        u_old = u_new;
        u_new = u_old + h * ODU(t, u_old); // Сдадийный процесс

        // Запись шага в файл
        data << i * h << " ";
        for (int elem = 0; elem < RANG; elem++){
            data << u_new[elem] << " ";
        }
        data << endl;

    }

    // Закрытие файла для записи
    data.close();
}


/* Метод Эйлера явный */

void Method_Euler_implicit(const double t, const vector<double> u0, const vector<double>& diapazon, double h){

    // Открытие файла для записи
    ofstream data("data/data1.txt");

    int RANG = u0.size();
    int n = (diapazon[1] - diapazon[0]) / h; // Количество разбиений отрезка
    vector<double> u_old = u0, u_new = u0;


    // Цикл по шагу
    for (int i = 0; i < n; i++){
//        u_old = u_new;
//        u_new = h * ODU(t, u_old) + u_old; // Стадийный процесс
//
//        // Запись шага в файл
//        data << i * h << " ";
//        for (int elem = 0; elem < RANG; elem++){
//            data << u_new[elem] << " ";
//        }
//        data << endl;



    }

    // Закрытие файла для записи
    data.close();
}



/* Метод Рунге-Кутты 4-го порядка точности */

/* Метод Адамса-Башформа 4-го порядка точности */

/* Метод "Прогноз-Коррекция" 4-го порядка точности */