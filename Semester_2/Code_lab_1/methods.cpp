#include "methods.h"
#include "tests.cpp"
#include "algebra.h"

using namespace std;


/* Метод Эйлера явный */

void Method_Euler(vector<double> (*ODU)(const vector<double>&), const vector<double> u0, const vector<double>& diapazon, double h){

    // Открытие файла для записи
    ofstream data("data/data1.txt");

    int RANG = u0.size();
    int n = (diapazon[1] - diapazon[0]) / h; // Количество разбиений отрезка
    vector<double> u_old = u0, u_new = u0;


    // Цикл по шагу
    for (int i = 0; i < n; i++){
        u_old = u_new;
        u_new = u_old + h * ODU(u_old); // Сдадийный процесс

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

vector<double> Method_Newton_for_Euler(vector<double> (*F)(const vector<double>&), vector<double> y_n, double h){

    // Объявление нелинейного уравнения
    auto equation = [y_n, h, F](const vector<double>& z){
        return z - y_n - h * F(z);
    };

    // Метод Ньютона

    double eps = h * 1e-5;
    int iterations = 0;
    int N = y_n.size(); // Размерность задачи
    vector<double> z0(N, 0); // Начальный вектор
    vector<double> z_old(N, 0), z(N, 1);

    do {
        //cout << equation(z) << endl;
        iterations += 1;
        // Вычисляем градиент функции
        vector<double> grad(N, 0);
        for (int i = 0; i < N; i++){
            vector<double> left_point(z);
            left_point[i] -= eps;
            vector<double> right_point(z);
            right_point[i] += eps;

            grad[i] = (equation(right_point)[i] - equation(left_point)[i]) / (2 * eps);
        }
        //cout << grad << endl;

        // Обновляем компоненты
        z_old = z;
        z = z - equation(z) / grad;
        //cout << sqr(z_old, z)<< endl;
        // Проверяем условие остановки
    } while ((abs(sqr(z_old, z)) > eps) and iterations < 50);
    //cout << iterations << endl;
    return z;
}



/* Метод Эйлера неявный */

void Method_Euler_implicit(vector<double> (*ODU)(const vector<double>&), const vector<double> u0, const vector<double>& diapazon, double h){

    // Открытие файла для записи
    ofstream data("data/data1.txt");

    int RANG = u0.size();
    int n = (diapazon[1] - diapazon[0]) / h; // Количество разбиений отрезка
    vector<double> u_old = u0, u_new = u0;


    // Цикл по шагу
    for (int i = 0; i < n; i++){
        u_old = u_new;
        u_new = Method_Newton_for_Euler(ODU, u_old, h);

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



/* Метод Рунге-Кутты 4-го порядка точности */

/* Метод Адамса-Башформа 4-го порядка точности */

/* Метод "Прогноз-Коррекция" 4-го порядка точности */