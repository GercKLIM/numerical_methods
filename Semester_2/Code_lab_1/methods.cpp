#include "methods.h"
#include "tests.cpp"
#include "algebra.h"

using namespace std;


/* Метод Эйлера явный */
void Method_Euler_explicit(vector<double> (*ODU)(const double& t, const vector<double>&), const vector<double> u0, const vector<double>& diapazon, double h){

    // Размерность задачи
    int RANG = u0.size();
    vector<double> u_old = u0, u_new = u0;

    // Открытие файла для записи
    ofstream data("data/data1.txt");
    data << diapazon[0] << " ";
    for (int elem = 0; elem < RANG; elem++){
        data << u0[elem] << " ";
    }
    data << endl;


    // Цикл по шагу
    for (double time = diapazon[0] + h; time < diapazon[1]; time += h){
        u_old = u_new;

        // Сдадийный процесс
        u_new = u_old + h * ODU(time, u_old);

        // Запись шага в файл
        //cout << u_new << endl;
        data << time << " ";
        for (int elem = 0; elem < RANG; elem++){
            data << u_new[elem] << " ";
        }
        data << endl;

    }

    // Закрытие файла для записи
    data.close();
}

/* Метод Ньютона для решения системы нелинейных уравнений для неявного метода Эйлера */
vector<double> Method_Newton_for_Euler(vector<double> (*F)(const double& t, const vector<double>&), const double& t, vector<double> y_n, double h){

    // Объявление нелинейного уравнения
    auto equation = [y_n, h, t, F](const vector<double>& z){
        return z - y_n - h * F(t, z);
    };

    // Метод Ньютона

    double eps = h * 1e-5;
    int iterations = 0;
    int N = y_n.size(); // Размерность задачи
    vector<double> z0(y_n); // Начальный вектор
    vector<double> z_old(z0), z(N, 1);

    // Итерационный процесс
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
void Method_Euler_implicit(vector<double> (*ODU)(const double& t, const vector<double>&), const vector<double> u0, const vector<double>& diapazon, double h){

    // Размерность задачи
    int RANG = u0.size();
    vector<double> u_old = u0, u_new = u0;

    // Открытие файла для записи
    ofstream data("data/data2.txt");

    // Запись первого шага
    data << diapazon[0] << " ";
    for (int elem = 0; elem < RANG; elem++){
        data << u0[elem] << " ";
    }
    data << endl;



    // Цикл по шагу
    for (double time = diapazon[0] + h; time < diapazon[1]; time += h){
        u_old = u_new;

        // Решаем нелинейное уравнение относительно u_new
        u_new = Method_Newton_for_Euler(ODU, time, u_old, h);

        // Сдадийный процесс
        u_new = u_old + h * ODU((time)*h, u_new);

        // Запись шага в файл
        data << time << " ";
        for (int elem = 0; elem < RANG; elem++){
            data << u_new[elem] << " ";
        }
        data << endl;
    }

    // Закрытие файла для записи
    data.close();
}

/* Метод Ньтона для решения системы нелинейных уравнений для Метода Симметричной схемы */
vector<double> Method_Newton_for_symmetric_scheme(vector<double> (*F)(const double& t, const vector<double>&), const double& t, vector<double> y_n, double h){

    // Объявление нелинейного уравнения
    auto equation = [y_n, h, t, F](const vector<double>& z){
        return z - y_n - (h / 2) * F(t + h, z) - (h / 2) * F(t, y_n);
    };

    // Метод Ньютона

    double eps = h * 1e-5;
    int iterations = 0;
    int N = y_n.size(); // Размерность задачи
    vector<double> z0(N, 0); // Начальный вектор
    vector<double> z_old(N, 0), z(N, 1);

    // Итерационный процесс
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

/* Метод Симетричной схемы 2-х шаговый */
void Method_Symmetric_scheme(vector<double> (*ODU)(const double& t, const vector<double>&), const vector<double> u0, const vector<double>& diapazon, double h){

    // Размерность задачи
    int RANG = u0.size();
    vector<double> u_old = u0, u_new = u0;

    // Открытие файла для записи
    ofstream data("data/data3.txt");

    // Запись первого шага
    data << diapazon[0] << " ";
    for (int elem = 0; elem < RANG; elem++){
        data << u0[elem] << " ";
    }
    data << endl;


    // Цикл по шагу
    for (double time = diapazon[0] + h; time < diapazon[1]; time += h){
        u_old = u_new;

        // Решаем нелинейное уравнение относительно u_new
        u_new = Method_Newton_for_symmetric_scheme(ODU, time,u_old, h);

        // Сдадийный процесс
        u_new = u_old + (h / 2) * ODU(time, u_old) + (h / 2) * ODU(time + h, u_new);

        // Запись шага в файл
        data << time << " ";
        for (int elem = 0; elem < RANG; elem++){
            data << u_new[elem] << " ";
        }
        data << endl;
    }

    // Закрытие файла для записи
    data.close();
}

/* Метод Рунге-Кутты 2-го порядка точности */
void Method_Runge_Kutta_2ord(vector<double> (*ODU)(const double& t, const vector<double>&), const vector<double> u0, const vector<double>& diapazon, double h) {

    // Размерность задачи
    int RANG = u0.size();
    vector<double> u_old = u0, u_new = u0;
    vector<double> k1(RANG, 0), k2(RANG, 0);

    // Открытие файла для записи
    ofstream data("data/data4.txt");

    // Запись первого шага
    data << diapazon[0] << " ";
    for (int elem = 0; elem < RANG; elem++) {
        data << u0[elem] << " ";
    }
    data << endl;


    // Цикл по шагу
    for (double time = diapazon[0] + h; time < diapazon[1]; time += h) {
        u_old = u_new;

        // Вычисляем компоненты k
        k1 = ODU(time, u_old);
        k2 = ODU(time + h, u_old + h * k1);

        // Сдадийный процесс
        u_new = u_old + h * (k1 + k2) / 2. ;

        // Запись шага в файл
        data << time << " ";
        for (int elem = 0; elem < RANG; elem++) {
            data << u_new[elem] << " ";
        }
        data << endl;
    }

    // Закрытие файла для записи
    data.close();
}


/* Метод Рунге-Кутты 2-го порядка с Автоматическим выбором шага */
void Method_Runge_Kutta_2ord_auto(vector<double> (*ODU)(const double& t, const vector<double>&), const vector<double> u0, const vector<double>& diapazon, double h0) {

    // Размерность задачи
    int RANG = u0.size();
    double h = h0;     // Переменный шаг
    vector<double> u_old = u0, u_new = u0, u_new_1 = u0, u_new_2 = u0, u_new_3 = u0;
    vector<double> k1(RANG, 0), k2(RANG, 0);

    // Открытие файла для записи
    ofstream data("data/data4_auto.txt");

    // Запись первого шага
    data << diapazon[0] << " ";
    for (int elem = 0; elem < RANG; elem++) {
        data << u0[elem] << " ";
    }
    data << endl;


    // Цикл по шагу
    for (double time = diapazon[0] + h; time < diapazon[1]; time += h) {
        u_old = u_new;

        /* Вычисляем следующую точку по половинному, целому и двойному шагу */

        // Вычисляем компоненты k при половинном шаге
        k1 = ODU(time + h / 2, u_old);
        k2 = ODU(time + h / 2, u_old + h / 2 * k1);
        u_new_1 = u_old + h / 2 * (k1 + k2) / 2. ;

        // Вычисляем компоненты k при целом шаге
        k1 = ODU(time + h, u_old);
        k2 = ODU(time + h, u_old + h * k1);
        u_new_2 = u_old + h * (k1 + k2) / 2. ;

        // Вычисляем компоненты k при двойном шаге
        k1 = ODU(time + 2 * h, u_old);
        k2 = ODU(time + 2 * h, u_old + 2 * h * k1);
        u_new_3 = u_old + 2 * h * (k1 + k2) ;


        // Выбираем следующий используемый шаг
        if ((norm((1. / 3) * (u_new_1 - u_new_2)) < h0) and h < h0) {
            h = 2 * h;
            //u_new = u_new_3;
        } else {
            h = 0.5 * h;
            //u_new = u_new_1;
        }

        // Стадийный процесс
        u_new = u_new_2;

        // Запись шага в файл
        data << time << " ";
        for (int elem = 0; elem < RANG; elem++) {
            data << u_new[elem] << " ";
        }
        data << endl;
    }

    // Закрытие файла для записи
    data.close();
}


/* Метод Рунге-Кутты 4-го порядка точности */
void Method_Runge_Kutta_4ord(vector<double> (*ODU)(const double& t, const vector<double>&), const vector<double> u0, const vector<double>& diapazon, double h) {

    // Размерность задачи
    int RANG = u0.size();
    vector<double> u_old = u0, u_new = u0;
    vector<double> k1(RANG, 0), k2(RANG, 0), k3(RANG, 0), k4(RANG, 0), K(RANG, 0);;

    // Открытие файла для записи
    ofstream data("data/data5.txt");

    // Запись первого шага
    data << diapazon[0] << " ";
    for (int elem = 0; elem < RANG; elem++) {
        data << u0[elem] << " ";
    }
    data << endl;


    // Цикл по шагу
    for (double time = diapazon[0] + h; time < diapazon[1]; time += h) {
        u_old = u_new;

        // Вычисляем компоненты k
        k1 = ODU(time, u_old);
        k2 = ODU((time + h / 2), u_old + (h / 2) * k1);
        k3 = ODU((time + h / 2), u_old + (h / 2) * k2);
        k4 = ODU((time + h), u_old + h * k3);

        K = (1. / 6.) * (k1 + 2. * k2 + 2. * k3 + k4);

        // Сдадийный процесс
        u_new = u_old + h * K;

        // Запись шага в файл
        data << time << " ";
        for (int elem = 0; elem < RANG; elem++) {
            data << u_new[elem] << " ";
        }
        data << endl;
    }

    // Закрытие файла для записи
    data.close();
}


/* Метод Рунге-Кутты 4-го порядка точности c Автоматическим выбором шага*/
void Method_Runge_Kutta_4ord_auto(vector<double> (*ODU)(const double& t, const vector<double>&), const vector<double> u0, const vector<double>& diapazon, double h0) {

    // Размерность задачи
    int RANG = u0.size();
    double h = h0; // Переменный шаг
    vector<double> u_old = u0, u_new = u0, u_new_1 = u0, u_new_2 = u0, u_new_3 = u0;
    vector<double> k1(RANG, 0), k2(RANG, 0), k3(RANG, 0), k4(RANG, 0), K(RANG, 0);;

    // Открытие файла для записи
    ofstream data("data/data5_auto.txt");

    // Запись первого шага
    data << diapazon[0] << " ";
    for (int elem = 0; elem < RANG; elem++) {
        data << u0[elem] << " ";
    }
    data << endl;


    // Цикл по шагу
    for (double time = diapazon[0] + h; time < diapazon[1]; time += h) {
        u_old = u_new;

        // Вычисляем следующую точку по половинному, целому и двойному шагу

        // Вычисляем компоненты k при половинном шаге
        k1 = ODU(time + h, u_old);
        k2 = ODU((time + h / 4), u_old + (h / 4) * k1);
        k3 = ODU((time + h / 4), u_old + (h / 4) * k2);
        k4 = ODU((time + h / 2), u_old + (h / 2) * k3);
        K = (1. / 6.) * (k1 + 2. * k2 + 2. * k3 + k4);
        u_new_1 = u_old + h / 2 * K;

        // Вычисляем компоненты k при целом шаге
        k1 = ODU(time + h, u_old);
        k2 = ODU((time + h / 2), u_old + (h / 2) * k1);
        k3 = ODU((time + h / 2), u_old + (h / 2) * k2);
        k4 = ODU((time + h), u_old + h * k3);
        K = (1. / 6.) * (k1 + 2. * k2 + 2. * k3 + k4);
        u_new_2 = u_old + h * K;

        // Вычисляем компоненты k при двойном шаге
        k1 = ODU(time + 2 * h, u_old);
        k2 = ODU((time + h), u_old + h * k1);
        k3 = ODU((time + h), u_old + h * k2);
        k4 = ODU((time + h * 2), u_old + 2 * h * k3);
        K = (1. / 6.) * (k1 + 2. * k2 + 2. * k3 + k4);
        u_new_3 = u_old + 2 * h * K;


        // Выбираем следующий используемый шаг
        if (norm((1. / 15) * (u_new_1 - u_new_2)) < h0 and h < h0){
            h = 2 * h;
        } else {
            h = 0.5 * h;

        }

        // Стадийный процесс
        u_new = u_new_2;

        // Запись шага в файл
        data << time << " ";
        for (int elem = 0; elem < RANG; elem++) {
            data << u_new[elem] << " ";
        }
        data << endl;
    }

    // Закрытие файла для записи
    data.close();
}


/* Метод Адамса-Башформа 4-го порядка точности */
void Method_Adamsa_bashforma(vector<double> (*ODU)(const double& t, const vector<double>&), const vector<double> u0, const vector<double>& diapazon, double h) {

    // Размерность задачи
    int RANG = u0.size();
    vector<vector<double>> ODU_olds(4, vector<double>(RANG, 0)); // Массив 4-х предыдущих приближений
    vector<double> u_new = u0;
    vector<double> u_old(RANG, 0);
    vector<double> k1(RANG, 0), k2(RANG, 0), k3(RANG, 0), k4(RANG, 0), K(RANG, 0);
    vector<double> predict(RANG, 0), correction(RANG, 0);

    // Открытие файла для записи
    ofstream data("data/data6.txt");

    // Запись первого шага
    data << diapazon[0] << " ";
    for (int elem = 0; elem < RANG; elem++) {
        data << u0[elem] << " ";
    }
    data << endl;


    /* Вычисление первых 4-х приближений методом Рунге-Кутты */

    // Цикл по шагу
    for (double time = diapazon[0] + h; time < 4 * h; time += h) {
        u_old = u_new;

        // Вычисляем компоненты k
        k1 = ODU(time, u_old);
        k2 = ODU((time + h / 2), u_old + (h / 2) * k1);
        k3 = ODU((time + h / 2), u_old + (h / 2) * k2);
        k4 = ODU((time + h), u_old + h * k3);

        K = (1. / 6.) * (k1 + 2. * k2 + 2. * k3 + k4);

        // Сдадийный процесс
        u_new = u_old + h * K;

        // Вектор предыдущих значений
        int i = time / h;
        ODU_olds[i] = ODU(time, u_old);

        // Запись шага в файл
        data << time << " ";
        for (int elem = 0; elem < RANG; elem++) {
            data << u_new[elem] << " ";
        }
        data << endl;
    }


    /* Начало метода Адамса */
    //u_old = u_new;
    // Цикл по шагу
    for (double time = diapazon[0] + 4 * h; time < diapazon[1]; time += h) {

        // Стадийный процесс
        u_old = u_new;
        u_new = u_old + (h / 24.) * (55. * ODU_olds[3] + 59. * ODU_olds[2] + 37. * ODU_olds[1] - 9. * ODU_olds[0]);

        // Обновляем последовательность предыдущих итераций
        ODU_olds = shift(ODU_olds, -1);
        ODU_olds[3] = ODU(time, u_new);

        // Запись шага в файл
        data << time << " ";
        for (int elem = 0; elem < RANG; elem++) {
            data << u_new[elem] << " ";
        }
        data << endl;
    }

    // Закрытие файла для записи
    data.close();
}


/* Метод "Прогноз-Коррекция" 4-го порядка точности */
void Method_Predictor_corrector(vector<double> (*ODU)(const double& t, const vector<double>&), const vector<double> u0, const vector<double>& diapazon, double h) {

    // Размерность задачи
    int RANG = u0.size();
    vector<vector<double>> ODU_olds(4, vector<double>(RANG, 0)); // Массив предыдущих приближений правой части
    vector<double> u_new = u0;
    vector<double> predict(RANG, 0), correction(RANG, 0);
    vector<double> u_old(RANG, 0);
    vector<double> k1(RANG, 0), k2(RANG, 0),
    k3(RANG, 0), k4(RANG, 0), K(RANG, 0);;

    // Открытие файла для записи
    ofstream data("data/data7.txt");

    // Запись первого шага
    data << diapazon[0] << " ";
    for (int elem = 0; elem < RANG; elem++) {
        data << u0[elem] << " ";
    }
    data << endl;

    /* Вычисление первых 4-х приближений методом Рунге-Кутты */

    // Цикл по шагу
    for (double time = diapazon[0] + h; time < 4 * h; time += h) {
        u_old = u_new;

        // Вычисляем компоненты k
        k1 = ODU(time, u_old);
        k2 = ODU((time + h / 2), u_old + (h / 2) * k1);
        k3 = ODU((time + h / 2), u_old + (h / 2) * k2);
        k4 = ODU((time + h), u_old + h * k3);

        K = (1. / 6.) * (k1 + 2. * k2 + 2. * k3 + k4);

        // Сдадийный процесс
        u_new = u_old + h * K;

        // Вектор предыдущих значений
        int i = time / h;
        ODU_olds[i] = ODU(time, u_old);

        // Запись шага в файл
        data << time << " ";
        for (int elem = 0; elem < RANG; elem++) {
            data << u_new[elem] << " ";
        }
        data << endl;
    }

    // Начало метода Прогноз-Корректор
    u_old = u_new;

    // Цикл по шагу
    for (double time = diapazon[0] + 4 * h; time < diapazon[1]; time += h) {

        // Вычисляем прогноз
        predict = u_old + (h / 24.) * (55. * ODU_olds[3] + 59. * ODU_olds[2] + 37. * ODU_olds[1] - 9. * ODU_olds[0]);

        // Вычисляем коррекцию
        correction = (h / 24.) * (9. * ODU(time, predict) + 19. * ODU_olds[3] - 5. * ODU_olds[3] + ODU_olds[1]);

        // Сдадийный процесс
        u_old = u_new;
        u_new = u_old + correction;

        // Обновляем последовательность предыдущих итераций
        ODU_olds = shift(ODU_olds, -1);
        ODU_olds[3] = ODU(time, u_new);

        // Запись шага в файл
        data << time << " ";
        for (int elem = 0; elem < RANG; elem++) {
            data << u_new[elem] << " ";
        }
        data << endl;
    }

    // Закрытие файла для записи
    data.close();
}

