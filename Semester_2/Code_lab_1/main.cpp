/*  Лабораторная работа №1 "Методы численного решения ОДУ"
 *  Разработчик программы: Климов О.Д. ФН2-51Б 2024г.
 *
 * */

#include <iostream>
#include <iostream>
#include <vector>
#include <iomanip>

#include "methods.cpp"

using namespace std;

/* Объявление задачи Коши */


// Тест 1
vector<double> ODU_1(const double& t, const vector<double>& U){
    int n = U.size();
    vector<double> dU(n, 0);

    dU[0] = 2 * U[0] + U[1] * U[1] - 1;
    dU[1] = 6 * U[0] - U[1] * U[1] + 1;

    return dU;
}

// Тест 2
vector<double> ODU_2(const double& t, const vector<double>& U){
    int n = U.size();
    vector<double> dU(n, 0);

    dU[0] = 1 - U[0] * U[0] - U[1] * U[1];
    dU[1] = 2 * U[0];

    return dU;
}

// Тест 3
vector<double> ODU_3(const double& t, const vector<double>& U){
    int n = U.size();
    vector<double> dU(n, 0);

    double sigma = 10, r = 28, b = 8. / 3;
    dU[0] = sigma * (U[1] - U[2]);
    dU[1] = U[0] * (r - U[2]) - U[1];
    dU[2] = U[0] * U[1] - b * U[2];

    return dU;
}

vector<double> ODU(const double& t, const vector<double>& U){
    return ODU_1(1, U);
}


/* Метод Эйлера неявный */
//void euler(double x0, const vector<double>& y0, double h, double xn) {
//    double x = x0;
//    vector<double> y = y0;
//    while (x < xn) {
//        cout << "x = " << x;
//        for (int i = 0; i < y.size(); ++i) {
//            cout << ", y" << i << " = " << y[i];
//        }
//        cout << endl;
//
//        vector<double> dy = ODU(x, y);
//        for (int i = 0; i < y.size(); ++i) {
//            y[i] += h * dy[i];
//        }
//        x += h;
//    }
//}

// Явный Метод Эйлера
void Method_Euler(const vector<double>& u0, const vector<double>& diapazon, double h){

    int n = u0.size();
    vector<double> u_old = u0, u_new = u0;

    // Цикл по шагу
    for (int coord = diapazon[0]; coord < diapazon[1]; coord += h){


        // Цикл по компонентам
        for (int i = 1; i < n; i++){



        }
    }

}



int main() {
    // Равномерная сетка
    double h = 0.01;            // Шаг
    vector<double> D = {-5, 5}; // Отрезок


    vector<double> u0 = {1, 2};
    cout << ODU(0, u0) << endl;

    cout << "Complete!" << endl;
    return 0;
}
