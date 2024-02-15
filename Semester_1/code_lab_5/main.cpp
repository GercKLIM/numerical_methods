#include <iostream>
#include <vector>
#include <iomanip>

#include"source.cpp"

using namespace std;


/* Тестирование уравнений */
void test1(){
    int KEY = 3;          // Выбор уравнения
    double EPS = 1e-06;  // Точность
    vector<double> line;

    if (KEY == 1) {
        line = {0,  1};
    }
    if (KEY == 2) {
        line = {-1, 10};
    }
    if (KEY == 3) {
        line = {0, 1};
    }

    /* 1) Локализация корней на отрезке */
    cout << endl << "1) Locale roots" << endl;
    //vector<vector<double> locale_line = {{-1, 10}};
    vector<vector<double>> locale_line = locale_roots(KEY, line[0], line[1],  13);
    cout << "Number of line: " << locale_line.size() << endl;
    print(locale_line);

    /* 2) Решение методом Бисекции */
    cout << endl << "2) Method Bisection" << endl;
    for (int i = 0; i < locale_line.size(); i++){
        double solve = Method_Bisection(KEY, locale_line[i], EPS);
        cout << "x = " << solve << endl;
        cout << "f(x) = " << f(solve, KEY) << endl;
    }

    /* 3) Решение методом Ньютона */
    cout << endl << "3) Method Newton" << endl << endl;

    for (int i = 0; i < locale_line.size(); i++){
        // Начальное приближение
        //double x0 = 0;
        //print(locale_line[i]);
        double X0 = 0; //(locale_line[i][1] + locale_line[i][0]) / 2 ;
        cout << "X0 = " << 8 << endl;
        double solve = Method_Newton(KEY, line, EPS, X0, true, true);
        cout << "x = " << solve << endl;
        cout << "f(x) = " << f(solve, KEY) << endl << endl;

        //cout << "Root №" << i << ": x = " << solve << ", f(x) = " << f(solve, KEY) << ", iterations = " << iteration << endl;
    }

    /* 4) Решение методом модифицированного Ньютона */
    cout << endl << "4) Method Newton modificate " << endl;
    for (int i = 0; i < locale_line.size(); i++){
        double X0 = 0;
        cout << "X0 = " << X0 << endl;
        double solve = Method_Newton_mod2(KEY, line, X0, EPS,  1, true, true);
        cout << "x = " << solve << endl;
        cout << /*fixed << setprecision(100) <<*/ "f(x) = " << f(solve, KEY) << endl << endl;
    }
}


/* Тестирование систем уравнений */
void test2(){
    int KEY1 = 1;          // Выбор уравнения
    int KEY2 = 2;          // Выбор уравнения
    double EPS = 1e-06;  // Точность

    /* Решение методом Ньютона */
    vector<double> x0 = {6.0, 2.0};

    vector<double> solve = Method_Newton_sys(KEY1, KEY2, x0, 10.0, 10.0, EPS, true);
    cout << "(x, y) = " ;
    print(solve);
    cout << "f1(x, y) = " << g(solve[0], solve[1], KEY1) << endl;
    cout << "f2(x, y) = " << g(solve[0], solve[1], KEY2) << endl;
}

void test3(){
    int KEY1 = 1;          // Выбор уравнения
    int KEY2 = 2;          // Выбор уравнения
    double EPS = 1e-06;    // Точности
    vector<double> solve2 = Method_Newton_sys2(KEY1, KEY2, 10.0, 10.0, EPS, true, true);
    //print(solve2);
}


int main(){
    test3();
    cout << "Complete!" << endl;
    return 0;
}
