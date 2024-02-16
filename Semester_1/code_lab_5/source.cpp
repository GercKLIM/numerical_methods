#include "header.h"

using namespace std;

/* Функция вывода матрицы на экран */
template <typename T>
void print(const vector<vector<T>>& matrix) {
    for (vector<T> row : matrix) {
        for (T value: row) {
            cout << value << ' ';
        }
        cout << '\n';
    }
    cout << endl;
}

/* Функция вывода вектора на экран */
template <typename T>
void print(const vector<T>& vec) {
    for (T value : vec) {
        cout << value << ' ';
    }
    cout << endl;
}

/* Функция численного вычисления производной функции в точке */
double df(double x, int key, double eps){
    double x_eps = x + eps;
    return (f(x_eps, key) - f(x, key)) / eps;
}

/* Функция локализации корней */
vector<vector<double>> locale_roots(int key, double l, double r, int div) {
    vector<vector<double>> line;
    double a;
    double b = l;
    double h = (r - l) / div;
    int iter = 0;
    for (int i = 0; i < div + 1; i++) {
        a = b;
        b = a + h;
        if (f(a, key) * f(b, key) < 0) {
            iter += 1;
            line.push_back({a, b});
        }
    }
    if (line.size() == 0) {
        line = {{ r, l }};
    }
    //cout << iter << endl;
    return line;
}


/* Функция решения уравнения методом Бисекции */
double Method_Bisection(int key, vector<double> line, double eps) {
    double b = line[0];
    double a = line[1];
    double mid = (b + a) / 2;
    int iteration = 1;
    //eps = eps / 10;

    // Проверка, что корень в начале отрезка


    double f1 = f(a, key);
    double f2 = f(b, key);
    double fmid;

    // Проверка кратности корня
    if (f1 * f2 > 0) {
        cout << "Error line in bisection" << endl;
        exit(1);
    }

    // Алгоритм бисекции

    double x = mid, x_old = 0.0;

    while (  /*(abs((b - a) / 2)) >= eps*/  abs(x - x_old) > eps) {
        fmid = f(mid, key);

        // Проверка границ отрезка на корень
        if (abs(f(a, key)) < eps) {
            x = a;
            break;
        }
        // Проверка, что корень в конце отрезка
        if (abs(f(b, key)) < eps) {
            x = b;
            break;
        }

        if (fmid * f2 < 0) {
            a = mid;
            f1 = fmid;
        } else {
            b = mid;
            f2 = fmid;
        }
        x_old = mid;
        mid = (b + a) / 2;
        x = mid;
        iteration += 1;
    }


    cout << "Iterations: " << iteration << endl;
    return x;
}


/* Функция для начального приближения решения методом Хорд */
double Method_Hord(int key, double xprev, double eps) {
    double a = xprev - eps;
    double b = xprev + eps;
    double fa = f(a, key);
    double fb = f(b, key);
    return (fa * b - fb * a) / (fa - fb);
}


/* Функция решения уравнения методом Ньютона */

double Method_Newton(int key, vector<double> line, double eps, double x0, bool flag, bool flag2) {

    double b = line[0];
    double a = line[1];
    double f1 = f(a, key);
    double f2 = f(b, key);
    double xcur;


    xcur = x0;

    double xtmp = xcur + 1;
    int iter = 0;

    cout << "Start:" << endl;
    while (true) {


        iter += 1;


        xtmp = xcur;
        if (flag) {
            xcur = xcur - f(xcur, key) / df(xcur, key);
        }
        else {
            xcur = xcur - f(xcur, key) / df(xcur, key, eps);
        }

        if ((xcur < line[0]) or (line[1] < xcur)){
            cout << "X IN OUT LINE: x_old = " << xtmp << " xcur = " << xcur << " f(xcur) = " << f(xcur, key) << " df(xcur) = " << df(xcur, key);

        }

        // Вывод итераций
        cout << xcur << endl;

        // Критерий остановы
        if (flag2) {
            if (abs(xcur - xtmp) < eps) {
                break;
            }
        }

        else {
            if (abs(f(xcur, key)) < eps) {
                break;
            }
        }

        // Ограничение итераций
        if (iter > 30) {
            iter = 30;
            xcur = xtmp;
            break;
        }
    }
    cout << "end." << endl;
    cout << "Iterations: " << iter << endl;
    return xcur;
}

//double Method_Newton2(int key, vector<double> line, double eps, double alpha, double x0, bool flag, bool flag2) {
//
//    double b = line[0];
//    double a = line[1];
//    double f1 = f(a, key);
//    double f2 = f(b, key);
//    double xcur;
//
//    if (x0 == -100) {
//        xcur = (f1 * b - f2 * a) / (f1 - f2);
//    }
//    else {
//        xcur = x0;
//    }
//    double xtmp = xcur + 1;
//    int iter = 0;
//    double tmp1 = 1;
//    double tmp2 = 1;
//
//    double  h = 1e-7;
//    while (true) {
//        if (iter > 30) {
//            iter = 30;
//            xcur = xtmp;
//            break;
//        }
//
//        //tmp1 = xcur - xtmp;
//        iter += 1;
//        xtmp = xcur;
//        if (flag) {
//            /*if (df(xcur, key) < 0.0001 * eps * eps) {
//                break;
//            }*/
//            xcur = xcur - alpha * f(xcur, key) / df(xcur, key);
//        }
//        else {
//            xcur = xcur - alpha * 2 * h * f(xcur, key) / (f(xcur + h, key) - f(xcur - h, key));
//        }
//
//        //tmp2 = xcur - xtmp;
//        //ot.push_back(tmp2 / tmp1);
//        if (flag2) {
//            if (abs(xcur - xtmp) < eps) {
//                break;
//            }
//            //cout << xcur << endl;
//        }
//
//        else {
//            if (abs(f(xcur, key)) < eps) {
//                //cout << xcur << " " << f(xcur, key) << endl;
//                break;
//            }
//            //cout << xcur << endl;
//        }
//        // Вывод итераций
//        //cout << xcur << endl;
//
//    }
//    cout << "Iterations: " << iter << endl;
//    return xcur;
//}

double Method_Newton_mod(int key, vector<double> line, double x0, double eps, double alpha, bool flag, bool flag2) {

    double b = line[0];
    double a = line[1];
    //cout << "line = (" << line[0] << ", " << line[1] << ")" << endl;
    double f1 = f(a, key);
    double f2 = f(b, key);
    double xcur = x0;

    double x_old = xcur + 1;
    int iter = 0;


    cout << "Start:" << endl;
    while (true) {
        if (iter > 30) {
            iter = 30;
            xcur = x_old;
            break;
        }

        iter += 1;
        //xtmp = xcur;


        x_old = xcur;
        // Выбор дифференцирования
        if (flag){
            xcur = xcur - alpha * f(xcur, key) / df(xcur, key);
        } else {
            xcur = xcur - alpha * f(xcur, key) / df(xcur, key, eps);
        }



        if ((xcur < line[0]) or (line[1] < xcur) or xcur == 1){
            // Запускаем метод хорд
            cout << "METHOD HORD: x_old = " << xcur;

            // Проверка начальной точки
            if ((x0 < line[0]) or (x0 < xcur)){

               //a = x_old - eps;
               //b = x_old + eps;

                if (xcur > x_old) {
                    a = x_old;
                    b = line[1];
                } else {
                    a = line[0];
                    b = x_old;
                }
                xcur = (f(a, key) * b - f(b, key) * a) / (f(a, key) - f(b, key));
                cout << " ->  x_new = " << xcur << endl;
            }
        } else {
            // Вывод итераций
            cout << xcur << endl;
        }


        // Критерий останова
        if (flag2) {
            if (abs(xcur - x_old) < eps) {
                break;
            }
        }

        else {
            if (abs(f(xcur, key)) < eps) {
                //cout << xcur << " " << f(xcur, key) << endl;
                break;
            }
        }

    }
    cout << "end." << endl;
    cout << "Iterations: " << iter << endl;
    return xcur;
}


bool in(vector<double> elems, double elem){
    int n = elems.size();

    for (int i = 0; i < n; i++){
        if (elems[i] == elem)
            return true;
    }
    return false;
}

double Method_Newton_mod2(int key, vector<double> line, double x0, double eps, double alpha, bool flag, bool flag2) {

    double b = line[0];
    double a = line[1];
    //cout << "line = (" << line[0] << ", " << line[1] << ")" << endl;
    double f1 = f(a, key);
    double f2 = f(b, key);
    double xcur = x0;

    double x_old = xcur + 1;
    int iter = 0;

    vector<double> all_xk = {xcur};

    cout << "Start:" << endl;
    while (true) {
        if (iter > 30) {
            iter = 30;
            xcur = x_old;
            break;
        }

        iter += 1;
        //xtmp = xcur;


        x_old = xcur;
        // Выбор дифференцирования
        if (flag){
            xcur = xcur - alpha * f(xcur, key) / df(xcur, key);
        } else {
            xcur = xcur - alpha * f(xcur, key) / df(xcur, key, eps);
        }



        if ((xcur < line[0]) or (line[1] < xcur) or in(all_xk, xcur)) {
            // Запускаем метод хорд
            cout << "METHOD HORD: x_old = " << xcur;

            if (xcur > x_old) {
                a = x_old;
                b = line[1];
            } else {
                a = line[0];
                b = x_old;
            }

            // Проверка начальной точки
            if ((x0 < line[0]) or (x0 < xcur)) {
                if (xcur > x0) {
                    a = x0;
                    b = line[1];
                } else {
                    a = line[0];
                    b = x0;
                }
                x0 = xcur;
            }

            xcur = (f(a, key) * b - f(b, key) * a) / (f(a, key) - f(b, key));
            cout << " ->  x_new = " << xcur << endl;

        } else {
            // Вывод итераций
            cout << xcur << endl;
        }

        all_xk.push_back(xcur);
        // Критерий останова
        if (flag2) {
            if (abs(xcur - x_old) < eps) {
                break;
            }
        }

        else {
            if (abs(f(xcur, key)) < eps) {
                //cout << xcur << " " << f(xcur, key) << endl;
                break;
            }
        }

    }
    cout << "end." << endl;
    cout << "Iterations: " << iter << endl;
    return xcur;
}


/* Функции для систем уравнений */


/* Функция, вычисляющая норму разности векторов */
double sqr(vector<double> vec1, vector<double> vec2) {
    int m = vec1.size();
    double sum;
    for (int i = 0; i < m; i++) {
        sum = (vec1[i] - vec2[i]) * (vec1[i] - vec2[i]);
    }
    return sum;
}

//vector<double> Method_Newton_sys(int key, vector<double> line, double eps, bool flag) {
//
//    double F11, F12, F21, F22;
//    double h = 1e-7;
//    vector<double> xtmp(2, 0), xcur(2, 0);
//    xtmp = line;
//    int iter = 0;
//    while (sqr(xcur, xtmp) > eps) {
//        iter++;
//        if (iter > 30) {
//            iter = 30;
//            break;
//        }
//        xtmp = xcur;
//        if (flag) {
//            vector<double> tmp1 = dg(xtmp[0], xtmp[1], key)[0];
//            F11 = tmp1[0];
//            F12 = tmp1[1];
//            vector<double> tmp2 = dg(xtmp[0], xtmp[1], key)[1];
//            F21 = tmp2[0];
//            F22 = tmp2[1];
//        } else {
//            F11 = ((g(xtmp[0] + h, xtmp[1], key))[0] - (g(xtmp[0] - h, xtmp[1], key))[0]) / (2 * h);
//            F12 = (g(xtmp[0], xtmp[1] + h, key)[0] - g(xtmp[0], xtmp[1] - h, key)[0]) / (2 * h);
//            F21 = (g(xtmp[0] + h, xtmp[1], key)[1] - g(xtmp[0] - h, xtmp[1], key)[1]) / (2 * h);
//            F22 = (g(xtmp[0], xtmp[1] + h, key)[1] - g(xtmp[0], xtmp[1] - h, key)[1]) / (2 * h);
//        }
//        double detf = F11 * F22 - F12 * F21;
//        if (abs(detf) < (10 * h) * (10 * h)) {
//            iter = 0;
//            break;
//        }
//        xcur[0] = xtmp[0] - (F22 * g(xtmp[0], xtmp[1], key)[0] - F12 * g(xtmp[0], xtmp[1], key)[1]) / detf;
//        xcur[1] = xtmp[1] - (-F21 * g(xtmp[0], xtmp[1], key)[0] + F11 * g(xtmp[0], xtmp[1], key)[1]) / detf;
//    }
//    cout << iter << endl;
//    cout << xcur[0] << " " << xcur[1] << endl;
//    return xcur;
//}


vector<double> Method_Newton_sys(int key1, int key2, vector<double> x0, double L1, double L2, double eps, bool flag) {

    double F11, F12, F21, F22;
    vector<double> xtmp(2, 0), xcur(2, 0), xprev(2, 0);
    vector<double> a(2, 0), b(2, 0);
    xcur = x0;
    xtmp = {0, 0};
    int iter = 0;

    cout << "START" << endl;
    while (abs(sqr(xcur, xtmp)) > eps / 10 ) {

        cout << "(" << xcur[0] << ", " << xcur[1] << ")" << endl;

        iter++;
        if (iter > 30) {
            iter = 30;
            break;
        }
        xtmp = xcur;

        // Выбор дифференцирования
        if (flag) {
            vector<double> tmp1 = dg(xtmp[0], xtmp[1], key1);
            F11 = tmp1[0];
            F12 = tmp1[1];
            vector<double> tmp2 = dg(xtmp[0], xtmp[1], key2);
            F21 = tmp2[0];
            F22 = tmp2[1];
        }
        else {
            F11 = (g(xtmp[0] + eps, xtmp[1], key1) - g(xtmp[0] - eps, xtmp[1], key1)) / (2 * eps);
            F12 = (g(xtmp[0], xtmp[1] + eps, key1) - g(xtmp[0], xtmp[1] - eps, key1)) / (2 * eps);
            F21 = (g(xtmp[0] + eps, xtmp[1], key2) - g(xtmp[0] - eps, xtmp[1], key2)) / (2 * eps);
            F22 = (g(xtmp[0], xtmp[1] + eps, key2) - g(xtmp[0], xtmp[1] - eps, key2)) / (2 * eps);
        }

        double detf = F11 * F22 - F12 * F21;
        if (detf < eps * eps) {
            cout << "det(Jacobian) = 0" << endl;
            iter = 0;
            return xcur;
        }
        xprev = xcur;

        double del1 = (F22 * g(xtmp[0], xtmp[1], key1) - F12 * g(xtmp[0], xtmp[1], key2)) / detf;
        double del2 = (-F21 * g(xtmp[0], xtmp[1], key1) + F11 * g(xtmp[0], xtmp[1], key2)) / detf;
        xcur[0] = xtmp[0] - del1;
        xcur[1] = xtmp[1] - del2;


        // Метод хорд, если вышли за допустимую область
        if ((abs(xcur[0]) >  L1) or (abs(xcur[1]) >  L2)) {

            cout << "Method Hord start: x_old =" << "(" << xcur[0] << ", " << xcur[1] << ")";
            double a1 = xprev[0] - eps;
            double b1 = xprev[0] + eps;

            double a2 = xprev[1] - eps;
            double b2 = xprev[1] + eps;

            double f1a = g(a1, xprev[1], key1);
            double f1b = g(xprev[0], b1, key1);

            double f2a = g(a2, xprev[1], key2);
            double f2b = g(xprev[0], b2, key2);

            xcur[0] = (f1a * b1 - f1b * a1) / (f1a - f1b);
            xcur[1] = (f2a * b2 - f2b * a2) / (f2a - f2b);
            cout << " -> x_new = (" << xcur[0] << ", " << xcur[1] << ")" << endl;

        }


    }
    cout << "END." << endl;
    cout << "Iterations: " << iter << endl;
    return xcur;
}

//vector<double> Method_Newton_sys3(int key1, int key2, vector<double> x0, double L1, double L2, double eps, bool flag) {

//    double J11, J12, J21, J22;
//
//    double xy_deltax, xy_deltay;
//    double alpha = 1;
//
//    int iter = 0;
//    double delta = 1e-10;
//    vector<double> temp1;
//    vector<double> temp2;
//
//    vector<double> xy_cur = { x0[0], x0[1] };
//    vector<double> temp = { 0, 0 };
//
//    while (sqr(xy_cur, temp) >= eps) {
//        if (iter > 30) {
//            iter = 30;
//            xy_cur = temp;
//            break;
//        }
//
//        temp = xy_cur;
//        J11 = (g(xy_cur[0] + delta, xy_cur[1], key1) - g(xy_cur[0] - delta, xy_cur[1], key1)) / (2 * delta);
//        J12 = (g(xy_cur[0], xy_cur[1] + delta, key1) - g(xy_cur[0], xy_cur[1] - delta, key1)) / (2 * delta);
//        J21 = (g(xy_cur[0] + delta, xy_cur[1], key2) - g(xy_cur[0] - delta, xy_cur[1], key2)) / (2 * delta);
//        J22 = (g(xy_cur[0], xy_cur[1] + delta, key2) - g(xy_cur[0], xy_cur[1] - delta, key2)) / (2 * delta);
//
//        iter = iter + 1;
//
//        double Jcb = J11 * J22 - J12 * J21;
//        cout << Jcb << endl;
//        if (abs(Jcb) < eps * eps * eps) {
//            cout << "break";
//            iter = 0;
//            break;
//        }
//
//        xy_deltax = (J22 * g(xy_cur[0], xy_cur[1], key1) - J12 * g(xy_cur[0], xy_cur[1], key2)) / Jcb;
//        xy_deltay = (-J21 * g(xy_cur[0], xy_cur[1], key1) + J11 * g(xy_cur[0], xy_cur[1], key2)) / Jcb;
//
//        xy_cur[0] -= alpha * xy_deltax;
//        xy_cur[1] -= alpha * xy_deltay;
//        //обработка случая выхода за пределы области
//        /*if ((abs(xy_cur[0]) > bx) or (abs(xy_cur[1]) > by))
//        {
//            xy_cur[0] += kf * xy_deltax;
//            xy_cur[1] += kf * xy_deltay;
//            kf /= 2;
//            xy_cur[0] -= kf * xy_deltax;
//            xy_cur[1] -= kf * xy_deltay;
//        }*/
//    }
//
//    cout << "Iterations: " << iter << endl;
//    return xy_cur;
//}



/* Функция метода Ньютона для рисовки бассейнов Ньютона */
vector<double> Method_Newton_sys2(int key1, int key2, double L1, double L2, double eps, bool flag, bool flag2) {

    double F11, F12, F21, F22;
    double h = 1e-6;
    int N = 200;
    vector<vector<int>> iteration(N, vector<int>(N, 0));

    /*double L1 = 10;
    double L2 = 10;*/
    double hx = 2 * L1 / N;
    double hy = 2 * L2 / N;

    vector<double> xtmp(2, 0), xcur(2, 0);
    vector<double> tmp1, tmp2;

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            xcur[0] = -L1 + i * hx;
            xcur[1] = -L2 + j * hy;
            xtmp = {0, 0};
            int iter = 0;
            while ((xcur[0] - xtmp[0])* (xcur[0] - xtmp[0]) + (xcur[1] - xtmp[1])* (xcur[1] - xtmp[1]) > eps * eps) {
                iter++;
                if (iter > 30) {
                    iter = 31;
                    xcur = xtmp;
                    break;
                }
                xtmp = xcur;
                if (flag) {
                    tmp1 = dg(xtmp[0], xtmp[1], key1);
                    F11 = tmp1[0];
                    F12 = tmp1[1];
                    tmp2 = dg(xtmp[0], xtmp[1], key2);
                    F21 = tmp2[0];
                    F22 = tmp2[1];
                }
                else {
                    F11 = (g(xtmp[0] + h, xtmp[1], key1) - g(xtmp[0] - h, xtmp[1], key1)) / (2 * h);
                    F12 = (g(xtmp[0], xtmp[1] + h, key1) - g(xtmp[0], xtmp[1] - h, key1)) / (2 * h);
                    F21 = (g(xtmp[0] + h, xtmp[1], key2) - g(xtmp[0] - h, xtmp[1], key2)) / (2 * h);
                    F22 = (g(xtmp[0], xtmp[1] + h, key2) - g(xtmp[0], xtmp[1] - h, key2)) / (2 * h);
                }
                double detf = F11 * F22 - F12 * F21;
                if (abs(detf) < 1e-20) {
                    iter = 0;
                    xcur = xtmp;
                    break;
                }

                vector<double> xprev = xcur;

                xcur[0] = xtmp[0] - (F22 * g(xtmp[0], xtmp[1], key1) - F12 * g(xtmp[0], xtmp[1], key2)) / detf;
                xcur[1] = xtmp[1] - (-F21 * g(xtmp[0], xtmp[1], key1) + F11 * g(xtmp[0], xtmp[1], key2)) / detf;

                if (abs(detf) < 1e-20 /*h * h*/) {
                    cout << "Method Hord start";
                    double a1 = xprev[0] - eps;
                    double b1 = xprev[0] + eps;

                    double a2 = xprev[1] - eps;
                    double b2 = xprev[1] + eps;

                    double f1a = g(a1, xprev[1], key1);
                    double f1b = g(xprev[0], b1, key1);

                    double f2a = g(a2, xprev[1], key2);
                    double f2b = g(xprev[0], b2, key2);

                    xcur[0] = (f1a * b1 - f1b * a1) / (f1a - f1b);
                    xcur[1] = (f2a * b2 - f2b * a2) / (f2a - f2b);
                }
            }
            iteration[i][j] = iter;
            /* if (((xcur[0]-4)*(xcur[0]-4)+(xcur[1]+1)*(xcur[1]+1) > eps*eps*eps) || ((xcur[0] + 4) * (xcur[0] + 4) + (xcur[1] - 1) * (xcur[1] - 1) > eps * eps*eps))
             {
                 cout << iter << endl;
                 cout << xtmp[0] << " " << xtmp[1] << endl;
                 cout << g(xcur[0], xcur[1], key) << endl;
             }*/
            //if ((xcur[0] - 4) + (xcur[1] + 1) < 2 * eps) {
                //cout << -L1 + i * hx << " " << -L2 + j * hy << endl;
                //cout << xcur[0] << " " << xcur[1] << endl;
            //}
        }
    }
    if (flag2) {
        ofstream file("data.txt");
        //file << xcur[0] << " " << xcur[1] << " " << N << endl;
        file << L1 << " " << L2 << " " << N << endl;
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < N; ++j) {
                file << iteration[j][i] << " ";
            }
            file << endl;
        }
        file.close();
    }
    //cout << xcur[0] << " " << xcur[1] << endl;
    return xcur;
}