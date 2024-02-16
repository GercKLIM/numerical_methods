#include <vector>

using namespace std;



/* ### Объявление тестов ### */



// Тест 0 (Маятник)
vector<double> ODU_0(const double& t, const vector<double>& U){
    int n = U.size();
    vector<double> dU(n, 0);

    double k = 1, m = 1;

    dU[0] = U[1];
    dU[1] = (-1 * k / m) * U[0];

    return dU;
}

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
    return ODU_0(1, U);
}
