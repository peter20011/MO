#define _USE_MATH_DEFINES 
#include <cmath>
#include <iomanip>
#include <iostream>
#include <fstream>

// ZADANE FUNKCJE
#define P 1.0
#define R -4.0

#define ALPHA 0.0
#define BETA 1.0
#define GAMMA -1.0

// ZADANE WSPOLCZYNNIKI
#define FI 0.0
#define PSI 1.0
#define THETA 0.0

// OBSZAR ROZWIAZANIA
#define START 0.0
#define STOP 1.0

using namespace std;

double* allocateVector(int n) {
    return new double[n];
}

double vectorNormMax(double* v, int n) {
    double currentMax = fabs(v[0]);

    for (int i = 1; i < n; i++)
        if (fabs(v[i]) > currentMax)
            currentMax = fabs(v[i]);

    return currentMax;
}

double analytical(double x) {
    return (exp(2.0 - 2.0 * x) - 4.0 * exp(4.0 - 2.0 * x) + 4.0 * exp(2.0 * x) - exp(2.0 + 2.0 * x) - x + x * exp(4.0)) / (4.0 - 4.0 * exp(4.0));
}

void Thomas(double* L, double* D, double* U, double* b, double* x, int n) {
    for (int i = 1; i < n; i++) 
        D[i] -= (U[i-1] * L[i]) / D[i - 1];

    for (int i = 1; i < n; i++)
        b[i] -= (L[i] * b[i - 1]) / D[i - 1];

    x[n - 1] = b[n - 1] / D[n -1];
    for (int i = n - 2; i >= 0; i--)
        x[i] = (b[i] - U[i] * x[i + 1]) / D[i];
}

void fillConventionalMatrix(double* L, double* D, double* U, double* b, double h, int n) {
    L[0] = 0.0;
    D[0] = BETA - (ALPHA / h);
    U[0] = ALPHA / h;
    b[0] = (-1) * GAMMA;

    for (int i = 1; i < n - 1; i++) {
        L[i] = P / (h * h);
        D[i] = R + (-2.0 * P) / (h * h);
        U[i] = P / (h * h);
        b[i] = i * h;
    }

    L[n - 1] = -FI / h;
    D[n - 1] = -FI / h + PSI;
    U[n - 1] = 0.0;
    b[n - 1] = (-1) * THETA;
}

void fillNumerovMatrix(double* L, double* D, double* U, double* b, double h, int n) {
    L[0] = 0.0;
    D[0] = BETA - (ALPHA / h);
    U[0] = ALPHA / h;
    b[0] = (-1) * GAMMA;

    for (int i = 1; i < n - 1; i++) {
        L[i] = P / (h * h) + R / 12.0;
        D[i] = (-2.0 * P) / (h * h) + R * (10.0 / 12.0);
        U[i] = P / (h * h) + R / 12.0;
        b[i] = (i * h - h) / 12.0 + (10.0 / 12.0) * (i * h) + (i * h + h) / 12.0;
    }

    L[n - 1] = -FI / h;
    D[n - 1] = -FI / h + PSI;
    U[n - 1] = 0.0;
    b[n - 1] = (-1) * THETA;
}

void getError(double* x, double* errors, double h, int n) {
    double x_n = START;

    for (int i = 0; i < n; i++) {
        errors[i] = fabs(x[i] - analytical(x_n));
        x_n += h;
    }
}

void saveToFile(string fileName, double* x, double h, int n) {
    ofstream file;
    file.open(fileName, ios::out);
    if (!file.good())
        exit(EXIT_FAILURE);

    double x_n = START;
    for (int i = 0; i < n; i++) {
        file << x_n << " " << x[i] << " " << analytical(x_n) << endl;
        x_n += h;
    }
}

double NumerowMethod(double h, int n) {
    double* L = allocateVector(n);
    double* D = allocateVector(n);
    double* U = allocateVector(n);
    double* b = allocateVector(n);
    double* x = allocateVector(n);
    double* errors = allocateVector(n);

    fillNumerovMatrix(L, D , U, b, h, n);
    Thomas(L, D, U, b, x, n);
    getError(x, errors, h, n);


    if (n == 210)
        saveToFile("Numerov.txt", x, h, n);

    double maxError = vectorNormMax(errors, n);

    delete[] L;
    delete[] D;
    delete[] U;
    delete[] b;
    delete[] x;
    delete[] errors;

    return maxError;
}

double conventionalMethod(double h, int n) {
    double* L = allocateVector(n);
    double* D = allocateVector(n);
    double* U = allocateVector(n);
    double* b = allocateVector(n);
    double* x = allocateVector(n);
    double* errors = allocateVector(n);

    fillConventionalMatrix(L, D, U, b, h, n);
    Thomas(L, D, U, b, x, n);
    getError(x, errors, h, n);

    if (n == 210)
        saveToFile("Conventional.txt", x, h, n);

    double maxError = vectorNormMax(errors, n);

    delete[] L;
    delete[] D;
    delete[] U;
    delete[] b;
    delete[] x;
    delete[] errors;

    return maxError;
}

int main()
{
    ofstream errors;
    errors.open("errors.txt", ios::out);
    errors << scientific;

    double h;

    for (int i = 10; i < 200000; i += 100) {
        h = (STOP - START) / (i - 1);
        errors << setw(16) << log10(h) << " " << setw(16) << log10(conventionalMethod(h, i)) << " " << setw(16) << log10(NumerowMethod(h, i)) << endl;
    }

    errors.close();
    return 0;
}

