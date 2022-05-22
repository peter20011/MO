#include <cmath>
#include <iomanip>
#include <iostream>
#include <fstream>
#define START 0.0
#define STOP 1.0
using namespace std;


double* allocateVector(int n) {
    return new double[n];
}

void printVector(double* vector, int n) {
    for (int i = 0; i < n; i++) {
            cout << setw(10) << vector[i] << endl;
    }
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
        D[i] -= (U[i] * L[i - 1]) / D[i - 1];

    for (int i = 1; i < n; i++)
        b[i] -= (U[i] * b[i - 1]) / D[i - 1];

    x[n - 1] = b[n - 1] / D[n -1];
    for (int i = n - 2; i >= 0; i--)
        x[i] = (b[i] - L[i] * x[i + 1]) / D[i];
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
void getError(double* x, double* errors, double h, int n) {
    double x_n = START;

    for (int i = 0; i < n; i++) {
        errors[i] = fabs(x[i] - analytical(x_n));
        x_n += h;
    }
}
