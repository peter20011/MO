#include<iostream>
#include<fstream>
#include<cmath>
#include<string>
#include<fstream>
#include <algorithm>
const double D = 1.0;
const double TAU = 0.1;
const double TMIN = 0.0;
const double TMAX = 2.0;
const double A = 6*sqrt(D*(TAU+TMAX))+0.1;
const double x_start =-A ;
const double x_end = A;
const double LAMBDA_BEZPOSREDNIE = 0.4;
const double LAMBDA_POSREDNIE = 1.0;
const double h = 0.025; 

using namespace std;

//Zapisywanie do pliku na wygląd siatki czasowo-przestrzennej
void zapis_do_pliku(const string &nazwa, double **macierz, int n, int m) {
    ofstream out;
    out.open("../" + nazwa);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            out << macierz[i][j] << " ";
        }
        out << endl;
    }
    out.close();
}

//alokowanie macierzy
double **tworzenie_macierzy(int n, int m) {
    auto **macierz = new double *[n];
    for (int i = 0; i < n; i++) {
        macierz[i] = new double[m];
    }
    return macierz;
}

//dealokowanie macierzy
void deleteMatrix(double **macierz, int n) {
    for (int i = 0; i < n; i++) {
        delete[] macierz[i];
    }
    delete[] macierz;
}

//wyliczanie maxError dla każdego t
double *maxError(double **error, int n, int m) {
    double maxE;
    auto *errors = new double[n];
    for (int i = 0; i < n; i++) {
        maxE = fabs(error[i][0]);
        for (int j = 0; j < m; j++) {
            if (maxE < fabs(error[i][j])) {
                maxE = fabs(error[i][j]);
            }
        }
        errors[i] = maxE;
    }
    return errors;
}

//obliczanie błędu bezwzględnego dla rozwiązania numerycznego
void countErrors(double **error, double **analytical, double **numerical, int n, int m) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            error[i][j] = fabs(numerical[i][j] - analytical[i][j]);
        }
    }
}

//wyliczenie kroku czasowego dla podanej wartości lambda oraz kroku przestrzennego
double getDelta_t(double lambda) {
    return (lambda * h * h) / D;
}

bool error(double *xn, double *xn_1, int m) {
    int count = 0;

    for (int i = 0; i < m; i++) {
        if (fabs(xn_1[i] - xn[i]) < 1e-10) count++;
    }
    return count == m;
}

bool residue(double **A, double *b, const double *xn, int m) {
    int count = 0;
    for (int i = 0; i < m; i++) {
        double sum = 0.0;
        for (int j = 0; j < m; j++) {
            sum += A[i][j] * xn[j];
        }
        if (fabs(sum - b[i]) < 1e-10) count++;
    }
    return count == m;
}