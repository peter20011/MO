#include<iostream>
#include<fstream>
#include<cmath>
#include<string>
#include<fstream>
#include <algorithm>
const double D = 1.0;
const double TAU = 0.1;

const double t_min = 0.0;
const double t_max = 2.0;
const double a = 6*sqrt(D*(TAU+t_max))+0.1;
const double x_start =-a ;
const double x_end = a;

const double lambda_kmb = 0.4;
const double lambda_laasonen = 1.0;

const double h = 0.025; 

using namespace std;

//Zapisywanie do pliku na wygląd siatki czasowo-przestrzennej
void saveToFile(const string &fileName, double **matrix, int n, int m) {
    ofstream out;
    out.open("../" + fileName);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            out << matrix[i][j] << " ";
        }
        out << endl;
    }
    out.close();
}

//alokowanie macierzy
double **allocMatrix(int n, int m) {
    auto **matrix = new double *[n];
    for (int i = 0; i < n; i++) {
        matrix[i] = new double[m];
    }
    return matrix;
}

//dealokowanie macierzy
void deleteMatrix(double **matrix, int n) {
    for (int i = 0; i < n; i++) {
        delete[] matrix[i];
    }
    delete[] matrix;
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