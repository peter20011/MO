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
const double A = 6*sqrt(D*(TAU+TMIN))+0.1;
const double XMIN = -A;
const double XMAX = A;
const double LAMBDA_BEZPOSREDNIE = 0.4;
const double LAMBDA_POSREDNIE = 1.0;
const double h = 0.10;


using namespace std;

//Zapisywanie do pliku na wygląd siatki czasowo-przestrzennej
void zapisz_do_pliku(const string &nazwa_pliku, double **macierz, int n, int m) {
    ofstream out;
    out.open("../" + nazwa_pliku);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            out << macierz[i][j] << " ";
        }
        out << "\n";
    }
    out.close();
}

//alokowanie macierzy
double **alokowanie_macierzy(int n, int m) {
    auto **macierz = new double *[n];
    for (int i = 0; i < n; i++) {
        macierz[i] = new double[m];
    }
    return macierz;
}

//dealokowanie macierzy
void  usuwanie_macierzy(double **macierz, int n) {
    for (int i = 0; i < n; i++) {
        delete[] macierz[i];
    }
    delete[] macierz;
}

//wyliczanie maxError dla każdego t
double *najwiekszy_blad(double **blad, int n, int m) {
    double najwiekszy_blad;
    auto *bledy = new double[n];
    for (int i = 0; i < n; i++) {
        najwiekszy_blad = fabs(blad[i][0]);
        for (int j = 0; j < m; j++) {
            if (najwiekszy_blad < fabs(blad[i][j])) {
                najwiekszy_blad = fabs(blad[i][j]);
            }
        }
        bledy[i] = najwiekszy_blad;
    }
    return bledy;
}

//obliczanie błędu bezwzględnego dla rozwiązania numerycznego
void obliczanie_bledu(double **blad, double **analityczne, double **numerycznie, int n, int m) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            blad[i][j] = fabs(numerycznie[i][j] - analityczne[i][j]);
        }
    }
}

//wyliczenie kroku czasowego dla podanej wartości lambda oraz kroku przestrzennego
double oblicz_dt(double lambda) {
    return (lambda * h * h) / D;
}