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
const double XMIN =-A ;
const double XMAX = A;
const double LAMBDA_BEZPOSREDNIE = 0.4;
const double LAMBDA_POSREDNIE = 1.0;
const double h = 0.025; 

using namespace std;

//Zapisywanie do pliku macierzy
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

//usuwanie macierzy
void usuwanie_macierzy(double **macierz, int n) {
    for (int i = 0; i < n; i++) {
        delete[] macierz[i];
    }
    delete[] macierz;
}

//wyliczanie bledu dla każdego przedziału czasowego
double *macierzBlad(double **blad, int n, int m) {
    double maxBlad;
    auto *bledy = new double[n];
    for (int i = 0; i < n; i++) {
        maxBlad = fabs(blad[i][0]);
        for (int j = 0; j < m; j++) {
            if (maxBlad < fabs(blad[i][j])) {
                maxBlad = fabs(blad[i][j]);
            }
        }
        bledy[i] = maxBlad;
    }
    return bledy;
}

//obliczanie błędu bezwzględnego dla rozwiązania numerycznego
void oblicz_blad(double **blad, double **anaalityczne, double **numerycczne, int n, int m) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            blad[i][j] = fabs(numerycczne[i][j] - anaalityczne[i][j]);
        }
    }
}

//wyliczenie kroku czasowego dla podanej wartości lambda 
double obliczdt(double lambda) {
    return (lambda * h * h) / D;
}
// do sora sprawdzenie warunku bledu
bool blad(double *x, double *x1, int m) {
    int count = 0;

    for (int i = 0; i < m; i++) {
        if (fabs(x1[i] - x[i]) < 1e-10) count++;
    }
    return count == m;
}
//do sora sprawdzenie  warunku residuum
bool residuum(double **macierz, double *b, const double *x, int m) {
    int licznik = 0;
    for (int i = 0; i < m; i++) {
        double sum = 0.0;
        for (int j = 0; j < m; j++) {
            sum += macierz[i][j] * x[j];
        }
        if (fabs(sum - b[i]) < 1e-16) licznik++;
    }
    return licznik == m;
}

void Thomas(double *L, double *D, double *U, double *b, double *x, const int m) {
    double *l = new double[m]; 
    for (int i = 1; i < m; i++) {
        l[i] = L[i] * (1 / D[i - 1]); 
        D[i] -= l[i] * U[i - 1]; 
    }
    for (int i = 1; i < m; i++) {
        b[i] -= l[i] * b[i - 1]; 
    }
    x[m - 1] = (1 / D[m - 1] * b[m - 1]); 
    for (int i = m - 2; i >= 0; i--) {
        x[i] = (1 / D[i]) * (b[i] - U[i] * x[i + 1]); 
    }
    delete[] l;
}



void sor(double **macierz, double *b, double *x, int m) {
     double omega = 1.0; //dla omega = 1.0 otrzymałem najszybszą zbierzność
    double *pomocniczy = new double[m]; //wektor przechowujący wyniki prawej strony wzoru operacyjnego
    double *przyblizenie = new double[m]; //wektor przyblizen poczatkowych
    double *x1 = new double[m]; //wektor kolejnego przyblizenia
    for (int i = 0; i < m; i++) {
        pomocniczy[i] = 0.0; 
        przyblizenie[i] = 1.0; 
    }

    for (int iter = 0; iter < 100; iter++) { //arbitralne ograniczenie na liczbe iteracji
        for (int i = 0; i < m; i++) {
            double tmp = ((1.0 - 1.0 / omega) * macierz[i][i]) * przyblizenie[i]; 

            for (int j = i + 1; j < m; j++) {
                tmp += macierz[i][j] * przyblizenie[j]; 
            }
            pomocniczy[i] = -tmp + b[i]; 
        }

        for (int i = 0; i < m; i++) {
            x1[i] = 0.0; //zerowanie nastepnego przyblizenia
        }

        for (int i = 0; i < m; i++) {
            double tmp = 0.0;
            for (int j = 0; j <= i; j++) {
                tmp += x1[j] * macierz[i][j]; // L*x_n
            }
            x1[i] = (pomocniczy[i] - tmp) / ((1.0 / omega) * macierz[i][i]); //wyliczenie x_n
        }

        bool blad_2 = blad(x1, przyblizenie, m); //kryterium dokładności wyznaczenia
        bool resuduum_2 = residuum(macierz, b, przyblizenie, m); // kryterium wiarygodności x_n
        if ( resuduum_2==true && blad_2==true) break;
        swap(x1, przyblizenie);
    }
    for (int i = 0; i < m; i++) {
        x[i] = x1[i];
    }
}