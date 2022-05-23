#include <iostream>
#include <cmath>
#include <fstream>
#include <cfloat>

using namespace std;
int numer = 0;

inline double rel_err(double calculated, double exact) {//Funkcja oblicza blad wzgledny - aby zapobiec dzieleniu przez 0 dgy wartość dokladna = 0 zwraca 0
    if (exact == 0) return 0;
    else return (exact - calculated) / exact;
}

int print_err(double **anal, const char arg0[], int imax, int jmax, double** arg1, double xpocz, double h, double dt) {
    // static int numer;
    fstream out;
    char nn [25];
    nn[0] = numer + 48;
    numer++;
    nn[1] = '.';
    int i;
    for (i = 2; i < 20; i++) {
        nn[i] = arg0[i - 2];
        if (!arg0[i - 2]) break;
    }
    nn[i] = '.';
    nn[i + 1] = 'e';
    nn[i + 2] = 'r';
    nn[i + 3] = 'r';
    nn[i + 4] = 0;
    out.open(nn, fstream::out);
    if (out.bad()) {
        return 0;
    }
    out << arg0 << "-bledy \th=" << h << "\tdt=" << dt << "\n x/t \t";
    for (int j = 0; j < imax; j++) {
        out << j * dt << "\t";
    }
    for (int i = 0; i < jmax; i++) {
        out << "\n";
        out << xpocz + h * i << "\t";
        for (int j = 0; j < imax; j++) {
            out << fabs(rel_err(arg1[j][i], anal[i][j])) << "\t";
        }
    }
    out << "\n\n" << arg0 << " bledy bezwgledne\th=" << h << "\tdt=" << dt << "\n x/t \t";
    for (int j = 0; j < imax; j++) {
        out << j * dt << "\t";
    }
    for (int i = 0; i < jmax; i++) {
        out << "\n";
        out << xpocz + h * i << "\t";
        for (int j = 0; j < imax; j++) {
            out << fabs(arg1[j][i] - anal[i][j]) << "\t";
        }
    }
    return 1;
}

int print_err1(double **anal, const char arg0[], double** arg1, int imax, int jmax, double xpocz, double h, double dt) {
    // static int numer;
    fstream out;
    char nn [25];
    nn[0] = numer + 48;
    numer++;
    nn[1] = '.';
    int i;
    for (i = 2; i < 20; i++) {
        nn[i] = arg0[i - 2];
        if (!arg0[i - 2]) break;
    }
    nn[i] = '.';
    nn[i + 1] = 'e';
    nn[i + 2] = 'r';
    nn[i + 3] = 'r';
    nn[i + 4] = 0;
    out.open(nn, fstream::out);
    if (out.bad()) {
        return 0;
    }
    out << arg0 << "-bledy \th=" << h << "\tdt=" << dt << "\n x/t \t";
    for (int j = 0; j < jmax; j++) {
        out << j * dt << "\t";
    }
    for (int i = 0; i < imax; i++) {
        out << "\n";
        out << xpocz + h * i << "\t";
        for (int j = 0; j < jmax; j++) {
            out << fabs(rel_err(arg1[i][j], anal[i][j])) << "\t";
        }
    }
    out << "\n\n" << arg0 << " bledy bezwgledne\th=" << h << "\tdt=" << dt << "\n x/t \t";
    for (int j = 0; j < jmax; j++) {
        out << j * dt << "\t";
    }
    for (int i = 0; i < imax; i++) {
        out << "\n";
        out << xpocz + h * i << "\t";
        for (int j = 0; j < jmax; j++) {
            out << fabs(arg1[i][j] - anal[i][j]) << "\t";
        }
    }
    return 1;
}


int print_err22(double **anal, const char arg0[], double** arg1, int imax, int jmax, double xpocz, double h, double dt) {
    double max=-DBL_MAX;
    for (int i = 0; i < imax-1; i++) {
        for (int j = 0; j < jmax-1; j++) {
            max=fabs(arg1[i][j] - anal[i][j]);
        }
    }
    cout<<max<<", ";
}

void print_err2(double **anal, const char arg0[], int imax, int jmax, double** arg1, double xpocz, double h, double dt) {
    double max=-DBL_MAX;
    for (int i = 0; i < jmax-1; i++) {
        for (int j = 0; j < imax-1; j++) {
            if( fabs(arg1[j][i] - anal[i][j])>max){
                max=fabs(arg1[j][i] - anal[i][j]);
            }
        }
    }
    cout<<max<<", ";
}