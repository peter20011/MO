#include<iostream>
#include<cmath>

using namespace std;

double funkcja_pierwsza(double x) {
    return pow(sin(x / 4.0), 2) - x;
}

double pochodna_funkcja_pierwsza(double x) {
    return (1.0 / 4.0) * (sin(x / 2.0)) - 1.0;
}

double fi_pierwsza(double x) {
    return pow(sin(x / 4.0), 2);
}

double pochodna_fi_pierwsza(double x) {
    return (1.0 / 4.0) * sin(x / 2.0);
}

double funkcja_druga(double x) {
    return tan(2.0 * x) - x - 1.0;
}

double pochodna_funkcja_druga(double x) {
    return -1.0 + 2.0 / pow(cos(2.0 * x), 2);
}

double fi_druga(double x) {
    return tan(2.0 * x) - 1.0;
}

double pochodna_fi_druga(double x) {
    return 2.0 / pow(cos(2.0 * x), 2);
}