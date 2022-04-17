#include <iostream>
#include <math.h>


using namespace std;

//kreyteria zakończenia iteracji
double TOL_F = 1e-16;
double TOL_X = 1e-16;
double MAX_ITERACJE = 30;


double max(double a, double b, double c){
    double najwieksza=a;

    if(b>najwieksza){
        najwieksza=b;
    }

     if(c>najwieksza){
        najwieksza=c;
    }

    return najwieksza;
}

//wzory z zad pierwszego

double funkcja_pierwsza(double x, double y, double z){
    return x*x+y*y+z*z-2.0;
}
double funkcja_druga(double x, double y,double z){
    return x*x+y*y-1;
}

double funkcja_trzecia(double x, double y,double z){
    return x*x-y;
}


double oblicz_a(double x, double y, double z) {
    return (pow(x, 2) - pow(y, 2) - 1.0 + 2.0 * pow(x, 2) * y) / (2.0 * x * (1.0 + 2.0 * y));
}

double oblicz_b(double x, double y, double z) {
    return (pow(y, 2) + y - 1.0) / (1.0 + 2.0 * y);
}

double oblicz_c(double x, double y, double z) {
    return (pow(z, 2) + 2.0 * pow(z, 2) * y - 2.0 * y - 1.0) / (2.0 * z * (1.0 + 2.0 * y));
}


bool czy_nalezy_do_dziedziny(double x, double y, double z){
    if (x == 0) {
        cout << "Nieprawidlowy x" << endl;
        return false;
    } else if (z == 0) {
        cout << "Nieprawidlowy z" << endl;
        return false;
    } else if (y == -1 / 2 || 4 * y * x + 2 * x == 0) {
        cout << "Nieprawidlowy x lub y" << endl;
        return false;
    } else
        return true;
}