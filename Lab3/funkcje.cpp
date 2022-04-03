#include "/mnt/c/Users/User/Desktop/Zajęcia/MO/Lab3/funkcje.h"
#include<iostream>
#include<cmath>

using namespace std;

double funkcja_pierwsza( double x)
{
    double tmp=1.0; //zmienna pomocnicza

    tmp=sin((double)x/4.0);

    return tmp*tmp -x ;
}

double pochodna_funkcja_pierwsza(double x)
{
    return 0.25*sin(x/2.0) -1.0;
}

double fi_z_pierwszej(double x){
    double tmp=sin(x/4.0);
    return tmp*tmp;
}

double pochodna_fi_z_pierwszej(double x){
    return 0.25*sin(x/2.0);
}


double funkcja_druga(double x)
{
    return tan(2.0*x) - x -1.0;
}

double pochodna_funkcja_druga(double x){
    return -1+2.0/(cos(2.0*x)*cos(2.0*x)) ;
}

double fi_z_drugiej(double x){
    return tan(2.0*x)-1.0;
}

double pochodna_fi_z_drugiej(double x){
     return 2.0 /(cos(2.0 * x)*cos(2.0 * x));
}