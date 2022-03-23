#include<iostream>
#include<cmath>

using namespace std;

double funkcja_pierwsza( double x)
{
    double tmp=1.0; //zmienna pomocnicza

    tmp=sin(x/4.0);

    return tmp*tmp -x ;
}


double funkcja_druga(double x)
{
    return tan(2.0*x) - x -1.0;
}


double pochodna_funkcja_pierwsza(double x)
{
    return 0.25*sin(x)
}