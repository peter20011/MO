#include<iostream>
#include<cmath>

typedef double (*funkcja)(double);

double Picard(funkcja funkcja_poczatkowa, funkcja fi, funkcja pochoda_fi,double x, double max_iteracji, double toleracja_bledu, double tolerancja_residuum);
double bisekjca(funkcja funkcja_poczatkowa, double a, double b, double max_iteracji, double tolerancja_residuum, double tolerancja_bledu);
double Newton(funkcja funkcja_poczatkowa,funkcja funkcja_pochodna,double x, double max_iteracji, double tolerancja_residuum,double tolerancja_bledu);
double siecznych(funkcja funkcja_poczatkowa, double x0, double x1, double max_iteracji, double tolerancja_residuum, double tolerancja_bledu);