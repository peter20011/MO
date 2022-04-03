#include "/mnt/c/Users/User/Desktop/Zajęcia/MO/Lab3/metody.h"
#include<iostream>
#include<cmath>

using namespace std;

typedef double (*funkcja)(double);

double Picard(funkcja funkcja_poczatkowa, funkcja fi, funkcja pochoda_fi,double x, double max_iteracji, double toleracja_bledu, double tolerancja_residuum){

    cout<<"Metoda Picara"<<endl;

    if(fabs(pochoda_fi(x))>=1.0){
        cout<<"Funkcja rozbieżna"<<endl;
        return 0;
    }


    double estymator = 0.0, residuum = 0.0, przyblizenie = x;
    int iteracja = 1;
    bool kontynuuj = true;





    
}