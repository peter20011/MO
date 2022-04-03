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

    cout << " ITERACJA - WARTOSC PRZYBLIZONA - RESIDUUM - ESTYMATOR BLEDU" << endl;

    while(kontynuuj){
        cout<< iteracja <<" - ";

        przyblizenie=fi(przyblizenie);
        cout << przyblizenie << " - ";

        estymator=fabs(przyblizenie-x);

        x=przyblizenie;

        residuum=funkcja_poczatkowa(x);
        cout << residuum << " - " << estymator << endl;

        if((fabs(residuum)<=tolerancja_residuum) && (estymator <= tolerancja_bledu) || (iteracja)>= max_iteracji){
            kontynuuj=false;
        }

        iteracja+=1;

    
    }

        return 0;

}


double bisekjca(funkcja funkcja_poczatkowa, double a, doubleb b, double max_iteracji, double tolerancja_residuum, double tolerancja_bledu)
{
    cout << "METODA BISEKCJI" << endl;

    if((funkcja_poczatkowa(a)>0 && funkcja_poczatkowa(b)>0) || (funkcja_poczatkowa(a)<0 && funkcja_poczatkowa(b)<0)){
        cout<<"Funkcja nie zmienia znaku dla danego przedzialu"<<endl;
        return 0;
    }

    double estymator = 0.0, residuum = 0.0, x = 0.0;
    int iteracja = 1;
    bool kontynuuj = true;

    cout << "ITERACJA - A - B - X- RESIDUUM - ESTYMATOR BLEDU" << endl;

    while(kontynuuj){
        x=(a+b)/2.0;
        estymator=fabs((b-a)/2.0)
        residuum=funkcja_poczatkowa(x);

        cout << iteracja << " - " << a << " - " << b << " - " << x << " - " << residuum << " - " << estymator << endl;
        
        if((funkcja_poczatkowa(a)>0 && funkcja_poczatkowa(b)>0) || (funkcja_poczatkowa(a)>0 && funkcja_poczatkowa(b)<0)){
            b=x;
        }
        else{
            a=x;
        }

        if((fabs(residuum)<=tolerancja_residuum) && (estymator <=tolerancja_residuum) || (iteracja >=max_iteracji)){
            kontynuuj=false;
        }

        iteracja+=1;
    }

        return 0;
}

double Newton(funkcja funkcja_poczatkowa,funkcja funkcja_pochodna,double x, double max_iteracji, double tolerancja_residuum,double tolerancja_bledu){

   cout << "Metoda Newtona" << endl;

    double x0 = x, x1 = 0.0, estymator = 0.0, residuum = 0.0;
    int iteracja = 1;
    bool kontynuuj = true; 

    while(kontynuuj){
        x1=x0-(funkcja_poczatkowa(x0)/funkcja_pochodna(x0));
        estymator=fabs(x0-x1);
        x0=x1;
        residuum=funkcja_poczatkowa(x0);

        cout << iteracja << " - " << x1 << " - " << residuum << " - " << estymator << endl;

        if((fabs(residuum)<=tolerancja_residuum) && (estymator <= tolerancja_bledu) || (iteracja >=max_iteracji)){
            kontynuuj=true;
        }

        iter_swap+=1;
    }

    return 0;
}

double siecznych(funkcja funkcja_poczatkowa, double x0, double x1, double max_iteracji, double tolerancja_residuum, double tolerancja_bledu){

    cout << "METODA SIECZNYCH" << endl;

    double x2 = 0.0, estymator = 0.0, residuum = 0.0;
    int iteracja = 0;
    bool kontynuuj = true;

    cout << "Iteracja - x1 - residuum - estymator bledu" << endl;

    while(kontynuuj){
        x2=x1-funkcja_poczatkowa(x1)/((funkcja_poczatkowa(x1)-funkcja_poczatkowa(x0))/(x1-x0));
        estymator=fabs(x2-x1);
        residuum=funkcja_poczatkowa(x2);
    
        cout << iteracja << " - " << x1 << " - " << residuum << " - " << estymator << endl;

        x0 = x1;
        x1 = x2;

        if((fabs(residuum)<=tolerancja_residuum) && (estymator<= tolerancja_bledu) || (iteracja>=max_iteracji)){
            kontynuuj=true;
        }
        iteracja+=1;
    }

    return
}