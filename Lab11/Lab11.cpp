#include<iostream>
#include<fstream>
#include<cmath>
#include<string>
#include "/mnt/c/Users/User/Desktop/Zajęcia/MO/Lab11/fun.cpp"
#define TOLX 10e-16
#define TOLF 10e-16
using namespace std;

const double D=1.0;
const double TAU=0.1;
const double TMAX=2.0;
const double TMIN=0.0;
const double A=6*pow(D*(TAU+TMAX),1./2)+0.01;  // wychodzi 8,6948 , iborę trochę większe dodaję 0.01
const double LAMBDA_POSREDNIE=1.0;
const double LAMBDA_BEZPOSREDNIE=0.4;
const double XMIN=-A;
const double XMAX=A;

double U_anal(double x, double t){
    return exp((-(x*x))/(4.0*D*(TAU+t)))*1./(2*pow(M_PI*D*(TAU+t),1./2));
}

double warunek_poczatkowy(double x){
    return U_anal(x,0);
}


void warunek_brzegowy(double** macierz, int n, int m){
    for(int i=0;i<n;i++){
        macierz[i][0]=0.0;
        macierz[i][m-1]=0.0;
    }
}

void warunek_poczatkowy(double **macierz, int n, int m){

}