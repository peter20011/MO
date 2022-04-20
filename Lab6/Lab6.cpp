#include <iostream>
#include <cmath>
#include <stdio.h>
#define N 6
using namespace std;

//Alokacja i wypelnienie macierzy A
double **matrix(){
	double **a = new double*[N];
	for(int i = 0; i < N; i++)
		a[i] = new double[3];

	a[0][0] = 0.0;			a[0][1] = 10.0; 	a[0][2] = 1.0 / 2.0;
	a[1][0] = 1.0 / 3.0; 	a[1][1] = 20.0;		a[1][2] = 1.0 / 4.0;
	a[2][0] = 1.0 / 5.0;	a[2][1] = 30.0;		a[2][2] = 1.0 / 6.0;
	a[3][0] = 1.0 / 7.0;	a[3][1] = 30.0;		a[3][2] = 1.0 / 8.0;
	a[4][0] = 1.0 / 9.0;	a[4][1] = 20.0;		a[4][2] = 1.0 / 10.0;
	a[5][0] = 1.0 / 11.0;	a[5][1] = 10.0;		a[5][2] = 0.0;
	
	return a;
}


//Alokacja i wypelnienie wektora b
double *vector(){
	double *b = new double[N];
	
	b[0] = 31.0;
	b[1] = 165.0 / 4.0;
	b[2] = 917.0 / 30.0;
	b[3] = 851.0 / 28.0;
	b[4] = 3637.0 / 90.0;
	b[5] = 332.0 / 11.0;
	
	return b;
}


//Wektor  eta -> nowe wartośc na diagonalnej 
double *vector_l_eta(){
	double *id = new double[N];
	for (int i = 0; i < N; i++) 
		id[i] = 0.0;
	return id;
}


//Obliczanie ety, v - wektor przechowujacy  iloczyn l[i]*eta^-1
void calculate_eta(double **A, double *v){
	for(int i = 1; i < N; i++){
		v[i] = A[i][0]  / A[i - 1][1];
		A[i][1] -= v[i] * A[i - 1][2];
	}
}


//Obliczanie r w wektorze b , nowe wartosci 
void calculate_r(double *b, double *v){
	for(int i = 1; i < N; i++)
		b[i] -= v[i] * b[i - 1];
}


//Obliczanie rozwiazania
void calculate_x(double** A, double *b, double* x){
	x[N - 1] =  b[N - 1] / A[N - 1][1];
	for(int i = N - 2; i >= 0; i--){
		x[i] = (b[i] - A[i][2] * x[i + 1]) / A[i][1];
	}
}


//Wypisanie wektora
void result(double* x){
	cout << "x = [ ";
	for(int i = 0; i < N; i++)
		printf("%.8lf ", x[i]);
	cout << "]" << endl;
}


//Algorytm Thomasa
void thomas(double **A, double *b, double*x, double *l_eta){
	calculate_eta(A, l_eta);
	calculate_r(b, l_eta);
	calculate_x(A, b, x);
	result(x);
}


//Usuwanie z pamieci zaalokowanych elementow
void delete_all(double** a, double *b, double *x, double *l_eta){
	for(int i = 0; i < N; i++)
		delete [] a[i];
	
	delete [] a;
	delete [] b;
	delete [] x;
	delete [] l_eta;
}

int main(){
	double **A = matrix();
	double *b = vector();
	double *l_eta = vector_l_eta();
	double *x = vector_l_eta();
	
	thomas(A, b, x, l_eta);
	
	delete_all(A, b, x, l_eta);
	return 0;
}
