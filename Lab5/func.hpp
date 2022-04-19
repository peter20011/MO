#include <iostream>
#include <stdio.h>
using namespace std;

//Alokacja i uzupe�nienie macierzy A
double **allocate_matrix_A(int N){
	double **m = new double*[N];
	for(int i = 0; i < N; i++)
		m[i] = new double[N];
		
	m[0][0] = 1.0;		m[0][1] = -20.0;	m[0][2] = 30.0;		m[0][3] = -4.0;
	m[1][0] = 2.0;		m[1][1] = -40.0;	m[1][2] = -6.0;		m[1][3] = 50.0;
	m[2][0] = 9.0;		m[2][1] = -180.0;	m[2][2] = 11.0;		m[2][3] = -12.0;
	m[3][0] = -16.0;	m[3][1] = 15.0;		m[3][2] = -140.0;	m[3][3] = 13.0;

	return m;
}


//Alokacja i uzupe�nienie wektora b
double *allocate_vector_b(int N){
	double* b = new double[N];
	
	b[0] = 35.0;
	b[1] = 104.0;
	b[2] = -366.0;
	b[3] = -354.0; 
	
	return b;
}


//Alokacja i wype�nienie wektora indeks�w
int *index_vector(int N){
	int *index = new int[N];
	
	for(int i = 0; i < N; i++)
		index[i] = i;

	return index;
}


//Alokacja wektora 
double *allocate_vector(int N){
	double* vector = new double[N];
}


//Funkcje usuwaj�ce macierze i wektory
void delete_matrix(double **m, int N){
	for(int i = 0; i < N; i++)
		delete [] m[i];
	delete [] m;
}

void delete_vector(double* v){
	delete [] v;
}

void delete_vector(int *v){
	delete [] v;
}


//Funkcja tworzy i wype�nie 1 macierz L
double** matrix_L(int N){
	double **l = new double*[N];
	for(int i = 0; i < N; i++)
		l[i] = new double[N];
		
	for(int i = 0; i < N; i++)
		for(int  j = 0; j < N; j++)
			if(i == j)
				l[i][j] = 1;
			else
				l[i][j] = 0;
				
	return l;
}


//Funkcja tworzy i wype�nia kolejno macierz U
double** matrix_U(int N, double **A, int* index){
	double **u = new double*[N];
	for(int i = 0; i < N; i++)
		u[i] = new double[N];
		
	for(int i = 0; i < N; i++)
		for(int j = 0; j < N; j++)
			u[i][j] = A[index[i]][j];
				
	return u;
}


//Funkcja sprawdzaj�ca rozwi�zanie
void check(double *x, int N){
	cout << "Sprawdzenie rozwiazania" << endl;
	double **A = allocate_matrix_A(N);
	double *b = allocate_vector_b(N);
	
	for(int i = 0; i < N; i++){
		double sum = 0;
		for(int j = 0; j < N; j++)
			sum += A[i][j] * x[j];
		sum -= b[i];
		cout << sum << " ";
	}
	cout << endl;
	delete_matrix(A, N);
	delete_vector(b);
}


//Funkcje wy�wietlaj�ce macierz oraz wektor
void show_matrix(char a, double** matrix, int* id, int N){
	cout << a << endl;
	for(int i = 0; i < N; i++){
		for(int j = 0; j < N; j++)
			printf("%.8lf\t", matrix[i][j]); 
		cout << endl;
	}
	cout << endl;
}

void show_vector(char x, double *vector, int N){
	cout << x << " = [ ";
	for(int i = 0; i < N; i++)
		printf("%.8f ", vector[i]);
	cout << "]" << endl;
}
