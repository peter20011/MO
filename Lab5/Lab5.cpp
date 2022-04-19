#include <iostream>
#include <cmath>
#include "func.hpp"

using namespace std;

#define N 4 // rozmiar macierzy

//Zamiana wektora indeksowego czeswiowy wybor element podstawowego
void swap_rows(double** A, int k, int* index){
	if(A[index[k]][index[k]] == 0){
		int max_row = index[k + 1];
		double max_val = fabs(A[index[k + 1]][k]);
		
		for(int i = k + 1; i < N; i++){
			double val = A[index[i]][k];
			if(fabs(val) > max_val){
				max_row = index[i];
				max_val = val;
			}
		}
		
		swap(index[max_row], index[k]);
	}	
}


//Redukcja kolumn i uzupelnienie macierzy L
void calculate(double** A, double **L,int *index, int k){
	for(int i = k; i < N - 1; i++){
		double reducer = (A[index[i + 1]][k] / A[index[k]][k]);
		for(int j = k; j < N; j++){
			A[index[i + 1]][j] -= A[index[k]][j] * reducer;
		}
		L[index[N - (i + 1)]][k] = reducer;
	}
}


//Obliczenie Ly = b
double *calculate_y(double **M, double *b, int* index_vector) {
    double *y = allocate_vector(N);
    
    for(int i = 0; i < N; i++){
    	y[i] = 0;
	}

    for (int i = 0; i < N; i++) {

        double sum = 0;
        for (int j = 0; j < i; j++) {
            sum += y[j] * M[i][j];
        }

        y[i] = b[index_vector[i]] - sum;
    }
    
    show_vector('y', y, N);
    return y;
}


//Obliczenie Ux = y
double *calculate_x(double **M, double *y) {
    double *x = allocate_vector(N);
    
    for(int i = 0; i < N; i++)
    	x[i] = 0;

    for (int i = N - 1; i >= 0; i--) {
        double sum = 0;
        for (int j = N - 1; j > i; j--) {
            sum += x[j] * M[i][j];
        }

        x[i] = (y[i] - sum) / M[i][i];
    }

	show_vector('x', x, N);
    return x;
}


//Rozklad LU
void decomposition(double** A,double *b,  int* index){
	double** L = matrix_L(N);
	
	for(int i = 0; i < N - 1; i++){
		swap_rows(A, i, index);
		calculate(A, L, index, i);
	}
	
	double** U = matrix_U(N, A, index);
	show_matrix('U', U, index, N);
	show_matrix('L', L, index, N);
	
	double *y = calculate_y(L, b, index);
	double *x = calculate_x(U, y);
	
	check(x, N);
	delete_vector(x);
	delete_vector(y);
	delete_matrix(L, N);
	delete_matrix(U, N);
}

int main()
{
	double **A = allocate_matrix_A(N);
	double *b = allocate_vector_b(N);
	int *index = index_vector(N);
	
	show_matrix('A', A,index,  N);
	
	decomposition(A, b, index);
	
	delete_matrix(A, N);
	delete_vector(b);
	delete_vector(index);
	return 0;
}
