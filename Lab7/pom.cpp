#include<iostream>
#include <cmath>
using namespace std;

int size=4;


double* multipyMatrixByVector(double **matrix, const double *vector_in){

double *vector_out=new double[size];

for(int i=0;i<size;i++){
    double sum=0.0;
    for(int j=0;j<size;j++){
        sum+=matrix[i][j]*vector_in[j];
    }
    vector_out[i]=sum;
}

    return vector_out;
}

double norm_max(double *vector){
    double biggest=fabs(vector[0]);
    for(int i=0;i<3;i++){
        if(fabs(vector[i]>biggest)){
            biggest=fabs(vector[i]);
        }
    }

    return biggest;
}

double *calculateVectorL(double **L, const double* vectorB, double omega){
    double *temp=new double[size];
    temp[0] = vectorB[0]/(L[0][0]*(1./omega));
    for(int i=1;i<size;i++){
        double sum=0;
        for(int j=0;j<i;j++){
            sum+=L[i][j]*temp[j];
        }

       temp[i] = (vectorB[i] - sum)/(L[i][i]*((1./omega)));
    }

    return temp;
}

void resetVariables(double *x_0, double **A, double *b){
    //matrixElement
    A[0][0] = 100.0; A[0][1] = -1.0; A[0][2] = 2.0; A[0][3] = -3.0;
    A[1][0] = 1.0; A[1][1] = 200.0; A[1][2] = -4.0; A[1][3] = 5.0;
    A[2][0] = -2.0; A[2][1] = 4.0; A[2][2] = 300.0; A[2][3] = -6.0;
    A[3][0] = 3.0; A[3][1] = -5.0; A[3][2] = 6.0; A[3][3] = 400.0;
    //vectorElement
    b[0] = 116.0;
    b[1] = -226.0;
    b[2] = 912.0;
    b[3] = -1174.0;
    //vectorElement
    x_0[0] = 2.0;
    x_0[1] = 2.0;
    x_0[2] = 2.0;
    x_0[3] = 2.0;
}

double **allocate_matrix_A(int N){
	double **m = new double*[N];
	for(int i = 0; i < N; i++)
		m[i] = new double[N];
		
	m[0][0] = 100.0;	m[0][1] = -1.0;	    m[0][2] = 2.0;		m[0][3] = -3.0;
	m[1][0] = 1.0;		m[1][1] = 200.0;	m[1][2] = -4.0;		m[1][3] = 5.0;
	m[2][0] = -2.0;		m[2][1] = 4.0;	    m[2][2] = 300.0;	m[2][3] = -6.0;
	m[3][0] = 3.0;	    m[3][1] = -5.0;	    m[3][2] = 6.0;	    m[3][3] = 400.0;

	return m;
}

double *allocate_vector_b(int N){
	double* b = new double[N];
	
	b[0] = 116.0;
	b[1] = -226.0;
	b[2] = 912.0;
	b[3] = -1174.0; 
	
	return b;
}

double *allocate_vector_x0(int N){
	double* c = new double[N];
	
	c[0] = 2.0;
	c[1] = 2.0;
	c[2] = 2.0;
	c[3] = 2.0; 
	
	return c;
}