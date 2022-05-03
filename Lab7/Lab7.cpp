#include<iostream>
#include<cmath>
#include "/mnt/c/Users/User/Desktop/Zajęcia/MO/Lab7/pom.cpp"
using namespace std;

int ITERACJE=50;
double TOLX=1e-10;
double  TOLF=1e-10;

void jacobi (double**A, double *x_0, const double *b, double *x_n){
    double *e=new  double[size];
    double *residual =new double[size];

    for(int n=0;n<ITERACJE;n++){
        double *vectorOut=new double[size];

        for(int i=0;i<size;i++){
            double sum=0.0;
            for(int j=0;j<size;j++){
                double matrixElement;
                if(i==j){
                    matrixElement=0;
                }
                else{
                    matrixElement=A[i][j];
                }
                sum+=-(1./A[i][j])*matrixElement*x_0[j];
            }

            vectorOut[i]=sum;
        }

        for(int z=0;z<size;z++){
            x_n[z] = vectorOut[z] + (1./A[z][z])*b[z];
            e[z] = x_n[z] - x_0[z];
            residual[z] = fabs(A[z][0]*x_n[0] + A[z][1]*x_n[1] + A[z][2]*x_n[2] + A[z][3]*x_n[3] - b[z]);
            x_0[z] = x_n[z];
        }
        cout<<"x_n: ";
        for(int i=0;i<size;i++){
            cout<<x_n[i]<<" ";
        }
        cout<<endl;

        cout << "e = " << norm_max(e) << "\t\t|f(x)| = " << norm_max(residual);
        cout<<endl;

        if(norm_max(e)<=TOLF && norm_max(residual)<=TOLF){
            cout<<"Zakończono na podstawie kreyterium numer 2 i 3"<<endl;
            return;
        }

        cout<<endl;
        cout<<endl;
    }
}


void gauss_seidel(double **A, double *x_0, const double *b, double *x_n){
    double *e=new  double[size];
    double *residual =new double[size];

    for(int n=0;n<ITERACJE;n++){
        double *vectorUxn=new double[size];

        for(int i=0;i<size;i++){
            double sum=0.0;
            for(int j=0;j<size;j++){
                double U_element;
                if(i<j){
                    U_element=A[i][j];
                }
                else{
                    U_element=0;
                }
                sum+=U_element*x_0[j];
            }
            vectorUxn[i]=sum;
        }

        double *y=new double[size];
        for(int z=0;z<size;z++){
            y[z]=-vectorUxn[z]+b[z];
        }
        x_n=calculateVectorL(A,y,1);

        for(int l=0;l<size;l++){
            e[l]=x_n[l]-x_n[l];
            residual[l]= fabs(A[l][0]*x_n[0] + A[l][1]*x_n[1] + A[l][2]*x_n[2] + A[l][3]*x_n[3] - b[l]);
            x_0[l]= x_n[l];
        }

        cout<<"x_n: ";
        for(int i=0;i<size;i++){
            cout<<x_n[i]<<" ";
        }
        cout<<endl;
        cout << "e = " << norm_max(e) << "\t\t|f(x)| = " << norm_max(residual);
        cout<<endl;

        if(norm_max(e)<=TOLX && norm_max(residual)<=TOLF){
            cout<<"Zakończono na podstawie kreyterium numer 2 i 3"<<endl;
            return;
        }

        cout<<endl;
        cout<<endl;
    }
} 

void SOR(double **A, double *x_0, const double *b, double *x_n) {
    double omega = 1./2.;
    double *e = new double[size];
    double *residual = new double[size];

    for(int n = 0; n < ITERACJE; n++) {

        double *vectorDUxn = new double[size];

        for(int i = 0; i < size; i++) {
            double sum = 0.0;
            for(int j = 0; j < size; j++) {
                double U_Element = 0;
                if(i < j) {
                    U_Element = A[i][j];
                } else if (i == j) {
                    U_Element = (1 - 1./omega) * A[i][j];
                }
                sum += U_Element * x_0[j];
            }
            vectorDUxn[i] = sum;
        }

        double *y = new double[size];

        for(int l = 0; l < size; l++) {
            y[l] = -vectorDUxn[l] + b[l];
        }

        x_n = calculateVectorL(A, y, omega);

        for(int l = 0; l < size; l++) {
            e[l] = x_n[l] - x_0[l];
            residual[l] = fabs(A[l][0]*x_n[0] + A[l][1]*x_n[1] + A[l][2]*x_n[2] + A[l][3]*x_n[3] - b[l]);
            x_0[l] = x_n[l];
        }

        std::cout << "x_n: ";
        for(int i = 0; i < size; i++) {
            std::cout << x_n[i] << " ";
        }
        std::cout << std::endl;

        cout << "e = " << norm_max(e) << "\t\t|f(x)| = " << norm_max(residual);
        cout << endl;

        if(norm_max(e) <= TOLX && norm_max(residual) <= TOLF) {
             cout<<"Zakończono na podstawie kreyterium numer 2 i 3"<<endl;
            return;
        }

        cout << endl;
        cout << endl;
    }
}

int main(){
    double **A=allocate_matrix_A(size);
    double *b=allocate_vector_b(size);
    double *x_0=allocate_vector_x0(size);

    double *x_n1=new double [size];

    cout<<"Jacobi: "<<endl;
    jacobi(A,x_0,b,x_n1);
    resetVariables(x_0,A,b);

    double *x_n2=new double [size];
    cout<<"Gauss-Seidel: "<<endl;
    gauss_seidel(A,x_0,b,x_n2);
    resetVariables(x_0,A,b);

    double *x_n3=new double[size];
    cout<<"SOR: "<<endl;
    SOR(A,x_0,b,x_n3);
    cout<<endl;

    return 0;

}