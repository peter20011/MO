#define _USE_MATH_DEFINES 
#include <iostream>
#include <cmath>
#include <fstream>
#include <time.h>
#define D 1
#define TETA 0.1
#define TMAX 2
#define ITERACJE 100
#define EPSILON 1e-16
using namespace std;

//rozwiąznie analityczne pomijamy D
double U_anal(double x ,double t){
    return exp(-(x*x)/(4*D*(t+TETA)))/(2.0*pow(M_PI*D*(t+TETA),(double)1/2));
}
//warunek poczatkowy
double inline initial_condition(double x){
    return U_anal(x,0);
}

// algorytm Thomasa

void Thomas_dist(double *l, double *d, double *u,int n){
        for(int i=1;i<n;i++){
        d[i]=d[i]-(u[i-1]*l[i])/d[i-1];
    }
}

double* Thomas_solve( double *l, double *d, double *u, double *b, int n){
    double *x=new double [n];
    

    for(int i=1;i<n;i++){
        b[i]=b[i]-(l[i]*b[i-1])/d[i-1];
    }

    x[n-1]=b[n - 1] /d[n - 1];

    for(int i=n-2;i>=0;i--){
        x[i]=(b[i]-u[i]*x[i+1])/d[i];
    }

    return x;
}

double** anal_solution(double x_min, double x_max, double tmax, double h, double dt){
    int count_x=(int)((x_max-x_max)/h),count_t=(int)(tmax/dt),i,k;

    double **u=new double*[count_x];
    for(int i=0;i<count_x; i++){
        u[i]=new double[count_t];
    }

    for(k=0;k<count_t;k++){
        u[0][k] = 0;
        u[0][count_x - 1] = 0;
        for (i = 1; i < count_x - 1; i++) {
            u[i][k] = U_anal(x_min + i*h, k * dt); //wyliczamy konkretne wartośći
        }
    }
    
    //dopisać do pliku !!!!!!!!!!!
    return u;
}
//klasyczna metoda bezpośrednia
double KMB(double x_min, double x_max, double tmax, double h, double dt){
      int count_x = (int) ((x_max - x_min) / h),count_t = (int) (tmax / dt),ex_time = clock(),i, k;
    double **u = new double*[count_x],
            lambda = D * dt / (h * h);
    for (i = 0; i < count_x; i++) {
        u[i] = new double[count_t];
    }
    u[0][0] = 0; //warunki brzegowe
    u[0][count_x - 1] = 0;
    for (i = 1; i < count_x - 1; i++) {
        u[i][0] = initial_condition(x_min + i * h);
    }

    for (k = 0; k < count_t - 1; k++) {
        u[0][k] = 0; //warunki brzegowe
        u[count_x - 1][k] = 0;

        for (i = 1; i < count_x - 1; i++) {
            u[i][k + 1] = u[i][k] + lambda * (u[i - 1][k] - 2 * u[i][k] + u[i + 1][k]); //obliczanie nowych warosci (metoda bezposrednia)
        }
    }
    //zapis do pliku
    double **tmp=anal_solution(x_min,x_min,tmax,h,dt);
    //zapis do pliku

    for(i=0;i<count_x; i++){
        delete u[i],tmp[i];
    }
    delete u,tmp;

    ex_time= clock()-ex_time;

    return ex_time;
}

double Laasonen_Tom(double x_min, double x_max, double tmax, double h, double dt){
        int count_x = (int) ((x_max - x_min) / h),count_t = (int) (tmax / dt),ex_time = clock(),i;
    double **result = new double*[tcount],
            *l = new double[count_x - 1],
            *d = new double[count_x],
            *b = new double[count_x],
            *u = new double[count_x - 1],
            lambda = D * dt / (h * h);
    result[0] = new double[count_x];
    for (i = 0; i < count_x; i++) {
        result[0][i] = initial_condition(x_min + i * h);
    }

    double lambda2 = -(1 + 2 * lambda);
    for (i = 1; i < count_x - 1; i++) { //zapełnienie ldu
        u[i] = lambda;
        l[i] = lambda;
        d[i] = lambda2;
    }
    u[0] = 0;
    l[0] = lambda;
    l[count_x - 2] = 0;
    d[0] = 1;
    d[count_x - 1] = 1;

    Thomas_dist(l, d, u, count_x); //1-wszy etap alg Thomasa
    for (i = 1; i < tcount; i++) {
        for (int j = 1; j < count_x - 1; j++) {
            b[j] = -result[i - 1][j];
        }
        b[0] = 0; //warunki brzegowe
        b[count_x - 1] = 0;

        result[i] = Thomas_solve(l, d, u, b, count_x); //wyliczenie nowych wartosci - metoda posrednia
    }

    delete l, d, u, b;

    //zapis do pliku
    double **tmp = anal_solution(x_min, x_max, tmax, h, dt);
    //zapis do pliku 

    for (i = 0; i < tcount; i++) {
        delete result[i], tmp[i];
    }
    delete result, tmp;
    ex_time = clock() - ex_time;
    return ex_time;
}

void SOR(double **A, double *x_0, const double *b, double *x_n) {
    double omega = 1./2.;
    double *e = new double[size];
    double *residual = new double[size];

    for(int n = 0; n < ITERACJE; n++) { //tworzymy pętle główną, której ograniczeniem jest liczba iteracji	

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
            residual[l] = fabs(A[l][0]*x_n[0] + A[l][1]*x_n[1] + A[l][2]*x_n[2] + A[l][3]*x_n[3] - b[l]); // przemnazamy wiersze razy wektor xx i odejmujemy wektor b
            x_0[l] = x_n[l];
        }

        cout << "x_n: ";
        for(int i = 0; i < size; i++) {
            cout << x_n[i] << " ";
        }
        cout << endl;

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

double norm_max(double *vector){
    double biggest=fabs(vector[0]);
    for(int i=0;i<3;i++){
        if(fabs(vector[i]>biggest)){
            biggest=fabs(vector[i]);
        }
    }

    return biggest;
}

double Laasonen_SOR(double x_min, double x_max, double tmax, double h, double dt){
    int count_x = (int) ((x_max - x_min) / h),tcount = (int) (tmax / dt),ex_time = clock(),i;
}


int main(){
    double x_max=6*pow(D*(TETA+TMAX),(double)1/2);
    double h=1.0, dt=0.4*h*h/D; //lambda 0.4-> KMB jest wtedy stabilna numerycznie 
    cout << "Rozwiazywanie rownania rozniczkowego czastkowego (typu Dyfuzji)" << endl;
    cout<<"KMB"<<endl;
    cout<<"h\terr"<<endl;
    while (h>0.009)
    {
        cout << h<<"\t";
        KMB(-x_min,x_max,TMAX,h,dt);
        cout<<endl;
        h=h*0.7;
        dt=0.4*h*h/D;
    }

    h=1.0;
    dt=h*h/D; // lambda 1 dla metody pośredniej Lassonen
    cout<<"Lassonen+Thomas"<<endl;
    while(h>0.009){
        cout<<h<<"\t";
        Laasonen_Tom(-xmax, xmax, TMAX, h, dt);
        cout<<endl;
        h=h*0.7;
        dt = h * h / D;
    }

    return 0;

}