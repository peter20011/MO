#define _USE_MATH_DEFINES 
#include <iostream>
#include <cmath>
#include <fstream>
#include <time.h>
#include<cmath>
#include "/mnt/c/Users/User/Desktop/Zajęcia/MO/Lab11/pom.cpp"
#define D 1
#define TETA 0.1
#define TMAX 2
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
    print_err22(tmp, "Klasyczna bezposrednia",u,count_x,count_t,x_min,h, dt);

    for(i=0;i<count_x; i++){
        delete u[i],tmp[i];
    }
    delete u,tmp;

    ex_time= clock()-ex_time;

    return ex_time;
}

double Laasonen_Tom(double x_min, double x_max, double tmax, double h, double dt){
    int count_x = (int) ((x_max - x_min) / h);
    int count_t = (int) (tmax / dt),ex_time = clock(),i;
    double **result = new double*[count_t],
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
    for (i = 1; i < count_t; i++) {
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

    for (i = 0; i < count_t; i++) {
        delete result[i], tmp[i];
    }
    delete result, tmp;
    ex_time = clock() - ex_time;
    return ex_time;
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
        KMB(-x_max,x_max,TMAX,h,dt);
        cout<<endl;
        h=h*0.7;
        dt=0.4*h*h/D;
    }

    h=1.0;
    dt=h*h/D; // lambda 1 dla metody pośredniej Lassonen
    cout<<"Lassonen+Thomas"<<endl;
    while(h>0.009){
        cout<<h<<"\t";
        Laasonen_Tom(-x_max, x_max, TMAX, h, dt);
        cout<<endl;
        h=h*0.7;
        dt = h * h / D;
    }

    return 0;

}