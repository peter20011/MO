#include <iostream>
#include <cmath>
#include <iomanip>

using namespace std;

double **nowa_macierz(int rozmiar){
    double **A=new double *[rozmiar];
    for(int i=0;i<rozmiar;i++){
        A[i]=new double[rozmiar];
    }
    return A;
}

void zapisz_do_macierzy(double **A, int wymiar){
    double dane_z_zadana[4][4]=
    {{1.0, -20.0, 30.0, -4.0},
         {2.0, -40.0, -6.0, 50.0},
         {9.0, -180.0, 11.0, -12.0},
         {-16.0, 15.0, -140.0, 13.0}};

    for (int i = 0; i < wymiar; i++)
        for (int j = 0; j < wymiar; j++)
            A[i][j] = dane_z_zadana[i][j];
}

void wypisz_maciarz(double **A, int rozmiar, int *indeks){
        for (int i = 0; i < rozmiar; i++) {
            for (int j = 0; j < rozmiar; j++)
                cout << setw(8) << A[indeks[i]][j] << "\t";
        cout << "\n";
}
}

void usun_macierz(double **A, int rozmiar) {
    for (int i = 0; i < rozmiar; i++)
        delete[] A[i];

    delete[] A;
}


void wypisz_wektor(double *wek, int rozmiar){
    for (int i = 0; i < rozmiar; i++){
        cout << wek[i] << endl;
    }
}


void wypisz_U(double **A, int rozmiar, int *indeks){
for(int i=0;i<rozmiar;++i){
    for(int j=0;j< rozmiar;++j){
        if(j>=i){
            cout<<setprecision(4)<<A[indeks[i]][j]<<"\t";
        }
        else{
            cout<<"0\t";
        }
    }

    cout<<"\n";
} 
}

void wypisz_L(double **A, int rozmiar, int *indeks){
    for (int i = 0; i < rozmiar; i++) {
            for (int j = 0; j < rozmiar; j++){
                if(j>i){
                    cout<<"0\t";
                }
                else{
                    if(j==i){
                        cout<<"1\t";
                    }
                    else{
                    cout<<setprecision(4)<<A[indeks[i]][j]<<"\t";
                    }
                }
            }
        cout << "\n";
}
}
