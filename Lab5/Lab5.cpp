#include <iostream>
#include <cmath>
#include <iomanip>
#include "/mnt/c/Users/User/Desktop/Zajęcia/MO/Lab5/funkcje.cpp"
using namespace std;

// funkcje główne do zadanaia 
 int czescwioy_wybor_elementu_podstawowego(double **A, int rozmiar, int pozycja, int *indeks){
         int numer_wiersza;
    // Szuka wartosci najwiekszej w kolumnie ponizej komorki [j][j]
    for (int i = pozycja; i < rozmiar; i++) {
        if (fabs(A[indeks[i]][pozycja]) < fabs(A[indeks[i + 1]][pozycja]))
            numer_wiersza = indeks[i + 1];
        else
            numer_wiersza = indeks[i];
    }

    return numer_wiersza;
 }

void Gauss(double **A, int *indeks){
    int wiersz;
    double v;   // Element ponizej diagonali
    // poruszanie sie po diagonali
    for (int k = 0; k < 3; k++) {
        if (A[indeks[k]][indeks[k]] == 0.0) {
            wiersz = czescwioy_wybor_elementu_podstawowego(A, 3, indeks[k], indeks);
            // Zapisanie zmian w tablicy indeksow
            indeks[wiersz] = indeks[k];
            indeks[k] = wiersz;
        }
        // Zejscie w dol po kolumnie
        for (int i = k + 1; i < 4; i++) {
            v = A[indeks[i]][k];
            // Poruszanie sie w prawo po wierszu
            for (int j = k + 1; j < 4; j++) {
                // Obliczanie wartosci dla macierzy U
                A[indeks[i]][j] = A[indeks[i]][j] - A[indeks[k]][j] * (v / A[indeks[k]][k]);
                // Zapisanie macierzy L
                // wspolczynniki przez ktore mnozymy kolejne wiersze (A[i][k] / A[k][k])
                A[indeks[i]][k] = v / A[indeks[k]][k];
            }
        }
    }
}


void wyznaczY(double **macierzL, double *wektorB, int *index, int rozmiar) {
    double suma = 0.0;
    // Poruszanie sie po kolumnie, zaczyna od lewego gornego
    for (int i = 0; i <= rozmiar; i++) {
        // Poruszanie sie wierszu w prawo
        for (int j = 0; j < i; j++)
            suma += macierzL[index[i]][j] * wektorB[index[j]];

        wektorB[index[i]] = (wektorB[index[i]] - suma) / 1.0;     //Na glownej przekatnej sa 1

        suma = 0.0;
    }
}

void wyznaczX(double **macierzU, double *wektorB, int *index, int rozmiar) {
    double suma = 0.0;
    // Zaczyna od prawego dolnego rogu
    for (int i = rozmiar; i >= 0; i--) {
        for (int j = i + 1; j <= rozmiar; j++)
            suma += macierzU[index[i]][j] * wektorB[index[j]];

        wektorB[index[i]] = (wektorB[index[i]] - suma) / (macierzU[index[i]][i]);     //UWAGA MOŻLIWE DZIELENIE PRZEZ 0

        suma = 0.0;
    }
}

int main() {
    double **M = nowa_macierz(4);
    int indeks[4] = {0, 1, 2, 3};

    double b[4] = {35.0, 104.0, -366.0, -354.0};
    zapisz_do_macierzy(M, 4);

    cout << "Macierz A: " << endl;
    wypisz_maciarz(M, 4, indeks);

    Gauss(M, indeks);

    cout<<endl;
    cout<<"Macierz L:"<<endl;
    wypisz_L(M,4, indeks);
    cout<<"Macierz U:"<<endl;
    wypisz_U(M,4, indeks);
    cout<<endl;
 
    cout << "\nRoziwiazanie ukladu rownan Ax = b: " << endl;
    wyznaczY(M, b, indeks, 4 - 1);
    cout << "\nWektor y: " << endl;
    wypisz_wektor(b, 4);
    cout << endl;

    wyznaczX(M, b, indeks, 4 - 1);
    cout << "\nWektor x (Ux = y): " << endl;
    wypisz_wektor(b, 4);
    usun_macierz(M, 4);
    return 0;
}