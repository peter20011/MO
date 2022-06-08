#include <iostream>
#include <fstream>
#include <ctime>
#include <cmath>
#include <algorithm>
#include "/mnt/c/Users/User/Desktop/Zajęcia/MO/Lab11/fun.cpp"
using namespace std;

//rozwiazanie analityczne
double U_anal(double x,double t){
    return (1./(sqrt(M_PI*D*(TAU+t))))*exp(-(x*x)/(4*D*(t+TAU)));
}


//wyliczanie rozwiązań analitycznych
double **analityczna(int n, int m, double delta_t) {
    double **wyniki =tworzenie_macierzy(n, m);

    double x = XMIN;
    double t = XMAX;

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            wyniki[i][j] =  U_anal(x,t);
            x += h;
        }
        x = XMIN;
        t += delta_t;
    }
    return wyniki;
}


double **KMB(int n, int m, double delta_t) {
    double **wyniki = tworzenie_macierzy(n, m); 
    double x = XMIN;
    for (int i = 0; i < m; i++) {
        wyniki[0][i] = U_anal(x,0);  //warunek początkowy
        x+=h;
    }

    for (int i = 0; i < n; i++) {
        wyniki[i][0] = 0.0;  //warunek brzegowy
        wyniki[i][m - 1] = 0.0; //warunek brzegowy 
    }

     x = XMIN;
    for (int k = 1; k < n; k++) {
        for (int i = 1; i < m - 1; i++) {
            wyniki[k][i] = wyniki[k - 1][i]+ LAMBDA_BEZPOSREDNIE * (wyniki[k - 1][i - 1] - (2 *  wyniki[k - 1][i]) +  wyniki[k - 1][i + 1]);
            x += h;
        }
        x =XMIN;
    }
    return wyniki;
}


double **LaasonenThomas(int n, int m, double dt) {

    double *l = new double[m];    //elementy pod przekątną
    double *d = new double[m]; //elementy na przekątnej
    double *u = new double[m]; //elementy nad przekątną
    double *b = new double[m]; // wektor b
    double *wynik = new double[m]; //wektor x z niewiadomymi

    double **macierz = tworzenie_macierzy(n, m); //macierz reprezentująca siatkę czasowo-przestrzenną
    double x = XMIN;
    for (int i = 0; i < m; i++) {
        macierz[0][i] = U_anal(x,0);    //warunek początkowy
        x += h;
    }


    for (int i = 0; i < n; i++) {
        macierz[i][0] = 0.0; //warunek brzegowy
        macierz[i][m - 1] = 0.0;//warunek brzegowy 
    }


    x = XMIN;
    for (int k = 1; k < n; k++) {
        l[0] = 0.0;  
        d[0] = 1.0; 
        u[0] = 0.0; 
        b[0] = macierz[k - 1][0]; 
        for (int i = 1; i < m - 1; i++) {
            l[i] = LAMBDA_POSREDNIE; 
            d[i] = -(1.0 + 2.0 * LAMBDA_POSREDNIE );
            u[i] = LAMBDA_POSREDNIE;
            b[i] = -macierz[k - 1][i]; 
            x += h;
        }
        x = XMIN;
        l[m - 1] = 0.0; 
        d[m - 1] = 1.0; 
        u[m - 1] = 0.0; 
        b[m - 1] = macierz[k - 1][m - 1]; 

        Thomas(l, d, u, b, wynik, m); //algorytm thomasa do wyliczenia rozwiazan układu

        for (int i = 1; i < m; i++) { 
            macierz[k][i] = wynik[i];
        }
    }

    return macierz;
}

double **LaasonenSOR(int n, int m, double delta_t) {
    double *b = new double[m]; //wektor b
    double *results = new double[m];// wektro x

    double **macierz_main = tworzenie_macierzy(n, m); 
    double **macierz = tworzenie_macierzy(m, m); 

    for (int i = 0; i < m; i++) {
        for (int j = 0; j < m; j++) {
            macierz[i][j] = 0.0;  
        }
    }
    double x = XMIN;
    for (int i = 0; i < m; i++) {
        macierz_main[0][i] = U_anal(x,0); //warunek początkowy
        x+= h;
    }

    for (int i = 0; i < n; i++) {   //warunki brzegowe
        macierz_main[i][0] = 0.0;
        macierz_main[i][m - 1] = 0.0;
    }

    double t = TMIN;
    x = XMIN;
    for (int k = 1; k < n; k++) {
        macierz[0][0] = 1.0; 
        b[0] = 0.0; 
        for (int i = 1; i < m - 1; i++) {
            macierz[i][i - 1] = LAMBDA_POSREDNIE; 
            macierz[i][i] = -(1.0 + 2.0 * LAMBDA_POSREDNIE);
            macierz[i][i + 1] = LAMBDA_POSREDNIE ;
            b[i] = -macierz_main[k - 1][i]; 
            x += h;
        }
        x = XMIN;
        macierz[m - 1][m - 1] = 1.0; 
        b[m - 1] = 0; 
        t += delta_t;

        sor(macierz, b, results, m); 

        for (int i = 1; i < m; i++) {
            macierz_main[k][i] = results[i];
        }
    }

    return macierz_main;
}

//dane do wykresów dla KMB
void KMB_DANE() {
    double dt = obliczdt(LAMBDA_POSREDNIE);
    int n = (TMAX - TMIN ) / dt + 2; //liczba kroków czasowych
    int m = (XMAX - XMIN) / h; // liczba kroków przestrzennych

    double **analityczna_macierz = analityczna(n, m, dt);
    double **wynikowa_macierz = KMB(n, m, dt);
    double **bledy = tworzenie_macierzy(n, m);

    oblicz_blad(bledy, analityczna_macierz, wynikowa_macierz, n, m);

    auto *krokit = new double[n];
    double t = TMAX;
    for (int i = 0; i < n; i++) {
        krokit[i] = t;
        t += dt;
    }


    ofstream file1, file2, file3;
    file1.open("../kmb125.txt");
    file2.open("../kmb250.txt");
    file3.open("../kmb500.txt");


    double x = XMIN;
    cout << krokit[125] << endl;
    for (int i = 0; i < m; i++) {

        file1 << x << " " << analityczna_macierz[125][i + 1] << " " << wynikowa_macierz[125][i + 1] << endl; 
        x += h;
    }
    x = XMIN;
    cout << krokit[250] << endl;
    for (int i = 0; i < m; i++) {
        file2 << x << " " << analityczna_macierz[250][i + 1] << " " << wynikowa_macierz[250][i + 1] << endl;  
        x += h;
    }
    
    x = XMIN;
    cout << krokit[500] << endl;
    for (int i = 0; i < m; i++) {
        file3 << x << " " << analityczna_macierz[500][i + 1] << " " << wynikowa_macierz[500][i + 1] << endl;  
        x += h;
    }


    file1.close();
    file2.close();
    file3.close();

    double *maxBledy = macierzBlad(bledy, n, m);
    ofstream fileError;
    fileError.open("../kmbmaxError.txt");
    for (int i = 0; i < n - 1; i++) {
        fileError << krokit[i + 1] << " " << maxBledy[i + 1] << endl;
    }
    fileError.close();

    zapis_do_pliku("analytic.txt", analityczna_macierz, n, m);
    zapis_do_pliku("kmbNumerical.txt", wynikowa_macierz, n, m);
    zapis_do_pliku("kmbError.txt", bledy, n, m);


    delete[] maxBledy;
    delete[] krokit;
    usuwanie_macierzy(analityczna_macierz, n);
    usuwanie_macierzy(wynikowa_macierz, n);
    usuwanie_macierzy(bledy, n);
}

//dane do wykresów dla Laasonen- algorytm Thomasa
void LassonenThomasDANE() {
    double dt=obliczdt(LAMBDA_POSREDNIE);
    int n =(TMAX - TMIN) / dt + 2;
    int m = (XMAX - XMIN) / h;

    double **analityczna_macierz = analityczna(n, m, dt);
    double **wynikowa_macierz = LaasonenThomas(n, m, dt);
    double **bledy = tworzenie_macierzy(n, m);

     oblicz_blad(bledy, analityczna_macierz, wynikowa_macierz, n, m);

    auto *krokit = new double[n];
    double t = TMIN;
    for (int i = 0; i < n; i++) {
        krokit[i] = t;
        t += dt;
    }


    ofstream lt50, lt100, lt200;
    lt50.open("../lt50.txt");
    lt100.open("../lt100.txt");
    lt200.open("../lt200.txt");


    double x = XMIN;
    cout << krokit[50] << endl;
    for (int i = 0; i < m; i++) {
        lt50 << x << " " << analityczna_macierz[50][i + 1] << " " << wynikowa_macierz[50][i + 1] << endl;  
        x += h;
    }
    x = XMIN;
    cout << krokit[100] << endl;
    for (int i = 0; i < m; i++) {
        lt100 << x << " " << analityczna_macierz[100][i + 1] << " " << wynikowa_macierz[100][i + 1] << endl;  
        x += h;
    }
    x = XMIN;
    cout << krokit[200] << endl;
    for (int i = 0; i < m; i++) {
        lt200 << x << " " << analityczna_macierz[200][i + 1] << " " << wynikowa_macierz[200][i + 1] << endl; 
        x += h;
    }


    lt50.close();
    lt100.close();
    lt200.close();

    double *maxBledy = macierzBlad(bledy, n, m);
    ofstream maxerror;
    maxerror.open("../LassTmaxError.txt");
    for (int i = 0; i < n - 1; i++) {
        maxerror << krokit[i + 1] << " " << log10(maxBledy[i + 1]) << endl;
    }
    maxerror.close();
    zapis_do_pliku("LassonenTNumerical.txt", wynikowa_macierz, n, m);
    zapis_do_pliku("analytic.txt", analityczna_macierz, n, m);
    zapis_do_pliku("LassonenTError.txt", bledy, n, m);


    delete[] maxBledy;
    delete[] krokit;
    usuwanie_macierzy(analityczna_macierz, n);
    usuwanie_macierzy(wynikowa_macierz, n);
    usuwanie_macierzy(bledy, n);


}

//dane do wykresów Lassonen - metoda SOR
void LaasonenSORDANE() {
    double dt=obliczdt(LAMBDA_POSREDNIE);
    int n =(TMAX - TMIN) / dt + 2;
    int m = (XMAX - XMIN) / h;

    double **analityczna_macierz = analityczna(n, m, dt);
    double **wynikowa_macierz = LaasonenSOR(n, m, dt);
    double **bledy = tworzenie_macierzy(n, m);

   oblicz_blad(bledy, analityczna_macierz, wynikowa_macierz, n, m);

    auto *krokit = new double[n];
    double t = TMIN;
    for (int i = 0; i < n; i++) {
        krokit[i] = t;
        t += dt;
    }

    ofstream sor50, sor100, sor200;
    sor50.open("../sro50.txt");
    sor100.open("../sor100.txt");
    sor200.open("../sor200.txt");


    double x = XMIN;
    cout << krokit[50] << endl;
    for (int i = 0; i < m; i++) {
        sor50 << x << " " << analityczna_macierz[50][i + 1] << " " << wynikowa_macierz[50][i + 1] << endl;  // t = 0.3125;
        x += h;
    }
    x = XMIN;
    cout << krokit[100] << endl;
    for (int i = 0; i < m; i++) {
        sor100 << x << " " << analityczna_macierz[100][i + 1] << " " << wynikowa_macierz[100][i + 1] << endl;  // t = 0.75;
        x += h;
    }
    x = XMIN;
    cout << krokit[200] << endl;
    for (int i = 0; i < m; i++) {
        sor200 << x << " " << analityczna_macierz[200][i + 1] << " " << wynikowa_macierz[200][i + 1] << endl;  // t = 1.25;
        x += h;
    }


    sor50.close();
    sor100.close();
    sor200.close();



    double *maxBledy = macierzBlad(bledy, n, m);
    ofstream maxerror;
    maxerror.open("../LassSORmaxError.txt");
    for (int i = 0; i < n - 1; i++) {
        maxerror << krokit[i + 1] << " " << log10(maxBledy[i + 1]) << endl;
    }
    maxerror.close();

    zapis_do_pliku("analytic.txt", analityczna_macierz, n, m);
    zapis_do_pliku("LassonenSORNumerical.txt", wynikowa_macierz, n, m);
    zapis_do_pliku("LassonenSORError.txt", bledy, n, m);


    delete[] maxBledy;
    delete[] krokit;
    usuwanie_macierzy(analityczna_macierz, n);
    usuwanie_macierzy(wynikowa_macierz, n);
    usuwanie_macierzy(bledy, n);
}


int main() {
    time_t begin, end;
    time(&begin);
    KMB_DANE();
    LassonenThomasDANE();
    LaasonenSORDANE();
    time(&end);
    time_t czas = end - begin;
    printf("Czas: %ld sek.", czas);
    return 0;
}