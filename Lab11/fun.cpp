#include<iostream>
#include<fstream>
#include<cmath>
#include<string>
const double D=1.0;
const double TAU=0.1;
const double TMAX=2.0;
const double TMIN=0.0;
const double A=6*sqrt(D*(TAU+TMAX))+0.1; // wychodzi 8,6948 , biorę trochę większe dodaję 0.1
const double LAMBDA_POSREDNIE=1.0;
const double LAMBDA_BEZPOSREDNIE=0.4;
const double XMIN=-A;
const double XMAX=A;
const double H=0.1;
const double OMEGA=0.9;
const int ITERACJE=100;

using namespace std;

double *utworz_wektor(int n){
    return new double[n];
}

double** utworz_macierz(int n, int m){
    double **macierz=new double*[n];

    for(int i=0;i<n; i++){
        macierz[i]=new double[m];
    }

    return macierz;
}

void usun_wektor(double *wektor){
    delete [] wektor;
}

void usun_macierz(double** macierz, int n){
    for(int i=0;i<n; i++)
    {
        delete[] macierz[i];
    }

    delete[] macierz;
}

void zapiszWektor(double *wektor, int n, string nazwapliku) {
  fstream file(nazwapliku.c_str(), ios::out);

  if (file.is_open()) {
    for (int i = 0; i < n; i++) {
      file << wektor[i] << endl;
    }
  }
}

void zapiszWektor2(double *wektor_1,double* wektor_2 ,int n,string nazwapliku){
   fstream file(nazwapliku.c_str(), ios::out);

  if (file.is_open()) {
    for (int i = 0; i < n; i++) {
      file << wektor_1[i] << "\t" << wektor_2[i] << endl;
    }
  }
}

void zapiszMacierz(double **matrix, int n, int m, string nazwapliku) {
  fstream file(nazwapliku.c_str(), ios::out);

  if (file.is_open()) {
    for (int i = 0; i < n; ++i) {
      for (int j = 0; j < m; ++j) {
        file << matrix[i][j] << ";";
      }
      file << "\n";
    }
  }

  file.close();
}

double norma_max(double *wektor, int n) {
  double max = fabs(wektor[0]);
  for (int i = 1; i < n; i++) {
    if (max < fabs(wektor[i])) {
      max = wektor[i];
    }
  }
  return max;
}

double **obliczBlad(double **dokladne, double **przyblizenie, int n, int m) {
  double **blad = utworz_macierz(n, m);
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < m; j++) {
      blad[i][j] = fabs(przyblizenie[i][j] - dokladne[i][j]);
    }
  }

  return blad;
}


double *maxBlad(double **blad, int n, int m) {
  double *wynik = utworz_wektor(n);

  for (int i = 0; i < n; i++)
    wynik[i] = norma_max(blad[i], m);

  return wynik;
}

double obliczDT(double lambda, double h, double d) {
  return (lambda * h * h) / d;
}


double *obliczOdstepyH(double dt, int n, int m) {
  double *wynik = utworz_wektor(m);
  double x = XMIN;

  for (int i = 0; i < m; i++) {
    wynik[i] = x;
    x += H;
  }

  return wynik;
}


double *obliczOdstepyDT(double dt, int n, int m) {
  double *wynik = utworz_wektor(n);
  double t = TMIN;

  for (int i = 0; i < n; i++) {
    wynik[i] = t;
    t += dt;
  }

  return wynik;
}

void zapiszDwaWektory(double *wektor1, double *wektor2, int n, string nazwa_pliku) {
  fstream file(nazwa_pliku.c_str(), ios::out);

  if (file.is_open()) {
    for (int i = 0; i < n; i++) {
      file << wektor1[i] << "\t" << wektor2[i] << endl;
    }
  }
}

void zapiszRozwiazanie_zad2(double **macierz, double *wektorKroki, int rozmiar, int pozycja, std::string nazwaPliku) {
  double *temp = utworz_wektor(rozmiar);
  for (int i = 0; i < rozmiar; ++i) {
    temp[i] = macierz[pozycja][i];
  }
  zapiszDwaWektory(wektorKroki, temp, rozmiar, nazwaPliku);
}

double estymator(double *xNowe, double *xPoprzednie, int m) {
  double max = 0.0;
  double *p = new double[m];

  for (int i = 0; i < m; i++)
    p[i] = xNowe[i] - xPoprzednie[i];

  if (fabs(p[0]) > fabs(p[1]))
    max = fabs(p[0]);
  else
    max = fabs(p[1]);

  for (int i = 0; i < m; i++) {
    if (fabs(p[i]) > max)
      max = fabs(p[i]);
  }

  delete[] p;
  return max;
}

double *residuum(double **macierz, double *b, double *x, int m) {
  double sum = 0.0;
  double *wynik = new double[m];
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < m; j++) {
      sum += macierz[i][j] * x[j];
    }
    wynik[i] = sum - b[i];
    sum = 0.0;
  }
  return wynik;
}