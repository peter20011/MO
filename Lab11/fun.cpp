#include<iostream>
#include<fstream>
#include<cmath>
#include<string>

using namespace std;

double *utwórz_wektor(int n){
    return new double[n];
}

double** utwórz_macierz(int n, int m){
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