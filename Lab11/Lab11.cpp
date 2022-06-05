#include<iostream>
#include<fstream>
#include<cmath>
#include<string>
#include "/mnt/c/Users/User/Desktop/Zajęcia/MO/Lab11/fun.cpp"
#define TOLX 10e-16
#define TOLF 10e-16
using namespace std;

double U_anal(double x, double t){
    return exp((-(x*x))/(4.0*D*(TAU+t)))*1./(2*pow(M_PI*D*(TAU+t),1./2));
}

void warunek_brzegowy(double** macierz, int n, int m){
    for(int i=1;i<n;i++){
        macierz[i][0]=0.0;
        macierz[i][m-1]=0.0;
    }
}

void warunek_poczatkowy(double **macierz, int n, int m){
    double x=XMIN;
      for (int i = 0; i < m; i++) {
          macierz[0][i]=U_anal(x,0);
          x += H;
    }
    
  }


double ** rozwiazanie_analityczne(int n, int m,double h, double dt){
    double **rozwiazanie=utworz_macierz(n,m);
    double x=XMIN;
    double t=TMIN;

      for (int i = 0; i < n; i++) {
    for (int j = 0; j < m; j++) {
      rozwiazanie[i][j] = U_anal(x, t);
      x += h;
    }
    x = XMIN;
    t += dt;
  }
  return rozwiazanie;

}
void algorytmThomasa(double *l, double *d, double *u, double *b, double *x, int m) {
  for (int i = 1; i < m; ++i) {
    d[i] = d[i] - (l[i - 1] / d[i - 1]) * u[i - 1];
    b[i] = b[i] - (l[i - 1] / d[i - 1]) * b[i - 1];
  }
  x[m - 1] = b[m - 1] / d[m - 1];
  for (int i = m - 2; i >= 0; --i) {
    x[i] = (b[i] - u[i] * x[i + 1]) / d[i];
  }
}
// najwyżej wrócić do thomasa 
void mlThomas(double **macierz, int n, int m) {
  double L_LAMBDA = 1.0 + 2.0 * LAMBDA_POSREDNIE;
  double *l = utworz_wektor(m);
  double *d = utworz_wektor(m);
  double *u = utworz_wektor(m);
  double *b = utworz_wektor(m);
  double *x = utworz_wektor(m);

  for (int k = 1; k < n; k++) {
    l[0] = 0.0;
    d[0] = 1.0;
    u[0] = 0.0;
    b[0] = macierz[k - 1][0];

    for (int i = 1; i < m - 1; i++) {
      l[i] = LAMBDA_POSREDNIE;
      d[i] = -L_LAMBDA;
      u[i] = LAMBDA_POSREDNIE;
      b[i] = -macierz[k - 1][i]; // zmiana 
    }

    l[m - 1] = 0.0;
    d[m - 1] = 1.0;
    u[m - 1] = 0.0;
    b[m - 1] = 0.0;

    algorytmThomasa(l, d, u, b, x, m);

    for (int i = 1; i < m - 1; i++) {
      macierz[k][i] = x[i];
    }
  }
  usun_wektor(l);
  usun_wektor(d);
  usun_wektor(u);
  usun_wektor(b);
 usun_wektor(x);
}


void SOR( double **macierz, double *b,double *x,int n, int m){

  double TOL=1e-16;
  double *x_nowy_wektor=utworz_wektor(m);
  double *x_pomocniczy= utworz_wektor(m);
  double sum = 0;
  for(int i = 0; i < m; i++){
    x_pomocniczy[i]=0;
  }
    for(int iter = 0; iter < ITERACJE; iter++){

        for(int i = 0; i < m; i++){

            double temp=0;

            for(int j = 0; j < m; j++){
                if(j>i){
                    temp += -macierz[i][j] * x[j];
                }
                if(j==i){
                    temp += -( (1. - 1./OMEGA) * macierz[i][j]) * x[j];
                }
            }

            temp +=b[i];
            x_nowy_wektor[i] = temp;
        }

        for(int i = 0; i < m; i++){
            double temp=0;
            for(int j = 0; j < m; j++){
                if(j<i) {
                    temp += macierz[i][j] * x_pomocniczy[j];
                }
            }
            x_pomocniczy[i] = (x_nowy_wektor[i] - temp)/((1./OMEGA)*macierz[i][i]);
        }
        for (int i = 0; i < m; i++) {
            x[i] = x_pomocniczy[i];
        }
            if ((fabs((norma_max(residuum(macierz, b, x, m), m))) < TOL) && (fabs((estymator(x_pomocniczy, x, m))) < TOL))
            {
              break;
            }
}
}


void mlSor(double **macierz, int n, int m){
  double tempLabda = 1.0 + 2.0 * LAMBDA_POSREDNIE;
  double *b = new double[m];
  double *wynik = new double[m];
  double **tempMacierz = new double *[m]; // poprawa

  for (int i = 0; i < m; i++) {
    b[i] = 0;
    wynik[i] = 0.0;
  }

  for (int i = 0; i < m; i++){
    tempMacierz[i] = new double[m];
  }

    for (int i = 0; i < m; i++) {
    for (int j = 0; j < m; j++)
      tempMacierz[i][j] = 0.0;
  }

 for (int k = 1; k < n; k++) {
    tempMacierz[0][0] = 1.0;
    b[0] = 0;

    for (int i = 1; i < m - 1; i++) {
      tempMacierz[i][i] = -tempLabda;
      tempMacierz[i][i + 1] = LAMBDA_POSREDNIE;
      tempMacierz[i][i - 1] = LAMBDA_POSREDNIE;
      b[i] = -macierz[k - 1][i];
      
    }
    b[m - 1] = 0.0;
    tempMacierz[m - 1][m - 1] = 1.0;

    SOR(tempMacierz, b, wynik, m, m);

    for (int i = 1; i < m - 1; i++)
      macierz[k][i] = wynik[i];
  }


}


void kmb(double **macierz, int n, int m) {
  for (int i = 1; i < n; i++) {
    for (int j = 1; j < m - 1; j++) {
      macierz[i][j] =
          macierz[i - 1][j]
              + LAMBDA_BEZPOSREDNIE * (macierz[i - 1][j - 1] - (2 * macierz[i - 1][j]) + macierz[i - 1][j + 1]);
    }
  }
}


double **kmbRozwiazanie(int n, int m) {
  double **rozwiazanie = utworz_macierz(n, m);
  warunek_poczatkowy(rozwiazanie, n, m);
  warunek_brzegowy(rozwiazanie, n, m);
  kmb(rozwiazanie, n, m);
  return rozwiazanie;
}

double **mlThomasRozwiazanie(int n, int m) {
  double **laasonen_thomas = utworz_macierz(n, m);
  warunek_poczatkowy(laasonen_thomas, n, m);
  warunek_brzegowy(laasonen_thomas, n, m);
  mlThomas(laasonen_thomas, n, m);
  return laasonen_thomas;
}

double **mlSorRozwiazanie(int n, int m){
  double **rozwiazanie = utworz_macierz(n, m);
  warunek_poczatkowy(rozwiazanie, n, m);
  warunek_brzegowy(rozwiazanie, n, m);
  mlSor(rozwiazanie, n, m);
  return rozwiazanie;
}

int main(){
double dt = obliczDT(LAMBDA_BEZPOSREDNIE, H, D);
  int n = ((TMAX - TMIN) / dt);
  int m = ((XMAX - XMIN) / H);
  double **rozwiazanieAnalityczne;
  double **rozwiazanieKmb;
  double **rozwiazanieLaasonenThomas;
  double **rozwiazanieLaasonenSOR;
  double **macierzBledy;
  double *wektorBledy;
  double *odstepDT;
  double *odstepX;

  //analityczne

    rozwiazanieAnalityczne = rozwiazanie_analityczne(n, m, H, dt);
     zapiszMacierz(rozwiazanieAnalityczne, n, m, "rozwiazanieAnalityczne.csv");

   //KMB
   rozwiazanieKmb = kmbRozwiazanie(n, m);
  zapiszMacierz(rozwiazanieKmb, n, m, "rozwiazanieKmb.csv");
 macierzBledy = obliczBlad(rozwiazanieAnalityczne, rozwiazanieKmb, n, m);
  wektorBledy = maxBlad(macierzBledy, n, m);
  zapiszWektor(wektorBledy, n, "maxErrKmb.csv");
  zapiszMacierz(macierzBledy, n, m, "errMacierzKmb.csv");

  odstepX = obliczOdstepyH(dt, n, m);
  odstepDT = obliczOdstepyDT(dt, n, m);
  zapiszWektor(odstepDT, n, "odstepyCzasoweKmb.csv");
  zapiszWektor(odstepX, n, "odstepyXKmb.csv");
  zapiszRozwiazanie_zad2(rozwiazanieKmb, odstepX, m, 84, "1rozKMB.csv");
  zapiszRozwiazanie_zad2(rozwiazanieAnalityczne, odstepX, m, 84, "1rozAnalityczne.csv");
  
  // Metoda Posrednia Laasonem + Algorytm Thomasa
  dt = obliczDT(LAMBDA_POSREDNIE, H, D);
  n = ((TMAX - TMIN) / dt);
  m = ((XMAX - XMIN) / H);

  rozwiazanieLaasonenThomas = mlThomasRozwiazanie(n, m);
  zapiszMacierz(rozwiazanieLaasonenThomas, n, m, "laasonenThomasRozwiazanie.csv");

  macierzBledy = obliczBlad(rozwiazanieAnalityczne, rozwiazanieLaasonenThomas, n, m);
  wektorBledy = maxBlad(macierzBledy, n, m);
  zapiszWektor(wektorBledy, n, "maxErrLaasonenThomas.csv");
  zapiszMacierz(macierzBledy, n, m, "errLaasonenThomas.csv");

  odstepX = obliczOdstepyH(dt, n, m);
  odstepDT = obliczOdstepyDT(dt, n, m);
  zapiszWektor(odstepX, n, "odstepyXLaasonenThomas.csv");
  zapiszWektor(odstepDT, n, "odstepyCzasoweLaasonenThomas.csv");

  zapiszRozwiazanie_zad2(rozwiazanieLaasonenThomas, odstepX, m, 84, "1rozLT.csv");
  // Laasonem + SOR
  dt = obliczDT(LAMBDA_POSREDNIE, H, D);
  n = ((TMAX - TMIN) / dt);
  m = ((XMAX - XMIN) / H);
  rozwiazanieLaasonenSOR = mlSorRozwiazanie(n, m);
  zapiszMacierz(rozwiazanieLaasonenSOR , n, m, "laasonenSorRozwiazanie.csv");

  macierzBledy = obliczBlad(rozwiazanieAnalityczne, rozwiazanieLaasonenSOR , n, m);
  wektorBledy = maxBlad(macierzBledy, n, m);
  zapiszWektor(wektorBledy, n, "maxErrLaasonenSOR.csv");
  zapiszMacierz(macierzBledy, n, m, "errLaasonenSOR.csv");

  odstepX = obliczOdstepyH(dt, n, m);
  odstepDT = obliczOdstepyDT(dt, n, m);
  zapiszWektor(odstepX, n, "odstepyXLaasonenSOR.csv");
  zapiszWektor(odstepDT, n, "odstepyCzasoweLaasonenSOR.csv");

  zapiszRozwiazanie_zad2(rozwiazanieLaasonenSOR, odstepX, m, 84, "1rozLSOR.csv");

usun_macierz(rozwiazanieAnalityczne, n);
usun_macierz(rozwiazanieLaasonenThomas, n);
usun_macierz(rozwiazanieLaasonenSOR, n);
usun_macierz(macierzBledy, n);
usun_wektor(wektorBledy);
usun_wektor(odstepX);
usun_wektor(odstepDT);

}