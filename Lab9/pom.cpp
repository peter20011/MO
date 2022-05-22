#include<iostream>
#include<fstream>
#include<cmath>

using namespace std;


//rownanie analityczne

double rownanieAnalityczne(double x){
    return ( exp( 2.0-2.0*x ) - 4*exp( 4.0 - x * 2.0  )+ 4*exp( x * 2.0 ) - exp( 2.0 + 2.0*x ) - x + x*exp(4.0)) / ( 4.0 - 4*exp(4.0) );
}


//funkcja wykonująca algorytm Thomasa na zadanych wektorach, ponieważ zadana macierz jest trójdiagonalna

void Thomas(double *l, double *d, double *u, double *b, double *x, int N){
    double *dKopia= new double [N];
    double  *bKopia= new double [N];

    //zmienne pomocnicze
    dKopia[0]=d[0];
    bKopia[0]=b[0];

    // przekształcenie Algorytm THOOMASA
    for(int i=1;i<N;i++){
       dKopia[i] = d[i] - (l[i] * u[i - 1]) /dKopia[i - 1];
    }

    for (int i = 1; i < N; i++) 
	{
        bKopia[i] = b[i] - (l[i] * bKopia[i - 1] )/ dKopia[i - 1];
    }

    x[N-1]=bKopia[N-1]/dKopia[N-1];

    
    for (int i = N - 2; i >= 0; i--)
	{
         x[i] = (bKopia[i] - u[i] * x[i + 1]) / dKopia[i];
    }

    //Usuwamy niepotrzebne tymczasowe wektory
    delete[] bKopia; 
	delete[] dKopia;


}


int najwiekszyBlad( double *blad, int N ) //wyznaczanie indeksu największego błędu
{
    //pierwszym maksimum jest pierwszy element tablicy blad, potem porównyjemy
	double maksymalny = fabs( blad[0] ); 
	int naj;
	for ( int i = 0; i < N; i++ )
		if ( fabs( blad[i] ) > maksymalny )
		{
		maksymalny = fabs(blad[i]);
		}
	return maksymalny;
}

