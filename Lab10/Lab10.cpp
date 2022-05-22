#include <iostream>
#include <fstream>
#include <math.h>
using namespace std;

double wzorAnalityczny( double t ); //funkcja do rozwiazywania równania analitycznie
double bezposredniaEulera( double h, double tmax ); //funkcja do rozwiazywania równania za pomocą metody bezpośredniej Eulera
double bezposredniaEuleraBlad( double h, int b); //funkcja do obliczania błędu dla metody bezpośredniej Eulera
double posredniaEulera( double h, double t ); //funkcja do rozwiazywania równania za pomocą metody bpośredniej Eulera
double posredniaEuleraBlad( double h, int ilosc ); //funkcja do obliczania błędu dla metody pośredniej Eulera
double trapezow( double h, double t ); //funkcja do rozwiazywania równania za pomocą metody trapezów
double trapezowBlad( double h, int ilosc ); //funkcja do obliczania błędu dla metody trapezów

int main()
{
	fstream bledy, wyniki, wyniki2; //tworzymy pliki
	bledy.open( "bledy.txt", fstream::out ); //otwieramy je
	wyniki.open( "wyniki.txt", fstream::out );
	wyniki2.open( "wyniki2.txt", fstream::out );
	wyniki<< scientific;
	wyniki2<< scientific; //ustalamy precyzję zapisu na 10 cyfrsaveToFile("Conventional.txt", x, h, n);
	bledy<< scientific;
	cout.precision(10);
	double h, analitycznieWynik, besWynik, bensWynik, peWynik, tWynik;
	int N = 1000; //maksymalna ilosć iteracji
	h=0.1;
	for ( h ; h > 1e-16; h = h/2 ) //pętla do obliczania błędów z metod, od granicy stabilności jaką jest 0.1
	{
		besWynik = log10( bezposredniaEuleraBlad( h, N ) ); //stabilna bezpośrednia Eulera
		peWynik = log10( posredniaEuleraBlad( h, N ) );
		tWynik = log10( trapezowBlad( h, N ) );
		bledy << log10( h ) << " " << besWynik << " " << peWynik << " " << tWynik << endl;
	}

	h = 0.01;
	for ( double t = 0; t < 5; t += 0.01 ) //pętla wyników do rysunku zależności y od zmiennej niezależnej t
	{
		analitycznieWynik = wzorAnalityczny( t );
		peWynik = posredniaEulera( h, t );
		tWynik = trapezow( h, t );
		besWynik = bezposredniaEulera( h, t );
		wyniki << t << " " << analitycznieWynik << " " << peWynik << " " << tWynik << " " << besWynik << " "<< bensWynik<< endl;
	}
	h=0.15;
	for(double t = 0 ; t < 4 ; t+=0.15)
    {
		bensWynik = bezposredniaEulera(h, t ); //niestabilna bezpośrednia Eulera
		wyniki2 <<t<<" "<<wzorAnalityczny( t )<<" "<<bensWynik<<endl;
	}
	wyniki2.close();
	bledy.close( );
	wyniki.close( );
}

double wzorAnalityczny( double t )
{
	return 1 - exp( -10.0 * ( t + atan( t ) ) );
}

double bezposredniaEulera( double h, double t )
{
	double y = 0.0;	//warunek poczatkowy y(0) = 0
	for ( double i = 0.0; i < t; i += h )
	{
		y = y + h *( -( ( 10.0 *i * i + 20.0 ) / ( i * i + 1.0 ) )*( y - 1.0 ) ); //obliczamy kolejny wyraz metodą bezpośrednią Eulera
	}
	return y;
}

double bezposredniaEuleraBlad( double h, int ilosc)
{
	double y = 0.0;	// z warunku poczatkowego y(0) = 0
	double blad = 0.0;
	double t = h;
	double wAnalityczna;
	double roznica = 0.0;
	for ( int i = 0; i < ilosc; i++ )
	{
		wAnalityczna = wzorAnalityczny( t );
		y = y + h *( -( ( 10.0 *t * t + 20.0 ) / ( t * t + 1.0 ) )*( y - 1.0 ) );
		roznica= fabs(wAnalityczna-y);
		if ( roznica > blad ) //znajdujemy największy błąd
			blad = roznica;
		t += h;
	}
	return blad;
}

double posredniaEulera( double h, double t )
{
	double y = 0.0;	
	double wspolczynnik;
	for ( double i = 0.0; i < t; i += h )
	{
		wspolczynnik = ( 10.0 * ( i + h ) * ( i + h ) + 20.0 ) / ( ( i + h ) * ( i + h ) + 1.0 ); //dla skrócenia obliczeń oddzielnie obliczamy ułamek z wyrażenia
		y = ( y + h * wspolczynnik ) / ( 1 + h * wspolczynnik ); //przekształcamy wzór (yk1-yk)/h+(10*t^2+20)/(t^2+1)*(yk1-1)=0 do postaci yk1=...
	}
	return y;
}

double posredniaEuleraBlad( double h, int ilosc )
{
	double y = 0.0;	
	double blad = 0.0;
	double t = h;
	double wAnalityczna;
	double wspolczynnik;
	double roznica= 0.0;
	
	for ( int i = 0; i < ilosc; i++ )
	{
		wAnalityczna = wzorAnalityczny( t );
		wspolczynnik = ( 10.0 * ( t + h ) * ( t + h ) + 20.0 ) / ( ( t + h ) * ( t + h ) + 1.0 );
		y = ( y + h * wspolczynnik ) / ( 1 + h * wspolczynnik ); 
		roznica= fabs(wAnalityczna-y);
		if ( roznica > blad ) //znajdujemy największy błąd
			blad = roznica;
		t += h;
	}
	return blad;
}

double trapezow( double h, double t )
{
	double y = 0.0, wspolczynnik0, wspolczynnik1;

	for ( double i = 0.0; i < t; i += h )
	{
		wspolczynnik0 = ( ( 10.0 * i * i + 20.0 ) / ( i * i + 1.0 ) ); //część f(tk,yk)
		wspolczynnik1 = ( 10.0 * ( i + h ) * ( i + h ) + 20.0 ) / ( ( i + h ) * ( i + h) + 1.0 ); //część f(tk+1,yk+1)
		y = ( ( -h/ 2.0 ) * ( wspolczynnik0 * ( y - 1.0 ) - wspolczynnik1 ) + y ) / ( 1.0 + ( h / 2.0 ) * wspolczynnik1 );
	} //przekształcamy wzór (yk1-yk)/h-(1/2)*((10*t^2+20)/(t^2+1)*(yk-1)+(10*(t+h)^2+20)/((t+h)^2+1)*(yk1-1))= 0 do postaci yk1=...
	return y;
}

double trapezowBlad( double h, int ilosc )
{
	double blad = 0.0;
	double t = h;
	double y = 0.0;	
	double wAnalityczna;
	double wspolczynnik0, wspolczynnik1;
	double roznica = 0.0;
	
	for ( int i = 0; i < ilosc; i++ )
	{
		wAnalityczna = wzorAnalityczny( t );
		wspolczynnik0 = ( ( 10.0 * t * t + 20.0 ) / ( t * t + 1.0 ) );
		wspolczynnik1 = ( 10.0 * ( t + h ) * ( t + h ) + 20.0 ) / ( ( t + h ) * ( t + h) + 1.0 ); ;
		y = ( ( -h/ 2.0 ) * ( wspolczynnik0 * ( y - 1.0 ) - wspolczynnik1 ) + y ) / ( 1.0 + ( h / 2.0 ) * wspolczynnik1 );
		roznica= fabs(wAnalityczna-y);
		if ( roznica > blad ) //znajdujemy największy błąd
			blad = roznica;
		t += h;
	}
	return blad;
}