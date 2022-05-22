#include<iostream>
#include <cmath>
#include <fstream>
#include "/mnt/c/Users/User/Desktop/Zajęcia/MO/Lab9/pom.cpp"
using namespace std;

double p = 1.0, q = 0.0, r = -4.0; //współczynniki w równaniu
double alfa = 0.0, beta = 1.0, gamma1 = -1.0; //współczynniki z warunków brzegowych
double fi = 0.0, psi = 1.0, teta = 0.0; 
double poczatekPrzedzialu = 0.0, koniecPrzedzialu = 1.0; 

double Numerowa( double h, int N ) 
{
	double *l, *d, *u, *b, *x, *blad, xp1 = poczatekPrzedzialu, xp2 = poczatekPrzedzialu; 
	fstream numerow, analitycznie;
	numerow.open( "wynikiNumerow.txt", ios_base::app ); //otwieramy pliki do zapisu wyników wyliczeń analitycznych oraz z Numerowa
	analitycznie.open( "wynikiAnalitycznie.txt", ios_base::app);
	analitycznie << scientific;
	numerow << scientific;
	cout.precision(10);
	l = new double[N]; 
	d = new double[N];
	u = new double[N];
	b = new double[N];
	x = new double[N];
	blad = new double[N];

	u[0] = alfa / h; //zgodnie z wzorami uzupełniamy wartości pierwszych wyrazów, które są inne niż pozostałe
	d[0] = beta - alfa / h;
	b[0] = -gamma1;
	
	for(int i = 1; i < N - 1; i++) //w pętli uzupełniamy wyrazy środkowe
	{
		l[i - 1] = p / (h * h) + r / 12.0;
		d[i] = (-2.0 * p) / (h * h) + r * (10.0 / 12.0);
		u[i] = p / (h * h) + r / 12.0;
		b[i] = (xp1 + i * h - h) / 12.0 + (10.0 / 12.0) *  (xp1 + i * h) + (xp1 + i * h + h) / 12.0;
	}

	l[N - 2] = -fi / h; //końcowe wyrazy mają także specjalne dane
	d[N - 1] = -fi / h + psi;
	b[N - 1] = -teta;

	Thomas( l, d, u, b, x, N ); //wykonujemy algorytm Thomasa na wyliczonych danych

	for ( int i = 0; i < N; i++ ) //obliczamy błąd bezwzględny między wyliczeniami Numerowa, a rozwiązaniem analitycznym
	{
		blad[i] = fabs( x[i] - rownanieAnalityczne( xp2 ) );
		xp2 += h;
	}
	int naj = najwiekszyBlad(blad, N); //znajdujemy największy błąd
	
	if(N==162) //losowa liczba
 	{
 		for ( int i = 0; i < N; i++ )
		{
			numerow<<xp1<<"\t"<<x[i]<<"\t"<<endl;
			analitycznie<<xp1<<"\t"<<rownanieAnalityczne( xp1)<<endl;
			xp1 += h;
		}
    }


	delete[] l; //usuwamy niepotrzebne wektory
	delete[] d;
	delete[] u;
	delete[] x;
	delete[] b;
	analitycznie.close(); //zamykamy pliki
	numerow.close();
	return naj; //zwracamy największy błąd
}

double konwencjonalnaTrzypunktowa( double h, int N ) //funkcja realizująca dyskretyzację konwencjonalną trzypunktową
{
	double *l, *d, *u, *b, *x, *blad, xp1 = poczatekPrzedzialu, xp2 = poczatekPrzedzialu; //wykonujemy analogiczne operacje jak w Numerowa
	fstream konwencjonalnie;
	konwencjonalnie.open( "wynikiKonwencjonalnie.txt", ios_base::app);
	konwencjonalnie << scientific;
	cout.precision(10);
	l = new double[N];
	d = new double[N];
	u = new double[N];
	b = new double[N];
	x = new double[N];
	blad = new double[N];

	u[0] = alfa / h;
	d[0] = beta - alfa / h;
	b[0] = -gamma1;

	for ( int i = 1; i < N - 1; i++ )
	{
		l[i - 1] = p / ( h * h ) - r / ( 2.0 * h ); //inne wartości w środku
		d[i] = ( -2.0 * p ) / ( h *h ) + r;
		u[i] = p / ( h * h ) + q / ( 2.0 * h );
		b[i] = (xp1+i*h); //nasze s
	}

	l[N - 2] = -fi / h;
	d[N - 1] = -fi / h + psi;
	b[N - 1] = -teta;

	Thomas( l, d, u, b, x, N );

	for ( int i = 0; i < N; i++ )
	{
		blad[i] = fabs( x[i] -rownanieAnalityczne( xp2 ) );
		xp2 += h;
	}
	
 	int naj = najwiekszyBlad(blad, N);
 	if(N==162) //losowa liczba
 	{
 		for ( int i = 0; i < N; i++ )
		{
			konwencjonalnie<<xp1<<"\t"<<x[i]<<endl;
			xp1 += h;
		}
    }

	delete[] l;
	delete[] d;
	delete[] u;
	delete[] x;
	delete[] b;

	return naj;
}

int main() //funkcja główna
{
	double h, bladNum, bladKonw;
	int N; //ilość iteracji

	fstream bledy, rzedy;
	bledy.open( "bledy.txt", fstream::out );
	bledy << scientific; //ustawienie precyzji
	cout.precision(10);
	
	for ( N = 2 ; N < 5000; N += 80)
	{
		h = ( koniecPrzedzialu - poczatekPrzedzialu ) / ( N-1 );
		bladKonw = ( konwencjonalnaTrzypunktowa( h, N ) );
		bladNum =( Numerowa( h, N ) );
		bledy<<log10(h)<<"\t"<<log10(bladKonw)<<"\t"<<log10(bladNum)<<endl;
	}
	bledy.close();
	return 0;
}
