#include <iostream>
#include <iomanip>

using namespace std;

void wypelnienie( double *l, double *d, double *u, double *b )
{
	u[0] = 1.0 / 2.0; u[1] = 1.0 / 4.0; u[2] = 1.0 / 6.0; u[3] = 1.0 / 8.0; u[4] = 1.0 / 10.0;u[5]=0.0;
	d[0] = 10.0; d[1] = 20.0; d[2] = 30.0; d[3] = 30.0; d[4] = 20.0; d[5] = 10.0;
	l[0]=0.0;l[1] = 1.0 / 3.0; l[2] = 1.0 / 5.0; l[3] = 1.0 / 7.0; l[4] = 1.0 / 9.0; l[5] = 1.0 / 11.0;

	b[0] = 31.0; b[1] = 165.0/ 4.0; b[2] = 917.0 / 30.0; b[3] =  851.0 / 28.0; b[4] = 3637.0 / 90.0; b[5] = 332.0 / 11.0;
}

void wypiszWek( double *wek, int n)
{
	for ( int i = 0; i < n; i++ )
	{
		cout << "| ";
		cout.width( 8 );
		cout <<setprecision(6) << wek[i];
		cout.width( 8 );
		cout << "| ";
		cout << endl;
	}
}
void operacjem(double *l, double *d, double *u, double *eta,int n){
	eta[0] = d[0];
	for (int i = 1; i <= n; i++)
	{
	eta[i]=d[i] -= l[i] / d[i - 1] * u[i - 1];
	}
}
void operacjew(double *l, double *d,double *r,double *b ,int n){
	r[0] = b[0];
	for (int i = 1; i <= n; i++)
	{
	r[i]=b[i] -= l[i] / d[i - 1] * b[i - 1];
	}
}
void rozwiazanie_ukladu( double *l, double *d, double *u, double *b, double *x,double *eta,double *r, int n )
{
	x[n-1] = (1.0 / eta[n-1])*r[n-1];
	for (int i = n-2; i >= 0; i--)
	{
		x[i] = ((r[i] - (u[i] * x[i+1]))/eta[i]);
	}
	cout << endl << "Wektor eta " << endl;
	wypiszWek( eta, n );
	cout << endl << "Wektor r " << endl;
	wypiszWek( r,  n);
}

void rozwiazanie()
{
	double *l, *d, *u, *b, *x,*eta,*r;
	int n = 6;

	l = new double[n];
	d = new double[n ];
	u = new double[n];
	b = new double[n];
	x = new double[n];
	eta= new double[n];
	r = new double[n];
	
	wypelnienie( l, d, u, b );
	
	cout << endl << "Wektor l " << endl;
	wypiszWek( l, n);
	cout << endl << "Wektor d " << endl;
	wypiszWek( d, n  );
	cout << endl << "Wektor u " << endl;
	wypiszWek( u, n );
	cout << endl << "Wektor b " << endl;
	wypiszWek( b, n );
	cout << endl;

	operacjem(l,d,u,eta,n);
	operacjew(l,d,r,b,n);
	rozwiazanie_ukladu( l, d, u, b, x, eta, r,  n );

	cout << endl << "Wektor x " << endl;
	wypiszWek( x, n );
	cout << endl;

	delete[] l;
	delete[] d;
	delete[] u;
	delete[] x;
	delete[] b;

}


int main()
{
	rozwiazanie( );
	return 0;
}
