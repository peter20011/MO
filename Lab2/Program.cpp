#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include<iostream>

using namespace std;

double oBlad(double x, double wynik)
{
	return fabs((x - wynik) / wynik);
}
double ObliczSzereg(double x, int N)
{
	double licznik, mianownik, wynik, wart_wyl, spr;
	int j;
	
	licznik = 1.0;
	mianownik = 1.0;
	wynik = 1.0;
	wart_wyl = wynik;

	

	
	for ( j = 2; j < N; j++)
	{
		licznik *= x;
		mianownik *= j;
		wynik = licznik / mianownik;
		wart_wyl = wart_wyl + wynik;
	}

}



int main(void)
{
    	double x, blad, wd, ww;


	printf("x                 |wartosc funkcji exp(x)             | Suma szeregu taylora        |Wynik z danych         | blad         \n");
	printf("------------------------------------------------------------------------------------------------------\n");

	
		FILE *plik;
	plik = fopen("dane1.txt","r");
    fclose(plik);
    return 0;
}