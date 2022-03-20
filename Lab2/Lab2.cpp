
#include <iostream>
#include <math.h>
#include <fstream>

using namespace std;

double funkcja_z_zadania(double x) {
  return (1 - exp(-x)) / x;
}

double szereg(double x)
{
    double wynik=1;
    double znak=-1;
    double krok=1;

    for(int i=2;i<30;++i)
    {
        krok=krok*(x/i);
        wynik+=znak*krok;
        znak=-znak;
    }

    return wynik;
}

int main(){
    ifstream plikDane;
    ofstream daneWynikowe, danePoprawione;
    double wartosc_plik=0,wartosc_obliczona=0,wartosc_log_10=0,argument=0,blad=0,blad_log=0;

    daneWynikowe.open("daneWynikowe.txt");
    danePoprawione.open("danePoprawione.txt");
    plikDane.open("dane.txt");

    cout<<"|      Argument Funkcji      |     Blad Wzgledny     |    log10(blad)    | "
     " ~Blad Wzgledny    |   ~log10(blad)    |\n";

    while(!plikDane.eof())
    {
        plikDane>>wartosc_log_10;
        plikDane>>argument;
        plikDane>>wartosc_plik;
    

    //obliczanie wartości funkcji oraz błędu względengo
    wartosc_obliczona=funkcja_z_zadania(argument);
    blad=abs((wartosc_obliczona-wartosc_plik)/wartosc_plik);

    //zależnośc logarytmu 10 z bezwzględną wartością błędu względnego
    blad_log=log10(blad);

    //zapisane wyniku do pliku

    daneWynikowe<<wartosc_log_10<<" "<<blad_log<<"\n";
    
    cout << "|          ";
    cout << argument << "         |         ";
    cout << blad << "         |         ";
    cout << blad_log << "         |         ";

    //Obliczanie wartości funkcji oraz obliczanie błedu względnego z szerefu 
    wartosc_obliczona=szereg(argument);
    blad=abs((wartosc_obliczona-wartosc_plik)/wartosc_plik);

    //zależnośc logarytmu 10 z bezwzględną wartością błędu względnego
    blad_log=log10(blad);

    //zapisanie wynik do pliku

    danePoprawione<<wartosc_log_10<< " "<<blad_log<<"\n";

    cout << blad << "         |         ";
    cout << blad_log << "         |\n";
    }

    daneWynikowe.close();
    daneWynikowe.close();
    plikDane.close();
}