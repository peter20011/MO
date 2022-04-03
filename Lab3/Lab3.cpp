#include<iostream>
#include "/mnt/c/Users/User/Desktop/Zajęcia/MO/Lab3/funkcje.h"
#include "/mnt/c/Users/User/Desktop/Zajęcia/MO/Lab3/metody.h"
using namespace std;


int main(){
    int max_iteracji =100;
    double blad = 10e-15;

    cout<<"\tPierwsze rownanie\t"<<endl;
    Picard(funkcja_pierwsza,&fi_z_pierwszej,&pochodna_fi_z_pierwszej,0.3,max_iteracji,blad,blad);
    bisekjca(funkcja_pierwsza,-1.0,1.0,max_iteracji,blad,blad);
    Newton(funkcja_pierwsza,&pochodna_funkcja_pierwsza,0.5,max_iteracji,blad,blad);
    siecznych(funkcja_pierwsza,-1.0,1.0,max_iteracji,blad,blad);

    cout<<"\tDrugie  rownanie\t"<<endl;
    Picard(funkcja_druga,&fi_z_drugiej,pochodna_fi_z_drugiej,1.0,max_iteracji,blad,blad);
    bisekjca(funkcja_druga,-4.0,10.0,max_iteracji,blad,blad);
    Newton(funkcja_druga,&pochodna_funkcja_druga,0.4,max_iteracji,blad,blad);
    siecznych(funkcja_druga,0.1,0.6,max_iteracji,blad,blad);

}