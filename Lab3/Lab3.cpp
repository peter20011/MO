#include<iostream>
#include "/mnt/c/Users/User/Desktop/Zajęcia/MO/Lab3/funkcje.cpp"
#include "/mnt/c/Users/User/Desktop/Zajęcia/MO/Lab3/metody.cpp"
using namespace std;

int main(){
    int max_iteracji =100;
    double blad = 1e-10;

cout << "\t\tPierwsze Rownanie\t\t" << endl;

    Picarda(funkcja_pierwsza, fi_pierwsza, pochodna_fi_pierwsza, 0.3, max_iteracji, blad, blad);
    bisekcja(funkcja_pierwsza, -0.6, 1.0, max_iteracji, blad, blad);
    Newton(funkcja_pierwsza, pochodna_funkcja_pierwsza, 0.5, max_iteracji, blad, blad);
    siecznych(funkcja_pierwsza, -0.6, 1.0, max_iteracji, blad, blad);

    cout << "\t\tDrugie Rownanie\t\t" << endl;
    
    Picarda(funkcja_druga, fi_druga, pochodna_fi_druga, 1.0, max_iteracji, blad, blad);
    bisekcja(funkcja_druga, -4.0, 10.0, max_iteracji, blad, blad);
    Newton(funkcja_druga, pochodna_funkcja_druga, 0.4, max_iteracji, blad, blad);
    siecznych(funkcja_druga, 0.1, 0.6, max_iteracji, blad, blad);


}