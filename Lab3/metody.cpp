#include<iostream>
#include<cmath>

using namespace std;

typedef double (*funkcja)(double);
/**
 * @brief metoda Picarda do rozwiązania algebraicznych równań nieliniowych
 * 
 * @param funkcja_poczatkowa funkcja dla której szukamy rozwiązania
 * @param fi przekształcenie funkcji do postaci x=
 * @param pochodna_fi pochodna z funkcji fi
 * @param x punkt początkowy
 * @param max_iteracji górna granic iteracji
 * @param tolerancja_bledu  zadana tolerancja bledu
 * @param tolerancja_residuum  zadanda tolerancja residuum
 */

double Picarda(funkcja funkcja_poczatkowa, funkcja fi, funkcja pochodna_fi, double x, double max_iteracji, double tolerancja_bledu, double tolerancja_residuum) {

    cout << "Metoda Picarda" << endl;

    if (fabs(pochodna_fi(x)) >= 1.0) {
        cout << "Funkcja rozbiezna" << endl;
        return 0;
    }

    double estymator = 0.0, residuum = 0.0, przyblizenie = x;
    int iteracja = 1;
    bool kontynuuj = true;

    cout << " ITERACJA - WARTOSC PRZYBLIZONA - RESIDUUM - ESTYMATOR BLEDU" << endl;

    while (kontynuuj) {
        cout << iteracja << " - ";

        // Obliczanie przybliżenia
        przyblizenie = fi(przyblizenie);
        cout << przyblizenie << " - ";

        estymator = fabs(przyblizenie - x);

        x = przyblizenie;
        residuum = funkcja_poczatkowa(x);

        cout << residuum << " - " << estymator << endl;

        //warunek przerwania
        if ((fabs(residuum) <= tolerancja_residuum) && (estymator <= tolerancja_bledu) || (iteracja >= max_iteracji)) {
            kontynuuj = false;
        }

        iteracja += 1;
    }

    return 0;
}

/**
 * @brief metoda bisekcji do rozwiazywania algebraicznych równań nielinowych
 * 
 * @param funkcja_poczatkowa funkcja dla której szukamy rozwiazywania
 * @param a lewy kraniec przedzialu
 * @param b  prawy kraniec przedzialu
 * @param max_iteracji górna granica iteracji
 * @param tolerancja_bledu zadana tolerancja bledu
 * @param tolerancja_residuum  zadana tolerancja residuum
 */

double bisekcja(funkcja funkcja_poczatkowa, double a, double b, double max_iteracji, double tolerancja_bledu, double tolerancja_residuum) {

    cout << "Metoda bisekcji" << endl;

    if ((funkcja_poczatkowa(a) > 0 && funkcja_poczatkowa(b) > 0) ||
        (funkcja_poczatkowa(a) < 0 && funkcja_poczatkowa(b) < 0)) {
        cout << "Funkcja nie zmienia znaku dla podanego przedzialu" << endl;
        return 0;
    }

    double estymator = 0.0, residuum = 0.0, x = 0.0;
    int iteracja = 1;
    bool kontynuuj = true;

    cout << "ITERACJA - A - B - X- RESIDUUM - ESTYMATOR BLEDU" << endl;

    while (kontynuuj) {
        // obliczanie srodka przedziału
        x = (a + b) / 2.0;
        estymator = fabs((b - a) / 2.0);
        residuum = funkcja_poczatkowa(x);

        cout << iteracja << " - " << a << " - " << b << " - " << x << " - " << residuum << " - " << estymator << endl;

        //wybranie nowego punktu skrajnego
        if ((funkcja_poczatkowa(a) < 0 && funkcja_poczatkowa(x) > 0) ||
            (funkcja_poczatkowa(a) > 0 && funkcja_poczatkowa(x) < 0)) {
            b = x;
        } else {
            a = x;
        }

        if ((fabs(residuum) <= tolerancja_residuum) && (estymator <= tolerancja_bledu) || (iteracja >= max_iteracji)) {
            kontynuuj = false;
        }

        iteracja += 1;
    }

    return 0;
}

/**
 * @brief metoda bisekcji do rozwiazywania algebraicznych równań nielinowych
 * @param funkcja_poczatkowa funkcja dla której szukamy rozwiazywania
 * @param funkcja_pochodna pochodna funkcji
 * @param x punkt początkowy
 * @param max_iteracji górna granica iteracji
 * @param tolerancja_bledu zadana tolerancja bledu
 * @param tolerancja_residuum zadana tolerancja residuum
 */
double Newton(funkcja funkcja_poczatkowa, funkcja funkcja_pochodna, double x, double max_iteracji, double tolerancja_bledu, double tolerancja_residuum) {

    cout << "Metoda Newtona" << endl;

    double x0 = x, x1 = 0.0, estymator = 0.0, residuum = 0.0;
    int iteracja = 1;
    bool kontynuuj = true;

    cout << "ITERACJA - X1 - RESIDUUM - ESTYMATOR BLEDU" << endl;

    while (kontynuuj) {
        //obliczanie wartosci x1
        x1 = x0 - (funkcja_poczatkowa(x0) / funkcja_pochodna(x0));
        estymator = fabs(x0 - x1);
        x0 = x1;
        residuum = funkcja_poczatkowa(x0);

        cout << iteracja << " - " << x1 << " - " << residuum << " - " << estymator << endl;

        //warunki przerwania
        if ((fabs(residuum) <= tolerancja_residuum) && (estymator <= tolerancja_bledu) || (iteracja >= max_iteracji)) {
            kontynuuj = false;
        }

        iteracja++;
    }

    return 0;
}

/**
 * @brief metoda siecznych do rozwiazywania algebraicznych równań nielinowych
 * @param funkcja_poczatkowa funkcja dla której szukamy rozwiazywania
 * @param x0 punkt
 * @param x1 punkt x + delta x
 * @param max_iteracji górna granica iteracji
 * @param tolerancja_bledu  zadana tolerancja bledu
 * @param tolerancja_residuum  zadana tolerancja residuum
 */
double siecznych(funkcja funkcja_poczatkowa, double x0, double x1, double max_iteracji, double tolerancja_bledu, double tolerancja_residuum) {

    cout << "Metoda siecznych" << endl;

    double x2 = 0.0, estymator = 0.0, residuum = 0.0;
    int iteracja = 0;
    bool kontynuuj = true;

    cout << "Iteracja - x1 - residuum - estymator bledu" << endl;

    while (kontynuuj) {
        x2 = x1 - funkcja_poczatkowa(x1) / ((funkcja_poczatkowa(x1) - funkcja_poczatkowa(x0)) / (x1 - x0));
        estymator = fabs(x2 - x1);
        residuum = funkcja_poczatkowa(x2);

        cout << iteracja << " - " << x1 << " - " << residuum << " - " << estymator << endl;

        // ustawienie nowych punktów
        x0 = x1;
        x1 = x2;

        //warunki przerwania
        if ((fabs(residuum) <= tolerancja_residuum) && (estymator <= tolerancja_bledu) || (iteracja >= max_iteracji)) {
            kontynuuj = false;
        }

        iteracja++;
    }

    return 0;
}