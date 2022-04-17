#include <iostream>
#include <math.h>
#include "/mnt/c/Users/User/Desktop/Zajęcia/MO/Lab4/funkcje.cpp"


using namespace std;

double **Jakobian(double x, double y,double z){
    double **matrix=new double *[3];
    for(int i=0;i<3;i++)
    {
        matrix[i]=new double[3];
    }

    matrix[0][0] = 2.0 * x;
    matrix[0][1] = 2.0 * y;
    matrix[0][2] = 2.0 * z;
    matrix[1][0] = 2.0 * x;
    matrix[1][1] = 2.0 * y;
    matrix[1][2] = 0.0;
    matrix[2][0] = 2.0 * x;
    matrix[2][1] = -1.0;
    matrix[2][2] = 0.0;

    return matrix;
}

void Newton(double x, double y,double z){
    double **macierz=NULL;
    int i = 1;
    bool kontynuuj = true;
    double estymator = 0.0, residuum = 0.0, wartosci_funkcji[3], wektor_ABC[3];

    while(kontynuuj){
        //czyszczenie wektora przy ponownej iteracji
        for(int j=0;j<3;j++){
            wektor_ABC[j]=0;
        }

        if(!czy_nalezy_do_dziedziny(x,y,z)){
            break;
        }

        //Jakobian 

        macierz=Jakobian(x,y,z);

        //wyznaczenie wektor abc (wektor delta)

        wektor_ABC[0] = oblicz_a(x, y, z);
        wektor_ABC[1] = oblicz_b(x, y, z);
        wektor_ABC[2] = oblicz_c(x, y, z);

        // nowe przybliżenie 

        x = x - wektor_ABC[0];
        y = y - wektor_ABC[1];
        z = z - wektor_ABC[2];

        //przypisujemy do funkcji nowe wartości

        wartosci_funkcji[0]=funkcja_pierwsza(x,y,z);
        wartosci_funkcji[1]=funkcja_druga(x,y,z);
        wartosci_funkcji[2]=funkcja_trzecia(x,y,z);

        //estymator
        estymator=max(fabs(wektor_ABC[0]),fabs(wektor_ABC[1]),fabs(wektor_ABC[2]));

        //residuum
        residuum=max(fabs(wartosci_funkcji[0]), fabs(wartosci_funkcji[1]), fabs(wartosci_funkcji[2]));

        cout << "Iteracja " << i << "\nx: " << x << "\ny: " << y << "\nz: " << z << "\nestymator: " << estymator<< "\nresiduum: " << residuum << "\n\n";

        //warunki przerwania

        if ((fabs(residuum) <= TOL_F) || (estymator <= TOL_X) || (i >= MAX_ITERACJE)) {
            kontynuuj = false;
        }
        i++;
    }

    

    //usuwanie pamięci

    if(macierz!=NULL){
        for(int j=0;j<3;j++){
            delete [] macierz[j];
        }
        
        delete macierz;
        macierz=NULL;
    }
}

int main(){
    Newton(5.0,3.0,5.0);
}

