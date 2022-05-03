#include<iostream>
#include<cmath>
#include<iomanip>
#include<fstream>
#define M_PI 3.14159265358979323846
using namespace std;
double TOLH=1e-16;

template <typename T>
T forwardDifference2(T x, T h) {
    return (sin(x + h) - sin(x))/h;
}

template <typename T>
T forwardDifference3(T x, T h) {
    return ((T)-1.5*sin(x) + (T)2.0*sin(x + h) - (T)0.5*sin(x + (T)2.0*h))/h;
}

template <typename T>
T backwardDifference2(T x, T h) {
    return (sin(x) - sin(x - h))/h;
}

template <typename T>
T backwardDifference3(T x, T h) {
    return ((T)0.5*sin(x - (T)2.0*h) - (T)2.0*sin(x - h) + (T)1.5*sin(x))/h;
}

template <typename T>
T centralDifference2(T x, T h) {
    return (sin(x+h)-sin(x-h))/((T)2.0*h);
}

template <typename T>
void calculateDerivatives(const string& filename) {

    ofstream wyniki;
    wyniki.open(filename);

    T x1 = 0.0;
    T x2 = M_PI/4;
    T x3 = M_PI/2;

    T x1_wart_dokladna = 1.0;
    T x2_wart_dokladna = sqrt(2.0)/2.0;
    T x3_wart_dokladna = 0.0;

    T *rzad = new T[9];
    T h1;
    T h2;
    T *w1 = new T[9];
    T *w2 = new T[9];

    T h = 0.1;
    int i = 0;
    while (h  > TOLH) {
        T roz_prog_2_x1 = forwardDifference2(x1, h);
        T roz_prog_3_x1 = forwardDifference3(x1, h);
        T roz_wstecz_2_x2 = backwardDifference2(x2, h);
        T roz_wstecz_3_x2 = backwardDifference3(x2, h);
        T roz_central_2_x2 = centralDifference2(x2, h);
        T roz_prog_2_x2 = forwardDifference2(x2, h);
        T roz_prog_3_x2 = forwardDifference3(x2, h);
        T roz_wstecz_2_x3 = backwardDifference2(x3, h);
        T roz_wstecz_3_x3 = backwardDifference3(x3, h);
/*
        std::cout << std::setw(12) << h << std::setw(12) << roz_prog_2_x1 << std::setw(12) << roz_prog_3_x1 << std::setw(12) << roz_wstecz_2_x2 <<
        std::setw(12) << roz_wstecz_3_x2 << std::setw(12) << roz_central_2_x2 << std::setw(12) << roz_prog_2_x2 << std::setw(12) << roz_prog_3_x2 <<
        std::setw(12) << roz_wstecz_2_x3 << std::setw(12) << roz_wstecz_3_x3 << std::endl;
*/
        cout << log10(h) << " " << log10(fabs(roz_prog_2_x1 - x1_wart_dokladna)) << " " << log10(fabs(roz_prog_3_x1 - x1_wart_dokladna)) << " "
               << log10(fabs(roz_wstecz_2_x2 - x2_wart_dokladna)) << " " << log10(fabs(roz_wstecz_3_x2 - x2_wart_dokladna)) << " "
               <<  log10(fabs(roz_central_2_x2 - x2_wart_dokladna)) << " " <<log10(fabs(roz_prog_2_x2 - x2_wart_dokladna))  << " "
               << log10(fabs(roz_prog_3_x2 - x2_wart_dokladna)) << " " << log10(fabs(roz_wstecz_2_x3 - x3_wart_dokladna)) << " "
               <<  log10(fabs(roz_wstecz_3_x3 - x3_wart_dokladna)) << endl;


        wyniki << log10(h) << " " << log10(fabs(roz_prog_2_x1 - x1_wart_dokladna)) << " " << log10(fabs(roz_prog_3_x1 - x1_wart_dokladna)) << " "
        << log10(fabs(roz_wstecz_2_x2 - x2_wart_dokladna)) << " " << log10(fabs(roz_wstecz_3_x2 - x2_wart_dokladna)) << " "
        <<  log10(fabs(roz_central_2_x2 - x2_wart_dokladna)) << " " <<log10(fabs(roz_prog_2_x2 - x2_wart_dokladna))  << " "
        << log10(fabs(roz_prog_3_x2 - x2_wart_dokladna)) << " " << log10(fabs(roz_wstecz_2_x3 - x3_wart_dokladna)) << " "
        <<  log10(fabs(roz_wstecz_3_x3 - x3_wart_dokladna)) <<endl;

        if(i == 1) {
            h1 = log10(h);
            w1[0] = log10(fabs(roz_prog_2_x1 - x1_wart_dokladna));
            w1[1] = log10(fabs(roz_prog_3_x1 - x1_wart_dokladna));
            w1[2] = log10(fabs(roz_wstecz_2_x2 - x2_wart_dokladna));
            w1[3] = log10(fabs(roz_wstecz_3_x2 - x2_wart_dokladna));
            w1[4] = log10(fabs(roz_central_2_x2 - x2_wart_dokladna));
            w1[5] = log10(fabs(roz_prog_2_x2 - x2_wart_dokladna));
            w1[6] = log10(fabs(roz_prog_3_x2 - x2_wart_dokladna));
            w1[7] = log10(fabs(roz_wstecz_2_x3 - x3_wart_dokladna));
            w1[8] = log10(fabs(roz_wstecz_3_x3 - x3_wart_dokladna));
        }
        if(i == 2) {
            h2 = log10(h);
            w2[0] = log10(fabs(roz_prog_2_x1 - x1_wart_dokladna));
            w2[1] = log10(fabs(roz_prog_3_x1 - x1_wart_dokladna));
            w2[2] = log10(fabs(roz_wstecz_2_x2 - x2_wart_dokladna));
            w2[3] = log10(fabs(roz_wstecz_3_x2 - x2_wart_dokladna));
            w2[4] = log10(fabs(roz_central_2_x2 - x2_wart_dokladna));
            w2[5] = log10(fabs(roz_prog_2_x2 - x2_wart_dokladna));
            w2[6] = log10(fabs(roz_prog_3_x2 - x2_wart_dokladna));
            w2[7] = log10(fabs(roz_wstecz_2_x3 - x3_wart_dokladna));
            w2[8] = log10(fabs(roz_wstecz_3_x3 - x3_wart_dokladna));
        }

        h /= 2;
        i++;
    }

    for(int l = 0; l < 9; l++) {
        rzad[i] = fabs(w2[l]-w1[l])/fabs(h2-h1);
        cout << rzad[i] <<endl;
    }

    wyniki.close();
}

int main() {

    cout << "DOUBLE: " << endl;
    calculateDerivatives<double>("double.txt");
    cout << endl;
    cout << "FLOAT: " << endl;
    calculateDerivatives<float>("float.txt");
    cout << endl;


    return 0;
}