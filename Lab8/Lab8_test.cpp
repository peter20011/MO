#include<iostream>
#include<cmath>
#include<iomanip>
#include<fstream>
#define M_PI 3.14159265358979323846
using namespace std;
double TOLH=1e-16;
using namespace std;

template <typename T>
T roznica_progresywna_2(T x, T h) {
    return (sin(x + h) - sin(x))/h;
}

template <typename T>
T  roznica_progresywna_3(T x, T h) {
    return ((T)-1.5*sin(x) + (T)2.0*sin(x + h) - (T)0.5*sin(x + (T)2.0*h))/h;
}

template <typename T>
T roznica_wsteczna_2(T x, T h) {
    return (sin(x) - sin(x - h))/h;
}

template <typename T>
T roznica_wsteczna_3(T x, T h) {
    return ((T)0.5*sin(x - (T)2.0*h) - (T)2.0*sin(x - h) + (T)1.5*sin(x))/h;
}


template <typename T>
T roznica_centralna_2(T x, T h) {
    return (sin(x+h)-sin(x-h))/((T)2.0*h);
}

template <template T> void oblicz_roznice(const string& file)