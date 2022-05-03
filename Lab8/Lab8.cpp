#include <cmath>
#include <iostream>
#include <stdio.h>
#include <fstream>
#include <stdlib.h>

#define METHODS 9
#define ITER 20

using namespace std;

//Wz�r funkcji
template <typename T>
T function(T x){
	return sin(x);
}


//Pochodna rzeczywista
template <typename T>
T derivative(T x){
	return cos(x);
}


//R�nica wsteczna dwupunktowa
template <typename T>
T backward_diff_2(T x, T h){
	return (function(x) - function(x - h)) / h;
}


//R�nica progresywna dwupunktowa
template <typename T>
T forward_diff_2(T x, T h){
	return (function(x + h) - function(x)) / h;
}


//R�nica centralna dwupunktowa
template <typename T>
T central_diff_2(T x, T h){
	return (function(x + h) - function(x - h)) / ((T)(2.0) * h);
}


//R�nica wsteczna 3-punktowa
template <typename T>
T backward_diff_3(T x, T h){
	return (((T)(1.0 / 2.0) * function(x - (T)2.0 * h)) - ((T)(2.0) * function(x - h)) + ((T)(3.0 / 2.0) * function(x))) / h;
}


//R�nica progresywna 3-punktowa
template <typename T>
T forward_diff_3(T x, T h){
	return (((T)(-3.0 / 2.0) * function(x)) + ((T)(2.0) * function(x + h)) - ((T)(1.0 / 2.0) * function(x + (T)2.0 * h))) / h;
}


//Macierz
template <typename T>
T **allocate(){
	T **er = new T*[ITER];
	
	for(int i = 0; i <= ITER; i++)
		er[i] = new T[METHODS];
		
	return er;
}


//Zapisanie do pliku
template <typename T>
void create_file(T** err, const char* filename){
	ofstream output;
	output.open(filename);
	
	if(!output.good()){
		cout << "Blad otwarcia pliku" << endl;
		exit(1);
	}
	
	for(int i = 0; i < ITER; i++){
		for(int j = 0; j <= METHODS; j++)
			output << scientific << log10(err[i][j]) << " ";
		output << endl;
	}
	output.close();
}


//Obliczanie b��du
template <typename T>
void calculate_err(const char* filename){
	T** err = allocate<T>();
	T h = 0.1;
	T start = 0;
	T end = M_PI / 2.0;
	T center = (start + end) / 2.0;
	
	for(int i = 0; i < ITER; i++){
		err[i][0] = h;
		err[i][1] = fabs(forward_diff_2(start, h) - derivative(start));
		err[i][2] = fabs(forward_diff_3(start, h) - derivative(start));
        err[i][3] = fabs(backward_diff_2(center,h)- derivative(center));
        err[i][4]= fabs(backward_diff_3(center,h)- derivative(center));
		err[i][5] = fabs(central_diff_2(center, h) - derivative(center));
        err[i][6] =fabsl(forward_diff_2(center,h)- derivative(center));
        err[i][7]  =fabs (forward_diff_3(center,h)-derivative(center));
		err[i][8] = fabs(backward_diff_2(end, h) - derivative(end));
		err[i][9] = fabs(backward_diff_3(end, h) - derivative(end));
		
		h *= (T)0.2;
		for(int j = 0; j < METHODS; j++)
			printf("%.6e\t", err[i][j]);
		cout << endl;
	}
	
	for(int i = 1; i < METHODS; i++)
		cout << "Rzad dokladnosci: " << (log10(err[3][i]) - log10(err[2][i])) / (log10(err[3][0]) - log10(err[2][0])) << endl;
	
	cout << endl;
	create_file(err, filename);
}

int main(){
	calculate_err<float>("float.txt");
	calculate_err<double>("double.txt");
	return 0;
}
