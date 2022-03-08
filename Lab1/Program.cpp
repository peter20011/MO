#include<iostream>
#include<math.h>

using namespace std;


int main(){

    //double
	double a,i,b,c; // a -prezycja , i-terator(licznik bitów), b -prezycja zmiennej, 
	float aa,ii,bb,cc;
	a=1.0; // ostatnia która ma 1 w ->cecha 
	i=0;
	
	do{
		i++;
		a/=2.0; // dzielmy przez 2 by dokonać 1 bitowego przemieszczenia po matysie 
		b=a+1.0; // kiedy mantysa dojdzie do postaci 0.000 dojdzie do przerwania pętli
		c=b-1.0;// zmienna po to by wykonał się while
	}while(c);
	

    //float       
	aa=1.0;
	ii=0;
	do{
		ii++;
		aa/=2.0;
		bb=aa+1.0;
		cc=bb-1.0f;
	}while(cc);
	
	cout<<"Epsylon (iloczyn maszynowy) dla flout wynosi: "<<aa*2.0<<endl; // iloczyn maszynowy = precyzja * 2
	cout<<"Liczba bitów mantysy dla flout wynosi: "<< ii-1<<endl;
	
	cout<<"Epsylon (iloczyn maszynowy) dla double wynosi: "<<a*2.0<<endl;  // iloczyn maszynowy = precyzja * 2
	cout<<"Liczba  bitów mantysy dla double wynosi: "<< i-1<<endl; // nie liczymy dla ostatniego przejścia pętli
}