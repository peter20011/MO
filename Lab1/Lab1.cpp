#include<iostream>


using namespace std;

int main()
{
float licznikbitowfloat=0.0;
double licznikbitowdouble=0.0;

float precyzjafloat=1.0f, pomocnicza_float=1.1f;

while(pomocnicza_float >1.0f){
    precyzjafloat/=2.0f;
    pomocnicza_float=1.0f+precyzjafloat;
    if(pomocnicza_float!=1.0f)
        licznikbitowfloat++;
}

double precyzjadouble=1.0f, pomocnicza_double=1.1f;

while(pomocnicza_double >1.0f){
    precyzjadouble/=2.0f;
    pomocnicza_double=1.0f+precyzjadouble;
    if(pomocnicza_double!=1.0f)
        licznikbitowdouble++;
}


cout<<"Liczba bitow mantysy  dla float wynosi:"<<licznikbitowfloat<<endl;
cout<<"Epsilon maszynosy dla float wynosi:"<<2*precyzjafloat<<endl;
cout<<endl;
cout<<"Liczba bitow mantysy  dla double wynosi:"<<licznikbitowdouble<<endl;
cout<<"Epsilon maszynosy dla dobule wynosi:"<<2*precyzjadouble<<endl;
}