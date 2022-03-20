#include<iostream>


using namespace std;

int main()
{
float licznikbitowfloat=0.0;
double licznikbitowdouble=0.0;

// dla float
float precyzjafloat=1.0f, pomocnicza_float=1.1f; //pomocznicza_float większa od 1 by while się wykonywał

while(pomocnicza_float >1.0f){
    precyzjafloat/=2.0f; // przechodzimy odpowiednio po bitach mantysy o 1
    pomocnicza_float=1.0f+precyzjafloat;
    if(pomocnicza_float!=1.0f) // nie zwiększamy licznika dla ostatniego przebiegu pętli
        licznikbitowfloat++;
}

//dla double

double precyzjadouble=1.0f, pomocnicza_double=1.1f;//pomocznicza_fdouble większa od 1 by while się wykonywał

while(pomocnicza_double >1.0f){
    precyzjadouble/=2.0f; // przechodzimy odpowiednio po bitach mantysy o 1
    pomocnicza_double=1.0f+precyzjadouble;
    if(pomocnicza_double!=1.0f) // nie zwiększamy licznika dla ostatniego przebiegu pętli
        licznikbitowdouble++;
}


cout<<"Liczba bitow mantysy  dla float wynosi:"<<licznikbitowfloat<<endl;
cout<<"Epsilon maszynosy dla float wynosi:"<<2*precyzjafloat<<endl; // epsilon maszynowy -> precyzja jest prawie 0 więc trzeba po wyjściu z pętli pomnożyć razy 2
cout<<endl;
cout<<"Liczba bitow mantysy  dla double wynosi:"<<licznikbitowdouble<<endl;
cout<<"Epsilon maszynosy dla dobule wynosi:"<<2*precyzjadouble<<endl; // epsilon maszynowy -> precyzja jest prawie 0 więc trzeba po wyjściu z pętli pomnożyć razy 2
}