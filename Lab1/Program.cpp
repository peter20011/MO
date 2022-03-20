#include<iostream>


using namespace std;

int main()
{
    float bity_flo=0.0f;
    float precyzja_flo, epsilon_,zmienna_pomocnicza;
    precyzja_flo=1.0f;
    do
    {   
        bity_flo++;
        precyzja_flo/=2.0f;
        epsilon_=1.0f+precyzja_flo;
        epsilon_=epsilon_-1.0f
        zmienna_pomocnicza=epsilon_+1.0f
    } while (zmienna_pomocnicza);
    cout<<"Precyzja arytmetyki float wynosi: "<<2.0f*precyzja_flo<<endl;
    cout<<"Liczba bitów mantysy wynosi: "<< bity_flo-1.0f<<endl;
}