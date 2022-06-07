#include <iostream>
#include <fstream>
#include <ctime>
#include <cmath>
#include <algorithm>
#include "/mnt/c/Users/User/Desktop/Zajęcia/MO/Lab11/fun.cpp"

using namespace std;

double U_anal(double x,double t){
    return (1./(sqrt(M_PI*D*(TAU+t))))*exp(-(x*x)/(4*D*(t+TAU)));
}


//wyliczanie rozwiązań analitycznych
double **analytical(int n, int m, double delta_t) {
    double **results = allocMatrix(n, m);

    double x = x_start;
    double t = t_min;

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            results[i][j] = U_anal(x,t);
            x += h;
        }
        x = x_start;
        t += delta_t;
    }
    return results;
}

//wyliczanie rozwiązań numerycznych dla metody KMB
// n - liczba kroków czasowych
// m - liczba kroków przestrzennych
// delta_t - wyliczone z lambda*h^2/D dla lambda = 0.4, D = 1
double **KMB(int n, int m, double delta_t) {
    double **results = allocMatrix(n, m); //macierz reprezentująca siatke czasowo-przestrzenną
    double x = x_start;
    for (int i = 0; i < m; i++) {
        results[0][i] = U_anal(x,0); 
        x+= h; 
    }


    for (int i = 0; i < n; i++) {
        results[i][0] = 0.0;  //warunek brzegowy U(r,t)
        results[i][m - 1] = 0.0;
 
    }

    x = x_start;
    for (int k = 1; k < n; k++) {
        for (int i = 1; i < m - 1; i++) {
                results[k][i] = results[k - 1][i]+ lambda_kmb * (results[k - 1][i - 1] - (2 *  results[k - 1][i]) +  results[k - 1][i + 1]);
            x += h;
            x += h;
        }
        x = x_start;
    }
    return results;
}

void Thomas(double *Lower, double *Diagonal, double *Upper, double *b, double *x, const int m) {
    auto *vector_l = new double[m]; //elementy na przekątnej przechowywujący eta
    for (int i = 1; i < m; i++) {
        vector_l[i] = Lower[i] * (1 / Diagonal[i - 1]); //l_i*eta_i-1^-1
        Diagonal[i] -= vector_l[i] * Upper[i - 1]; //eta_i = d_i - l_i*eta_i-1^-1 * u_i-1
    }
    for (int i = 1; i < m; i++) {
        b[i] -= vector_l[i] * b[i - 1]; //r_i = b_i - l_i*eta_i-1^-1 * r_i-1
    }
    x[m - 1] = (1 / Diagonal[m - 1] * b[m - 1]); // x_n = eta_n^-1 * r_n
    for (int i = m - 2; i >= 0; i--) {
        x[i] = (1 / Diagonal[i]) * (b[i] - Upper[i] * x[i + 1]); //x_i = eta_i^-1(r_i - u_i * x_i+1)
    }
    delete[] vector_l;
}

//wyliczanie rozwiązań numerycznych dla metody Laasonen
// n - liczba kroków czasowych
// m - liczba kroków przestrzennych
// delta_t - wyliczone z lambda*h^2/D dla lambda = 1, D = 1
double **LaasonenThomas(int n, int m, double delta_t) {
    //rozwiazywanie okladu rownan za pomoca macierzy trójdiagonalnej
    auto *lower = new double[m];    //elementy pod przekątną
    auto *diagonal = new double[m]; //elementy na przekątnej
    auto *upper = new double[m]; //elementy nad przekątną
    auto *b = new double[m]; // wektor b
    auto *results = new double[m]; //wektor x z niewiadomymi

    double **matrixA = allocMatrix(n, m); //macierz reprezentująca siatkę czasowo-przestrzenną
    double x = x_start;
    for (int i = 0; i < m; i++) {
        matrixA[0][i] = U_anal(x,0);   
        x+=h; //warunek początkowy
    }


    for (int i = 0; i < n; i++) {
        matrixA[i][0] = 0.0; //warunek brzegowy U(r,t)
        matrixA[i][m - 1] = 0.0;//warunek brzegowy U(r+a,t)
    }


    x = x_start;
    for (int k = 1; k < n; k++) {
        lower[0] = 0.0;  //zawsze 0
        diagonal[0] = 1.0; //alfa=0, beta =1
        upper[0] = 0.0; //alfa = 0
        b[0] = matrixA[k - 1][0]; //gamma = 0, pierwszy warunek brzegowy
        for (int i = 1; i < m - 1; i++) {
            lower[i] = lambda_laasonen; // wyliczone w sprawozdaniu
            diagonal[i] = -(1.0 + 2.0 * lambda_laasonen );
            upper[i] = lambda_laasonen;
            b[i] = -matrixA[k - 1][i]; //b[i] = -U_(n-1,k)
            x += h;
        }
        x = x_start;
        lower[m - 1] = 0.0; //phi =0
        diagonal[m - 1] = 1.0; // phi=0 psi=1
        upper[m - 1] = 0.0; // zawsze 0
        b[m - 1] = matrixA[k - 1][m - 1]; // theta = 1-(r/(r+a))erfc(a/(2*sqrt(D*t))

        Thomas(lower, diagonal, upper, b, results, m); //algorytm thomasa do wyliczenia rozwiazan układu

        for (int i = 1; i < m; i++) { //dodawanie przybliżeń do siatki dla kolejnych t
            matrixA[k][i] = results[i];
        }
    }

    return matrixA;
}

bool error(double *xn, double *xn_1, int m) {
    int count = 0;

    for (int i = 0; i < m; i++) {
        if (fabs(xn_1[i] - xn[i]) < 1e-10) count++;
    }
    return count == m;
}

bool residue(double **A, double *b, const double *xn, int m) {
    int count = 0;
    for (int i = 0; i < m; i++) {
        double sum = 0.0;
        for (int j = 0; j < m; j++) {
            sum += A[i][j] * xn[j];
        }
        if (fabs(sum - b[i]) < 1e-10) count++;
    }
    return count == m;
}

//matrix - macierz A
//b - wektor b
//x - wektor x
//m - rozmiar macierzy mxm
void sor(double **matrix, double *b, double *x, int m) {
    const double omega = 1.0; //dla omega = 1.0 otrzymałem najszybszą zbierzność
    auto *right = new double[m]; //wektor przechowujący wyniki prawej strony wzoru operacyjnego
    auto *approx = new double[m]; //wektor przyblizen poczatkowych
    auto *xn_1 = new double[m]; //wektor kolejnego przyblizenia
    for (int i = 0; i < m; i++) {
        right[i] = 0.0; //na początku 0
        approx[i] = 1.0; //ponieważ w poprzednich metodach wychodziło 1 to dałem przybliżenie początkowe 1
    }

    for (int limit = 0; limit < 10000; limit++) { //arbitralne ograniczenie na liczbe iteracji
        for (int i = 0; i < m; i++) {
            double tmp = ((1.0 - 1.0 / omega) * matrix[i][i]) * approx[i]; // ((1- 1/omega)D*x_n-1

            for (int j = i + 1; j < m; j++) {
                tmp += matrix[i][j] * approx[j]; //U*x_n-1
            }
            right[i] = -tmp + b[i]; // ((1- 1/omega)D+U)x_n-1
        }

        for (int i = 0; i < m; i++) {
            xn_1[i] = 0.0; //zerowanie nastepnego przyblizenia
        }

        for (int i = 0; i < m; i++) {
            double tmp = 0.0;
            for (int j = 0; j <= i; j++) {
                tmp += xn_1[j] * matrix[i][j]; // L*x_n
            }
            xn_1[i] = (right[i] - tmp) / ((1.0 / omega) * matrix[i][i]); //wyliczenie x_n
        }

        bool en = error(xn_1, approx, m); //kryterium dokładności wyznaczenia
        bool re = residue(matrix, b, approx, m); // kryterium wiarygodności x_n
        if (en && re) break;
        swap(xn_1, approx);
    }
    for (int i = 0; i < m; i++) {
        x[i] = xn_1[i];
    }
}


//wyliczanie rozwiązań numerycznych dla metody Laasonen
// n - liczba kroków czasowych
// m - liczba kroków przestrzennych
// delta_t - wyliczone z lambda*h^2/D dla lambda = 1, D = 1
double **LaasonenSOR(int n, int m, double delta_t) {
    auto *b = new double[m]; //wektor b
    auto *results = new double[m];// wektro x

    double **matrixA = allocMatrix(n, m); //macierz reprezentujaca siatke czasowo-przestrzenną
    double **matrix = allocMatrix(m, m); //macierz do obliczeń rozwiązań układu równań

    for (int i = 0; i < m; i++) {
        for (int j = 0; j < m; j++) {
            matrix[i][j] = 0.0;  //zerowanie elementów macierzy
        }
    }
    double x = x_start;
    for (int i = 0; i < m; i++) {
        matrixA[0][i] = U_anal(x,0);
        x += h;
    }

    for (int i = 0; i < n; i++) {   //warunki brzegowe
        matrixA[i][0] = 0.0;
        matrixA[i][m - 1] = 0.0;

    }

    double t = t_min;
    x = x_start;
    for (int k = 1; k < n; k++) {
        matrix[0][0] = 1.0; // d1: beta = 1
        b[0] = 0.0; // gamma = 0, 1 warunek brzegowy
        for (int i = 1; i < m - 1; i++) {
            matrix[i][i - 1] = lambda_laasonen; //wyliczone w sprawozdaniu
            matrix[i][i] = -(1.0 + 2.0 * lambda_laasonen);
            matrix[i][i + 1] = lambda_laasonen;
            b[i] = -matrixA[k - 1][i]; ////b[i] = -U_(n-1,k)
            x += h;
        }
        x = x_start;
        matrix[m - 1][m - 1] = 1.0; // dn: psi = 1
        b[m - 1] = 0.0; /// 2 warunek brzegowy
        t += delta_t;

        sor(matrix, b, results, m); // metoda iteracyjne sor do obliczena ukladu rownan

        for (int i = 1; i < m; i++) {
            matrixA[k][i] = results[i];//dodawanie przybliżeń do siatki dla kolejnych t
        }
    }

    return matrixA;
}

//dane do wykresów dla KMB
void KMB_save() {
    double delta_t = getDelta_t(lambda_kmb);
    int n = (int) ((t_max - t_min) / delta_t) + 2; //liczba kroków czasowych
    int m = (int) ((x_end - x_start) / h); // liczba kroków przestrzennych

    double **analytic = analytical(n, m, delta_t);
    double **results = KMB(n, m, delta_t);
    double **error = allocMatrix(n, m);

    countErrors(error, analytic, results, n, m);

    auto *steps_t = new double[n];
    double t = t_min;
    for (int i = 0; i < n; i++) {
        steps_t[i] = t;
        t += delta_t;
    }


    ofstream file1, file2, file3;
    file1.open("../kmb125.txt");
    file2.open("../kmb250.txt");
    file3.open("../kmb500.txt");


    double x = x_start;
    cout << steps_t[125] << endl;
    for (int i = 0; i < m; i++) {
        if (abs(results[125][i + 1]) < 1e-16)
            break;
        file1 << x << " " << analytic[125][i + 1] << " " << results[125][i + 1] << endl;  // t = 0.36;
        x += h;
    }
    x = x_start;
    cout << steps_t[250] << endl;
    for (int i = 0; i < m; i++) {
        if (abs(results[250][i + 1]) < 1e-16)
            break;
        file2 << x << " " << analytic[250][i + 1] << " " << results[250][i + 1] << endl;  // t = 1.08;
        x += h;
    }
    x = x_start;
    cout << steps_t[500] << endl;
    for (int i = 0; i < m; i++) {
        if (abs(results[500][i + 1]) < 1e-16)
            break;
        file3 << x << " " << analytic[500][i + 1] << " " << results[500][i + 1] << endl;  // t = 1.8;
        x += h;
    }


    file1.close();
    file2.close();
    file3.close();

    double *maxErrors = maxError(error, n, m);
    ofstream fileError;
    fileError.open("../kmbmaxError.txt");
    for (int i = 0; i < n - 1; i++) {
        if ((maxErrors[i + 1]) < 1e-16)
            break;
        fileError << steps_t[i + 1] << " " << maxErrors[i + 1] << endl;
    }
    fileError.close();

    saveToFile("analytic.txt", analytic, n, m);
    saveToFile("kmbNumerical.txt", results, n, m);
    saveToFile("kmbError.txt", error, n, m);


    delete[] maxErrors;
    delete[] steps_t;
    deleteMatrix(analytic, n);
    deleteMatrix(results, n);
    deleteMatrix(error, n);
}

//dane do wykresów dla Laasonen- algorytm Thomasa
void LassonenThomasSave() {
    double delta_t = getDelta_t(lambda_laasonen);
    int n = (int) ((t_max - t_min) / delta_t) + 2;
    int m = (int) ((x_end - x_start) / h);

    double **analytic = analytical(n, m, delta_t);
    double **results = LaasonenThomas(n, m, delta_t);
    double **error = allocMatrix(n, m);

    countErrors(error, analytic, results, n, m);

    auto *steps_t = new double[n];
    double t = t_min;
    for (int i = 0; i < n; i++) {
        steps_t[i] = t;
        t += delta_t;
    }

    ofstream lt100, lt300, lt400;
    lt100.open("../lt100.txt");
    lt300.open("../lt300.txt");
    lt400.open("../lt400.txt");


    double x = x_start;
    cout << steps_t[500] << endl;
    for (int i = 0; i < m; i++) {
        if (abs(results[500][i + 1]) < 1e-16)
            break;
        lt100 << x << " " << analytic[500][i + 1] << " " << results[500][i + 1] << endl;  // t = 0.3125;
        x += h;
    }
    x = x_start;
    cout << steps_t[1200] << endl;
    for (int i = 0; i < m; i++) {
        if (abs(results[1200][i + 1]) < 1e-16)
            break;
        lt300 << x << " " << analytic[1200][i + 1] << " " << results[1200][i + 1] << endl;  // t = 0.75;
        x += h;
    }
    x = x_start;
    cout << steps_t[2000] << endl;
    for (int i = 0; i < m; i++) {
        if (abs(results[2000][i + 1]) < 1e-16)
            break;
        lt400 << x << " " << analytic[2000][i + 1] << " " << results[2000][i + 1] << endl;  // t = 1.25;
        x += h;
    }


    lt100.close();
    lt300.close();
    lt400.close();

    double *maxErrors = maxError(error, n, m);
    ofstream maxerror;
    maxerror.open("../LassTmaxError.txt");
    for (int i = 0; i < n - 1; i++) {
        if ((maxErrors[i + 1]) < 1e-16)
            break;
        maxerror << steps_t[i + 1] << " " << maxErrors[i + 1] << endl;
    }
    maxerror.close();

    saveToFile("analytic.txt", analytic, n, m);
    saveToFile("LassonenTNumerical.txt", results, n, m);
    saveToFile("LassonenTError.txt", error, n, m);


    delete[] maxErrors;
    delete[] steps_t;
    deleteMatrix(analytic, n);
    deleteMatrix(results, n);
    deleteMatrix(error, n);


}

//dane do wykresów Lassonen - metoda SOR
void LaasonenSORSave() {
    double delta_t = getDelta_t(lambda_laasonen);
    int n = (int) ((t_max - t_min) / delta_t) + 2;
    int m = (int) ((x_end - x_start) / h);

    double **analytic = analytical(n, m, delta_t);
    double **results = LaasonenSOR(n, m, delta_t);
    double **error = allocMatrix(n, m);

    countErrors(error, analytic, results, n, m);

    auto *steps_t = new double[n];
    double t = t_min;
    for (int i = 0; i < n; i++) {
        steps_t[i] = t;
        t += delta_t;
    }

    ofstream sor100, sor300, sor400;
    sor100.open("../sro100.txt");
    sor300.open("../sor300.txt");
    sor400.open("../sor400.txt");


    double x = x_start;
    cout << steps_t[200] << endl;
    for (int i = 0; i < m; i++) {
        if (abs(results[200][i + 1]) < 1e-16)
            break;
        sor100 << x << " " << analytic[200][i + 1] << " " << results[200][i + 1] << endl;  // t = 0.3125;
        x += h;
    }
    x = x_start;
    cout << steps_t[700] << endl;
    for (int i = 0; i < m; i++) {
        if (abs(results[700][i + 1]) < 1e-16)
            break;
        sor300 << x << " " << analytic[700][i + 1] << " " << results[700][i + 1] << endl;  // t = 0.75;
        x += h;
    }
    x = x_start;
    cout << steps_t[1000] << endl;
    for (int i = 0; i < m; i++) {
        if (abs(results[1000][i + 1]) < 1e-16)
            break;
        sor400 << x << " " << analytic[1000][i + 1] << " " << results[1000][i + 1] << endl;  // t = 1.25;
        x += h;
    }


    sor100.close();
    sor300.close();
    sor400.close();


    double *maxErrors = maxError(error, n, m);
    ofstream maxerror;
    maxerror.open("../LassSORmaxError.txt");
    for (int i = 0; i < n - 1; i++) {
        if ((maxErrors[i + 1]) < 1e-16)
            break;
        maxerror << steps_t[i + 1] << " " << maxErrors[i + 1] << endl;
    }
    maxerror.close();

    saveToFile("analytic.txt", analytic, n, m);
    saveToFile("LassonenSORNumerical.txt", results, n, m);
    saveToFile("LassonenSORError.txt", error, n, m);


    delete[] maxErrors;
    delete[] steps_t;
    deleteMatrix(analytic, n);
    deleteMatrix(results, n);
    deleteMatrix(error, n);
}


int main() {
    time_t begin, end;
    time(&begin);
    KMB_save();
    //LassonenThomasSave();
    //LaasonenSORSave();
    time(&end);
    time_t czas = end - begin;
    printf("Czas: %ld sek.", czas);
    return 0;
}
