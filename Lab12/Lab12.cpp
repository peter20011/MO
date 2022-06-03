#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
using namespace std;

void exact_to_file(double, double, double);
void newton_to_file(vector<double>, vector<double>, double, double, double, string);
vector<double> get_Newton_c(vector<double>, vector<double>);
double horner_alg(vector<double>, vector<double>, double);
vector<double> get_keys(double, double, double);
vector <double> get_values(const vector<double>);
vector <double> generate_czebyszew(double, double, double);


int main() {
    double start = -1.0, finish = 1.0, interval = 0.2;
    exact_to_file(start, finish, 0.01);
 
    vector <double> keys = get_keys(start, finish, interval);
    vector <double> values = get_values(keys);
    vector<double> c = get_Newton_c(keys, values);
    newton_to_file(keys, c, start, finish, 0.01, "newton.dat");
    vector <double> keys_cheb = generate_czebyszew(start, finish, interval);
    vector <double> values_cheb = get_values(keys_cheb);
    vector <double> c1 = get_Newton_c(keys_cheb, values_cheb);
    newton_to_file(keys_cheb, c1, start, finish, 0.01, "czebyszew.dat");
    return 0;
}


vector<double> get_keys(double start, double finish, double interval) {
    vector<double> range;
    while(start <= finish){
        range.push_back(start);
        start += interval;
    }
    return range;
}

vector <double> get_values(const vector<double> keys) {
    vector <double> values;
    values.reserve(keys.size());
    for(double x: keys){
        values.push_back(1.0 / (1.0 + 10.0 * pow(x, 6)));
    }
    return values;
}

void exact_to_file(double start, double finish, double interval) {
    ofstream exact;
    exact.open("funkcja_dokladna.dat");
    vector<double> keys = get_keys(start, finish, interval);
    vector <double> values = get_values(keys);
    for(int i = 0; i < keys.size(); i++){
        exact << keys[i] << " " << values[i] << endl;
    }
}

vector<double> get_Newton_c(vector<double> keys, vector<double> values) {
    vector <double> c(values.size());
    c[0] = values[0];
    int iter = keys.size() - 1;
    int i = 1;
    while(i <= iter) {
        vector<double> second_vec(values.size() - 1);
        for(int j = 0; j < second_vec.size(); j++) {
            second_vec[j] = (values[j + 1] - values[j])/(keys[j+i] - keys[j]);
        }
        c[i] = second_vec[0];
        values = second_vec;
        i++;
    }
    return c;
}

void newton_to_file(vector<double> keys, vector<double> c, double start, double finish, double interval, string name) {
    ofstream newton;
    newton.open(name);
    while(start <= finish) {
        newton << start << " " << horner_alg(keys, c, start) << endl;
        start += interval;
    }
}

double horner_alg(vector<double> keys, vector<double> c, double x) {
    double cn = c[c.size() - 1];
    for(int i = 1; i < c.size(); i++) {
        cn = cn * (x - keys[keys.size() - 1 - i]) + c[c.size() - 1 - i];
    }
    return cn;
}

vector<double> generate_czebyszew(double start, double finish, double interval) {
    vector<double> keys;
    int size =int((finish - start) / interval) + 1;
    double i = 0.0;
    double c;
    while(i <= size) {
        c = (2.0 * i + 1.0) / (2.0 * size + 2.0);
        keys.push_back((finish + start) / 2.0 + (finish - start)*cos(c * M_PI) / 2.0);
        i += 1;
    }
    return keys;
}


