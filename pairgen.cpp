#include <bits/stdc++.h>
#include "unifcrown.hpp"
using namespace std;

double taking_input() {
    double value;
    cout<< endl << "Enter the value of R-: " ;
    cin >> value;
    
}
int main() {
    double R_min = taking_input();
    double R_max = taking_input();
   
    for (int i = 0; i < d/2; i++) {
        vector <double> wow = unifcrown(R_min,R_max);
        double theta_x = assign();
        double theta_y = assign();
        double x_re = cos(2*acos(0.0)*theta_x) * wow[0];
        double x_im = sin(2*acos(0.0)*theta_x) * wow[0];
        double y_re = cos(2*acos(0.0)*theta_y) * wow[1];
        double y_im = sin(2*acos(0.0)*theta_y) * wow[1];

        double comp1 = comp(x_re,x_im);
        double comp2 = comp(y_re,y_im);
        double phi_x = abs(wow[0]) *comp1;
        double phi_y = abs(wow[1]) * comp2;


    cout << "Hello world" << endl;
}
