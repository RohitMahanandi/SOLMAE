#include <bits/stdc++.h>
using namespace std;
double assign() {
    const long max_rand = 1000000L;
 
    double lower_bound = 0;
    double upper_bound = 1;
    double k = lower_bound+(upper_bound - lower_bound) * (random() % max_rand) / max_rand;
    //cout << "k: " << k << endl;
    return k;
}

vector <double> unifcrown() {
    double R_min, R_max;
    cout << "Enter the value of R-: " ;
    cin >> R_min;
    cout << endl << "Enter the value of R+: ";
    cin >> R_max;
    double mu_x,mu_y;
    double mu_row;

    srandom(time(NULL));
 
    mu_row = assign();
    mu_x = assign();
    mu_y = assign();
    //cout << "mu_row = " << mu_row << endl;
    //cout << "mu_x = " << mu_x << endl;
    //cout << "mu_y = " << mu_y << endl;

    
    double pi_value =  acos(0.0);
    //cout << "pi_value = " << pi_value << endl;
    
    double inside = pow(R_min,2) + mu_row*(pow(R_max,2)-pow(R_min,2));
    double rhow = sqrt(inside);

    //cout << "rhow: " << rhow << endl;
    double x = rhow*cos(pi_value*mu_x);
    double y = rhow*sin(pi_value*mu_y);
    
    //cout << "x: " << x << endl<< "y: " << y << endl;
    return {x,y};
    
    
   
    //cout << "Hello world" << endl;
}