#include <bits/stdc++.h>
#include <cmath>
#include<complex>
#include "fft.hpp"
using namespace std;
typedef complex<double> Complex;


//params start

const int SOLMAE_D = 512;
const int SOLMAE_Q = 12289;

// Logarithmic values
std::map<int, int> logn = {
    {2, 1},
    {4, 2},
    {8, 3},
    {16, 4},
    {32, 5},
    {64, 6},
    {128, 7},
    {256, 8},
    {512, 9},
    {1024, 10}
};

// Parameters
const int HEAD_LEN = 1;
const int SALT_LEN = 320;
const int SEED_LEN = 56;

// Parameter sets for Solmae
struct SolmaeParams {
    int d;
    double sigma;
    double sigmin;
    long long sig_bound;
    int sig_bytelen;
    double smoothing;
    double quality;
    double correction;
    double lower_radius;
    double upper_radius;
    double signature_width;
    double slack;
};

std::map<int, SolmaeParams> Params = {
    // SolmaeParam(512, 128)
    {512, {
        512,
        1.32,
        1.2778336969128337,
        33870790,
        666,
        1.338,
        1.17,
        0.065,
        101.95,
        122.49,
        173.54,
        1.04
    }},
    // SolmaeParam(1024, 256)
    {1024, {
        1024,
        1.32,
        1.298280334344292,
        134150669,
        1375,
        1.351,
        1.64,
        0.3,
        100.85,
        148.54,
        245.62,
        1.04
    }}
};

// params end



struct UnifcrownResult {
    double x;
    double y;
};

UnifcrownResult Unifcrown(double R_min, double R_max) {
    std::random_device rd;
    std::default_random_engine generator(rd());
    std::uniform_real_distribution<double> distribution(0.0, 1.0);

    double u_rho = distribution(generator);
    double u_theta = distribution(generator);
    
    double rho = sqrt(R_min * R_min + u_rho * (R_max * R_max - R_min * R_min));
    double x = rho * cos(M_PI / 2 * u_theta);
    double y = rho * sin(M_PI / 2 * u_theta);

    UnifcrownResult result;
    result.x = x;
    result.y = y;
    return result;
}



pair<vector<int>, vector<int>> pair_generation() {
    double R_min = Params[SOLMAE_D].lower_radius;
    double R_max = Params[SOLMAE_D].upper_radius;
    double quality = Params[SOLMAE_D].quality;

    while(true)
    {
        bool flag = true;
        vector<double> x_list;
    vector<double> y_list;

    vector<complex<double>> f_fft(SOLMAE_D/2);
    vector<complex<double>> g_fft(SOLMAE_D/2);


    // Call the Unifcrown function from UnifCrown.cpp
   for(int i = 0; i < SOLMAE_D/2; i++)
    {
        UnifcrownResult result = Unifcrown(R_min, R_max);
        x_list.push_back(result.x);
        y_list.push_back(result.y);
        //cout << "X: "<< result.x <<" Y:"<< result.y << endl ;
        UnifcrownResult theta = Unifcrown(0,1);
        double phi_f_i_r= (result.x) * cos(theta.x);
        double phi_f_i_i = (result.x) * sin(theta.x);

        f_fft[i] = complex<double>(phi_f_i_r, phi_f_i_i);

         double phi_g_i_r= (result.y) * cos(theta.y);
        double phi_g_i_i = (result.y) * sin(theta.y);

        g_fft[i] = complex<double>(phi_g_i_r, phi_g_i_i);
        

    }

    vector<int> fR, gR;
    vector<Complex> res_f_fft = fft(f_fft);
    vector<Complex> res_g_fft = fft(g_fft);
        for (int i = 0; i < SOLMAE_D / 2; i++) {
            double norm_sq = res_f_fft[i].real() * res_f_fft[i].real() + res_f_fft[i].imag() * res_f_fft[i].imag() +
                             res_g_fft[i].real() * res_g_fft[i].real() + res_g_fft[i].imag() * res_g_fft[i].imag();

            double quality_sq = quality * quality;
            if (norm_sq < SOLMAE_Q / quality_sq || norm_sq > SOLMAE_Q * quality_sq) {
                flag = false;
                continue;
            }

            // Convert the complex FFT results to integers and append to fR and gR
            fR.push_back(static_cast<int>(round(res_f_fft[i].real())));
            gR.push_back(static_cast<int>(round(res_g_fft[i].real())));
        }
        cout << "fR vector: ";
        for (int i = 0; i < fR.size(); i++) {
            cout << fR[i] << " ";
        }
        cout << endl;

        // Print the gR vector
        cout << "gR vector: ";
        for (int i = 0; i < gR.size(); i++) {
            cout << gR[i] << " ";
        }
        cout << endl;
        if (flag) {
            return make_pair(fR, gR);
        }

    }
    ;
   // cout << f.size() << " " << g.size();n
}



int main() {
    // ... Your other code ...

    pair<vector<int>, vector<int>> result = pair_generation();

    // Print the fR vector
    cout << "fR vector: ";
    for (int i = 0; i < result.first.size(); i++) {
        cout << result.first[i] << " ";
    }
    cout << endl;

    // Print the gR vector
    cout << "gR vector: ";
    for (int i = 0; i < result.second.size(); i++) {
        cout << result.second[i] << " ";
    }
    cout << endl;

    return 0;
}
