

#include <iostream>
#include <cmath>
#include <random>


// Use high-quality randomness
std::random_device rd;
std::mt19937 gen(rd());

// Upper bound on all the values of sigma
const double MAX_SIGMA = 1.8205;
const double INV_2SIGMA2 = 1 / (2 * pow(MAX_SIGMA, 2));

// Precision of RCDT
const int RCDT_PREC = 72;

// ln(2) and 1 / ln(2), with ln the natural logarithm
const double LN2 = 0.69314718056;
const double ILN2 = 1.44269504089;

// RCDT is the reverse cumulative distribution table of a distribution that
// is very close to a half-Gaussian of parameter MAX_SIGMA.
unsigned long long RCDT[] = {
    3024686241123004913666,
    1564742784480091954050,
    636254429462080897535,
    199560484645026482916,
    47667343854657281903,
    8595902006365044063,
    1163297957344668388,
    117656387352093658,
    8867391802663976,
    496969357462633,
    20680885154299,
    638331848991,
    14602316184,
    247426747,
    3104126,
    28824,
    198,
    1};

// C contains the coefficients of a polynomial that approximates exp(-x)
// More precisely, the value:
// (2 ** -63) * sum(C[12 - i] * (x ** i) for i in range(i))
// Should be very close to exp(-x).
// This polynomial is lifted from FACCT: https://doi.org/10.1109/TC.2019.2940949
unsigned long long C[] = {
    19127174051,
    233346759686,
    2542029181962,
    25415798087749,
    228754078003076,
    1830034511206115,
    0x002D82D8305B0FEA,
    0x011111110E066FD0,
    0x0555555555070F00,
    0x155555555581FF00,
    0x400000000002B400,
    0x7FFFFFFFFFFF4800,
    9223372036854775808};

// Sample z0 in {0, 1, ..., 18} with a distribution
// very close to the half-Gaussian D_{Z+, 0, MAX_SIGMA}.
int basesampler() {
    std::uniform_int_distribution<> dis(0, RCDT_PREC >> 3);
    int u = dis(gen);

    int z0 = 0;
    for (auto elt : RCDT) {
        z0 += int(u < elt);
    }
    return z0;
}

long long approxexp(double x, double ccs) {
    unsigned long long y = C[0];
    unsigned long long z = int(x * (1 << 63));
    for (int i = 1; i < sizeof(C) / sizeof(C[0]); i++) {
        y = C[i] - ((z * y) >> 63);
    }
    z = int(ccs * (1 << 63)) << 1;
    y = (z * y) >> 63;
    return y;
}

bool berexp(double x, double ccs) {
    const double ILN2 = 1.4426950408889634;
    const double LN2 = 0.6931471805599453;
    
    int s = int(x * ILN2);
    double r = x - s * LN2;
    s = std::min(s, 63);
    double z = (std::exp(r) - 1) / std::pow(2, s);
    
    std::uniform_int_distribution<int> dist(0, 255);
    std::mt19937 gen(rd());
    
    for (int i = 56; i >= -8; i -= 8) {
        int p = dist(gen);
        int w = p - ((int)(z / std::pow(2, i)) & 0xFF);
        if (w != 0) {
            return w < 0;
        }
    }
    
    return false;
}

int samplerz(double mu, double sigma, double sigmin) {
    int s = int(std::floor(mu));
    double r = mu - s;
    double dss = 1 / (2 * sigma * sigma);
    double ccs = sigmin / sigma;

    while (true) {
        // Sampler z0 from a Half-Gaussian
        int z0 = basesampler();
        // Convert z0 into a pseudo-Gaussian sample z
        std::uniform_int_distribution<> dis(0, 1);
        int b = dis(gen);
        b &= 1;
        int z = b + (2 * b - 1) * z0;
        // Rejection sampling to obtain a true Gaussian sample
        double x = pow((z - r), 2) * dss;
        x -= pow(z0, 2) * INV_2SIGMA2;
        if (berexp(x, ccs)) {
            return z + s;
        }
    }
}