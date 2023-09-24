#include <bits/stdc++.h>
#include <vector>
#include <cmath>
#include <complex>
#include "fft.hpp"
#include "samplerz.hpp"
#include "ntt.hpp"


//  int q = 12 * 1024 + 1;

typedef std::vector<int> Polynomial;
typedef std::vector<std::complex<double>> ComplexPolynomial;

// FFT, IFFT, and other FFT-related functions go here (you'll need to implement them).

Polynomial karatsuba(const Polynomial &a, const Polynomial &b, int n) {
    if (n == 1) {
        return {a[0] * b[0], 0};
    } else {
        int n2 = n / 2;
        Polynomial a0(a.begin(), a.begin() + n2);
        Polynomial a1(a.begin() + n2, a.end());
        Polynomial b0(b.begin(), b.begin() + n2);
        Polynomial b1(b.begin() + n2, b.end());

        Polynomial ax(n2);
        Polynomial bx(n2);
        for (int i = 0; i < n2; i++) {
            ax[i] = a0[i] + a1[i];
            bx[i] = b0[i] + b1[i];
        }

        Polynomial a0b0 = karatsuba(a0, b0, n2);
        Polynomial a1b1 = karatsuba(a1, b1, n2);
        Polynomial axbx = karatsuba(ax, bx, n2);

        for (int i = 0; i < n; i++) {
            axbx[i] -= (a0b0[i] + a1b1[i]);
        }

        Polynomial ab(2 * n, 0);
        for (int i = 0; i < n; i++) {
            ab[i] += a0b0[i];
            ab[i + n] += a1b1[i];
            ab[i + n2] += axbx[i];
        }

        return ab;
    }
}

Polynomial karamul(const Polynomial &a, const Polynomial &b) {
    int n = a.size();
    Polynomial ab = karatsuba(a, b, n);
    Polynomial abr(n);
    for (int i = 0; i < n; i++) {
        abr[i] = ab[i] - ab[i + n];
    }
    return abr;
}

Polynomial galois_conjugate(const Polynomial &a) {
    int n = a.size();
    Polynomial result(n);
    for (int i = 0; i < n; i++) {
        result[i] = ((i % 2 == 0) ? 1 : -1) * a[i];
    }
    return result;
}

Polynomial field_norm(const Polynomial &a) {
    int n2 = a.size() / 2;
    Polynomial ae(n2);
    Polynomial ao(n2);
    for (int i = 0; i < n2; i++) {
        ae[i] = a[2 * i];
        ao[i] = a[2 * i + 1];
    }

    Polynomial ae_squared = karamul(ae, ae);
    Polynomial ao_squared = karamul(ao, ao);
    Polynomial res = ae_squared;

    for (int i = 0; i < n2 - 1; i++) {
        res[i + 1] -= ao_squared[i];
    }
    res[0] += ao_squared[n2 - 1];

    return res;
}

Polynomial lift(const Polynomial &a) {
    int n = a.size();
    Polynomial res(2 * n, 0);
    for (int i = 0; i < n; i++) {
        res[2 * i] = a[i];
    }
    return res;
}

int bitsize(int a) {
    int val = std::abs(a);
    int res = 0;
    while (val) {
        res += 8;
        val >>= 8;
    }
    return res;
}

// Implement FFT, IFFT, and other FFT-related functions here.

Polynomial reduce(const Polynomial &f, const Polynomial &g, Polynomial F, Polynomial G) {
    int n = f.size();
    int size = std::max(53, std::max({bitsize(*std::min_element(f.begin(), f.end())),
                                      bitsize(*std::max_element(f.begin(), f.end())),
                                      bitsize(*std::min_element(g.begin(), g.end())),
                                      bitsize(*std::max_element(g.begin(), g.end()))}));

    ComplexPolynomial f_adjust(n);
    ComplexPolynomial g_adjust(n);
    for (int i = 0; i < n; i++) {
        f_adjust[i] = f[i] >> (size - 53);
        g_adjust[i] = g[i] >> (size - 53);
    }

    ComplexPolynomial fa_fft = fft(f_adjust);
    ComplexPolynomial ga_fft = fft(g_adjust);

    while (true) {
        int Size = std::max(53, std::max({bitsize(*std::min_element(F.begin(), F.end())),
                                          bitsize(*std::max_element(F.begin(), F.end())),
                                          bitsize(*std::min_element(G.begin(), G.end())),
                                          bitsize(*std::max_element(G.begin(), G.end()))}));
        if (Size < size) {
            break;
        }

        ComplexPolynomial F_adjust(n);
        ComplexPolynomial G_adjust(n);
        for (int i = 0; i < n; i++) {
            F_adjust[i] = F[i] >> (Size - 53);
            G_adjust[i] = G[i] >> (Size - 53);
        }

        ComplexPolynomial Fa_fft = fft(F_adjust);
        ComplexPolynomial Ga_fft = fft(G_adjust);

        ComplexPolynomial den_fft = add_fft(mul_fft(fa_fft, adj_fft(fa_fft)), mul_fft(ga_fft, adj_fft(ga_fft)));
        ComplexPolynomial num_fft = add_fft(mul_fft(Fa_fft, adj_fft(fa_fft)), mul_fft(Ga_fft, adj_fft(ga_fft)));

        ComplexPolynomial k_fft = div_fft(num_fft, den_fft);
        ComplexPolynomial k = ifft(k_fft);
        Polynomial k2;
        for (int i = 0; i < n; i++) {
            k2[i] = static_cast<int>(round(real(k[i])));
        }


        if (std::all_of(k.begin(), k.end(), [](int elt) { return elt == 0; })) {
            break;
        }

        Polynomial fk = karamul(f, k2);
        Polynomial gk = karamul(g, k2);

        for (int i = 0; i < n; i++) {
            F[i] -= fk[i] << (Size - size);
            G[i] -= gk[i] << (Size - size);
        }
    }
    return F;
}

std::tuple<int, int, int> xgcd(int b, int n) {
    int x0 = 1, x1 = 0, y0 = 0, y1 = 1;
    while (n != 0) {
        int q = b / n;
        std::tie(b, n) = std::make_tuple(n, b % n);
        std::tie(x0, x1) = std::make_tuple(x1, x0 - q * x1);
        std::tie(y0, y1) = std::make_tuple(y1, y0 - q * y1);
    }
    return std::make_tuple(b, x0, y0);
}

Polynomial ntru_solve(const Polynomial &f, const Polynomial &g) {
    int n = f.size();
    if (n == 1) {
        int f0 = f[0];
        int g0 = g[0];
        auto [d, u, v] = xgcd(f0, g0);
        if (d != 1) {
            throw std::runtime_error("xgcd failed");
        } else {
            return {-q * v};
        }
    } else {
        Polynomial fp = field_norm(f);
        Polynomial gp = field_norm(g);
        Polynomial Fp = ntru_solve(fp, gp);
        Polynomial F = karamul(lift(Fp), galois_conjugate(g));
        Polynomial Gp = ntru_solve(gp, fp);
        Polynomial G = karamul(lift(Gp), galois_conjugate(f));
        F = reduce(f, g, F, G);
        return F;
    }
}

double sqnorm(const Polynomial &poly) {
    double result = 0.0;
    for (int coeff : poly) {
        result += coeff * coeff;
    }
    return result;
}

double gs_norm(const std::vector<Polynomial> &matrix) {
    double result = 0.0;
    for (const Polynomial &poly : matrix) {
        result += sqnorm(poly);
    }
    return result;
}

Polynomial gen_poly(int n) {
    double sigma = 1.43300980528773;
    assert(n < 4096);
    std::vector<double> f0(4096);

    for (int i = 0; i < 4096; i++) {
        f0[i] = samplerz(0, sigma, sigma - 0.001);  // Implement samplerz function
    }

    Polynomial f(n, 0);
    int k = 4096 / n;

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < k; j++) {
            f[i] += f0[i * k + j];
        }
    }

    return f;
}

Polynomial ntru_gen(int n) {
    while (true) {
        Polynomial f = gen_poly(n);
        Polynomial g = gen_poly(n);
        std:: vector <Polynomial> matrix;
        matrix.push_back(f);
        matrix.push_back(g);
        //matrix.push_back(q);
        vector<int>f_ntt = ntt(f);
        if (gs_norm(matrix) > (1.17 * 1.17 * q)) {
            continue;
        }

        // Implement ntt function here for f_ntt calculation

        if (std::any_of(f_ntt.begin(), f_ntt.end(), [](const ComplexPolynomial &elem) {
                return std::any_of(elem.begin(), elem.end(), [](const std::complex<double> &c) {
                    return std::abs(c) < 1e-9;
                });
            })) {
            continue;
        }

        try {
            Polynomial F = ntru_solve(f, g);
            // Implement galois_conjugate function here
            Polynomial G = ntru_solve(galois_conjugate(f), galois_conjugate(g));
            F = reduce(f, g, F, G);
            return F;
        } catch (const std::runtime_error &e) {
            continue;
        }
    }
}

int main() {
    int n = 1024;  // Adjust n as needed
    Polynomial f = ntru_gen(n);
    for (int coeff : f) {
        std::cout << coeff << " ";
    }
    return 0;
}
