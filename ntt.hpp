#include <iostream>
#include <vector>
#include "common.hpp"
#include "ntt_constants.hpp"

const int i2 = 6145;
const int sqr1 = roots_dict_Zq[2][0];
int q = 12 * 1024 + 1;

std::vector<vector<int>> split_ntt(std::vector<int> f_ntt) {
    int n = f_ntt.size();
    std::vector<int> w = roots_dict_Zq[n];
    std::vector<int> f0_ntt(n / 2, 0);
    std::vector<int> f1_ntt(n / 2, 0);
    for (int i = 0; i < n / 2; i++) {
        f0_ntt[i] = (i2 * (f_ntt[2 * i] + f_ntt[2 * i + 1])) % q;
        f1_ntt[i] = (i2 * (f_ntt[2 * i] - f_ntt[2 * i + 1]) * inv_mod_q[w[2 * i]]) % q;
    }
    return {f0_ntt, f1_ntt};
}

std::vector<int> merge_ntt(std::vector<std::vector<int>> f_list_ntt) {
    std::vector<int> f0_ntt = f_list_ntt[0];
    std::vector<int> f1_ntt = f_list_ntt[1];
    int n = 2 * f0_ntt.size();
    std::vector<int> w = roots_dict_Zq[n];
    std::vector<int> f_ntt(n, 0);
    for (int i = 0; i < n / 2; i++) {
        f_ntt[2 * i + 0] = (f0_ntt[i] + w[2 * i] * f1_ntt[i]) % q;
        f_ntt[2 * i + 1] = (f0_ntt[i] - w[2 * i] * f1_ntt[i]) % q;
    }
    return f_ntt;
}

std::vector<int> ntt(std::vector<int> f) {
    int n = f.size();
    if (n > 2) {
        std::vector<int> f0 = std::vector<int>(f.begin(), f.begin() + n / 2);
        std::vector<int> f1 = std::vector<int>(f.begin() + n / 2, f.end());
        std::vector<int> f0_ntt = ntt(f0);
        std::vector<int> f1_ntt = ntt(f1);
        return merge_ntt({f0_ntt, f1_ntt});
    } else if (n == 2) {
        std::vector<int> f_ntt(2, 0);
        f_ntt[0] = (f[0] + sqr1 * f[1]) % q;
        f_ntt[1] = (f[0] - sqr1 * f[1]) % q;
        return f_ntt;
    }
    return f;
}

std::vector<int> intt(std::vector<int> f_ntt) {
    int n = f_ntt.size();
    if (n > 2) {
        std::vector<int> f0_ntt, f1_ntt;
        std: vector< vector<int>> hello = split(f_ntt);
        f0_ntt = hello[0];
        f1_ntt = hello[1];
        std::vector<int> f0 = intt(f0_ntt);
        std::vector<int> f1 = intt(f1_ntt);
        return merge({f0, f1});
    } else if (n == 2) {
        std::vector<int> f(2, 0);
        f[0] = (i2 * (f_ntt[0] + f_ntt[1])) % q;
        f[1] = (i2 * inv_mod_q[1479] * (f_ntt[0] - f_ntt[1])) % q;
        return f;
    }
    return f_ntt;
}

std::vector<int> add_zq(std::vector<int> f, std::vector<int> g) {
    assert(f.size() == g.size());
    int deg = f.size();
    std::vector<int> result(deg);
    for (int i = 0; i < deg; i++) {
        result[i] = (f[i] + g[i]) % q;
    }
    return result;
}

std::vector<int> neg_zq(std::vector<int> f) {
    int deg = f.size();
    std::vector<int> result(deg);
    for (int i = 0; i < deg; i++) {
        result[i] = (-f[i]) % q;
    }
    return result;
}

std::vector<int> sub_zq(std::vector<int> f, std::vector<int> g) {
    return add_zq(f, neg_zq(g));
}



std::vector<int> add_ntt(std::vector<int> f_ntt, std::vector<int> g_ntt) {
    return add_zq(f_ntt, g_ntt);
}

std::vector<int> sub_ntt(std::vector<int> f_ntt, std::vector<int> g_ntt) {
    return sub_zq(f_ntt, g_ntt);
}

std::vector<int> mul_ntt(std::vector<int> f_ntt, std::vector<int> g_ntt) {
    assert(f_ntt.size() == g_ntt.size());
    int deg = f_ntt.size();
    std::vector<int> result(deg);
    for (int i = 0; i < deg; i++) {
        result[i] = (f_ntt[i] * g_ntt[i]) % q;
    }
    return result;
}

std::vector<int> div_ntt(std::vector<int> f_ntt, std::vector<int> g_ntt) {
    assert(f_ntt.size() == g_ntt.size());
    int deg = f_ntt.size();
    if (std::any_of(g_ntt.begin(), g_ntt.end(), [](int elt) { return elt == 0; })) {
        throw std::runtime_error("Division by zero");
    }
    std::vector<int> result(deg);
    for (int i = 0; i < deg; i++) {
        result[i] = (f_ntt[i] * inv_mod_q[g_ntt[i]]) % q;
    }
    return result;
}
std::vector<int> mul_zq(std::vector<int> f, std::vector<int> g) {
    return intt(mul_ntt(ntt(f), ntt(g)));
}


std::vector<int> div_zq(std::vector<int> f, std::vector<int> g) {
    try {
        return intt(div_ntt(ntt(f), ntt(g)));
    } catch (std::exception& e) {
        throw;
    }
}

double ntt_ratio = 1;
