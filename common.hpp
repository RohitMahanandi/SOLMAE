#include <vector>

// q is the integer modulus which is used in Falcon.
//const int q = 12 * 1024 + 1;

std::vector<std::vector<int>> split(const std::vector<int>& f) {
    // Split a polynomial f into two polynomials.

    int n = f.size();
    std::vector<int> f0(n / 2);
    std::vector<int> f1(n / 2);
    for (int i = 0; i < n / 2; i++) {
        f0[i] = f[2 * i];
        f1[i] = f[2 * i + 1];
    }
    return {f0, f1};
}

std::vector<int> merge(const std::vector<std::vector<int>>& f_list) {
    // Merge two polynomials into a single polynomial f.

    const std::vector<int>& f0 = f_list[0];
    const std::vector<int>& f1 = f_list[1];
    int n = 2 * f0.size();
    std::vector<int> f(n);
    for (int i = 0; i < n / 2; i++) {
        f[2 * i] = f0[i];
        f[2 * i + 1] = f1[i];
    }
    return f;
}

int sqnorm(const std::vector<std::vector<int>>& v) {
    // Compute the square euclidean norm of the vector v.

    int res = 0;
    for (const std::vector<int>& elt : v) {
        for (int coef : elt) {
            res += coef * coef;
        }
    }
    return res;
}

