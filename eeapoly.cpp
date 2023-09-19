#include <iostream>
#include <vector>

using namespace std;

// Function to perform polynomial division
vector<int> polyDivide(vector<int> A, vector<int> B) {
    vector<int> quotient;
    int n = A.size();
    int m = B.size();

    while (n >= m) {
        int factor = A[n - 1] / B[m - 1];
        quotient.push_back(factor);

        for (int i = 0; i < m; i++) {
            A[n - i - 1] -= factor * B[m - i - 1];
        }

        n--;
    }

    return quotient;
}

// Function to calculate the GCD of two polynomials using the EEA
vector<int> polyGCD(vector<int> A, vector<int> B) {
    if (B.empty()) {
        return A;
    }

    vector<int> remainder = polyDivide(A, B);
    return polyGCD(B, remainder);
}

int main() {
    // Input polynomials represented as vectors of coefficients
    vector<int> A = {2,-3,0,1,5}; // f1(x) = 2x^4-3x^3+0x^2+1x^1+5x^0
    vector<int> B = {1, 2};       // B(x) = 1x^1+2x^0

    // Calculate the GCD of A and B
    vector<int> gcd = polyGCD(A, B);

    
    cout << "GCD of A and B is: ";
    int result_degree = gcd.size();
    for (int i = 0; i < gcd.size(); i++) {
        cout <<"(" <<gcd[i] << "x^"<< result_degree-1 << ")"<< " ";
        result_degree -=1;
    }
    cout << endl;

    return 0;
}
