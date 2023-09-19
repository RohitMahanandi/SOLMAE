#include <iostream>
using namespace std;

int EEA(int x, int y,int t1,int t2) {
    int q = x/y;
    int r = x%y;
    int t = t1-(t2*q);
    int answer;
    if (r == 0) return t2;
    else return EEA(y,r,t2,t);

}
int main() {
    int a,b;
    cout << "Remember that a>b" << endl;
    cout << "Enter the value of a: ";
    cin >> a;
    cout << "Enter the value of b: ";
    cin >> b;
    int answer = EEA(a,b,0,1);
    if (answer < 0) answer = answer+a;
    cout<< endl << "multiplicative inverse of a mod b: "<< answer << endl;



}