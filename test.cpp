#include <iostream>
#include <limits>

int main()
{   
     std::cout << "ull\t"
              << std::numeric_limits<unsigned long long>::lowest() << '\t'
              << std::numeric_limits<unsigned long long>::max() << std::endl ;
}