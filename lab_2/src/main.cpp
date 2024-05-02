#include <iostream>
#include "integral.hpp"

int main()
{
    double a = 0, b = 0, error_rate;
    int processes = 0;
    std::cin >> processes >> a >> b >> error_rate;
    std::cout << my_integrate(processes, a, b, error_rate) << std::endl;
}