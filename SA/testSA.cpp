#include <iostream>
#include <cstdlib>
#include "SimulatedAnnealing.hpp"


double func(const double * const x) {
    return /*(x[0] - 10)*(x[0] - 10)*/ /*5 * sin(2 * x[0]) + x[0] * x[0]*/ (x[0]*x[0] - 1) / (x[0]*x[0] + 1);
}

int main() {
    double num1[1], num2[1];
    num1[0] = -700;
    num2[0] = 700;
    double res[1];
    simulated_Annealing<double>(1, num1, num2, res, func, 0.1);
    std::cout << res[0] << "|-----|" << func(res);
    return 0;

}
