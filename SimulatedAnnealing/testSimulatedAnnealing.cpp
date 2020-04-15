/*
* File: testSimulatedAnnealing.cpp
* Author: Maksim Galynchik
*/

#include <cstdlib>
#include <iostream>
#include <cmath>
#include "SimulatedAnnealing.hpp"

using namespace std;

//input function(for example)
double func(const double* const x) {
    return /*(x[0] - 10)*(x[0] - 10)*/ /*5 * sin(2 * x[0]) + x[0] * x[0]*/ (x[0] * x[0] - 1) / (x[0] * x[0] + 1);
}

int main() {
    double num1[1], num2[1];
    num1[0] = -30;
    num2[0] = 70;
    double num[1];
    Option<double> opt;
    opt.lower_bound = num1, opt.upper_bound = num2, opt.E = func;
    simulated_Annealing<double>(opt, num);
    std::cout << about(opt);
    std::cout << "Found " << func(num) << " at " << num[0] << "\n";
    return 0;

}
