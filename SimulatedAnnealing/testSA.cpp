/*
* File: testSimulatedAnnealing.cpp
* Author: Maksim Galynchik
*/

#include <cstdlib>
#include <iostream>
#include <cmath>
#include "SimulatedAnnealing.hpp"
#include "SAComponents.cpp"

using namespace std;

//input function(for example)
double func(const double* const x) {
    return /*(x[0] - 10)*(x[0] - 10)*/ /*5 * sin(2 * x[0]) + x[0] * x[0]*/ (x[0] * x[0] - 1) / (x[0] * x[0] + 1);
}

using namespace panther;

int main() {
    RandomCandidate<double> D;
    Metropolis<double> A;
    Cooling<double> Temp;
    double num1[1], num2[1];
    num1[0] = -30;
    num2[0] = 70;
    double num[1];
    std::fill(num, num + 1, 20);
    SimulatedAnnealing<double> SA(D, A, Temp);
    std::cout << SA.about();
    std::cout << SA.search(1, num, num1, num2, func);
    return 0;
}
