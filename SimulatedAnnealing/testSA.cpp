/*
* File: testSimulatedAnnealing.cpp
* Author: Maksim Galynchik
*/

#include <cstdlib>
#include <iostream>
#include <cmath>
#include <ctime>
//#include "SAstandart.hpp"
//#include "SAparallel1.hpp"
#include "SAparallel2.hpp"
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
    StandartCooling<double> Temp;
    StandartStoping<double> Stop;
    double num1[1], num2[1];
    num1[0] = -30;
    num2[0] = 70;
    double num[1];
    std::fill(num, num + 1, 2);
    srand(time(0));
    SimulatedAnnealing<double> SA(D, A, Temp, Stop);
    std::cout << SA.about();
    std::cout << SA.search(1, num, num1, num2, func) << std::endl;
    std::cout << "runtime = " << clock() / 1000.0 << std::endl;
    std::cout << num[0] << std::endl;
    return 0;

}
