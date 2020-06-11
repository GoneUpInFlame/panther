/*
 * File: testSimulatedAnnealing.cpp
 * Author: Maksim Galynchik
 */

#include <cmath>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include "sannestand.hpp"
#include "sannecomponents.hpp"

using namespace std;

//input function(for example)
double func(const double* const x) {
    return /*(x[0] - 10)*(x[0] - 10)*/ 5 * sin(2 * x[0]) + x[0] * x[0] /*(x[0] * x[0] - 1) / (x[0] * x[0] + 1)*/;
}

using namespace panther;

int main() {
    const int n = 1;
    double lowerBound[n], upperBound[n], startPoint[n];
    std::fill(lowerBound, lowerBound + n, -3);
    std::fill(upperBound, upperBound + n, 7);
    std::fill(startPoint, startPoint + n, 2);
    srand(time(0));

    double delta = 0.025;
    RandomCandidate<double> D(delta);

    Metropolis<double> A;

    StandartCooling<double> Temp;

    unsigned int maxIter = 3000;
    double accuracy = 0.01;
    unsigned int stopingIter = 10;
    StandartStoping<double> Stop(maxIter, accuracy, stopingIter);

    StandartSimulatedAnnealing<double> SA(D, A, Temp, Stop);
    std::cout << SA.about();
    std::cout << SA.search(1, startPoint, lowerBound, upperBound, func) << std::endl;
    std::cout << Stop.aboutStoping();

    std::cout << "runtime = " << clock() / 1000.0 << std::endl;
    std::cout << startPoint[0] << std::endl;
    return 0;

}
