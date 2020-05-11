/*
* File: testSimulatedAnnealing.cpp
* Author: Maksim Galynchik
*/

#include <cstdlib>
#include <iostream>
#include <cmath>
#include "SimulatedAnnealing.hpp"

using namespace std;

template<class T> class RandomCandidate : public panther::NextCandidateDistribution<T> {
public:
    void nextCandidate(int n, T* point, const T* lower_bound, const T* upper_bound, T delta) {
        for (int i = 0; i < n; i++) {
            T num;
            do {
                num = point[i] + (2 * ((double)rand() / (RAND_MAX)) - 1) * (upper_bound[i] - lower_bound[i]) * delta;
            } while (num >= upper_bound[i] || num <= lower_bound[i]);
            point[i] = num;
        }
    }
};

template<class T> class Metropolis : public panther::AcceptanceFunction<T> {
public:
    bool acceptance(const T state1, const T state2, const size_t t) {
        return (state1 > state2) ? (1 >= (double)rand() / (RAND_MAX)) :
               (exp((state1 - state2) * pow(t, -1)) >= (double)rand() / (RAND_MAX));
    }
};

template<class T> class Cooling : public panther::CoolingSchedule<T> {
public:
    T coolingSchedule(const size_t max_k, size_t iteration) {
        return max_k - iteration;
    }
};

