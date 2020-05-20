/*
* File: testSimulatedAnnealing.cpp
* Author: Maksim Galynchik
*/

#include <cstdlib>
#include <iostream>
#include <cmath>
#include "SAstandart.hpp"

using namespace std;

template<class T> class RandomCandidate : public panther::NextCandidateDistribution<T> {
public:
    void nextCandidate(int n, T* point, const T* lowerBound, const T* upperBound, T delta) {
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_int_distribution<> dis(-1, 1);
        for (int i = 0; i < n; i++) {
            T num;
            do {
                num = point[i] + dis(gen) * (upperBound[i] - lowerBound[i]) * delta;
            } while (num >= upperBound[i] || num <= lowerBound[i]);
            point[i] = num;
        }
    }
};

template<class T> class Metropolis : public panther::AcceptanceFunction<T> {
public:
    bool acceptance(T state1, T state2, size_t t) const {
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<> dis(0.0, 1.0);
        return (state1 > state2) ? (1 >= dis(gen)) :
            (exp((state1 - state2) * pow(t, -1)) >= dis(gen));
    }
};

template<class T> class StandartCooling : public panther::CoolingSchedule<T> {
public:
    T coolingSchedule(size_t iteration) {
        return 1 / (exp(iteration));
    }
};

template<class T> class StandartStoping : public panther::StopingCriterion<T> {
public:
    bool stoping(size_t maxIter, size_t iter, T fOldPoint, T fNewPoint, T delta, size_t stopingIter = 1) {
        if (fOldPoint - fNewPoint <= delta && fOldPoint - fNewPoint >= 0) {
            nowIter += 1;
            if (stopingIter == nowIter) {
                return 0;
            }
        }
        else {
            nowIter = 0;
        }
        if (iter == maxIter) {
            std::cout << "Stoping with iteration limit" << std::endl;
            return 0;
        }
        return 1;
    }
private:
    size_t nowIter = 0;
};

