/*
 * File: sannecomponents.hpp
 * Author: Maksim Galynchik
 */

#ifndef SANNECOMPONENTS_HPP
#define SANNECOMPONENTS_HPP

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <random>
#include "sannestand.hpp"
#include "sanneutilities.hpp"
#include <string>

using namespace std;

template<class T> class RandomCandidate : public panther::NextCandidateDistribution<T> {
public:
    RandomCandidate(T step = 0.025): delta(step) {}
    void nextCandidate(int n, T* point, const T* lowerBound, const T* upperBound) override {
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_int_distribution<> dis(-1., 1.);
        for (unsigned int i = 0; i < n; i++) {
            T num;
            do {
                num = point[i] + dis(gen) * (upperBound[i] - lowerBound[i]) * delta;
            } while (num >= upperBound[i] || num <= lowerBound[i]);
            point[i] = num;
        }
    }

    std::string about() override {
        std::ostringstream options;
        options << "NextCandidateDistribution : RandomCandidate\n";
        options << "options:\n";
        options << "Step of the NextCandidateDistribution proportional to the function domain " << delta << "\n";
        return options.str();
    }
private:
   /**
    * Step of the NextCandidateDistribution proportional to the function domain
    */
    T delta;
};

template<class T> class Metropolis : public panther::AcceptanceFunction<T> {
public:
    bool acceptance(T state1, T state2, unsigned int t) const override {
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<> dis(0.0, 1.0);
        return (state1 > state2) ? true:
            (exp((state1 - state2) * pow(t, -1)) >= dis(gen));
    }

    std::string about() override {
        std::ostringstream options;
        options << "AcceptanceFunction : Metropolis\n";
        return options.str();
    }
};

template<class T> class StandartCooling : public panther::CoolingSchedule<T> {
public:
    T coolingSchedule(unsigned int iteration) const override {
        return 1. / (exp(iteration));
    }

    std::string about() override {
        std::ostringstream options;
        options << "CoolingSchedule : StandartCooling\n";
        return options.str();
    }
};

template<class T> class StandartStoping : public panther::StopingCriterion<T> {
public:
    StandartStoping(unsigned int maxI = 3000, T acc = 0.01, unsigned int stopIt = 1) : maxIter(maxI), accuracy(acc), stopingIter(stopIt) {}
    bool stoping(unsigned int iter, T fOldPoint, T fNewPoint) override {
        if (fOldPoint - fNewPoint <= accuracy && fOldPoint - fNewPoint >= 0) {
            nowIter += 1;
            if (stopingIter == nowIter) {
                stopInfo += "Stoping with low value difference on iteration " + std::to_string(iter) + "\n";
                return 0;
            }
        }
        else {
            nowIter = 0;
        }
        if (iter == maxIter) {
            stopInfo += "Stoping with iteration limit.\n";
            return 0;
        }
        return 1;
    }

    std::string about() override {
        std::ostringstream options;
        options << "StopingCriterion : StandartStoping\n";
        options << "options:\n";
        options << "number of steps = " << maxIter << "\n";
        options << "Minimum difference between the value of the previous point and the new " << accuracy << "\n";
        options << "Maximum number of iterations at which the difference between points is less than accuracy " << stopingIter << "\n";
        return options.str();
    }

    std::string aboutStoping() override {
        std::ostringstream options;
        options << stopInfo;
        return options.str();
    }
private:
    /**
     * Curent iteration number with difference lower accuracy
     */
    unsigned int nowIter = 0;
    /**
     * Maximum number of iterations of Simulated Annealing
     */
    unsigned int maxIter;
    /**
     * Minimum difference between the value of the previous point and the new
     */
    T accuracy;
    /**
     * Maximum number of iterations at which the difference between points is less than accuracy
     */
    unsigned int stopingIter;
    std::string stopInfo;
};

#endif
