/*
 * File:   SimulatedAnnealing.hpp
 * Author: Maksim Galynchik
 */

#ifndef SIMULATEDANNEALING_HPP
#define  SIMULATEDANNEALING_HPP

#include <sstream>
#include <common/bbsolver.hpp>
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <omp.h>
#include <functional>
#include <algorithm>
#include <random>


namespace panther {
    template<class T> class NextCandidateDistribution {
    public:
        virtual void nextCandidate(int n, T* point, const T* lowerBound, const T* upperBound, T delta) = 0;
    };

    template<class T> class CoolingSchedule {
    public:
        virtual T coolingSchedule(size_t iteration) = 0;
    };

    template<class T> class StopingCriterion {
    public:
        virtual bool stoping(size_t maxIter, size_t iter, T fOldPoint, T fNewPoint, T delta, size_t stopingIter = 1) = 0;
    };

    template<class T> class AcceptanceFunction {
    public:
        virtual bool acceptance(T state1, T state2, size_t t) const = 0;
    };

    template<class T> class SimulatedAnnealing : public BlackBoxSolver<T> {
    public:
        SimulatedAnnealing(NextCandidateDistribution<T>& next, AcceptanceFunction<T>& accept,
            CoolingSchedule<T>& t, StopingCriterion<T>& stop) : mD(next), mA(accept), mTemp(t), mStop(stop) {}
        struct Option {

            //number of iterations of simulated_Annealing
            size_t maxIter = 3000;

            //f(oldPoint) - f(newPoint) < accurancy? //stoping rule
            T accuracy = 5;

            //step from NextCandidateDistribution proportional to the function domain
            T delta = 0.0065;

            size_t stopingIter = 4;

        };

        //Standart SA
        T search(int n, T* currentPoint, const T* lowerBound, const T* upperBound, const std::function<T(const T*)>& f) override {
            T* newPoint = new T[n];
            std::copy(currentPoint, currentPoint + n, newPoint);
            T fOldPoint = f(currentPoint);
            size_t i = 0;
            do {
                mD.nextCandidate(n, newPoint, lowerBound, upperBound, mOpt.delta);
                if (mA.acceptance(f(currentPoint), f(newPoint), mTemp.coolingSchedule(i))) {
                    fOldPoint = f(currentPoint);
                    std::copy(newPoint, newPoint + n, currentPoint);
                }
                i++;
            } while (mStop.stoping(mOpt.maxIter, i - 1, fOldPoint, f(newPoint), mOpt.accuracy, mOpt.stopingIter));
            delete[] newPoint;
            return f(currentPoint);
        }

        std::string about() {
            std::ostringstream options;
            options << "Simulated Annealing\n";
            options << "options:\n";
            //options << "space dimencion = " << opt.n << "\n";
            options << "number of steps = " << mOpt.maxIter << "\n";
            options << "accuracy-neighbourhood of the descent (after T=0) = " << mOpt.accuracy << "\n";
            options << "delta-neighbourhood of the temperature jump = " << mOpt.delta << "\n";
            return options.str();
        }

        Option& getOptions() {
            return mOpt;
        }

    private:
        NextCandidateDistribution<T>& mD;
        AcceptanceFunction<T>& mA;
        CoolingSchedule<T>& mTemp;
        StopingCriterion<T>& mStop;
        Option mOpt;
    };
}

#endif