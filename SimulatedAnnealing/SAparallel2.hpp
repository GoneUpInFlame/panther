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

        //parallel SA2
        T search(int n, T* currentPoint, const T* lowerBound, const T* upperBound, const std::function<T(const T*)>& f) override {
            int np = omp_get_num_procs();
            T* domain = new T[n];
            for (int i = 0; i < n; i++) {
                domain[i] = (upperBound[i] - lowerBound[i]) / np;
            }
            omp_set_num_threads(np);
            T* points = new T[np * n];
            T* bounds = new T[(np + 1) * n];
            std::copy(lowerBound, lowerBound + n, bounds);
            for (size_t i = 1; i < np + 1; i++) {
                std::copy(bounds + (i - 1) * n, bounds + i * n, bounds + i * n);
                for (size_t j = 0; j < n; j++) {
                    bounds[i * n + j] += domain[j];
                }
                std::copy(bounds + i * n, bounds + (i + 1) * n, points + (i - 1) * n);
                for (size_t j = 0; j < n; j++) {
                    points[(i - 1) * n + j] -= domain[j] / 2;
                }
            }
#pragma omp parallel num_threads(np)
            {
                size_t tid = omp_get_thread_num();
                T fOldPoint = f(points + tid * n);
                T* newPoint = new T[n];
                std::copy(points + tid * n, points + (tid + 1) * n, newPoint);
                size_t i = 0;
                do {
                    mD.nextCandidate(n, newPoint, bounds + tid * n, bounds + (tid + 1) * n, mOpt.delta);
                    if (mA.acceptance(f(points + tid * n), f(newPoint), mTemp.coolingSchedule(i))) {
                        fOldPoint = f(points + tid * n);
                        std::copy(newPoint, newPoint + n, points + tid * n);
                    }
                    i++;
                } while (mStop.stoping(mOpt.maxIter, i - 1, fOldPoint, f(newPoint), mOpt.accuracy, mOpt.stopingIter));

                std::cout << "Thread " << tid << " find " << f(points + tid * n) << " in " << points[tid * n] << std::endl;
                delete[] newPoint;
            }

            size_t minId = 0;
            for (size_t i = 1; i < np; i++) {
                if (f(points + minId * n) > f(points + i * n)) {
                    minId = i;
                }
            }
            std::copy(points + minId * n, points + (minId + 1) * n, currentPoint);

            delete[] bounds;
            delete[] points;
            delete[] domain;

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