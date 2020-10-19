/*
 * File:   parallel.hpp
 * Author: Maksim Galynchik, Romanova Karina
 */

#ifndef PARALLEL_HPP
#define PARALLEL_HPP

#include <sstream>
#include <common/bbsolver.hpp>
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <functional>
#include <algorithm>
#include <random>
#include <omp.h>
#include "parallelutilities.hpp"


namespace panther {
    /**
     * Universal parallel framework
     */
    template<class T> class ParallelSolver {
    public:
        /**
         * Construct Parallel Solver class
         * @param bbsolver searching method class
         * @param genp class that generate start points
         * @param pointsNum number of start points
         * @param pl determines whether class parallel or no
         */
        ParallelSolver(BlackBoxSolver<T>& bbsolver, GenStartPoints<T>& genp, unsigned int pointsNum = 1, bool pl = false) : method(bbsolver), genPoints(genp) {
            mOpt.pointsNum = pointsNum;
            mOpt.pl = pl;

        }
        struct Options {
            /*
             * Number of start points
             */
            unsigned int pointsNum;
            /*
             * Is parallel
             */
            bool pl;

        };


        /**
	       * Performs search from different start points
	       * @param n space dimension
	       * @param currentPoint result
	       * @param lowerBound lower bound of domain
         * @param upperBound upper bound of domain
	       * @param f target function
	       * @return value of the result point
	       */
        T search(int n, T* currentPoint, const T* lowerBound, const T* upperBound, const std::function<T(const T*)>& f) {
            T* points = new T [mOpt.pointsNum * n];
            genPoints.genStartPoints(n, mOpt.pointsNum, points, lowerBound, upperBound);
            T* fresult = new T[mOpt.pointsNum];

#pragma omp parallel for if(mOpt.pl) num_threads(omp_get_num_procs())
            for (int tid = 0; tid < mOpt.pointsNum; tid++) {
                fresult[tid] = method.search(n, points + tid * n, lowerBound, upperBound, f);
            }

            size_t minId = 0;
            for (size_t i = 1; i < mOpt.pointsNum; i++) {
                if (fresult[minId] > fresult[i]) {
                    minId = i;
                }
            }
            std::copy(points + minId * n, points + (minId + 1) * n, currentPoint);
            delete[] points;
            delete[] fresult;
            return f(currentPoint);

        }


        std::string about() {
            std::ostringstream options;
            options << "Parallel\n";
            options << "options:\n";
            options << "Number of start points " << mOpt.pointsNum << "\n";
            if (mOpt.pl) {
                options << "Is parallel\n";
            } else {
                options << "Is consequent\n";
            }
            return options.str();
        }

        Options& getOptions() {
            return mOpt;
        }

    private:
        BlackBoxSolver<T>& method;
        GenStartPoints<T>& genPoints;
        Options mOpt;
    };
}

#endif
