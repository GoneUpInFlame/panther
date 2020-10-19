/*
 * File: parallelutilities.cpp
 * Author: Maksim Galynchik, Romanova Karina
 */

#ifndef PARALLELUTILITIES_HPP
#define PARALLELUTILITIES_HPP

using namespace std;

namespace panther {
    /**
     * Generate Start Points class
     */
    template<class T> class GenStartPoints {
    public:
        /**
         * Performs generation of different start points
         * @param numPoints number of strart points
         * @param n space dimension
         * @param points result
         * @param lowerBound lower bound of domain
         * @param upperBound upper bound of domain
         */
        virtual void genStartPoints(int n, int numPoints, T* points, const T* lowerBound, const T* upperBound) const = 0;

    };
}

#endif
