/*
 * File: sanneutilities.cpp
 * Author: Maksim Galynchik
 */

#ifndef SANNEUTILITIES_HPP
#define SANNEUTILITIES_HPP

using namespace std;

namespace panther {
    /*
     * Next Candidate Distribution
     */
    template<class T> class NextCandidateDistribution {
    public:
	/*
         * Performs search for a candidate point
         * @param n space dimension
         * @param point start point and result
         * @param lowerBound lower bound of target function domain
         * @param upperBound upper bound of target function domain
         * @param delta step proportional to the function domain
         */
        virtual void nextCandidate(int n, T* point, const T* lowerBound, const T* upperBound) = 0;
        virtual std::string about() = 0;
    };

    /*
     * Cooling Schedule
     */
    template<class T> class CoolingSchedule {
    public:
        /*
         * Calculates current temperature
         * @param iteration current iteration number
         * @return current temperature
         */
        virtual T coolingSchedule(unsigned int iteration) const = 0;
        virtual std::string about() = 0;
    };
    /*
     * Stoping Criterion
     */
    template<class T> class StopingCriterion {
    public:
	/*
         * Attempt to stop Simulated Annealing
         * @param iter current iteration of Simulated Annealing
         * @param fOldPoint function value in previous point
         * @param fNewPoint function value in new point
         * @return true if the work of Simulated Annealing continues
         */
        virtual bool stoping(unsigned int iter, T fOldPoint, T fNewPoint) = 0;
        virtual std::string about() = 0;
        virtual std::pair<unsigned int, unsigned int> aboutStoping() = 0;
    };

    /*
     * Acceptance Function
     */
    template<class T> class AcceptanceFunction {
    public:
        /*
         * Attempt to accept candidate point
         * @param state1 function value at current point
         * @param state2 function value at candidate point
         * @param t temperature value
         * @return true if candidate point needs to be accepted
         */
        virtual bool acceptance(T state1, T state2, T t) const = 0;
        virtual std::string about() = 0;
    };
}

#endif
