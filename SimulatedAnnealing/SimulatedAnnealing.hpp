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
#include <functional>
#include <algorithm>


namespace panther {
    template<class T> class NextCandidateDistribution {
    public:
        virtual void nextCandidate(int n, T* point, const T* lower_bound, const T* upper_bound, T delta) = 0;
    };

    template<class T> class CoolingSchedule {
    public:
        virtual T coolingSchedule(const size_t max_k, size_t iteration) = 0;
    };

    template<class T> class StopingCriterion {
    public:
        virtual bool stoping() = 0;
    };

    template<class T> class AcceptanceFunction {
    public:
        virtual bool acceptance(const T state1, const T state2, const size_t t) = 0;
    };

    template<class T> class SimulatedAnnealing : public BlackBoxSolver<T> {
    public:
        SimulatedAnnealing(NextCandidateDistribution<T>& next, AcceptanceFunction<T>& accept, CoolingSchedule<T>& t) : D(next), A(accept), Temp(t) {}
        struct Option {

            //number of iterations of simulated_Annealing(create a stoping rule!!!!!)
            size_t max_k = 300;

            //DownHill step per coordinate(used in DownHill->neighbour)
            T accuracy = 0.001;

            //delta-neighbourhood from TemperatureJump proportional to the function domain
            T delta = 0.25;

        };

        //search for a less-energy neighbour in accurancy-neighbourhood
        void neighbour(int n, T* current_point, T* new_point, const T* lower_bound, const T* upper_bound, const std::function<T(const T*)>& f) {
            T num, num2, num3;
            for (int i = 0; i < n && f(current_point) <= f(new_point); i++) {
                std::copy(current_point, current_point + n, new_point);
                num = (upper_bound[i] - lower_bound[i]) * opt.accuracy;
                num2 = std::min(num, new_point[i] - lower_bound[i]);
                num3 = std::min(num, upper_bound[i] - new_point[i]);
                new_point[i] -= num2;
                if (f((new_point)) >= f(current_point))
                    new_point[i] += num2;
                new_point[i] += num3;
                if (f(new_point) >= f(current_point))
                    new_point[i] -= num3;
            }
        }

        //descent to a local minimum
        void downHill(int n, T* current_point, const T* lower_bound, const T* upper_bound, const std::function<T(const T*)>& f) {
            T* new_point = new T[n];
            std::copy(current_point, current_point + n, new_point);
            neighbour(n, current_point, new_point, lower_bound, upper_bound, f);
            while (f(current_point) > f(new_point)) {
                std::copy(new_point, new_point + n, current_point);
                neighbour(n, current_point, new_point, lower_bound, upper_bound, f);
            }
            delete[] new_point;
        }



        T search(int n, T* current_point, const T* lower_bound, const T* upper_bound, const std::function<T(const T*)>& f) override {
            T* new_point = new T[n];
            std::copy(current_point, current_point + n, new_point);
            for (size_t i = 0; i < opt.max_k; i++) { //stoping criterion
                D.nextCandidate(n, new_point, lower_bound, upper_bound, opt.delta);
                if (A.acceptance(f(current_point), f(new_point), Temp.coolingSchedule(opt.max_k, i)))
                    std::copy(new_point, new_point + n, current_point);
            }
            downHill(n, current_point, lower_bound, upper_bound, f);
            delete[] new_point;
            return f(current_point);
        }

        std::string about() {
            std::ostringstream options;
            options << "Simulated Annealing\n";
            options << "options:\n";
            options << "number of steps = " << opt.max_k << "\n";
            options << "accuracy-neighbourhood of the descent (after T=0) = " << opt.accuracy << "\n";
            options << "delta-neighbourhood of the temperature jump = " << opt.delta << "\n";
            return options.str();
        }

        Option& getOptions() {
            return opt;
        }

    private:
        NextCandidateDistribution<T>& D;
        AcceptanceFunction<T>& A;
        CoolingSchedule<T>& Temp;
        Option opt;
    };
}

#endif