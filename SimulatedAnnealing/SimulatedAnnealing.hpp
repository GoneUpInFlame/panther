/* 
 * File:   SimulatedAnnealing.hpp
 * Author: Maksim Galynchik 
 */

#include <sstream>
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <functional>
#include <algorithm>

template<class T>
struct Option {

    //space dimencion
    size_t n = 1;

    //number of iterations of simulated_Annealing
    size_t max_k = 300;

    //DownHill step per coordinate(used in DownHill->neighbour)
    T accuracy = 0.001;

    //delta-neighbourhood from TemperatureJump proportional to the function domain
    T delta = 0.25;

    //domain bounds
    T * lower_bound;
    T * upper_bound;

    //Energy function
    std::function<T(const T* const)> E;

    //Temperature change function
    std::function<size_t(const size_t, const size_t)> Temp = 
        [](const size_t max_k, size_t iteration) {return max_k - iteration; };

    //Probability function
    std::function<bool(const T, const T, const size_t)> P =
        [](const T state1, const T state2, const size_t t) {return (state1 > state2) ? (1 >= (double)rand() / (RAND_MAX)) :
                                                            (exp((state1 - state2) * pow(t, -1)) >= (double)rand() / (RAND_MAX)); };

};

template<class T>
void TemperatureJump(struct Option<T> opt, T* state, T* new_state) {
    for (int i = 0; i < opt.n; i++) {
        T num;
        do {
            num = new_state[i] + (2 * ((double)rand() / (RAND_MAX)) - 1) * (opt.upper_bound[i] - opt.lower_bound[i]) * opt.delta;
        } while (num >= opt.upper_bound[i] || num <= opt.lower_bound[i]);
        new_state[i] = num;
    }
}

//sezrch for a less-energy neighbour in accurancy-neighbourhood
template<class T>
void neighbour(const struct Option<T>& opt, T* state, T* new_state) {
    T num, num2, num3;
    for (int i = 0; i < opt.n && opt.E(state) <= opt.E(new_state); i++) {
        std::copy(state, state + opt.n, new_state);
        num = (opt.upper_bound[i] - opt.lower_bound[i]) * opt.accuracy;
        num2 = std::min(num, new_state[i] - opt.lower_bound[i]);
        num3 = std::min(num, opt.upper_bound[i] - new_state[i]);
        new_state[i] -= num2;
        if (opt.E((new_state)) >= opt.E(state)) 
            new_state[i] += num2;
        new_state[i] += num3;
        if (opt.E(new_state) >= opt.E(state)) 
            new_state[i] -= num3;
    }
}

//descent to a local minimum
template<class T>
void DownHill(const struct Option<T>& opt, T* state) {
    T* new_state = new T[opt.n];
    std::copy(state, state + opt.n, new_state);
    neighbour(opt, state, new_state);
    while (opt.E(state) > opt.E(new_state)) {
        std::copy(new_state, new_state + opt.n, state);
        neighbour(opt, state, new_state);
    }
    delete[] new_state;
}

template<class T>
T* simulated_Annealing(const struct Option<T>& opt, T* state) {
    std::copy(opt.lower_bound, opt.lower_bound + opt.n, state);
    T* new_state = new T[opt.n];
    std::copy(state, state + opt.n, new_state);
    for (size_t i = 0; i < opt.max_k; i++) {
        TemperatureJump(opt, state, new_state);
        if (opt.P(opt.E(state), opt.E(new_state), opt.Temp(opt.max_k, i)))
            std::copy(new_state, new_state + opt.n, state);
    }
    DownHill(opt, state);
    delete[] new_state;
    return state;
}

template<class T>
std::string about(const struct Option<T>& opt) {
    std::ostringstream options;
            options << "Simulated Annealing\n";
            options << "options:\n";
            options << "space dimencion = " << opt.n << "\n";
            options << "number of steps = " << opt.max_k << "\n";
            options << "accuracy-neighbourhood of the descent (after T=0) = " << opt.accuracy << "\n";
            options << "delta-neighbourhood of the temperature jump = " << opt.delta << "\n";
            options << "lower_bound = [ ";
            for(size_t i = 0; i < opt.n; i++) 
                options << " " << opt.lower_bound[i];
            options << " ]\n";
            options << "upper_bound = [ ";
            for(size_t i = 0; i < opt.n; i++) 
                options << " " << opt.upper_bound[i];
            options << " ]\n";
            return options.str();
}
