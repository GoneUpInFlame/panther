#include <iostream>
#include <cmath>
#include <cstdlib>
#include <functional>
#include <algorithm>



template<class T>
double P(const T state1, const T state2, size_t t) {
    if (state1 > state2) return 1;
    return exp((state1 - state2) * pow(t, -1));
}

template<class T>
void a_sum(int n, T* lower_bound, T* upper_bound, T*state) {
    for (size_t i = 0; i < n; i++) {
        state[i] = (lower_bound[i] + upper_bound[i]) / 2;
    }
}

template<class T>
void RandomState(int n, T* lower_bound, T* upper_bound, T* state, T* new_state, T delta) {
    for (int i = 0; i < n; i++) {
        T num;
        do {
            num = state[i] + (2 * ((T)rand() / (RAND_MAX)) - 1) * delta;
        } while (num >= upper_bound[i] || num <= lower_bound[i]);
        new_state[i] = num;
    }
}

template<class T>
void neighbour(int n, T* lower_bound, T* upper_bound, T* state, const std::function<T(const T* const)>& E, T delta, size_t t) {
    /*T num;
    for (int i = 0; i < n && E(state) <= E(new_state); i++) {
        std::copy(state, state + n, new_state);
        num = (upper_bound[i] - lower_bound[i]) / (T)1000;
        new_state[i] -= num;
        if (E((new_state)) >= E(state)) new_state[i] += 2 * num;
        if(E(new_state) >= E(state)) new_state[i] -= num;
    }*/
    T*new_state = new T[n];
    while(true) {
        RandomState(n, lower_bound, upper_bound, state, new_state, delta);
        if (P(E(state), E(new_state), t) >= ((double)rand() / RAND_MAX) ) {
            break;
        }
    }
    std::copy(new_state, new_state + n, state);
    delete[] new_state;
}



/*template<class T>
void DownHill(int n, T* lower_bound, T* upper_bound, const std::function<T(const T* const)>& E, T* state) { //нету ограничений на область определения
    T* new_state = new T[n];
    std::copy(state, state + n, new_state);
    neighbour(state, n, lower_bound, upper_bound, new_state, E);
    while (E(state) > E(new_state)) {
        std::copy(new_state, new_state + n, state);
        neighbour(state, n, lower_bound, upper_bound, new_state, E);
    }
}*/

template<class T>
void DownHill(int n, T* lower_bound, T* upper_bound, T* state, const std::function<T(const T* const)>& E, T delta) {
    T d = delta / 100;
    T * state_n = new T[n];
    std::copy(state, state + n, state_n);
    while(true) {
        std::copy(state_n, state_n + n, state);
        for (int i = 0; i < n; i++) {
            T d1 = std::min(d, upper_bound[i] - state_n[i]);
            state_n[i] += d1;
            if (E(state) > E(state_n)) {}
            else {
                state_n[i] -= d1;
                T d2 = std::min(d, state_n[i] - lower_bound[i]);
                state_n[i] -= d2;
                if (E(state) > E(state_n)) {}
                else {state_n[i] += d2;}
            }
        }
        if (E(state) == E(state_n)) {
            break;
        }
    }
    delete[] state_n;
}

template<class T>
void simulated_Annealing(int n, T* lower_bound, T* upper_bound, T* res, const std::function<T(const T* const)>& E, T delta, size_t max_k=100) {

    T * state = new T[n];
    a_sum(n, lower_bound, upper_bound, state);

    //std::vector<T> new_state(state.begin(), state.end());
    DownHill(n, lower_bound, upper_bound, state,  E, delta);

    for (size_t k = 0; k < max_k; k++) {
        //RandomState(state, lower_boundv, upper_boundv, new_state, delta);
        //size_t state_t = 10 * max_k - (i + 1);
        //if (P(E(state), E(new_state), state_t) >= (double)rand() / (RAND_MAX)) std::copy(new_state, new_state + n, state);
        //DownHill(n, lower_bound, upper_bound, E, state);
        neighbour(n, lower_bound, upper_bound, state, E, delta, max_k - k);
        DownHill(n, lower_bound, upper_bound, state, E, delta);
    }
    //DownHill(n, lower_bound, upper_bound, state, E, delta);
    std::copy(state, state + n, res);
    delete[] state;
    //DownHill(n, lower_bound, upper_bound, E, state);
}
