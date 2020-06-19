/*
 * File: parallelcomponents.hpp
 * Author: Maksim Galynchik
 */

#ifndef PARALLELCOMPONENTS_HPP
#define PARALLELCOMPONENTS_HPP

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <random>
#include "parallelutilities.hpp"

using namespace std;

template<class T> class UniformStartPoints : public panther::GenStartPoints<T> {
public:
    void genStartPoints(int n, int numPoints, T* points, const T* lowerBound, const T* upperBound) const override {
        std::mt19937 gen;
        T* point = new T[n];
        for (unsigned int i = 0; i < numPoints; i++) {
            for (unsigned int j = 0; j < n; j++) {
                std::uniform_real_distribution<> dis(lowerBound[j], upperBound[j]);

                point[j] = dis(gen);
            }
            std::copy(point, point + n, points + i * n);
        }
        delete[] point;
    }
};

#endif
