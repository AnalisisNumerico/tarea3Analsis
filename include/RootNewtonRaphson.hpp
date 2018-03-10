/**
 * Copyright (C) 2018
 * Área Académica de Ingeniería en Computadoras, ITCR, Costa Rica
 *
 * This file is part of the numerical analysis lecture CE3102 at TEC
 *
 * @Author: Pablo Alvarado
 * @Date  : 10.02.2018
 */

#include <cmath>
#include <limits>
#include <functional>

#include "Exception.hpp"

#ifndef ANPI_NEWTON_RAPHSON_HPP
#define ANPI_NEWTON_RAPHSON_HPP

namespace anpi {

    /** Returns the evaluated derivate of a given function
     * by means of forward aproximation
     *
     *
     * @param funct a functor of the form "T funct(T x)"
     * @param value
     * @return evalueted derivate aproximation
     */
    template<typename T>
    T primeraDerivada(const std::function<T(T)>& funct, T x,T eps) {
        const T h = std::abs(eps) / T(2);
        return ((funct(x+h) - funct(x))/h);
    }

    /**
     * Find the roots of the function funct looking by means of the
     * Newton-Raphson method
     *
     * @param funct a functor of the form "T funct(T x)"
     * @param xi initial root guess
     *
     * @return root found, or NaN if none could be found.
     *
     * @throws anpi::Exception if inteval is reversed or both extremes
     *         have same sign.
     */
    template<typename T>
    T rootNewtonRaphson(const std::function<T(T)>& funct,T xi,const T eps) {

        int const MAX_ITERATIONS = 20;

        T x = xi;
        T dx;

        for(int i = 0; i < MAX_ITERATIONS; i++) {
            dx = ((funct(x))/(primeraDerivada(funct,x,eps)));
            x = x - dx;
            if(std::abs(dx) < eps) {
                return x;
            }
        }

        // Return NaN if no root was found
        return std::numeric_limits<T>::quiet_NaN();
    }

}

#endif