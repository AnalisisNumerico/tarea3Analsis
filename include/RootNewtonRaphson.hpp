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

    template<typename T>
    T primeraDerivada(const std::function<T(T)>& funct, T x) {
        const T h = 0.001;
        return ((funct(x+h) - funct(x-h))/2*h);
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

        T xp = xi;
        T xa;
        T h;
        T epsilon = eps/10;

        for(int i = 0; i < MAX_ITERATIONS; i++) {
            h = (funct(xp))/(primeraDerivada(funct,xp));
            xa = xp - h;
            if(abs(xa - xp) <= epsilon) {
                return xa;
            }
            xp = xa;
        }

        // Return NaN if no root was found
        return std::numeric_limits<T>::quiet_NaN();
    }

}


#endif