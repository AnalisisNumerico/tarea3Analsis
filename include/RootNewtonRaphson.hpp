#include <cmath>
#include <limits>
#include <functional>

#include "Exception.hpp"

#ifndef ANPI_NEWTON_RAPHSON_HPP
#define ANPI_NEWTON_RAPHSON_HPP

namespace anpi {

    template<typename T>
    T primeraDerivada(const std::function<T(T)>& funct, T x,T eps) {
        const T h = eps/2;
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

        T eps = eps/10;

        double err, x1;
        int it, maxit=100;
        it=0;
        err=eps+1;
        while( err > eps && it < maxit )
        {
            x1 = xi - funct(xi)/(primeraDerivada(funct,xi,eps));
            err=fabs(x1-xi);
            xi=x1;
            it++;
        }
        if( err <= eps ){
            return xi;
        }

        // Return NaN if no root was found
        return std::numeric_limits<T>::quiet_NaN();
    }

}


#endif