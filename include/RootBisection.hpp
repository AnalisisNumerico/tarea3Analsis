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

#ifndef ANPI_ROOT_BISECTION_HPP
#define ANPI_ROOT_BISECTION_HPP

namespace anpi {

    /**
     * Find the roots of the function funct looking for it in the
     * interval [xl,xu], using the bisection method.
     *
     * @param funct a std::function of the form "T funct(T x)"
     * @param xl lower interval limit
     * @param xu upper interval limit
     *
     * @return root found, or NaN if none could be found.
     *
     * @throws anpi::Exception if inteval is reversed or both extremes
     *         have same sign.
     */
    template<typename T>
    T rootBisection(const std::function<T(T)>& funct,T xl,T xu,const T eps) {

        T dx,f,fmid,xmid,rtb;
        f=funct(xl);
        fmid=funct(xu);
        if(xl>xu || f * fmid >= 0){
            throw anpi::Exception("received invalid values");
        }        rtb = f < 0.0 ? (dx=xu-xl,xl) : (dx=xl-xu,xu);
        for (int j=1;j<=std::numeric_limits<T>::digits;j++) {
            fmid=funct(xmid=rtb+(dx *= 0.5));
            if (fmid <= 0.0) rtb=xmid;
            if (fabs(dx) < eps || fmid == 0.0) return rtb;
        }

        // Return NaN if no root was found
        return std::numeric_limits<T>::quiet_NaN();
    }

}
#endif