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

#ifndef ANPI_ROOT_INTERPOLATION_HPP
#define ANPI_ROOT_INTERPOLATION_HPP


namespace anpi {

    /**
     * Find the roots of the function funct looking for it in the
     * interval [xl,xu], by means of the interpolation method.
     *
     * @param funct a functor of the form "T funct(T x)"
     * @param xl lower interval limit
     * @param xu upper interval limit
     *
     * @return root found, or NaN if none could be found.
     *
     * @throws anpi::Exception if inteval is reversed or both extremes
     *         have same sign.
     */
    template<typename T>
    T rootInterpolation(const std::function<T(T)>& funct,T xl,T xu,const T eps) {

        T xr=xl;    //hay que iniciar con algo válido
        T fl = funct(xl);
        T fu = funct(xu);
        T ea=T ( );
        int iu (0), il (0); //contadores para detectar estancamientos

        if(xl>xu || fl * fu > 0){
            throw anpi::Exception("received invalid values");
        }


        for(int i =std::numeric_limits<T>::digits; i > 0; --i) {
            T xrold(xr);   //para cálculo de error
            xr=xu - fu * (xl-xu)/(fl - fu);
            T fr =funct (xr);
            //para evitar división por cero
            if ( std::abs(xr) > eps ) {
                ea   =   std::abs( ( xr-xrold) / xr) * T(100);
            }
            T cond= fl * fr; //cual intervalo tiene la raíz
            if(cond < T(0)) {   //el lado izquierdo tiene la raíz
                xu=xr;
                fu= fr;
                iu =0;
                il ++;
                if ( il >= 2) {
                    fl   /=  T ( 2 ) ;
                }
            } else if ( cond > T (0) ) {//el lado derecho tiene la raíz
                xl = xr;
                fl = fr;
                il =0;
                iu ++;
                if (iu >=2) {
                    fu /= 2;
                }
            } else {
                ea = T(0);  //no hay error
                xr = (fl == T(0)) ? xl : xu;
            }if ( ea < eps ){   //si se alcanzó precisión, termine
                return xr;
            }
        }
        return std::numeric_limits<T>::quiet_NaN();
    }
}
#endif