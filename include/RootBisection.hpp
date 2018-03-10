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

#include <iostream>

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
    T rootBisection(const std::function<T(T)>& funct, //puntero a función
                    T xl, //límite inferior de intervalo
                    T xu, //límite superior de intervalo
                    const T eps) {

        T xr= xl; //hay que iniciar con algo válido
        T fl =funct(xl); //sombra para ahorrar evaluaciones de f()
        T fu =funct(xu);
        T ea=T ( ); //error aproximado

        if(xl>xu || fl * fu > 0){
            throw anpi::Exception("received invalid values");
        }
        for(int i =std::numeric_limits<T>::digits; i > 0;--i){
            T   xrold(xr);  //para cálculo de error
            xr =( xl +xu)/T(2);  //nueva estimación de raíz, centrada
            T   fr =funct(xr) ;  //sombra de f en el centro

            //para evitar división por cero
            if ( std::abs(xr) > eps ){

                ea=std::abs ( ( xr - xrold )/ xr ) * T (100); //nuevo error aprox.
            }
            T   cond= fl * fr;  //esto es negativo si extremo inferior y el
                                // nuevo centro tienen signos distintos
            if( cond < T (0) ) {
                xu=xr;           //Es negativo siga con el lado izquierdo
            } else if ( cond > T(0)) {
                xl =xr;            //Es positivo siga con el lado derecho
                fl = fr;
            } else {
                ea = T(0);          //No hay error implica que algún borde es cero
                xr = (std::abs(fl) < eps) ? xl : xr; //fl==0
            }if ( ea < eps ){       //si se alcanzó precisión, termine
                return xr;
            }

        }

        return xr;
    }

}
#endif