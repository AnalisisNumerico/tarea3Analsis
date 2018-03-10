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

#ifndef ANPI_ROOT_SECANT_HPP
#define ANPI_ROOT_SECANT_HPP

namespace anpi {
  
  /**
   * Find a root of the function funct looking for it starting at xi
   * by means of the secant method.
   *
   * @param funct a functor of the form "T funct(T x)"
   * @param xi initial position
   * @param xii second initial position 
   *
   * @return root found, or NaN if no root could be found
   */
  template<typename T>
  T rootSecant(const std::function<T(T)>& funct,T xi,T xii,const T eps) {
      T fl,f,dx,swap,xl,rts;
      fl=funct(xi);
      f=funct(xii);
      if (fabs(fl) < fabs(f)) { //Pick  the  bound  with  the  smaller  function  value  as
                                // the  most  recent  guess.
        rts=xi;
        xl=xii;
        swap=fl;
        fl=f;
        f=swap;
      } else {
        xl=xi;
        rts=xii;
      }
      for (int j=1;j<=std::numeric_limits<T>::digits;j++) { //Secant  loop.
                                                            // Increment  with  respect  to  latest  value.
                dx=(xl-rts)*f/(f-fl);
                xl=rts;
        fl=f;
        rts += dx;
        f=funct(rts);
        if (fabs(dx) < eps || f == 0.0){
          return rts;                                       //Convergencia
        }
      }
    
    // Return NaN if no root was found
    return std::numeric_limits<T>::quiet_NaN();
  }

}
  
#endif

