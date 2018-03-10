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

#ifndef ANPI_ROOT_BRENT_HPP
#define ANPI_ROOT_BRENT_HPP

namespace anpi {

    /**
     * Find the roots of the function funct looking for it in the
     * interval [xl,xu], using the Brent's method.
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
    T rootBrent(const std::function<T(T)>& funct,T xl,T xu,const T eps) {

      T raizA = xl;
      T raizB = xu;

      T fA = funct(raizA);
      T fB = funct(raizB);
      if(fA * fB >= 0) {
        return std::numeric_limits<T>::quiet_NaN();
      }
      if(abs(fA) < abs(fB)) {
        T temp = raizA;
        raizA = raizB;
        raizB = temp;
      }
      T raizC = raizA;
      T raizD;
      T raizS;
      bool mflag = true;
      while((funct(raizB) != 0) && (funct(raizS) !=0) && (abs(raizB - raizA) >= eps)) {
        if((funct(raizA) != funct(raizC)) && (funct(raizB) != funct(raizC))) {
          raizS = ((raizA * funct(raizB) * funct(raizC))/((funct(raizA) - funct(raizB)) * (funct(raizA) - funct(raizC)))) +
                  ((raizB * funct(raizA) * funct(raizC))/((funct(raizB) - funct(raizA)) * (funct(raizB) - funct(raizC)))) +
                  ((raizC * funct(raizA) * funct(raizB))/((funct(raizC) - funct(raizA)) * (funct(raizC) - funct(raizB))));
        }
        else {
          raizS = raizB - (funct(raizB) * ((raizB - raizA)/(funct(raizB) - funct(raizA))));
        }
        if(!((((3*raizA + raizB)/(4)) <= raizS) && (raizS <= raizB)) ||
           (mflag && (abs(raizS - raizB) >= abs(raizB - raizC)/2)) ||
           (!mflag && (abs(raizS - raizB) >= abs(raizC - raizD)/2)) ||
           (mflag && (abs(raizB - raizC) < abs(eps))) ||
           (!mflag && (abs(raizC - raizD) < abs(eps)))) {
          raizS = (raizA + raizB) / 2;
          mflag = true;
        }
        else {
          mflag = false;
        }
        raizD = raizC;
        raizC  = raizB;
        if((funct(raizA) * funct(raizS)) < 0) {
          raizB = raizS;
        }
        else {
          raizA = raizS;
        }
        if(abs(fA) < abs(fB)) {
          T temp = raizA;
          raizA = raizB;
          raizB = temp;
        }
      }

      // Return NaN if no root was found
      return raizB;
    }
}

#endif
