


#include <boost/test/unit_test.hpp>
#include <functional>
#include <iostream>
#include <exception>
#include <cstdlib>
#include <complex>
#include <python2.7/Python.h>
#include <PlotPy.hpp>

#include <cmath>
#include <RootSecant.hpp>
#include <RootInterpolation.hpp>
#include <RootBisection.hpp>
#include <RootBrent.hpp>
#include <RootNewtonRaphson.hpp>

/**
 *
 * Clase que encapsula a cualquier funcion std::function <T(T)>.
 * Puede ser utilizada con los solucionadores implementados, es decir, puede
 * ser utilizada a su vez como std::function<T(T)>.
 * Intercepta llamadas a la funcion original e incrementa un contador,
 * que dicha clase permite acceder para saber cuantas veces se ha llamado la funcion
 * encapsulada.
 *
 * @author gabriel
 * @date   06.03.2018
 */
template<typename T>
class Encapsuladora {

private:
    // contador de llamadas
    mutable int _contador;
    // función
    std::function<T(const T)> _funct;

public:
    // Constructor
    Encapsuladora();
    Encapsuladora(const std::function<T(const T)>);
    // Sobrecarga del operador () para
    T operator ()(const T pEntrada){
        ++_contador;
        return _funct(pEntrada);
    }
    int getCuenta();
};

template<typename T>
Encapsuladora<T>::Encapsuladora() {
    _contador = 0;  //inicia el contador en cero
}

template<typename T>
Encapsuladora<T>::Encapsuladora(const std::function<T(const T)> pfunct) {
    _contador = 0;  //inicia el contador en cero
    _funct = pfunct;    //iguala la función de entrada con la del atributo para encapsularla
}

template<typename T>
int Encapsuladora<T>::getCuenta() {
    return _contador;   //get del contador para poder ver la cuenta
}

namespace anpi {
    namespace bench {

        /// Square of a number
        template<typename T>
        inline T sqr(const T x) { return x*x; }

        /// First testing function for roots |x|=e^(-x)
        template<typename T>
        T t1(const T x)  { return std::abs(x)-std::exp(-x); }

        /// Second testing function for roots e^(-x²) = e^(-(x-3)²/3 )
        template<typename T>
        T t2(const T x) { return std::exp(-x*x) - std::exp(-sqr(x-T(3))/T(3)); }

        /// Third testing function for roots x² = atan(x)
        template<typename T>
        T t3(const T x)  { return x*x-std::atan(x); }


        /**
         *
         * Método para graficar los vectores donde se guarda el epsilon y
         * las raíces según cada una de las función
         *
         * @param pMetodo es el título que llevará a gráfica
         * @param pEjeError vector con los epsilons utilizados
         * @param pEje1 vector con las raíces calculadas de la función 1
         * @param pEje2 vector con las raíces calculadas de la función 2
         * @param pEje3 vector con las raíces calculadas de la función 3
         *
         * @author gabriel
         * @date   06.03.2018
         */
        template<typename T>
        void grafica(std::string pMetodo, std::vector<T> pEjeError, std::vector<T> pEje1, std::vector<T> pEje2, std::vector<T> pEje3){
            anpi::Plot2d<T> plot2d;
            plot2d.initialize(1);
            plot2d.setTitle(pMetodo);
            plot2d.setYLabel("Llamadas a funcion");
            plot2d.setXLabel("Error");
            plot2d.setGridSize(10);
            plot2d.setYRange(0,400);
            plot2d.setXRange(0,8);


            plot2d.plot(pEjeError,pEje1,"1","g");
            plot2d.plot(pEjeError,pEje2,"2","r");
            plot2d.plot(pEjeError,pEje3,"3","b");

            plot2d.show();

        }


        /// measure the given closed root finder
        template<typename T>
        void benchTest(const std::function<T(const std::function<T(T)>&,
                                            T,
                                            T,
                                            const T)>& solver, std::string pMetodo) {

            std::vector<T> _error;
            Encapsuladora<T> EncapsuladoraF1(t1<T>);
            std::vector<T> _F1Llamadas;
            std::function<T(const T)> capsulaF1 = EncapsuladoraF1;
            Encapsuladora<T> EncapsuladoraF2(t2<T>);
            std::vector<T> _F2Llamadas;
            std::function<T(const T)> capsulaF2 = EncapsuladoraF2;
            Encapsuladora<T> EncapsuladoraF3(t3<T>);
            std::vector<T> _F3Llamadas;
            std::function<T(const T)> capsulaF3 = EncapsuladoraF3;

            for (T eps=T(1)/T(10); eps>static_cast<T>(1.0e-7); eps/=T(2)) {
                T sol = solver(capsulaF1, T(0), T(2), eps);
                T llamadas = T(capsulaF1.template target<Encapsuladora<T>>()->getCuenta());
                BOOST_CHECK(0 < llamadas);
                _F1Llamadas.push_back(llamadas);

                sol = solver(capsulaF2, T(0), T(2), eps);
                llamadas = T(capsulaF2.template target<Encapsuladora<T>>()->getCuenta());
                BOOST_CHECK(0 < llamadas);
                _F2Llamadas.push_back(llamadas);

                sol = solver(capsulaF3,T(0),T(0.5),eps);
                llamadas = T(capsulaF3.template target<Encapsuladora<T>>()->getCuenta());
                BOOST_CHECK(0 < llamadas);
                _F3Llamadas.push_back(llamadas);

                _error.push_back(eps*100);
            }
            bench::grafica(pMetodo, _error, _F1Llamadas,_F2Llamadas,_F3Llamadas);
        }


        /// measure the given open root finder
        template<typename T>
        void benchTest(const std::function<T(const std::function<T(T)>&,
                                             T,
                                            const T)>& solver, std::string pMetodo) {

            std::vector<T> _error;

            Encapsuladora<T> EncapsuladoraF1(t1<T>);
            std::vector<T> _F1Llamadas;
            std::function<T(const T)> capsulaF1 = EncapsuladoraF1;
            Encapsuladora<T> EncapsuladoraF2(t2<T>);
            std::vector<T> _F2Llamadas;
            std::function<T(const T)> capsulaF2 = EncapsuladoraF2;
            Encapsuladora<T> EncapsuladoraF3(t3<T>);
            std::vector<T> _F3Llamadas;
            std::function<T(const T)> capsulaF3 = EncapsuladoraF3;

            for (T eps=T(1)/T(10); eps>static_cast<T>(1.0e-7); eps/=T(2)) {
                T sol = solver(capsulaF1,T(0),eps);
                T llamadas = T(capsulaF1.template target<Encapsuladora<T>>()->getCuenta());
                BOOST_CHECK(0 < llamadas);
                _F1Llamadas.push_back(llamadas);

                sol = solver(capsulaF2,T(2),eps);
                llamadas = T(capsulaF2.template target<Encapsuladora<T>>()->getCuenta());
                BOOST_CHECK(0 < llamadas);
                _F2Llamadas.push_back(llamadas);

                sol = solver(capsulaF3,T(0),eps);
                llamadas = T(capsulaF3.template target<Encapsuladora<T>>()->getCuenta());
                BOOST_CHECK(0 < llamadas);
                _F3Llamadas.push_back(llamadas);

                _error.push_back(eps*100);
            }
            bench::grafica(pMetodo, _error, _F1Llamadas,_F2Llamadas,_F3Llamadas);
        }
    } // bench
}  // anpi


BOOST_AUTO_TEST_SUITE( Bench )

    BOOST_AUTO_TEST_CASE(Bisection)
    {
        anpi::bench::benchTest<float>(anpi::rootBisection<float>, "Presicion simple Biseccion");
        anpi::bench::benchTest<double>(anpi::rootBisection<double>, "Presicion doble Biseccion");
    }

    BOOST_AUTO_TEST_CASE(Interpolation)
    {
        anpi::bench::benchTest<float>(anpi::rootInterpolation<float>, "Presicion simple Interpolacion");
        anpi::bench::benchTest<double>(anpi::rootInterpolation<double>, "Presicion doble Interpolacion");
    }

    BOOST_AUTO_TEST_CASE(Secant)
    {
        anpi::bench::benchTest<float>(anpi::rootSecant<float>, "Presicion simple Secante");
        anpi::bench::benchTest<double>(anpi::rootSecant<double>, "Presicion doble Secante");
    }

    BOOST_AUTO_TEST_CASE(NewtonRaphson)
    {
        anpi::bench::benchTest<float>(anpi::rootNewtonRaphson<float>, "Presicion simple Newton-Raphson");
        anpi::bench::benchTest<double>(anpi::rootNewtonRaphson<double>, "Presicion doble Newton-Raphson");
    }

    BOOST_AUTO_TEST_CASE(Brent)
    {
        anpi::bench::benchTest<float>(anpi::rootBrent<float>, "Presicion simple Brent");
        anpi::bench::benchTest<double>(anpi::rootBrent<double>, "Presicion doble Brent");
    }

BOOST_AUTO_TEST_SUITE_END()