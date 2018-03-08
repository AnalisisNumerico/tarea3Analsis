#include <functional>

//
// Created by gabriel on 07/03/18.
//

/**
 *
 * Clase que encapsula a cualquier funcion std::function <T(T)>.
 * Puede ser utilizada con los solucionadores implementados, es decir, debe
 * disenar su clase para que pueda ser utilizada a su vez como std::function<T(T)>.
 * Intercepta llamadas a la funcion original e incrementa un contador,
 * que dicha clase permita acceder para saber cuantas veces se ha llamado la funcion
 * encapsulada.
 *
 * @author gabriel
 * @date   06.03.2018
 */


template<typename T>
T benchmarkRootFinders(const std::function<T(T)>& funct) {
    mutable int i;
}
