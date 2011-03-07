#include "D.hpp"
#include "C.hpp"
#include "B.hpp"
#include "A.hpp"

template <typename T>
D<T>::D() {

}

template <typename T>
D<T>::~D() {

}

template <typename T>
void D<T>::f() {

		C<T>::set_vec1(&v1);
		C<T>::set_vec2(&v2);
}

template class D<A>;
