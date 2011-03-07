
#include "C.hpp"
#include "B.hpp"
#include "A.hpp"

template <typename T>
std::vector<T>* C<T>::vec1 = 0;

template <typename T>
std::vector<B<T> >* C<T>::vec2 = 0;

template <typename T>
void C<T>::set_vec1(std::vector<T>* vec) { vec1 = vec; }

template <typename T>
void C<T>::set_vec2(std::vector<B<T> >* vec) { vec2 = vec; }

template <typename T>
C<T>::C(int x) : x(x) { }

template <typename T>
C<T>::~C() { }

template <typename T>
C2<T>::C2() : C<T>(0), y(0) {

}

template <typename T>
void C2<T>::f() {

	C<T>::value();
}

//template class C<A>;
template class C2<A>;

