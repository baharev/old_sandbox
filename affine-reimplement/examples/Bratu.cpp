//==============================================================================
//
// This code is part of ASOL (nonlinear system solver using affine arithmetic)
//
// Copyright (C) 2010, 2011 Ali Baharev
// All rights reserved. E-mail: <my_first_name.my_last_name@gmail.com>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program. If not, see <http://www.gnu.org/licenses/>.
//
//==============================================================================

#include "Bratu.hpp"
#include "builder.hpp"

namespace asol {

template <typename T>
int Bratu<T>::number_of_variables() const {

	return SIZE;
}

template <typename T>
T* Bratu<T>::initial_box() const {

	T* x = new T[SIZE];

	for (int i=0; i<SIZE; ++i) {

		x[i] = T(-10.0, 10.0);
	}

	return x;
}

template <typename T>
void Bratu<T>::evaluate(const T x[]) const {

	const double c1 = -2.0, c2 = 1.04058273e-03;

	T eq_0 = x[1] + c1*x[0] + c2*exp(x[0]);

	eq_0.equals(0.0);

	for (int i=1; i<SIZE-1; ++i) {

		T eq_i = x[i+1] + c1*x[i] + c2*exp(x[i]) + x[i-1];

		eq_i.equals(0.0);
	}

	const int N = SIZE-1;

	T eq_N = x[N-1] + c1*x[N] + c2*exp(x[N]);

	eq_N.equals(0.0);
}

template class Bratu<builder>;

}
