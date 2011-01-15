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

#include "Challenge.hpp"
#include "builder.hpp"

namespace asol {

template <typename T>
int Example_challange<T>::number_of_variables() const {

	return SIZE;
}

template <typename T>
T* Example_challange<T>::initial_box() const {

	T* box = new T[SIZE];

	box[W] = T(-0.9, -0.6);
	box[X] = T(-0.1,  0.2);

	box[Y] = T( 0.3,  0.7);
	box[Z] = T(-0.2,  0.1);

	box[B] = T(-1, 1);
	box[C] = T(-1, 1);

	box[A] = T(7.0, 9.0);

	return box;
}

template <typename T>
void Example_challange<T>::evaluate(const T var[]) const {

	const T& w = var[W];
	const T& x = var[X];

	const T& y = var[Y];
	const T& z = var[Z];

	const T& b = var[B];
	const T& c = var[C];

	const T& a = var[A];

	const T u = sqr(w) + sqr(x);

	u.mark_as_common_subexpression();

	const T v = sqr(y) + sqr(z);

	v.mark_as_common_subexpression();

	const T d = b*(x*y - w*z) + c*(x*z + w*y);

	const T f = (a*(u-v)+2*d)/(u+v);

	f.equals(0);
}

template class Example_challange<builder>;

}
