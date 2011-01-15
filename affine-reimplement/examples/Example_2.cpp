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

#include "Example_2.hpp"
#include "builder.hpp"

namespace asol {

template <typename T>
int Example_2<T>::number_of_variables() const {

	return SIZE;
}

template <typename T>
T* Example_2<T>::initial_box() const {

	T* box = new T[SIZE];

	box[X] = T(1.0, 2.0);

	return box;
}

template <typename T>
void Example_2<T>::evaluate(const T v[]) const {

	const T& x = v[X];

	const T  y = (sqr(x)+x)/(16*x-9);
}

template class Example_2<builder>;

}
