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

#include "Hansen.hpp"
#include "builder.hpp"

namespace asol {

template <typename T>
int Hansen_example<T>::number_of_variables() const {

	return SIZE;
}

template <typename T>
T* Hansen_example<T>::initial_box() const {

	T* box = new T[SIZE];

	//box[X] = T(1.33073);
	//box[Y] = T(1);

	//box[X] = T(1);
	//box[Y] = T(10);

	box[X] = T(1, 10);
	box[Y] = T(1, 10);

	return box;
}

template <typename T>
void Hansen_example<T>::evaluate(const T v[]) const {

	const T& x = v[X];
	const T& y = v[Y];

	const T xy = x*y;

	xy.mark_as_common_subexpression();

	const T z = (5*x-4*sqr(y)+14*xy)/(sqr(x)+y+xy);

	z.equals(14.5);
}

template class Hansen_example<builder>;

}
