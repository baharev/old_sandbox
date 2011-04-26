//==============================================================================
//
// This code is part of ASOL (nonlinear system solver using affine arithmetic)
//
// Copyright (C) 2011 Ali Baharev
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

#ifndef EVALUATE_HPP_
#define EVALUATE_HPP_

#include "affine.hpp"
#include "builder.hpp"
#include "interval.hpp"

namespace asol {

template <typename T>
inline void add(T& z, const T& x, const T& y) {

	z.assign(x+y);
}

void affine_add_friend();

template <>
inline void add<affine>(affine& z, const affine& x, const affine& y) {

	affine_add_friend();
}

/*
inline void add(interval& z, const interval& x, const interval& y) {

	z.assign(x+y);
}

inline void add(builder& z, const builder& x, const builder& y) {

	z.assign(x+y);
}

inline void add(affine& z, const affine& x, const affine& y) {

}
*/

}

#endif // EVALUATE_HPP_
