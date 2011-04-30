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

template <typename T>
inline void sub(T& z, const T& x, const T& y) {

	z.assign(x-y);
}

template <typename T>
inline void mul(T& z, const T& x, const T& y) {

	z.assign(x*y);
}

template <typename T>
inline void div(T& z, const T& x, const T& y) {

	z.assign(x/y);
}

template <typename T>
inline void exp(T& z, const T& x) {

	z.assign(exp(x));
}

template <typename T>
inline void log(T& z, const T& x) {

	z.assign(log(x));
}

template <typename T>
inline void sqr(T& z, const T& x) {

	z.assign(sqr(x));
}

template <>
inline void add<affine>(affine& z, const affine& x, const affine& y) {

	aa_addition(z, x, y);
}

template <>
inline void sub<affine>(affine& z, const affine& x, const affine& y) {

	aa_substraction(z, x, y);
}

template <>
inline void mul<affine>(affine& z, const affine& x, const affine& y) {

	aa_multiplication(z, x, y);
}

template <>
inline void div<affine>(affine& z, const affine& x, const affine& y) {

	aa_division(z, x, y);
}

template <>
inline void exp<affine>(affine& z, const affine& x) {

	aa_exp(z, x);
}

template <>
inline void log<affine>(affine& z, const affine& x) {

	aa_log(z, x);
}

template <>
inline void sqr<affine>(affine& z, const affine& x) {

	aa_sqr(z, x);
}

}

#endif // EVALUATE_HPP_
