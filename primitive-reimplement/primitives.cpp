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

#include "primitives.hpp"
#include "interval.hpp"
#include "recorder.hpp"

namespace asol {

template <typename T>
std::vector<T>* primitive<T>::v(0);

template <typename T>
primitive<T>::primitive(int val, int arg1) : z(val), x(arg1) { }

template <typename T>
primitive<T>::~primitive() { }

template <typename T>
unary_primitive<T>::unary_primitive(int z, int x) : primitive<T>(z, x) { }

template <typename T>
bool unary_primitive<T>::common_subexpressions(const primitive<T>* p) const {

	bool ret_val = false;

	if (const unary_primitive<T>* other = downcast(p)) {

		ret_val = (this->x==other->x);
	}

	return ret_val;
}

template <typename T>
binary_primitive<T>::binary_primitive(int value, int arg1, int arg2)
: primitive<T>(value, arg1), y(arg2)
{

}

template <typename T>
bool binary_primitive<T>::common_subexpressions(const primitive<T>* p) const {

	bool ret_val = false;

	if (const binary_primitive<T>* other = downcast(p)) {

		ret_val = (this->x==other->x && y==other->y);
	}

	return ret_val;
}

template <typename T>
addition<T>::addition(int z, int x, int y) : binary_primitive<T>(z, x, y) { }

template <typename T>
void addition<T>::evaluate() const {

	this->val().assign( this->arg1() + this->arg2() );
}

template <typename T>
void addition<T>::revise() const {

	addition_inverse(this->val(), this->arg1(), this->arg2());
}

template <typename T>
const binary_primitive<T>* addition<T>::downcast(const primitive<T>* p) const {

	return dynamic_cast<const addition<T>*> (p);
}

template <typename T>
void addition<T>::record(recorder* rec) const {

	rec->addition(this->z, this->x, this->y);
}


template <typename T>
square<T>::square(int value, int arg) : unary_primitive<T>(value, arg) { }

template <typename T>
void square<T>::evaluate() const {

	this->val().assign( sqr(this->arg()) );
}

template <typename T>
void square<T>::revise() const {

	sqr_inverse(this->val(), this->arg());
}

template <typename T>
const unary_primitive<T>* square<T>::downcast(const primitive<T>* p) const {

	return dynamic_cast<const square<T>*> (p);
}

template <typename T>
void square<T>::record(recorder* rec) const {

	rec->square(this->z, this->x);
}

template class addition<interval>;
template class square<interval>;

}
