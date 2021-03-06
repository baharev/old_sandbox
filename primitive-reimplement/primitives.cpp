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
std::vector<T> primitive<T>::v = std::vector<T>();

template <typename T>
primitive<T>::primitive(int lhs) : z(lhs) { }

template <typename T>
primitive<T>::~primitive() { }

template <typename T>
unary_primitive<T>::unary_primitive(int value, int arg) :
primitive<T>(value), x(arg)
{

}

template <typename T>
bool unary_primitive<T>::common_subexpressions(const primitive<T>* p) const {

	bool ret_val = false;

	if (const unary_primitive<T>* other = downcast(p)) {

		ret_val = (x==other->x);
	}

	return ret_val;
}

template <typename T>
binary_primitive<T>::binary_primitive(int value, int arg1, int arg2)
: primitive<T>(value), x(arg1), y(arg2)
{

}

template <typename T>
bool binary_primitive<T>::common_subexpressions(const primitive<T>* p) const {

	bool ret_val = false;

	if (const binary_primitive<T>* other = downcast(p)) {

		ret_val = (x==other->x && y==other->y);
	}

	return ret_val;
}

template <typename T>
addition<T>::addition(int z, int x, int y) :
binary_primitive<T>(z, x, y)
{

}

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
substraction<T>::substraction(int z, int x, int y) :
binary_primitive<T>(z, x, y)
{

}

template <typename T>
void substraction<T>::evaluate() const {

	this->val().assign( this->arg1() - this->arg2() );
}

template <typename T>
void substraction<T>::revise() const {

	substraction_inverse(this->val(), this->arg1(), this->arg2());
}

template <typename T>
const binary_primitive<T>* substraction<T>::downcast(const primitive<T>* p) const {

	return dynamic_cast<const substraction<T>*> (p);
}

template <typename T>
void substraction<T>::record(recorder* rec) const {

	rec->substraction(this->z, this->x, this->y);
}

template <typename T>
multiplication<T>::multiplication(int z, int x, int y) :
binary_primitive<T>(z, x, y)
{

}

template <typename T>
void multiplication<T>::evaluate() const {

	this->val().assign( this->arg1() * (this->arg2()) );
}

template <typename T>
void multiplication<T>::revise() const {

	multiplication_inverse(this->val(), this->arg1(), this->arg2());
}

template <typename T>
const binary_primitive<T>* multiplication<T>::downcast(const primitive<T>* p) const {

	return dynamic_cast<const multiplication<T>*> (p);
}

template <typename T>
void multiplication<T>::record(recorder* rec) const {

	rec->multiplication(this->z, this->x, this->y);
}

template <typename T>
division<T>::division(int z, int x, int y) :
binary_primitive<T>(z, x, y)
{

}

template <typename T>
void division<T>::evaluate() const {
	// TODO Extended division with intersection
	this->val().assign( this->arg1() / this->arg2() );
}

template <typename T>
void division<T>::revise() const {

	division_inverse(this->val(), this->arg1(), this->arg2());
}

template <typename T>
const binary_primitive<T>* division<T>::downcast(const primitive<T>* p) const {

	return dynamic_cast<const division<T>*> (p);
}

template <typename T>
void division<T>::record(recorder* rec) const {

	rec->division(this->z, this->x, this->y);
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

template <typename T>
exponential<T>::exponential(int z, int x) : unary_primitive<T>(z, x) { }

template <typename T>
void exponential<T>::evaluate() const {

	this->val().assign( exp(this->arg()) );
}

template <typename T>
void exponential<T>::revise() const {

	exp_inverse(this->val(), this->arg());
}

template <typename T>
const unary_primitive<T>* exponential<T>::downcast(const primitive<T>* p) const {

	return dynamic_cast<const exponential*> (p);
}

template <typename T>
void exponential<T>::record(recorder* rec) const {

	rec->exponential(this->z, this->x);
}

template <typename T>
equality_constraint<T>::equality_constraint(int body, int index, double rhs) :
primitive<T>(body), x(index), rhs(rhs)
{

}

template <typename T>
void equality_constraint<T>::evaluate() const {

	this->val().equals(rhs);
}

template <typename T>
void equality_constraint<T>::revise() const {

	equality_constraint_inverse(this->val(), rhs);
}

template <typename T>
bool equality_constraint<T>::common_subexpressions(const primitive<T>* p) const {

	return false;
}

template <typename T>
void equality_constraint<T>::record(recorder* rec) const {

	rec->equality_constraint(this->z, this->x);
}

template class addition<interval>;
template class substraction<interval>;
template class multiplication<interval>;
template class division<interval>;
template class square<interval>;
template class exponential<interval>;
template class equality_constraint<interval>;

}
