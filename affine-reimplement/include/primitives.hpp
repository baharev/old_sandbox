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

#ifndef PRIMITIVES_HPP_
#define PRIMITIVES_HPP_

#include <vector>

namespace asol {

class recorder;
template <typename T> struct gap_info;

template <typename T>
class primitive {

public:

	virtual void evaluate() const = 0;

	virtual void revise() const = 0;

	// TODO Is there a way to do it without the downcast?
	virtual bool common_subexpressions(const primitive<T>* other) const = 0;

	virtual void record(recorder* rec) const = 0;

	virtual ~primitive();

	static void set_vector(std::vector<T>* vec) { v = vec; }

	static void set_gap_container(std::vector<gap_info<T> >* vec) { gaps = vec; }

protected:

	explicit primitive(int lhs);

	T& val() const { return v->at(z); }

	const int z;

	static std::vector<T>* v;

	static std::vector<gap_info<T> >* gaps;
};

template <typename T>
class unary_primitive : public primitive<T> {

protected:

	unary_primitive(int value, int arg);

	virtual bool common_subexpressions(const primitive<T>* other) const;

	virtual const unary_primitive<T>* downcast(const primitive<T>* other) const = 0;

	T& arg() const { return primitive<T>::v->at(x); }

	const int x;
};

template <typename T>
class binary_primitive : public primitive<T> {

protected:

	binary_primitive(int value, int arg1, int arg2);

	virtual bool common_subexpressions(const primitive<T>* other) const;

	virtual const binary_primitive<T>* downcast(const primitive<T>* other) const = 0;

	T& arg1() const { return primitive<T>::v->at(x); }

	T& arg2() const { return primitive<T>::v->at(y); }

	const int x;

	const int y;
};

template <typename T>
class addition : public binary_primitive<T> {

public:

	addition(int value, int arg1, int arg2);

private:

	virtual void evaluate() const;

	virtual void revise() const;

	virtual const binary_primitive<T>* downcast(const primitive<T>* p) const;

	virtual void record(recorder* rec) const;
};

template <typename T>
class substraction : public binary_primitive<T> {

public:

	substraction(int value, int arg1, int arg2);

private:

	virtual void evaluate() const;

	virtual void revise() const;

	virtual const binary_primitive<T>* downcast(const primitive<T>* p) const;

	virtual void record(recorder* rec) const;
};

template <typename T>
class multiplication : public binary_primitive<T> {

public:

	multiplication(int value, int arg1, int arg2);

private:

	virtual void evaluate() const;

	virtual void revise() const;

	virtual const binary_primitive<T>* downcast(const primitive<T>* p) const;

	virtual void record(recorder* rec) const;
};

template <typename T>
class division : public binary_primitive<T> {

public:

	division(int value, int arg1, int arg2);

private:

	virtual void evaluate() const;

	virtual void revise() const;

	virtual const binary_primitive<T>* downcast(const primitive<T>* p) const;

	virtual void record(recorder* rec) const;
};

template <typename T>
class square : public unary_primitive<T> {

public:

	square(int value, int arg);

private:

	virtual void evaluate() const;

	virtual void revise() const;

	virtual const unary_primitive<T>* downcast(const primitive<T>* other) const;

	virtual void record(recorder* rec) const;
};

template <typename T>
class exponential : public unary_primitive<T> {

public:

	exponential(int value, int arg);

private:

	virtual void evaluate() const;

	virtual void revise() const;

	virtual const unary_primitive<T>* downcast(const primitive<T>* other) const;

	virtual void record(recorder* rec) const;
};

template <typename T>
class logarithm : public unary_primitive<T> {

public:

	logarithm(int value, int arg);

private:

	virtual void evaluate() const;

	virtual void revise() const;

	virtual const unary_primitive<T>* downcast(const primitive<T>* other) const;

	virtual void record(recorder* rec) const;
};

template <typename T>
class equality_constraint : public primitive<T> {

public:

	equality_constraint(int body, int index, double rhs);

private:

	virtual void evaluate() const;

	virtual void revise() const;

	virtual bool common_subexpressions(const primitive<T>* p) const;

	virtual void record(recorder* rec) const;

	const int x;

	const double rhs;
};

template <typename T>
class common_subexpression : public primitive<T> {

public:

	common_subexpression(int index, int ordinal);

private:

	virtual void evaluate() const;

	virtual void revise() const;

	virtual bool common_subexpressions(const primitive<T>* p) const;

	virtual void record(recorder* rec) const;

	const int x;
};

}

#endif // PRIMITIVES_HPP_
