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

#include <stdexcept>
#include <vector>

namespace asol {

class recorder;

template <typename T>
class primitive {

public:

	virtual void evaluate() const = 0;

	virtual void revise() const = 0;

	// TODO Is there a way to do it without the cast?
	virtual bool common_subexpressions(const primitive<T>* other) const = 0;

	virtual void record(recorder* rec) const = 0;

	virtual ~primitive();

	static void set_vector(std::vector<T>* vec) { v = vec; }

protected:

	explicit primitive(int value_offset, int first_arg);

	const int z;

	const int x;

	static std::vector<T>* v;

private:

	// Unused in all subclasses, just to eliminate compiler warning C4512
	primitive& operator=(const primitive& );
};

template <typename T>
class unary_primitive : public primitive<T> {

protected:

	unary_primitive(int value, int arg);

	virtual bool common_subexpressions(const primitive<T>* other) const;

	virtual const unary_primitive<T>* downcast(const primitive<T>* other) const = 0;

	T& val() const { return primitive<T>::v->at(this->z); }

	T& arg() const { return primitive<T>::v->at(this->x); }

private:

	unary_primitive& operator=(const unary_primitive& );
};

template <typename T>
class binary_primitive : public primitive<T> {

protected:

	binary_primitive(int value, int arg1, int arg2);

	virtual bool common_subexpressions(const primitive<T>* other) const;

	virtual const binary_primitive<T>* downcast(const primitive<T>* other) const = 0;

	T& val()  const { return primitive<T>::v->at(this->z); }

	T& arg1() const { return primitive<T>::v->at(this->x); }

	T& arg2() const { return primitive<T>::v->at(this->y); }

	const int y;

private:

	binary_primitive& operator=(const binary_primitive& );
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

	addition& operator=(const addition& );
};
/*
class substraction : public binary_primitive {

public:

	substraction(int value, int arg1, int arg2)
	: binary_primitive(value, arg1, arg2) { }

private:

	virtual void evaluate(operations* op) const {
		op->substraction(z, x, y);
	}

	virtual void revise(operations* op) const {
		op->substraction_revise(z, x, y);
	}

	virtual bool common_subexpressions(const primitive* p) const {

		bool ret_val = false;

		if (const substraction* other=dynamic_cast<const substraction*>(p)) {
			ret_val = args_equal(other);
		}
		return ret_val;
	}

	virtual void record(recorder* rec) const {
		rec->record(this);
	}

	substraction& operator=(const substraction& );
};

class multiplication : public binary_primitive {

public:

	multiplication(int value, int arg1, int arg2)
	: binary_primitive(value, arg1, arg2) { }

private:

	virtual void evaluate(operations* op) const {
		op->multiplication(z, x, y);
	}

	virtual void revise(operations* op) const {
		op->multiplication_revise(z, x, y);
	}

	virtual bool common_subexpressions(const primitive* p) const {

		bool ret_val = false;

		if (const multiplication* other=dynamic_cast<const multiplication*>(p)){
			ret_val = args_equal(other);
		}
		return ret_val;
	}

	virtual void record(recorder* rec) const {
		rec->record(this);
	}

	multiplication& operator=(const multiplication& );
};

class division : public binary_primitive {

public:

	division(int value, int arg1, int arg2)
	: binary_primitive(value, arg1, arg2) { }

private:

	virtual void evaluate(operations* op) const {
		op->division(z, x, y);
	}

	virtual void revise(operations* op) const {
		op->division_revise(z, x, y);
	}

	virtual bool common_subexpressions(const primitive* p) const {

		bool ret_val = false;

		if (const division* other=dynamic_cast<const division*>(p)) {
			ret_val = args_equal(other);
		}
		return ret_val;
	}

	virtual void record(recorder* rec) const {
		rec->record(this);
	}

	division& operator=(const division& );
};
*/


template <typename T>
class square : public unary_primitive<T> {

public:

	square(int value, int arg);

private:

	virtual void evaluate() const;

	virtual void revise() const;

	virtual const unary_primitive<T>* downcast(const primitive<T>* other) const;

	virtual void record(recorder* rec) const;

	square& operator=(const square& );
};
/*
class exponential : public unary_primitive {

public:

	exponential(int value, int arg)
	: unary_primitive(value, arg) { }

private:

	virtual void evaluate(operations* op) const {
		op->exponential(z, x);
	}

	virtual void revise(operations* op) const {
		op->exponential_revise(z, x);
	}

	virtual bool common_subexpressions(const primitive* p) const {

		bool ret_val = false;

		if (const exponential* other=dynamic_cast<const exponential*>(p)) {
			ret_val = args_equal(other);
		}
		return ret_val;
	}

	virtual void record(recorder* rec) const {
		rec->record(this);
	}

	exponential& operator=(const exponential& );
};

class equality_constraint : public primitive {

public:

	equality_constraint(int body, int rhs)
	: primitive(body, rhs) { }

private:

	virtual void evaluate(operations* op) const {
		op->equality_constraint(z, x);
	}

	virtual void revise(operations* op) const {
		op->equality_constraint_revise(z, x);
	}

	virtual void record_indices(std::set<int>& ) const {
		// TODO Figure out what and how to do
		throw std::logic_error("record_indices called on a constraint");
	}

	virtual bool common_subexpressions(const primitive* p) const {

		bool ret_val = false;

		if (const equality_constraint* other=dynamic_cast<const equality_constraint*>(p)) {
			ret_val = (x == other->x);
		}
		return ret_val;
	}

	virtual void record(recorder* rec) const {
		rec->record(this);
	}

	equality_constraint& operator=(const equality_constraint& );
};

*/

}

#endif // PRIMITIVES_HPP_
