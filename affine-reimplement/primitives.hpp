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

#include <set>
#include <stdexcept>
#include "operations.hpp"
#include "recorder.hpp"

namespace asol {

class primitive {

public:

	virtual void evaluate(operations* op) const = 0;

	virtual void revise(operations* op) const = 0;

	// TODO Is there a way to do it without the cast?
	virtual bool common_subexpressions(const primitive* other) const = 0;

	virtual void record(recorder* rec) const = 0;

	virtual void record_indices(std::set<int>& index_set) const = 0;

	virtual ~primitive() { }

	const int z;

	const int x;

protected:

	explicit primitive(int value_offset, int first_arg)
	: z(value_offset), x(first_arg) { }

private:

	// Unused in all subclasses, just to eliminate compiler warning C4512
	primitive& operator=(const primitive& );
};

class unary_primitive : public primitive {

protected:

	unary_primitive(int value, int arg)
	: primitive(value, arg) { }

	virtual void record_indices(std::set<int>& index_set) const {
		index_set.insert(z);
		index_set.insert(x);
	}

	bool args_equal(const unary_primitive* other) const {
		return x==other->x;
	}

private:

	unary_primitive& operator=(const unary_primitive& );
};

class binary_primitive : public primitive {

public:

	const int y;

protected:

	binary_primitive(int value, int arg1, int arg2)
	: primitive(value, arg1), y(arg2) { }

	virtual void record_indices(std::set<int>& index_set) const {
		index_set.insert(z);
		index_set.insert(x);
		index_set.insert(y);
	}

	bool args_equal(const binary_primitive* other) const {
		return x==other->x && y==other->y;
	}

private:

	binary_primitive& operator=(const binary_primitive& );
};

class addition : public binary_primitive {

public:

	addition(int value, int arg1, int arg2)
	: binary_primitive(value, arg1, arg2) { }

private:

	virtual void evaluate(operations* op) const {
		op->addition(z, x, y);
	}

	virtual void revise(operations* op) const {
		op->addition_revise(z, x, y);
	}

	virtual bool common_subexpressions(const primitive* p) const {

		bool ret_val = false;

		if (const addition* other=dynamic_cast<const addition*>(p)) {
			ret_val = args_equal(other);
		}
		return ret_val;
	}

	virtual void record(recorder* rec) const {
		rec->record(this);
	}

	addition& operator=(const addition& );
};

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

class square : public unary_primitive {

public:

	square(int value, int arg)
	: unary_primitive(value, arg) { }

private:

	virtual void evaluate(operations* op) const {
		op->square(z, x);
	}

	virtual void revise(operations* op) const {
		op->square_revise(z, x);
	}

	virtual bool common_subexpressions(const primitive* p) const {

		bool ret_val = false;

		if (const square* other=dynamic_cast<const square*>(p)) {
			ret_val = args_equal(other);
		}
		return ret_val;
	}

	virtual void record(recorder* rec) const {
		rec->record(this);
	}

	square& operator=(const square& );
};

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

}

#endif // PRIMITIVES_HPP_
