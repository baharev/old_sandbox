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

#include "operations.hpp"

namespace asol {

class primitive {

public:

	virtual void evaluate(operations* op) const = 0;

	virtual ~primitive() { }

protected:

	explicit primitive(int value_offset) : z(value_offset) { }

	const int z;
};

class unary_primitive : public primitive {

protected:

	unary_primitive(int value, int arg)
	: primitive(value), x(arg) { }

	const int x;
};

class binary_primitive : public primitive {

protected:

	binary_primitive(int value, int arg1, int arg2)
	: primitive(value), x(arg1), y(arg2) { }

	const int x;
	const int y;
};

class addition : public binary_primitive {

public:

	addition(int value, int arg1, int arg2)
	: binary_primitive(value, arg1, arg2) { }

private:

	virtual void evaluate(operations* op) const {
		op->addition(z, x, y);
	}
};

class substraction : public binary_primitive {

public:

	substraction(int value, int arg1, int arg2)
	: binary_primitive(value, arg1, arg2) { }

private:

	virtual void evaluate(operations* op) const {
		op->substraction(z, x, y);
	}
};

class multiplication : public binary_primitive {

public:

	multiplication(int value, int arg1, int arg2)
	: binary_primitive(value, arg1, arg2) { }

private:

	virtual void evaluate(operations* op) const {
		op->multiplication(z, x, y);
	}
};

class division : public binary_primitive {

public:

	division(int value, int arg1, int arg2)
	: binary_primitive(value, arg1, arg2) { }

private:

	virtual void evaluate(operations* op) const {
		op->division(z, x, y);
	}
};

class square : public unary_primitive {

public:

	square(int value, int arg)
	: unary_primitive(value, arg) { }

private:

	virtual void evaluate(operations* op) const {
		op->square(z, x);
	}
};

class equality_constraint : public primitive {

public:

	equality_constraint(int body, int rhs)
	: primitive(body), x(rhs) { }

private:

	virtual void evaluate(operations* op) const {

		op->equality_constraint(z, x);
	}

	const int x;
};

}

#endif // PRIMITIVES_HPP_
