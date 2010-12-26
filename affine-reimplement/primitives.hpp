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

class addition : public primitive {

public:

	addition(int value, int arg1, int arg2)
	: primitive(value), x(arg1), y(arg2) { }

private:

	virtual void evaluate(operations* op) const {
		op->addition(z, x, y);
	}

	const int x;
	const int y;
};

class substraction : public primitive {

public:

	substraction(int value, int arg1, int arg2)
	: primitive(value), x(arg1), y(arg2) { }

private:

	virtual void evaluate(operations* op) const {
		op->substraction(z, x, y);
	}

	const int x;
	const int y;
};

class multiplication : public primitive {

public:

	multiplication(int value, int arg1, int arg2)
	: primitive(value), x(arg1), y(arg2) { }

private:

	virtual void evaluate(operations* op) const {
		op->multiplication(z, x, y);
	}

	const int x;
	const int y;
};

class division : public primitive {

public:

	division(int value, int arg1, int arg2)
	: primitive(value), x(arg1), y(arg2) { }

private:

	virtual void evaluate(operations* op) const {
		op->division(z, x, y);
	}

	const int x;
	const int y;
};

class square : public primitive {

public:

	square(int value, int arg)
	: primitive(value), x(arg) { }

private:

	virtual void evaluate(operations* op) const {
		op->square(z, x);
	}

	const int x;
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
