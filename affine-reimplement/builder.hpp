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

#ifndef BUILDER_HPP_
#define BUILDER_HPP_

#include <vector>

namespace asol {

class primitive;

class builder {

public:

	builder();

	builder(double value);

	builder(double lb, double ub);

	friend const builder operator+(const builder& x, const builder& y);

	friend const builder operator-(const builder& x, const builder& y);

	friend const builder operator*(const builder& x, const builder& y);

	friend const builder operator/(const builder& x, const builder& y);

	friend const builder sqr(const builder& x);

	static int number_of_variables();

	static const std::vector<primitive*>& get_primitives();

	static void reset();

	void dbg_consistency() const;

private:

	static int unused_index;

	static std::vector<primitive*> primitives;

	int index;
};

void dbg_consistency(const builder& x, const builder& y);

}

#endif // BUILDER_HPP_
