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
#include "typedefs.hpp"

namespace asol {

class primitive;

class builder {

public:

	builder();

	explicit builder(double value);

	builder(double lb, double ub);

	friend const builder operator+(const builder& x, const builder& y);

	friend const builder operator+(const builder& x, double y);

	friend const builder operator-(const builder& x, const builder& y);

	friend const builder operator-(double x, const builder& y);

	friend const builder operator*(const builder& x, const builder& y);

	friend const builder operator*(double x, const builder& y);

	friend const builder operator/(const builder& x, const builder& y);

	friend const builder sqr(const builder& x);

	void mark_as_common_subexpression() const;

	void equals(double value) const;

	static void record_occurence_info();

	static int number_of_arguments();

	static int number_of_variables();

	static const PrimVector& get_primitives();

	static const PairVector& get_numeric_constants();

	// TODO Figure out how to use this info?
	static const IntVector& get_common_subexpressions();

	static const PairVector& get_rhs_of_constraints();

	static const BoundVector& get_initial_box();

	static void reset();

	static void dbg_dump_type_of_primitives();

	static void dbg_show_info();

	void dbg_consistency() const;

private:

	static int last_constraint_offset();

	static void occurence_info_of_constraint(const int k);

	static int number_of_vars;

	static int unused_index;

	static PrimVector primitives;

	static PairVector numeric_constants;

	static IntVector common_subexpressions;

	static PairVector constraints_rhs;

	static BoundVector initial_box;

	int index;
};

const builder operator+(double x, const builder& y);

const builder operator-(const builder& x, double y);

void dbg_consistency(const builder& x, const builder& y);

}

#endif // BUILDER_HPP_
