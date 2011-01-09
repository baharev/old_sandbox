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

#include <iosfwd>
#include <map>
#include <vector>
#include "primitives.hpp"
#include "typedefs.hpp"

namespace asol {

class index_set;

class builder {

public:

	builder();

	builder(double lb, double ub);

	friend const builder operator+(const builder& x, const builder& y);

	friend const builder operator+(const builder& x, double y);

	friend const builder operator-(const builder& x, const builder& y);

	friend const builder operator-(double x, const builder& y);

	friend const builder operator*(const builder& x, const builder& y);

	friend const builder operator*(double x, const builder& y);

	friend const builder operator/(const builder& x, const builder& y);

	friend const builder operator/(double x, const builder& y);

	friend const builder sqr(const builder& x);

	friend const builder exp(const builder& x);

	void mark_as_common_subexpression() const;

	void equals(double value) const;

	void assign(const builder& ) const { }

	static int number_of_arguments();

	static int number_of_variables();

	static const std::vector<primitive<builder>*>& get_primitives();

	static const std::map<int,double>& get_numeric_constants();

	// TODO Figure out how to use this info?
	static const IntVector& get_common_subexpressions();

	static const PairVector& get_rhs_of_constraints();

	static const BoundVector& get_initial_box();

	static void reset();

	static void print_primitives(std::ostream& out);

	static void print_index_set(std::ostream& out);

	static void print_type1_common_subexpressions(std::ostream& out);

	static void print_type2_common_subexpressions(std::ostream& out);

	static void print_type3_common_subexpressions(std::ostream& out);

	static void print_info(std::ostream& out);

	void dbg_consistency() const;

private:

	explicit builder(double value);

	static int last_constraint_offset();

	static void common_subexpressions_type1(const int i, std::ostream& out);

	static void insert_numeric_constant(const int index, const double value);

	static index_set* record_index_set();

	static int primitives_size();

	static int number_of_vars;

	static int unused_index;

	static std::vector<primitive<builder>*> primitives;

	static std::map<int,double> numeric_constants;

	static IntVector common_subexpressions;

	static PairVector constraints_rhs;

	static BoundVector initial_box;

	int index;
};

template <typename T>
extern const std::vector<primitive<T>*> convert(const std::vector<primitive<builder>*>& v);

const builder operator+(double x, const builder& y);

const builder operator-(const builder& x, double y);

void dbg_consistency(const builder& x, const builder& y);

void addition_inverse(builder& , builder& , builder& );

void substraction_inverse(builder& , builder& , builder& );

void multiplication_inverse(builder& , builder& , builder& );

void division_inverse(builder& , builder& , builder& );

void sqr_inverse(builder& , builder& );

void exp_inverse(builder& , builder& );

void equality_constraint_inverse(builder& , double );

}

#endif // BUILDER_HPP_
