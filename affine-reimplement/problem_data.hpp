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

#ifndef PROBLEM_DATA_HPP_
#define PROBLEM_DATA_HPP_

#include <iosfwd>
#include <map>
#include <vector>
#include "builder.hpp"
#include "primitives.hpp"
#include "typedefs.hpp"

namespace asol {

class index_set;

class problem_data {

public:

	typedef primitive<builder> Primitive;

	problem_data();

	~problem_data();

	void increment_index_counter();

	void add_variable(double lb, double ub);

	void add_primitive(const Primitive* p);

	void add_numeric_constant(int index, double value);

	void add_common_subexpression(int index);

	void add_constraint_rhs(double value);

	int number_of_arguments() const;

	int number_of_variables() const;

	typedef std::vector<Primitive*> VecPrimitive;

	typedef std::map<int,double> Map;

	const VecPrimitive& get_primitives() const;

	const Map& get_numeric_constants() const;

	// TODO Figure out how to use this info?
	const IntVector& get_common_subexpressions() const;

	const PairVector& get_rhs_of_constraints() const;

	const BoundVector& get_initial_box() const;

	void print_primitives(std::ostream& out) const;

	void print_index_set(std::ostream& out) const;

	void print_type1_common_subexpressions(std::ostream& out) const;

	void print_type2_common_subexpressions(std::ostream& out) const;

	void print_type3_common_subexpressions(std::ostream& out) const;

	void print_info(std::ostream& out) const;

private:

	problem_data(const problem_data& );

	problem_data& operator=(const problem_data& );

	int last_constraint_offset() const;

	void common_subexpressions_type1(const int i, std::ostream& out) const;

	void insert_numeric_constant(const int index, const double value);

	index_set* record_index_set() const;

	int primitives_size() const;

	int number_of_vars;

	VecPrimitive primitives;

	Map numeric_constants;

	IntVector common_subexpressions;

	PairVector constraints_rhs;

	BoundVector initial_box;
};

}

#endif // PROBLEM_DATA_HPP_
