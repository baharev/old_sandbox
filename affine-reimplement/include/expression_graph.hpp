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

#ifndef EXPRESSION_GRAPH_HPP_
#define EXPRESSION_GRAPH_HPP_

#include <iosfwd>
#include <map>
#include <vector>
#include "typedefs.hpp"

namespace asol {

template <typename T> class primitive;
class problem_data;

template <typename T>
class expression_graph {

public:

	expression_graph(const problem_data* problem);

	void evaluate_all();

	void revise_all();

	void evaluate_constraint(int i); // TODO Make it private?

	void revise_constraint(int i); // TODO Make it private?

	void revise_all2(); // TODO Remove

	void evaluate_up_to(const int i);

	void revise_up_to(const int i);

	void probing();

	void show_variables(std::ostream& out) const;

	const T& last_value() const;

	void evaluate_primitive(int i);

	~expression_graph();

private:

	expression_graph(const expression_graph& );
	expression_graph& operator=(const expression_graph& );

	int v_size() const { return static_cast<int> (v.size()); }
	int constraints_size() const { return static_cast<int> (constraints.size()); }
	int primitives_size() const { return static_cast<int> (primitives.size()); }

	void set_variables();
	void set_non_variables(const int length);
	void set_numeric_consts(const int length);

	int constraint_begin(int i) const;
	int constraint_end(int i) const;

	typedef std::vector<std::vector<int> > IntArray2D;
	typedef std::vector<primitive<T>*> PrimVector;
	typedef std::map<int,double> Map;

	std::vector<T> v;
	const PrimVector primitives;
	const Map constants;
	const BoundVector initial_box;
	const IntArray2D index_sets;
	const IntVector constraints;
};

}

#endif // EXPRESSION_GRAPH_HPP_
