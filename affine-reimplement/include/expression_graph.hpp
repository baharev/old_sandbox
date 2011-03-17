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
#include <vector>
#include "containment.hpp"
#include "typedefs.hpp"

namespace asol {

template <typename T> class primitive;
class problem_data;
template <typename T> struct gap_info;

template <typename T>
class expression_graph {

public:

	expression_graph(const problem_data* problem);

	void set_box(const T* box, const int length);

	const T* get_box() const;

	containment<T> contains(const std::vector<double>& solution) const;

	void evaluate_all();

	void revise_all();

	void revise_all2(); // TODO Remove

	void iterative_revision();

	void iterative_revision_save_gaps();

	void probing();

	void show_variables(std::ostream& out) const;

	const T& last_value() const;

	~expression_graph();

private:

	expression_graph(const expression_graph& );
	expression_graph& operator=(const expression_graph& );

	int constraints_size() const { return static_cast<int> (constraints.size()); }

	void set_variables();
	void set_non_variables();
	void set_numeric_consts();
	void set_constant(const Pair p);

	int constraint_begin(int i) const;
	int constraint_end(int i) const;

	int first_not_strictly_contained(const std::vector<double>& solution) const;
	int first_not_easily_contained(const std::vector<double>& solution, int from) const;

	void evaluate_up_to(const int i);
	void evaluate_constraint(int i);

	void revise_up_to(const int i);
	void revise_constraint(int i);

	void probe_in_constraint(const int i);
	void revise_up_to_with_hull_saved(const int i);
	void save_current_as_orig();
	void set_orig_as_v();
	void save_hull();
	void write_hull_to_orig();
	bool intersect_hull_and_orig();
	bool compute_intersection();

	typedef std::vector<primitive<T>*> PrimVector;

	std::vector<T> v;
	const int n_vars;
	const PrimVector primitives;
	const PairVector constants;
	const BoundVector initial_box;
	const IntArray2D index_sets;
	const IntVector constraints;

	std::vector<T> orig;
	std::vector<T> hull;

	std::vector<gap_info<T> > gaps;
};

}

#endif // EXPRESSION_GRAPH_HPP_
