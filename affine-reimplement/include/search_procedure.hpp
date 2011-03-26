//==============================================================================
//
// This code is part of ASOL (nonlinear system solver using affine arithmetic)
//
// Copyright (C) 2011 Ali Baharev
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

#ifndef SEARCH_PROCEDURE_HPP_
#define SEARCH_PROCEDURE_HPP_

#include <deque>
#include <vector>

namespace asol {

template <typename T> class expression_graph;
template <typename T> class problem;
class builder;
class interval;
class problem_data;

class search_procedure {

public:

	search_procedure(const problem<builder>* problem_to_solve);

	void run();

	~search_procedure();

private:

	search_procedure(const search_procedure& );
	search_procedure& operator=(const search_procedure& );

	void build_problem_representation();
	void evaluate_with_builder() const;
	void push_initial_box_to_deque();

	bool has_more_boxes() const;
	void get_next_box();
	void process_box();
	void split_if_not_discarded();
	void print_statistics() const;

	void iteration_step();
	bool not_done_with_box() const;
	bool sufficient_progress();
	bool sufficient(const double max_progress) const;
	double compute_max_progress() const;
	void split();
	int select_index_to_split() const;

	void delete_box();
	void print_box() const;
	void contracting_step();
	void check_convergence();

	void dbg_check_infeasibilty() const;
	void dbg_solution_count();
	void dbg_initial_box_from_dump();

	const problem<builder>* prob;

	const int n_vars;

	const problem_data* representation;

	expression_graph<interval>* expr_graph;

	std::deque<interval*> pending_boxes;

	interval* box_orig;

	int solutions_found;

	int splits;

	int boxes_processed;
};

}

#endif // SEARCH_PROCEDURE_HPP_
