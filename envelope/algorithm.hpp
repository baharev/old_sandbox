//==============================================================================
//
// This code is part of ASOL (nonlinear system solver using affine arithmetic)
//
// Copyright (C) 2010 Ali Baharev
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

#ifndef ALGORITHM_HPP_
#define ALGORITHM_HPP_

#include <memory>
#include <deque>

namespace asol {

class var;
class interval;
class problem;

class algorithm {

public:

	explicit algorithm(const problem* const to_solve);

	void run();

	~algorithm();

private:

	algorithm(const algorithm& );
	algorithm& operator=(const algorithm& );

	void add_initial_box();
	void contracting_step();
	void delete_box();
	void evaluate();
	void increment_counters();
	void iteration_step();
	void lp_pruning();
	void prepare_to_repeat();
	void rollback();
	void set_current_box();
	void split();
	bool sufficient_progress();

	const int n;

	const std::auto_ptr<const problem> prob;

	var* const box;

	interval* box_orig;

	std::deque<interval*> pending;

	int depth;

	int boxes_processed;

};

}

#endif /* ALGORITHM_HPP_ */
