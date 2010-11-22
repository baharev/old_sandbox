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

#ifndef LP_PAIR_HPP_
#define LP_PAIR_HPP_

namespace lp_solver {

class lp_impl;

class lp_pair {

public:

	explicit lp_pair();

	int add_col_nonbasic(double lb, double ub);

	// x + y = c
	void add_sum_row(int x_index, int y_index, double c);

	// x + y - z = 0
	void add_add_row(int x_index, int y_index, int z_index);

	// x - z = -y
	void add_shift_row(int x_index, int z_index, double y);

	// x - y - z = 0
	void add_sub_row(int x_index, int y_index, int z_index);

	// c <= ax - z
	void add_lo_row(double a, int x_index, int z_index, double c);

	// ax - z <= c
	void add_up_row(double a, int x_index, int z_index, double c);

	// c <= ax + by - z
	int add_lo_row(double a, int x, double b, int y, int z, double c);

	// ax + by - z <= c
	int add_up_row(double a, int x, double b, int y, int z, double c);

	void add_mult_envelope(int x, double xL, double xU, int y, double yL, double yU, int z, bool reset);

	// c*x - z = 0
	void add_cx_row(double c, int x_index, int z_index);

	void fix_col(int index, double value);

	void set_bounds(int index, double lb, double ub);

	bool tighten_col(int index, double& lb, double& ub);

	bool col_type_db_or_fx(int index) const;

	void reset();

	void dump(const char* filename);

	~lp_pair();

	static void free_environment();

private:

	lp_pair(const lp_pair& );
	lp_pair& operator=(const lp_pair& );

	void fix_narrow_slight_infeas(int index, double& lb, double& ub);
	void save_row_status();
	void remove_mult_envelope();
	void restore_row_status();

	lp_solver::lp_impl* const lp_min;
	lp_solver::lp_impl* const lp_max;

	int rows[5];
	int stat_min[5];
	int stat_max[5];
};

}

#endif /* LP_PAIR_HPP_ */
