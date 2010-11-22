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

#include <string>
#include <assert.h>
#include "lp_pair.hpp"
#include "lp_impl.hpp"

namespace lp_solver {

lp_pair::lp_pair() : lp_min(new lp_impl), lp_max(new lp_impl) {

	for (int i=0; i<5; ++i) {
		    rows[i] = -1;
		stat_min[i] = -1;
		stat_max[i] = -1;
	}
}

lp_pair::~lp_pair() {

	delete lp_min;
	delete lp_max;
}

void lp_pair::dump(const char* file) {

	using std::string;

	string fmin(file);
	string fmax(file);

	fmin += "_min";
	fmax += "_max";

	lp_min->dump(fmin.c_str());
	lp_max->dump(fmax.c_str());
}

int lp_pair::add_col_nonbasic(double lb, double ub) {

	assert(lb <= ub);

	int index_min = lp_min->add_col_nonbasic_on_lb(lb, ub);
	int index_max = lp_max->add_col_nonbasic_on_ub(lb, ub);

	assert(index_min == index_max);

	return index_min;
}

// x + y = c
void lp_pair::add_sum_row(int x_index, int y_index, double c) {

	lp_min->add_sum_row(x_index, y_index, c);
	lp_max->add_sum_row(x_index, y_index, c);
}

// x + y - z = 0
void lp_pair::add_add_row(int x_index, int y_index, int z_index) {

	lp_min->add_add_row(x_index, y_index, z_index);
	lp_max->add_add_row(x_index, y_index, z_index);
}

// x - z = -y
void lp_pair::add_shift_row(int x_index, int z_index, double y) {

	lp_min->add_shift_row(x_index, z_index, y);
	lp_max->add_shift_row(x_index, z_index, y);
}

// x - y - z = 0
void lp_pair::add_sub_row(int x_index, int y_index, int z_index) {

	lp_min->add_sub_row(x_index, y_index, z_index);
	lp_max->add_sub_row(x_index, y_index, z_index);
}

// c <= ax - z
void lp_pair::add_lo_row(double a, int x_index, int z_index, double c) {

	lp_min->add_lo_row(a, x_index, z_index, c);
	lp_max->add_lo_row(a, x_index, z_index, c);
}

// ax - z <= c
void lp_pair::add_up_row(double a, int x_index, int z_index, double c) {

	lp_min->add_up_row(a, x_index, z_index, c);
	lp_max->add_up_row(a, x_index, z_index, c);
}

// c <= ax + by - z
int lp_pair::add_lo_row(double a, int x, double b, int y, int z, double c) {

	int i_min = lp_min->add_lo_row(a, x, b, y, z, c);
	int i_max = lp_max->add_lo_row(a, x, b, y, z, c);
	assert(i_min == i_max);
	return i_min;
}

// ax + by - z <= c
int lp_pair::add_up_row(double a, int x, double b, int y, int z, double c) {

	int i_min = lp_min->add_up_row(a, x, b, y, z, c);
	int i_max = lp_max->add_up_row(a, x, b, y, z, c);
	assert(i_min == i_max);
	return i_min;
}

void lp_pair::fix_col(int index, double value) {

	lp_min->fix_col(index, value);
	lp_max->fix_col(index, value);
}

void lp_pair::save_row_status() {

	lp_min->get_row_status(rows, stat_min);
	lp_max->get_row_status(rows, stat_max);
}

void lp_pair::remove_mult_envelope() {

	lp_min->remove_envelope(rows);
	lp_max->remove_envelope(rows);
}

void lp_pair::restore_row_status() {

	lp_min->set_row_status(rows, stat_min);
	lp_max->set_row_status(rows, stat_max);
}

void lp_pair::add_mult_envelope(
		int x,
		double xL,
		double xU,
		int y,
		double yL,
		double yU,
		int z,
		bool reset)
{

	if (reset) {
		save_row_status();
		remove_mult_envelope();
	}

	// yL*xU <= yL*x + xU*y - z
	rows[1] = add_lo_row(yL, x, xU, y, z, yL*xU);

	// yU*xL <= yU*x + xL*y - z
	rows[2] = add_lo_row(yU, x, xL, y, z, yU*xL);

	// yL*x + xL*y - z <= yL*xL
	rows[3] = add_up_row(yL, x, xL, y, z, yL*xL);

	// yU*x + xU*y - z <= yU*xU
	rows[4] = add_up_row(yU, x, xU, y, z, yU*xU);

	if (reset) {
		restore_row_status();
	}
}

// c*x - z = 0
void lp_pair::add_cx_row(double c, int x_index, int z_index) {

	lp_min->add_cx_row(c, x_index, z_index);
	lp_max->add_cx_row(c, x_index, z_index);
}

void lp_pair::set_bounds(int index, double lb, double ub) {

	lp_min->set_bounds(index, lb, ub);
	lp_max->set_bounds(index, lb, ub);
}

void lp_pair::fix_narrow_slight_infeas(int index, double& lb, double& ub) {

	if (col_should_be_fixed(lb, ub)) {
		double mid = (lb+ub)/2;
		lb = ub = mid;
	}
}

bool lp_pair::tighten_col(int index, double& lb, double& ub) {

	if (lp_min->is_fixed(index)) {
		assert(lp_max->is_fixed(index));
		return false;
	}

	bool min_improved = lp_min->tighten_col_lb(index, lb);
	bool max_improved = lp_max->tighten_col_ub(index, ub);

	bool improved = min_improved || max_improved;

	if (improved) {
		fix_narrow_slight_infeas(index, lb, ub);
		set_bounds(index, lb, ub);
	}

	return improved;
}

bool lp_pair::col_type_db_or_fx(int index) const {

	return lp_min->col_type_db_or_fx(index) && lp_max->col_type_db_or_fx(index);
}

void lp_pair::reset() {

	lp_min->reset();
	lp_max->reset();
}

void lp_pair::free_environment() {

	lp_impl::free_environment();
}

}
