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

#ifndef ENVELOPE_HPP_
#define ENVELOPE_HPP_

#include <iosfwd>
#include "interval.hpp"

namespace asol {

class lp_pair;
class dag;

class var {

public:

	var();

	var(double lb, double ub);

	void fix_at(double val);
	// FIXME Make it private - Only digression uses it
	bool tighten_bounds();
	// FIXME Make it private - Only digression uses it
	void intersect(double lb, double ub);

	double width() const;

	bool contains_zero() const;

	void check_consistency() const;
	// FIXME Remove it - Only digression uses it
	void copy_bounds(double& lb, double& ub) const;

	friend void copy_bounds(const var arr[], interval bounds[], int size);

	friend void init_variables(var x[], const interval var_bounds[], int size);

	friend int find_max_width(const var x[], int size); // TODO abs/rel width

	friend const var operator+(const var& x, const var& y);

	friend const var operator+(const var& x, double y);

	friend const var operator-(double x, const var& y);

	friend const var operator-(const var& x, const var& y);

	friend const var operator*(const var& x, const var& y);

	friend const var operator*(const double c, const var& x);

	friend const var operator/(const var& x, const var& y);

	friend const var lin_comb_n(const double c[], const var x[], int length);

	// TODO Implement!
	friend void lin_con_n(double val, const double c[], const var x[], int n);

	friend const var sqr(const var& x);

	friend const var y_eq(const var& x);

	friend const var H_Liq(const var& x);

	friend const var H_Vap(const var& x);

	friend std::ostream& operator<<(std::ostream& , const var& );

	static void tighten_all();

	static void tighten_up_to(int index);

	static void dump_lp(const char* file);

	static void reset(); // FIXME Make it private

	static void release_all();

private:

	var(const interval& range);

	bool intersect_in_dag(const interval& range);
	const interval bounds() const;
	const interval compute_bounds() const;
	const interval lp_tighten_col(bool& improved) const;
	friend const interval sum(const double c[], const var x[], int n);

	int index;

	static lp_pair* lp;
	static dag* ia_dag;
};

const var operator+(double x, const var& y);

const var operator-(const var& x, double y);

const var operator*(const var& x, const double c);

// sum c*x = value
template <int n>
void linear_constraint(const double (&c)[n], const var (&x)[n], double value) {
	return lin_con_n(value, c, x, n);
}

// z = sum c*x
template <int n>
const var linear_combination(const double (&c)[n], const var (&x)[n]) {
	return lin_comb_n(c, x, n);
}

}

#endif /* ENVELOPE_HPP_ */
