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

#ifndef AFFINE_HPP_
#define AFFINE_HPP_

#include <iosfwd>
#include <vector>
#include "interval.hpp"

namespace asol {

struct epsilon {

	epsilon() : index(-1), coeff(1.0e+300) { }

	epsilon(int i, double c) : index(i), coeff(c) { }

	int    index;

	double coeff;
};

class lp_solver;

class affine {

public:

	affine() : range_index(-1) { }

	void set_range_index(int i) { range_index = i; }

	void make_variable();

	void recompute_variable(int i);

	void make_numeric_constant();

	void check_if_numeric_constant();

	void set_range_with_epsilon_bounds(const double lb, const double ub);

	void renormalize(const double lb, const double ub);

	~affine();

	void equals(double value);

	void less_than_or_equal_to(affine& rhs);

	friend void aa_exp(affine& z, const affine& x);

	friend void aa_log(affine& z, const affine& x);

	friend void aa_sqr(affine& z, const affine& x);

	friend void aa_reciprocal(affine& z, const affine& x);

	friend void aa_multiplication(affine& z, const affine& x, const affine& y);

	friend void aa_division(affine& z, const affine& x, const affine& y);

	friend std::ostream& operator<<(std::ostream& , const affine& );

	void dbg_consistency() const;

	friend void dbg_consistency(const affine& z, const affine& x);

	static void set_vector(std::vector<interval>* vec);

	static void set_lp_solver(lp_solver* lp);

	static void reset_counter();

	friend class affine_pair_iterator;

	template <typename> friend class binary_operation;

	friend class lp_solver;

	const interval range() const { return get_range(); }

	const interval merge() const;

private:

	explicit affine(const interval& range) : range_index(-1), ia_range(range) { }

	int size() const { return static_cast<int>(noise_vars.size()); }

	interval& get_range() {	return range_index>=0 ? v->at(range_index) : ia_range; }

	const interval& get_range() const { return range_index>=0 ? v->at(range_index) : ia_range; }

	double  central_value() const { return noise_vars.at(0).coeff; }

	double& central_value()       { return noise_vars.at(0).coeff; }

	void force_intersection(double lb, double ub) { get_range().force_intersection(lb, ub); }

	void force_intersection(const interval& new_range) { get_range().force_intersection(new_range); }

	void unary_op(const affine& x, double alpha, double zeta, double delta);

	void add_noise_var(int index, double coeff);

	void set_var_range(int index, const interval& rng);

	void condense_last_two_noise_vars();

	std::vector<epsilon> noise_vars;

	int range_index;

	interval ia_range;

	static int max_used_index;

	static std::vector<interval>* v;

	static lp_solver* lp;

	static const double NARROW;
};

void aa_addition(affine& z, const affine& x, const affine& y);

void aa_substraction(affine& z, const affine& x, const affine& y);

void dbg_consistency(const affine& z, const affine& x, const affine& y);

}

#endif // AFFINE_HPP_
