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

class affine {

public:

	affine() : range_index(-1) { }

	void set_range_index(int i) { range_index = i; }

	void make_variable();

	void recompute_variable(int i);

	void make_numeric_constant();

	void check_if_numeric_constant();

	~affine();

	void assign(const affine& other) { operator=(other); }

	void equals(double value);

	void less_than_or_equal_to(affine& rhs);

	friend void aa_exp(affine& z, const affine& x);

	friend void aa_log(affine& z, const affine& x);

	friend void aa_sqr(affine& z, const affine& x);

	friend void aa_addition(affine& z, const affine& x, const affine& y);

	friend const affine operator-(const affine& x, const affine& y);

	friend const affine operator*(const affine& x, const affine& y);

	friend const affine operator/(const affine& x, const affine& y);

	friend std::ostream& operator<<(std::ostream& , const affine& );

	void dbg_consistency() const;

	friend void unary_op(affine& z, const affine& x, double alpha, double zeta, double delta);

	static void set_vector(std::vector<interval>* vec);

	static void reset_counter();

private:

	int size() const { return static_cast<int>(noise_vars.size()); }

	const interval range() const { return v->at(range_index); }

	void add_noise_var(int index, double coeff);

	void set_var_range(int index, const interval& rng);

	std::vector<epsilon> noise_vars;

	int range_index; // TODO Range index can be eliminated

	static int max_used_index;

	static std::vector<interval>* v;
};

void dbg_consistency(const affine& x, const affine& y);

}

#endif // AFFINE_HPP_
