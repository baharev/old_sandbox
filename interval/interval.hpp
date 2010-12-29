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

#ifndef INTERVAL_HPP_
#define INTERVAL_HPP_

#include <iosfwd>

namespace asol {

class interval {

public:

	interval();

	explicit interval(double value);

	interval(double lower_bound, double upper_bound);

	bool intersect(const double l, const double u);

	bool intersect(const interval& other);

	// Used by the expression graph
	void assign(const interval& other);
	void equals(double value);

	friend void copy_array(const interval src[], interval dstn[], int size);

	interval& operator+=(const interval& x);

	friend const interval operator+(const interval& x, const interval& y);

	friend const interval operator+(const interval& x, double y);

	friend const interval operator-(const interval& x);

	friend const interval operator-(const interval& x, const interval& y);

	friend const interval operator-(double x, const interval& y);

	friend const interval operator*(const interval& x, const interval& y);

	friend const interval operator*(double x, const interval& y);

	friend const interval operator/(const interval& x, const interval& y);

	friend const interval sqr(const interval& x);

	friend const interval sqrt(const interval& x);

	friend const interval exp(const interval& x);

	friend const interval log(const interval& x);

	bool degenerate() const;

	bool contains(double value) const;

	double midpoint() const;

	double diameter() const;

	double radius() const;

	double inf() const;

	double sup() const;

	friend std::ostream& operator<<(std::ostream& , const interval& );

private:

	double lb;

	double ub;
};

void addition_inverse(interval& z, interval& x, interval& y);

void substraction_inverse(interval& z, interval& x, interval& y);

void multiplication_inverse(interval& z, interval& x, interval& y);

void division_inverse(interval& z, interval& x, interval& y);

void sqr_inverse(interval& z, interval& x);

void exp_inverse(interval& z, interval& x);

void equality_constraint_inverse(interval& z, interval& x);

void propagate_mult(interval& z, interval& x, interval& y);

}

#endif
