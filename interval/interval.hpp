//==============================================================================
//
// This code is part of ASOL (nonlinear system solver using affine arithmetic)
//
// Copyright (C) 2010, 2011, Ali Baharev
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

	bool intersect(const interval& x) { return intersect(x.lb, x.ub); } // Profiler says...

	bool prechecked_intersection(const double l, const double u);

	// Used by the primitives
	void assign(const interval& other) { intersect(other); } // Profiler says...
	void equals(double value);
	void less_than_or_equal_to(interval& rhs);

	friend void copy_array(const interval src[], interval dstn[], int size);

	interval& operator+=(const interval& x);

	interval& operator+=(double x);

	friend const interval operator+(const interval& x, const interval& y);

	friend const interval operator+(const interval& x, double y);

	friend const interval operator-(const interval& x);

	friend const interval operator-(const interval& x, const interval& y);

	friend const interval operator-(double x, const interval& y);

	friend const interval operator*(const interval& x, const interval& y);

	friend const interval operator*(double x, const interval& y);

	friend const interval operator/(const interval& x, const interval& y);

	friend bool true_extended_division(const interval& x, const interval& y, interval& z, interval& gap);

	friend const interval sqr(const interval& x);

	friend const interval sqrt(const interval& x);

	friend const interval exp(const interval& x);

	friend const interval log(const interval& x);

	friend const interval hull_of(const interval& x, const interval& y);

	bool true_subset_of(const interval& other) const;

	bool subset_of(const interval& other) const;

	friend bool disjoint(const interval& x, const interval& y);

	friend bool lessByLb(const interval& x, const interval& y);

	bool is_narrow() const { return is_narrow(1.0e-4); } // FIXME Magic number NARROW

	bool is_narrow(const double TOLERANCE) const; // FIXME Only for testing

	bool degenerate() const; // Convex envelopes use it

	bool contains(double value) const;

	friend bool easy_containment(double x, const interval& y);

	double midpoint() const;

	double diameter() const;

	double radius() const;

	double inf() const;

	double sup() const;

	double unchecked_inf() const;

	double unchecked_sup() const;

	bool valid() const;

	friend std::ostream& operator<<(std::ostream& , const interval& );

	static const interval ANY_REAL() { return interval(-1.0e+150, 1.0e+150); }

private:

	double lb;

	double ub;
};

const interval operator*(const interval& x, double y);

bool extended_division(const interval& x, const interval& y, interval& z, interval& gap);

void addition_inverse(interval& z, interval& x, interval& y);

void substraction_inverse(interval& z, interval& x, interval& y);

bool division_inverse(interval& z, interval& x, interval& y, interval& gap);

bool sqr_inverse(interval& z, interval& x, interval& gap);

void exp_inverse(interval& z, interval& x);

void log_inverse(interval& z, interval& x);

void equality_constraint_inverse(interval& z, double rhs);

void propagate_mult(interval& z, interval& x, interval& y);

const interval hull_of(const interval& x, const interval& y);

bool lessByLb(const interval& x, const interval& y);

const interval intersection(const interval& x, const interval& y);

}

#endif
