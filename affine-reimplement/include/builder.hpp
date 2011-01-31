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

#ifndef BUILDER_HPP_
#define BUILDER_HPP_

namespace asol {

class problem_data;

class builder {

public:

	builder();

	builder(double lb, double ub);

	friend const builder operator+(const builder& x, const builder& y);

	friend const builder operator+(const builder& x, double y);

	friend const builder operator-(const builder& x, const builder& y);

	friend const builder operator-(double x, const builder& y);

	friend const builder operator*(const builder& x, const builder& y);

	friend const builder operator*(double x, const builder& y);

	friend const builder operator/(const builder& x, const builder& y);

	friend const builder operator/(double x, const builder& y);

	friend const builder sqr(const builder& x);

	friend const builder exp(const builder& x);

	void mark_as_common_subexpression() const;

	void equals(double value) const;

	void assign(const builder& ) const { }

	void dbg_consistency() const;

	static const problem_data* get_problem_data();

	static void add_solution(const double* sol, const int length);

	static void finished();

	static void reset();

	static void release();

private:

	explicit builder(double );

	static problem_data* problem;

	int index;
};

const builder operator+(double x, const builder& y);

const builder operator-(const builder& x, double y);

void dbg_consistency(const builder& x, const builder& y);

}

#endif // BUILDER_HPP_
