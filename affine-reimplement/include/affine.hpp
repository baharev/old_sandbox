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
#include <interval.hpp>

namespace asol {

struct epsilon {

	epsilon(int i, double c) : index(i), coeff(c) { }

	int    index;

	double coeff;
};

class affine {

public:

	affine();

	explicit affine(double value);

	affine(double lb, double ub);

	~affine();

	void assign(const affine& other);

	void equals(double value);

	void less_than_or_equal_to(affine& rhs);

	friend const affine exp(const affine& x);

	friend const affine log(const affine& x);

	friend const affine sqr(const affine& x);

	friend const affine operator+(const affine& x, const affine& y);

	friend const affine operator-(const affine& x, const affine& y);

	friend const affine operator*(const affine& x, const affine& y);

	friend const affine operator/(const affine& x, const affine& y);

	friend std::ostream& operator<<(std::ostream& , const affine& );

private:

	int size() const { return static_cast<int>(noise_vars.size()); }

	const interval range() const { return v->at(range_index); }

	std::vector<epsilon> noise_vars;

	int range_index;

	static int max_used_index;

	static std::vector<interval>* v;
};

}

#endif // AFFINE_HPP_
