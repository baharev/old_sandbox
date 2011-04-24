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

namespace asol {

class affine {

public:

	affine();

	explicit affine(double value);

	affine(double lb, double ub);

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

};

std::ostream& operator<<(std::ostream& , const affine& );

}

#endif // AFFINE_HPP_