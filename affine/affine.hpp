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

#ifndef AFFINE_HPP_
#define AFFINE_HPP_

#include <iosfwd>
#include "interval.hpp"

namespace asol {

class affine {

public:

	affine(double lower_bound, double upper_bound);

	affine(const affine& other);

	friend const affine operator+(const affine& x, const affine& y);

	friend const affine operator*(const affine& x, const affine& y);

	void check_consistency() const;

	static void reset_max_used_index();

	friend std::ostream& operator<<(std::ostream& , const affine& );

	~affine();

	//==========================================================================

	static void set_max_used_index(int idx);

	static int get_max_used_index();

private:

	affine& operator=(const affine& other);

	affine(bool, int size);

	int index_at(int i) const;

	friend bool next(const affine& x, const affine& y, int& idx, double& xi, double& yi);

	int* const index;

	double* const value;

	int n;

	interval range;

	static int max_used_index;

};

void check_consistency(const affine& x, const affine& y);

}

#endif
