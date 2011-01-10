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

#ifndef COMBINATION_HPP_
#define COMBINATION_HPP_

#include <vector>
#include "index_range.hpp"

namespace asol {

class combination {

public:

	typedef std::vector<index_range> Vector;

	combination(const Vector& index_bound, int equal_parts);

	~combination();

private:

	combination(const combination& );
	combination& operator=(const combination& );

	void copy_if_not_narrow(const index_range ir);

	int* const index;

	int* const max_part;

	int* const counter;

	interval* const part;

	int length;

	int pos;
};

}

#endif // COMBINATION_HPP_
