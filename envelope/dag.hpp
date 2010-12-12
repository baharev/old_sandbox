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

#ifndef DAG_HPP_
#define DAG_HPP_

#include <vector>

namespace asol {

class interval;

class dag {

public:

	dag();

	void add(int index, const interval& bounds);

	void add(int index, double lb, double ub);

	const interval bounds(int index) const;

	bool intersect(int index, const interval& other);

	void reset();

private:

	dag(const dag& );
	dag& operator=(const dag& );

	void init();

	std::vector<interval> variables;
};

}

#endif /* DAG_HPP_ */
