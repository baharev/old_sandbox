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

#ifndef SOL_TRACKER_HPP_
#define SOL_TRACKER_HPP_

#include <map>
#include "problem.hpp"

namespace asol {

class builder;
class interval;

class sol_tracker {

public:

	explicit sol_tracker(const problem<builder>* prob);

	void save_containment_info(const interval* box);

private:

	enum type { NOT, EASY, STRICT };

	const int n_vars;
	DoubleArray2D solutions;
	std::map<int,type> containment;
};

}

#endif // SOL_TRACKER_HPP_
