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

#ifndef LP_SOLVER_HPP_
#define LP_SOLVER_HPP_

namespace asol {

class affine;
class lp_impl;

class lp_solver {

public:

	lp_solver();

	~lp_solver();

	static void free_environment();

private:

	lp_solver(const lp_solver& );
	lp_solver& operator=(const lp_solver& );

	lp_impl* lp;

};

}

#endif // LP_SOLVER_HPP_
