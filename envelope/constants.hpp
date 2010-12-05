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

#ifndef CONSTANTS_HPP_
#define CONSTANTS_HPP_

namespace lp_solver {

	const double TOL_MIN_REL_DIAM = 1.0e-4; // TODO abs/rel width
}

namespace asol {

	const double TOL_SOLVED = 10*lp_solver::TOL_MIN_REL_DIAM;

	const double TOL_RANGE  = 1.0e-6;
}

#endif /* CONSTANTS_HPP_ */
