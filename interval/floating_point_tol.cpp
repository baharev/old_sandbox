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

#include <algorithm>
#include <cmath>

namespace asol {

double add_tol(const double x, const double TOLERANCE) {

	const double abs_tol = x + TOLERANCE;

	const double rel_tol = x + TOLERANCE*std::fabs(x);

	return std::max(abs_tol, rel_tol);
}

double sub_tol(const double x, const double TOLERANCE) {

	const double abs_tol = x - TOLERANCE;

	const double rel_tol = x - TOLERANCE*std::fabs(x);

	return std::min(abs_tol, rel_tol);
}

}
