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

#ifndef ENVELOPE_HPP_
#define ENVELOPE_HPP_

#include <iosfwd>

namespace lp_solver {

	class lp_impl;
}

namespace asol {

class var {

public:

	var(double lb, double ub);

	void fix_at(double val);

	friend const var operator+(const var& x, const var& y);

	friend std::ostream& operator<<(std::ostream& , const var& );

	static void dump_lp(const char* file);

private:

	int index;
	double lb;
	double ub;

	static lp_solver::lp_impl* lp;

};



}

#endif /* ENVELOPE_HPP_ */
