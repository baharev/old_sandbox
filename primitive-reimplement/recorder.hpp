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

#ifndef RECORDER_HPP_
#define RECORDER_HPP_

#include "primitives.hpp"

namespace asol {

class recorder {

public:

	virtual void addition(int z, int x, int y) = 0;
	virtual void substraction(int z, int x, int y) = 0;
	virtual void multiplication(int z, int x, int y) = 0;
	virtual void division(int z, int x, int y) = 0;
	virtual void square(int z, int x) = 0;
	virtual void exponential(int z, int x) = 0;
	virtual void equality_constraint(int z, int x) = 0;

	virtual ~recorder() { }
};

}

#endif // RECORDER_HPP_
