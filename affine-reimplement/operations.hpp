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

#ifndef OPERATIONS_HPP_
#define OPERATIONS_HPP_

namespace asol {

class operations {

public:

	virtual void addition(int z, int x, int y) = 0;

	virtual void addition_revise(int z, int x, int y) = 0;

	virtual void substraction(int z, int x, int y) = 0;

	virtual void substraction_revise(int z, int x, int y) = 0;

	virtual void multiplication(int z, int x, int y) = 0;

	virtual void multiplication_revise(int z, int x, int y) = 0;

	virtual void division(int z, int x, int y) = 0;

	virtual void division_revise(int z, int x, int y) = 0;

	virtual void square(int z, int x) = 0;

	virtual void square_revise(int z, int x) = 0;

	virtual void equality_constraint(int body, int rhs) = 0;

	virtual void equality_constraint_revise(int body, int rhs) = 0;

	virtual ~operations() { }
};

}

#endif // OPERATIONS_HPP_
