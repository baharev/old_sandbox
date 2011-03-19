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

#ifndef WILSON16_HPP_
#define WILSON16_HPP_

#include "problem.hpp"

namespace asol {

template <typename T>
class Wilson16 : public problem<T> {

private:

	virtual int number_of_variables() const;

	virtual T* initial_box() const;

	virtual void evaluate(const T x[]) const;

	virtual int number_of_stored_solutions() const;

	virtual const DoubleArray2D solutions() const;

	enum { X1, X2, X3, t,
		   S1, S2, S3,
		   T1, T2, T3,
		   U1, U2, U3,
		   LN_K1, LN_K2, LN_K3,
		   SIZE, SOLS = 7 };
	};
};

#endif // WILSON16_HPP_
