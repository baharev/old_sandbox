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

#ifndef CONTAINMENT_HPP_
#define CONTAINMENT_HPP_

namespace asol {

enum type { NOT_CONT, EASY_CONT, STRICT_CONT };

template <typename T>
class containment {

public:

	containment(type t, int index, const T& val = T()) : t(t), i(index), v(val) { }

	bool strict()   const { return t == STRICT_CONT; }
	bool easy()     const { return t == EASY_CONT; }
	bool no_sol()   const { return t == NOT_CONT; }
	int  index()    const { return i; }
	const T value() const { return v; }

private:

	type t;
	int  i;
	T    v;
};

}

#endif // CONTAINMENT_HPP_
