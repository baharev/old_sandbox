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

#ifndef BOX_GENERATOR_HPP_
#define BOX_GENERATOR_HPP_

#include <memory>
#include <vector>
#include "interval.hpp"

namespace asol {

class combination;

class box_generator {

public:

	typedef std::vector<interval> IVector;
	typedef std::vector<int> IntVector;

	box_generator(IVector& v, const IntVector& index_set, int equal_parts);

	bool empty() const;

	bool get_next();

	void set_box();

	~box_generator();

private:

	typedef std::vector<IVector> IArray2D;

	box_generator(const box_generator& );
	box_generator& operator=(const box_generator& );

	void reserve(int index_set_size);
	void generate_parts(int i);
	void cut_into_equal_parts(const double LB, const double UB);

	IVector& v;
	const int parts_to_generate;

	IntVector index;
	IArray2D parts;
	std::auto_ptr<combination> index_generator;
};

}

#endif // BOX_GENERATOR_HPP_
