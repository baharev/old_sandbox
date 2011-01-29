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

#ifndef COMBINATION_HPP_
#define COMBINATION_HPP_

#include <vector>

namespace asol {

class combination {

public:

	typedef std::vector<int> IntVector;

	combination();

	explicit combination(int size, int parts_to_generate);

	bool step_counters();

	const IntVector& counters() const;

private:

	combination(const combination& );
	combination& operator=(const combination& );

	bool has_more_counters() const;
	bool counter_at_max() const;
	void next(const bool overflow);
	void handle_overflow();

	IntVector counter;
	int position;
	int counter_max;
	int size;
};

}

#endif // COMBINATION_HPP_
