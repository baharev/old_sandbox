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

#include "sol_tracker.hpp"
#include "builder.hpp"
#include "diagnostics.hpp"
#include "interval.hpp"

using namespace std;

namespace asol {

sol_tracker::sol_tracker(const problem<builder>* prob)
: n_vars(prob->number_of_variables()), solutions(prob->solutions())
{
	// TODO Make number_of_stored_solutions() return size()
	const int size = static_cast<int>(solutions.size());
	ASSERT2(size==prob->number_of_stored_solutions(),"size: "<<size);
}

void sol_tracker::save_containment_info(const interval* box) {

	// TODO Move containment-related methods from expression_graph to here
}

}
