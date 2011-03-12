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

#include "expression_graph.hpp"
#include "interval.hpp"
#include "gap_probing.hpp"

namespace asol {

gap_probing::gap_probing(expression_graph<interval>& g, interval* initial_box, int length)
: graph(g), n_vars(length), box(0)
{
	pending.push_back(initial_box);
}

interval* gap_probing::contracted_box() {

	box = pending.back();

	graph.set_box(box, n_vars);

	graph.iterative_revision_save_gaps();

	return box;
}

}

