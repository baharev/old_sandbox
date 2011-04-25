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

#include "affine.hpp"
#include "builder.hpp"
#include "diagnostics.hpp"
#include "expression_graph.hpp"
#include "interval.hpp"
#include "problem.hpp"
#include "problem_data.hpp"

namespace asol {

// FIXME Code triplication
const problem_data* build_repr(const problem<builder>* prob) {

	builder::reset();

	builder* x = prob->initial_box();

	prob->evaluate(x);

	delete[] x;

	builder::finished();

	const problem_data* representation = builder::get_problem_data();

	ASSERT(representation->number_of_variables() == prob->number_of_variables());

	return representation;
}

void affine_expr_graph_test(const problem<builder>* prob) {

	const problem_data* representation = build_repr(prob);

	expression_graph<interval> ia_dag(representation, DoubleArray2D());

	ia_dag.evaluate_all();

	affine::set_vector(ia_dag.get_v());

	expression_graph<affine> aa_dag(representation, IntArray2D());

	aa_dag.evaluate_all();

	aa_dag.reset_vars();

	aa_dag.evaluate_all();

	delete prob;

	builder::reset();
}

}
