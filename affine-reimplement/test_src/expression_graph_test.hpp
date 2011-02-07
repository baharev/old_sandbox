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

#ifndef EXPRESSION_GRAPH_TEST_HPP_
#define EXPRESSION_GRAPH_TEST_HPP_

namespace asol {

template <typename T> class problem;
class builder;

void dag_test(const problem<builder>* prob);

void print_sparsity(const problem<builder>* prob);

void test_probing_on_initial_box(const problem<builder>* prob);

void test_solutions_revise(const problem<builder>* prob);

void test_solutions_revise2(const problem<builder>* prob);

void test_solutions_iterative_revise(const problem<builder>* prob);

void test_solutions_probing(const problem<builder>* prob);

}

#endif // EXPRESSION_GRAPH_TEST_HPP_
