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

#include <algorithm>
#include <iterator>
#include <ostream>
#include "index_set.hpp"
#include "diagnostics.hpp"
#include "primitives.hpp"

namespace asol {

index_set::index_set(const int num_of_vars, const Map& numeric_constants)
: number_of_variables(num_of_vars), numeric_const(numeric_constants)
{
	current = new Set;
}

void index_set::finished() {

	delete current;

	current = 0;
}

int index_set::number_of_constraints() const {

	return static_cast<int> (constraint_index_sets.size());;
}

index_set::~index_set() {

	const int n = number_of_constraints();

	for (int i=0; i<n; ++i) {

		delete constraint_index_sets.at(i);
	}

	delete current;
}

void index_set::record_primitive(const primitive* p) {

	current->insert(p->z);

	record_arg(p->x);
}

bool index_set::is_numeric_constant(const int index) const {

	Map::const_iterator i = numeric_const.find(index);

	return (i!=numeric_const.end()) ? true : false;
}

void index_set::record_arg(const int index) {

	if (!is_numeric_constant(index)) {

		current->insert(index);
	}
}

void index_set::record_binary_primitive(const binary_primitive* p) {

	record_primitive(p);

	record_arg(p->y);
}

void index_set::record(const addition* p) {

	record_binary_primitive(p);
}

void index_set::record(const substraction* p) {

	record_binary_primitive(p);
}

void index_set::record(const multiplication* p) {

	record_binary_primitive(p);
}

void index_set::record(const division* p) {

	record_binary_primitive(p);
}

void index_set::record(const square* p) {

	record_primitive(p);
}

void index_set::record(const exponential* p) {

	record_primitive(p);
}

void index_set::record(const equality_constraint* p) {

	constraint_index_sets.push_back(current);

	current = new Set;
}

void index_set::print_constraint(const int k, std::ostream& out) const {

	out << "Constraint " << k << '\n';

	const Set* const idx = constraint_index_sets.at(k);

	std::copy(idx->begin(), idx->end(), std::ostream_iterator<int>(out, "\n"));

	out << '\n' << std::flush;
}

void index_set::print(std::ostream& out) const {

	const int n = number_of_constraints();

	for (int i=0; i<n; ++i) {

		print_constraint(i, out);
	}
}

const std::set<int> index_set::non_variables(const int i) const {

	const Set* s_i = constraint_index_sets.at(i);

	return Set(s_i->lower_bound(number_of_variables), s_i->end());
}

void index_set::check_for_common_subexpressions(const int i) {

	const Set si(non_variables(i));

	for (int j=i+1; j<number_of_constraints(); ++j) {

		const Set sj(non_variables(j));

		Set tmp;

		std::set_intersection(si.begin(), si.end(),
							  sj.begin(), sj.end(),
							  std::inserter(tmp, tmp.begin()) );

		type2_cse.insert(tmp.begin(), tmp.end());
	}
}

void index_set::collect_type2_common_subexpressions() {

	ASSERT2(current==0,"recording not finished or not run");

	for (int i=0; i<number_of_constraints(); ++i) {

		check_for_common_subexpressions(i);
	}
}

const std::set<int>& index_set::type2_common_subexpressions() const {

	return type2_cse;
}

}
