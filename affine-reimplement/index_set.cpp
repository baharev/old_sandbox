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

#include <ostream>
#include "index_set.hpp"
#include "primitives.hpp"

namespace asol {

// FIXME Pass in numeric constants
index_set::index_set(int max_variable_index) : max_var_index(max_variable_index) {

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

void index_set::record_arg(const int index) {

	// FIXME Bug: also filters out common subexpressions!!! E.g. Eq15 v48
	if (index <= max_var_index) {

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

	out << "Constraint " << k << std::endl;

	const Set* const idx = constraint_index_sets.at(k);

	Set::const_iterator i = idx->begin();

	while (i!=idx->end()) {

		out << *i << std::endl;

		++i;
	}

	out << std::endl;
}

void index_set::print(std::ostream& out) const {

	const int n = number_of_constraints();

	for (int i=0; i<n; ++i) {

		print_constraint(i, out);
	}
}

}
