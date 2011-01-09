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

namespace asol {

index_set::index_set(const int num_of_vars,
		             const std::map<int,double>& numeric_constants)
: number_of_variables(num_of_vars), numeric_const(numeric_constants)
{

}

void index_set::finished() {

	if (!current.empty()) {

		push_back_current();
	}
}

void index_set::push_back_current() {

	Set* tmp = new Set;

	for (Map::const_iterator i=current.begin(); i!=current.end(); ++i) {

		const int index = i->first;
		const int count = i->second;

		tmp->insert(index);

		push_back_type3_cse(index, count);
	}

	ASSERT(!tmp->empty());

	constraint_index_sets.push_back(tmp);

	current.clear();
}

void index_set::push_back_type3_cse(const int index, const int count) {

	if (index>=number_of_variables && count>=2) {

		type3_cse.insert(index);
	}
}

int index_set::number_of_constraints() const {

	return static_cast<int> (constraint_index_sets.size());;
}

index_set::~index_set() {

	const int n = number_of_constraints();

	for (int i=0; i<n; ++i) {

		delete constraint_index_sets.at(i);
	}
}

void index_set::record_unary_primitive(int z, int x) {

	std::pair<Map::iterator,bool> res = current.insert(Pair(z, 0));

	ASSERT2(res.second, "index already inserted: "<<res.first->first);

	record_arg(x);
}

bool index_set::is_numeric_constant(const int index) const {

	std::map<int,double>::const_iterator i = numeric_const.find(index);

	return (i!=numeric_const.end()) ? true : false;
}

void index_set::record_arg(const int index) {

	if (is_numeric_constant(index)) {

		return;
	}

	std::pair<Map::iterator,bool> res = current.insert(Pair(index, 0));

	bool new_index = res.second;

	if (!new_index) {

		int& count = res.first->second;

		++count;
	}
}

void index_set::record_binary_primitive(int z, int x, int y) {

	record_unary_primitive(z, x);

	record_arg(y);
}

void index_set::addition(int z, int x, int y) {

	record_binary_primitive(z, x, y);
}

void index_set::substraction(int z, int x, int y) {

	record_binary_primitive(z, x, y);
}

void index_set::multiplication(int z, int x, int y) {

	record_binary_primitive(z, x, y);
}

void index_set::division(int z, int x, int y) {

	record_binary_primitive(z, x, y);
}

void index_set::square(int z, int x) {

	record_unary_primitive(z, x);
}

void index_set::exponential(int z, int x) {

	record_unary_primitive(z, x);
}

void index_set::equality_constraint(int z, int x, double val) {

	push_back_current();
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

	ASSERT2(current.empty(),"recording not finished or not run");
	ASSERT2(type2_cse.empty(),"this function has already been called");

	for (int i=0; i<number_of_constraints(); ++i) {

		check_for_common_subexpressions(i);
	}
}

const std::set<int>& index_set::type2_common_subexpressions() const {

	return type2_cse;
}

const std::set<int>& index_set::type3_common_subexpressions() const {

	return type3_cse;
}

}
