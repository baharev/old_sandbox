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
#include <functional>
#include <ostream>
#include <sstream>
#include "printer.hpp"
#include "primitives.hpp"
#include "diagnostics.hpp"

using std::endl;

namespace asol {

printer::printer(std::ostream& os, const PairVector& num_const)
: out(os), numeric_const(num_const)
{
	dbg_check_if_sorted();
}

int printer::numeric_const_size() const {

	return static_cast<int>(numeric_const.size());
}

void printer::dbg_check_if_sorted() {

	const int n = numeric_const_size();

	if (n<=1) {

		return;
	}

	for (int i=1; i<n; ++i) {

		ASSERT(numeric_const.at(i-1) < numeric_const.at(i));
	}
}

struct compareIndices : std::less<Pair> {

	bool operator() (const Pair& p1, const Pair& p2) {

		return p1.first < p2.first;
	}

} byIndex;

char printer::type(const int index) const {

	typedef PairVector::const_iterator itr;

	Pair p(index, 0.0);

	itr i = std::lower_bound(numeric_const.begin(), numeric_const.end(), p, byIndex);

	char ret_val = 'v';

	if (i!=numeric_const.end() && i->first == index) {

		ret_val = 'n';
	}

	return ret_val;
}

const std::string printer::arg(const int index) const {

	std::ostringstream os;

	os << type(index) << index << std::flush;

	return os.str();
}

void printer::record_unary_primitive(const unary_primitive* p, const char* op) {

	out << 'v' << p->z << " = " << op << '(' << arg(p->x) << ')' << endl;
}

void printer::record_binary_primitive(const binary_primitive* p, const char* op) {

	out << 'v' << p->z << " = " << arg(p->x) << ' ' << op << ' ' << arg(p->y) << endl;
}

void printer::record(const addition* p) {

	record_binary_primitive(p, "+");
}

void printer::record(const substraction* p) {

	record_binary_primitive(p, "-");
}

void printer::record(const multiplication* p) {

	record_binary_primitive(p, "*");
}

void printer::record(const division* p) {

	record_binary_primitive(p, "/");
}

void printer::record(const square* p) {

	record_unary_primitive(p, "sqr");
}

void printer::record(const exponential* p) {

	record_unary_primitive(p, "exp");
}

void printer::record(const equality_constraint* p) {

	out << "Eq: v" << p->z << " = rhs(" << p->x << ")" << endl << endl;
}

}
