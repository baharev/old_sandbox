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

#ifndef INDEX_SET_HPP_
#define INDEX_SET_HPP_

#include <iosfwd>
#include <vector>
#include <map>
#include <set>
#include "recorder.hpp"

namespace asol {

class index_set : public recorder {

public:

	index_set(const int number_of_variables,
			  const std::map<int,double>& numeric_constants);

	void print(std::ostream& out) const;

	void collect_type2_common_subexpressions();

	void collect_variable_set(const std::vector<int>& marked_common_subexpressions);

	const std::set<int>& type2_common_subexpressions() const;

	const std::set<int>& type3_common_subexpressions() const;

	const std::vector<std::vector<int> >& variable_set() const;

	void finished();

	~index_set();

private:

	typedef std::set<int> Set;

	index_set(const index_set& );
	index_set& operator=(const index_set& );

	virtual void addition      (int z, int x, int y);
	virtual void substraction  (int z, int x, int y);
	virtual void multiplication(int z, int x, int y);
	virtual void division      (int z, int x, int y);

	virtual void square     (int z, int x);
	virtual void exponential(int z, int x);

	virtual void equality_constraint(int z, int x, double val);

	virtual void common_subexpression(int z, int x);

	void record_unary_primitive(int z, int x);
	void record_arg(const int index);
	void record_binary_primitive(int z, int x, int y);

	void print_constraint(const int i, std::ostream& out) const;
	const Set non_variables(const int from_constraint_i) const;
	void check_for_common_subexpressions(const int i);

	int number_of_constraints() const;
	bool is_numeric_constant(const int index) const;
	void push_back_current();
	void push_back_type3_cse(const int index, const int count);

	bool not_variable(int index) const;
	void copy_vars(const Set* indices);
	void look_for_cse_mismatch();

	const int number_of_variables;

	const std::map<int,double>& numeric_const;

	std::vector<Set*> constraint_index_sets;

	typedef std::vector<std::vector<int> > IntArray2D;
	typedef std::map<int,int> Map;
	typedef std::pair<int,int> Pair;

	IntArray2D constraint_variable_set;

	Map current;

	Set marked_cse;

	Set type2_cse;

	Set type3_cse;
};

}

#endif // INDEX_SET_HPP_
