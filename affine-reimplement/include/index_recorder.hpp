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

#ifndef INDEX_RECORDER_HPP_
#define INDEX_RECORDER_HPP_

#include <map>
#include <set>
#include <vector>
#include "recorder.hpp"

namespace asol {

class problem_data;

class index_recorder : public recorder {

public:

	// TODO It should build the problem data
	index_recorder(const problem_data* prob);

	const std::vector<std::vector<int> >& constraint_index_sets() const;

	void dump() const;

private:

	virtual void addition      (int z, int x, int y);
	virtual void substraction  (int z, int x, int y);
	virtual void multiplication(int z, int x, int y);
	virtual void division      (int z, int x, int y);

	virtual void square     (int z, int x);
	virtual void exponential(int z, int x);
	virtual void logarithm  (int z, int x);

	virtual void equality_constraint(int z, int x, double val);
	virtual void common_subexpression(int z, int x);
	virtual void less_than_or_equal_to(int z, int x);

	void record_unary_primitive(int z, int x);
	void record_binary_primitive(int z, int x, int y);
	void record_arg(const int index);

	void resolve_def_var_dependecies();
	const std::set<int> extract_variables(const std::set<int>& index_set);
	void push_back_current();
	void compute_constraint_index_set();
	void merge_up_to(const int last_primitive_index);

	int  primitive_index() const;
	bool is_numeric_constant(const int index) const;
	bool is_defined_variable(const int index) const;
	bool is_variable(const int index) const;

	std::vector<std::set<int> > indices; // indices occurring in constraints AND def vars
	std::vector<int> boundary;           // last primitive index in the con OR def var
	std::map<int,int> def_var_indices;   // def_var index -> indices index to get dependencies

	std::vector<int> constraint_end;     // last primitive index in con ONLY
	std::vector<std::vector<int> > constraint_indices; // indices occurring in con ONLY

	std::map<int,double> numeric_const; // numeric constants
	int n_vars;                         // number of variables
	std::set<int> current;
	int pos;
	int idx;
};

}

#endif // INDEX_RECORDER_HPP_
