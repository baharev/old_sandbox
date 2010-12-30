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
#include <set>
#include "recorder.hpp"

namespace asol {

class index_set : public recorder {

public:

	index_set(int max_variable_index);

	void print(std::ostream& out) const;

	void finished();

	~index_set();

private:

	index_set(const index_set& );
	index_set& operator=(const index_set& );

	virtual void record(const addition*            );
	virtual void record(const substraction*        );
	virtual void record(const multiplication*      );
	virtual void record(const division*            );
	virtual void record(const square*              );
	virtual void record(const exponential*         );
	virtual void record(const equality_constraint* );

	void record_primitive(const primitive* p);
	void record_arg(const int index);
	void record_binary_primitive(const binary_primitive* p);

	void print_constraint(const int i, std::ostream& out) const;

	int number_of_constraints() const;

	typedef std::set<int> Set;

	const int max_var_index;

	std::vector<Set*> constraint_index_sets;

	Set* current;
};

}

#endif // INDEX_SET_HPP_
