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

#ifndef SPLITTING_STRATEGY_HPP_
#define SPLITTING_STRATEGY_HPP_

namespace asol {

class interval;

class splitting_strategy {

public:

	virtual int index_to_split(const interval* box) const = 0;

protected:

	splitting_strategy(const int n_vars) : n_vars(n_vars) { }

	const int n_vars;
};

class max_diam_selector : public splitting_strategy {

public:

	max_diam_selector(const int n_vars) : splitting_strategy(n_vars) { }

private:

	virtual int index_to_split(const interval* box) const;
};

class Jacobsen_x1_D : public splitting_strategy {

public:

	Jacobsen_x1_D(const int n_vars);

private:

	virtual int index_to_split(const interval* box) const;
};

class eco9_sparsity : public splitting_strategy {

public:

	eco9_sparsity(const int n_vars);

private:

	virtual int index_to_split(const interval* box) const;

	int widest_containing_zero(const interval* box) const;
};

}

#endif // SPLITTING_STRATEGY_HPP_
