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

#ifndef PORT_IMPL_HPP_
#define PORT_IMPL_HPP_

#include <vector>
#include <stdint.h>
#include "lp_impl.hpp"

namespace asol {

class port_impl : public lp_impl {

public:

	port_impl();

private:

	virtual ~port_impl();

	virtual void reset();

	virtual void add_cols(int n);

	// index[1] ... index[length]
	virtual void add_eq_row(const int index[], const double value[], int length, double lb, double ub);

	virtual void set_col_bounds(int index, double lb, double ub);

	virtual void run_simplex();

	virtual double tighten_col_lb(int i, const double lb);

	virtual double tighten_col_ub(int i, const double ub);

	virtual int num_cols() const;

	virtual int num_rows() const;

	virtual col_status col_stat(int i) const;

	virtual double col_val(int i) const;

	virtual double col_lb(int i) const;

	virtual double col_ub(int i) const;

	virtual double col_dual_val(int i) const;

	virtual bool is_fixed(int index) const;

	virtual void dump(const char* file) const;

	virtual void show_iteration_count() const;

	//===================================

	port_impl(const port_impl& );

	port_impl& operator=(const port_impl& );

	void init();

	void set_simple_bounds(int index, double lb, double ub);

	int find_index_position(int index) const;

	void solve_max_x(int index, const double obj_coeff);

	//===================================

	enum ERROR_CODES {
		NO_FEASIBLE_SOLUTION = 8
	};

	// Data chunk passed to the FORTRAN LP solver
	std::vector<double> A; // Columns-major coefficient matrix, A[j*M+i]
	int M;  // Number of rows, only equality constraints are implemented
	int N;  // Number of columns including slack variables <= 2*M
	int IA; // Number of rows of the coefficient matrix A
	std::vector<double> B; // Constraint RHS
	std::vector<double> C; // Objective coefficients
	std::vector<double> X; // Initial estimate / solution vector
	int MAXITR; // Doc recommends 3*N
	double CTX; // Objective, Ctranspose*X
//	int IS;     // Number of variable bounds (_S_imple bounds), 2*N
	std::vector<double> SIMP; // Variable bound, lower if ISIMP(I)>0, upper if ISIMP(I)<0
	std::vector<int>   ISIMP; // Index of var, SIMP(I) is lower boind if ISIMP(I)>0, upper if ISIMP(I)<0
	int IE;   // IE = M, The ﬁrst IE constraints in A are equality constraints
	int IERR; // Error code returned

	int rows_added;
	int slacks_added;

	uint64_t sum_itr_count;
};

}

#endif // PORT_IMPL_HPP_
