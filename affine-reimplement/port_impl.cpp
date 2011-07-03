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

#include <algorithm>
#include <iostream>
#include "port_impl.hpp"
#include "diagnostics.hpp"

extern "C" {

// Defined in LINPAS.f
void linpas_(double* A, int* M, int* N, int* IA, double* B, double* C, double* X, int* MAXITR, double* CTX, int* , double* SIMP, int* ISIMP, int* IE, int* ERRCOD);

// Defined at the bottom of the file
void cprint_(double* A, int* M, int* N, int* IA, double* B, double* C, double* X, double* CTX, int* IS, double* SIMP, int* ISIMP, int* IE, int* ITER, int* IPTG, int* IAG, int* IAS, double* U, int* DONE);

}

namespace {

int itr_count; // callback function cprint_ updates this at each iteration
}

namespace asol {

port_impl::port_impl() {

	init();

	sum_itr_count = 0;
}

port_impl::~port_impl() {

}

void port_impl::init() {

	M  = 0;
	N  = 0;
	IA = 0;

	MAXITR = 0;
	CTX    = 0;
	IS     = 0;

	IE     = 0;
	IERR   = 0;

	rows_added   = 0;
	slacks_added = 0;
}

void port_impl::reset() {

	sum_itr_count += itr_count;

	itr_count = 0;

	init();
}

void port_impl::add_cols(int n) {

	// FIXME Assumes n==number of constraints, well defined system of equations
	M = n;
	N = 2*M; // at most, equality constraints won't require slack variables

	A.assign(N*M, 0);
	B.assign(M, 0);
	C.assign(N, 0);
	X.assign(N, 0);

	IS = 2*N;

	SIMP.assign( IS, 0);
	ISIMP.assign(IS, 0);

	for (int i=0; i<n; ++i) {

		set_simple_bounds(i, -1, 1);
	}
}

void port_impl::set_simple_bounds(int i, double lb, double ub) {

	const int pos = 2*i;

	ASSERT2(ISIMP.at(pos)==0 && SIMP.at(pos)==0.0,"bounds are already set, i: "<<i);

	ISIMP.at(pos) =  i+1; // TODO replace with push_back and eliminate slacks_added?
	SIMP.at( pos) =  lb;

	ISIMP.at(pos+1) =-(i+1);
	SIMP.at( pos+1) =  ub;
}

// index[1] ... index[length]
void port_impl::add_eq_row(const int index[], const double value[], int length, double lb, double ub) {

	const int i = rows_added;

	for (int k=1; k<=length; ++k) {

		int j       = index[k] - 1;

		A.at(j*M+i) = value[k];
	}

	if (lb < ub) {

		int j = M+slacks_added;

		A.at(j*M+i) = -1;
		B.at(i) = 0;

		set_simple_bounds(j, lb, ub);

		++slacks_added;
	}
	else {

		B.at(i) = lb; // == ub, equality constraints do not require slacks
	}

	++rows_added;
}

void port_impl::set_col_bounds(int index, double lb, double ub) {
	ASSERT(1<=index && index<=M);
#warning Bounds for noise variables are hard-coded
}

void port_impl::run_simplex() {

// FIXME Only used by search_procedure to make a feasible basis before pruning
}

double port_impl::tighten_col_lb(int i, const double lb) {

	solve_max_x(i, -1);

	const double new_lb = -X.at(i-1);

	return (new_lb > lb) ? new_lb : lb;
}

double port_impl::tighten_col_ub(int i, const double ub) {

	solve_max_x(i,  1);

	const double new_ub =  X.at(i-1);

	return (new_ub < ub) ? new_ub : ub;
}

void port_impl::solve_max_x(int i, const double obj_coeff) {

	ASSERT(i<=M); // C.size()==N

	C.at(i-1) = obj_coeff;

	try {

		// FIXME Call FORTRAN solver
	}
	catch (...) {

		C.at(i-1) = 0.0;

		throw;
	}

	C.at(i-1) = 0.0;
}

int port_impl::num_cols() const {
	// FIXME Eliminate as only used for saving dual costs and
	// that should be changed too to save dual costs of rows
	return M; // Somewhat obscure, would expect N...
}

int port_impl::num_rows() const {

	return M;
}

col_status port_impl::col_stat(int i) const {

	ASSERT(false);
}

double port_impl::col_val(int i) const {

	ASSERT(i<=M); // X.size() == N

	double val = X.at(i-1);

	const double lb = col_lb(i);

	const double ub = col_ub(i);

	ASSERT2(lb<=ub,"lb, ub: "<<lb<<", "<<ub);

	if ((val+1.0e-4 < lb) || (val > ub+1.0e-4)) {

		ASSERT(false);
		//throw numerical_problems();
	}

	if (val < lb) {

		val = lb;
	}
	else if (val > ub) {

		val = ub;
	}

	return val;
}

double port_impl::col_lb(int i) const {

	return SIMP.at(find_index_position( i));
}

double port_impl::col_ub(int i) const {

	return SIMP.at(find_index_position(-i));
}

int port_impl::find_index_position(int i) const {

	std::vector<int>::const_iterator itr = std::find(ISIMP.begin(), ISIMP.end(), i);

	ASSERT2(itr!=ISIMP.end(),"index not found: "<<i);

	return itr - ISIMP.begin();
}

double port_impl::col_dual_val(int i) const {

	ASSERT(false);
}

bool port_impl::is_fixed(int index) const {

	ASSERT(false);
}

void port_impl::dump(const char* file) const {

	ASSERT(false);
}

void port_impl::show_iteration_count() const {

	uint64_t count = sum_itr_count + itr_count;

	std::cout << "Simplex iterations: " << count << std::endl;
}

}

void cprint_(
		double* A,
		int* M,
		int* N,
		int* IA,
		double* B,
		double* C,
		double* X,
		double* CTX,
		int* IS,
		double* SIMP,
		int* ISIMP,
		int* IE,
		int* ITER,
		int* IPTG,
		int* IAG,
		int* IAS,
		double* U,
		int* DONE)
{
	::itr_count = *ITER;
}
