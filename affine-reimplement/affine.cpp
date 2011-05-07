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

#include <cmath>
#include <ostream>
#include "affine.hpp"
#include "affine_pair_iterator.hpp"
#include "diagnostics.hpp"
#include "exceptions.hpp"
#include "lp_solver.hpp"

namespace asol {

int affine::max_used_index(0);

std::vector<interval>* affine::v(0);

lp_solver* affine::lp(new lp_solver);

const double affine::NARROW(1.0e-6);

void affine::set_vector(std::vector<interval>* vec) {

	v = vec;
}

void affine::reset_counter() {

	max_used_index = 0;
}

void affine::add_noise_var(int index, double coeff) {

	noise_vars.push_back(epsilon(index, coeff));
}

void affine::make_variable() {

	ASSERT(noise_vars.empty());

	const interval rng = range();

	ASSERT2(!rng.degenerate() , "index, range: " << range_index << ", " << rng);

	add_noise_var(               0, rng.midpoint());

	add_noise_var(++max_used_index, rng.radius()  );
}

void affine::recompute_variable(int i) {

	ASSERT(noise_vars.size()==2);

	ASSERT(range_index==i);

	ASSERT(max_used_index==i);

	const interval rng = range();

	ASSERT2(!rng.degenerate() , "index, range: " << range_index << ", " << rng);

	set_var_range(++max_used_index, rng);
}

void affine::set_var_range(int index, const interval& rng) {

	epsilon& e_0 = noise_vars.at(0);

	ASSERT(e_0.index==0);

	e_0.coeff = rng.midpoint();

	epsilon& e_i = noise_vars.at(1);

	ASSERT(e_i.index==index);

	e_i.coeff = rng.radius();
}

void affine::make_numeric_constant() {

	ASSERT(noise_vars.empty());

	const interval rng = range();

	ASSERT2(rng.degenerate() , "index, range: " << range_index << ", " << rng);

	add_noise_var(0, rng.inf());
}

void affine::check_if_numeric_constant() {

	ASSERT(noise_vars.size()==1);

	ASSERT2(range().degenerate() , "index, range: " << range_index << ", " << range());
}

std::ostream& operator<<(std::ostream& os, const affine& x) {

	for (int i=0; i<x.size(); ++i) {

		epsilon e = x.noise_vars.at(i);

		os << e.index << ": " << e.coeff << '\n';
	}

	os << "IA: " << x.range() << "\n\n";

	return os;
}

template <typename T>
void binary_op(const affine& x, const affine& y, T bin_op) {

	affine_pair_iterator itr(x, y);

	bin_op.init_z0(itr.x_i(), itr.y_i());

	while ( itr.increment() ) {

		bin_op.set_zi(itr.index(), itr.x_i(), itr.y_i());
	}

	bin_op.set_z_range(x.range(), y.range());
}

template <typename T>
class binary_operation {

public:

	binary_operation(affine& z) : z(z), rad(0), c(0) { z.noise_vars.clear(); }

	void init_z0(double x0, double y0);

	void set_zi(int index, double x_i, double y_i);

	void set_z_range(const interval& x, const interval& y);

private:

	void intersect_range() {

		const double mid = z.central_value();

		z.force_intersection(mid - rad, mid + rad);
	}

	void process_zi(int index, double tmp) {

		z.add_noise_var(index, tmp);

		rad += std::fabs(tmp);
	}

	affine& z;
	double rad;
	double c; // Only DivP uses this
};

struct Add { };

template <>
inline void binary_operation<Add>::init_z0(double x0, double y0) {

	z.add_noise_var(0, x0+y0);
}

template <>
inline void binary_operation<Add>::set_zi(int index, double x_i, double y_i) {

	process_zi(index, x_i+y_i);
}

template <>
inline void binary_operation<Add>::set_z_range(const interval& x, const interval& y) {

	z.force_intersection(x+y);

	intersect_range();
}

struct Sub { };

template <>
inline void binary_operation<Sub>::init_z0(double x0, double y0) {

	z.add_noise_var(0, x0-y0);
}

template <>
inline void binary_operation<Sub>::set_zi(int index, double x_i, double y_i) {

	process_zi(index, x_i-y_i);
}

template <>
inline void binary_operation<Sub>::set_z_range(const interval& x, const interval& y) {

	z.force_intersection(x-y);

	intersect_range();
}

void aa_addition(affine& z, const affine& x, const affine& y) {

	binary_op(x, y, binary_operation<Add>(z));
}

void aa_substraction(affine& z, const affine& x, const affine& y) {

	binary_op(x, y, binary_operation<Sub>(z));
}

void aa_multiplication(affine& z, const affine& x, const affine& y) {

	z.noise_vars.clear();
	z.noise_vars.reserve(x.size()+y.size()+1);

	affine_pair_iterator itr(x, y);

	const double x0 = itr.x_i();
	const double y0 = itr.y_i();

	double rad(0), rad_x(0), rad_y(0), c(0);

	z.add_noise_var(0, 0);

	while ( itr.increment() ) {

		const double x_i = itr.x_i();
		const double y_i = itr.y_i();

		c += x_i*y_i;

		const double tmp = x0*y_i+y0*x_i;

		z.add_noise_var(itr.index(), tmp);

	    rad   += std::fabs(tmp);
		rad_x += std::fabs(x_i);
		rad_y += std::fabs(y_i);
	}

	c /= 2.0;

	const double mid = x0*y0 + c;

	z.central_value() = mid;

	const double delta = rad_x*rad_y - std::fabs(c);

	ASSERT(delta >= 0.0);

	z.add_noise_var(++affine::max_used_index, delta);

	rad += delta;

	z.force_intersection(mid-rad, mid+rad);
	z.force_intersection(x.range()*y.range());
}

struct DivP { };

template <>
inline void binary_operation<DivP>::init_z0(double x0, double y0) {

	c = x0/y0; // y.central_value() != 0.0 has already been checked

	z.add_noise_var(0, 0);
}

template <>
inline void binary_operation<DivP>::set_zi(int index, double x_i, double y_i) {

	process_zi(index, x_i-c*y_i);
}

template <>
inline void binary_operation<DivP>::set_z_range(const interval& x, const interval& y) {

	z.force_intersection(x-c*y);

	intersect_range();
}

void aa_division(affine& z, const affine& x, const affine& y) {

	ASSERT(y.central_value()!=0.0);

	const int start_index = affine::max_used_index;

	affine P(interval::ANY_REAL());

	binary_op(x, y, binary_operation<DivP>(P));

	affine Q(interval::ANY_REAL());

	aa_reciprocal(Q, y);

	const interval z_old_range = z.range();

	z.get_range() = interval::ANY_REAL();

	aa_multiplication(z, P, Q);

	const double c = x.central_value() / y.central_value();

	z.central_value() += c;

	z.get_range() += c;

	z.force_intersection(z_old_range);

	z.force_intersection(x.range()/y.range());

	const int new_indices = affine::max_used_index - start_index;

	if (new_indices == 2) {

		z.condense_last_two_noise_vars();
	}
}

void affine::condense_last_two_noise_vars() {

	ASSERT(affine::max_used_index == noise_vars.back().index);

	--affine::max_used_index;

	epsilon& before_last = noise_vars.at(size()-2);

	ASSERT(before_last.index == affine::max_used_index);

	ASSERT(std::fabs(before_last.coeff) < 1.0e-12);

	const double last_value = noise_vars.back().coeff;

	ASSERT(last_value >= 0);

	before_last.coeff += last_value;

	noise_vars.pop_back();
}

// TODO Finish
void affine::equals(double value) {

	get_range().force_intersection(value, value);
}

// TODO Finish
void affine::less_than_or_equal_to(affine& rhs) {

	// TODO Not exactly forced
	get_range().less_than_or_equal_to(rhs.get_range());
}

void affine::unary_op(const affine& arg, double alpha, double zeta, double delta) {

	ASSERT2(delta > 0, "delta: " << delta);

	const std::vector<epsilon>& x = arg.noise_vars;

	const int n = arg.size();

	noise_vars.clear();

	noise_vars.reserve(n+1);

	for (int i=0; i<n; ++i) {

		const epsilon& e = x.at(i);

		add_noise_var(e.index, alpha*e.coeff);
	}

	add_noise_var(++affine::max_used_index, delta);

	central_value() += zeta;
}

void aa_exp(affine& z, const affine& x) {

	dbg_consistency(z, x);

	const interval x_rng = x.range();

	ASSERT(!x_rng.degenerate());

	// TODO Either use derivative instead of slope or use only the ia_dag in the search_procedure
	if (x_rng.diameter() < affine::NARROW) {

		throw numerical_problems();
	}

	const double a = x_rng.inf();

	const double b = x_rng.sup();

	const double e_a = std::exp(a);

	const double e_b = std::exp(b);

	const double alpha = (e_b-e_a)/(b-a);

	const double ln_alpha = std::log(alpha);

	const double zeta  = (-alpha*ln_alpha+e_b-alpha*b+alpha)/2.0;

	const double delta = ( alpha*ln_alpha+e_b-alpha*b-alpha)/2.0;

	z.force_intersection(e_a, e_b);

	z.unary_op(x, alpha, zeta, delta);
}

void aa_log(affine& z, const affine& x) {

	dbg_consistency(z, x);

	const interval x_rng = x.range();

	ASSERT(!x_rng.degenerate());

	if (x_rng.diameter() < affine::NARROW) {

		throw numerical_problems();
	}

	const double a = x_rng.inf();

	const double b = x_rng.sup();

	ASSERT2( a > 0, "invalid argument in ln(): " << x_rng);

	const double log_b =  std::log(b);

	const double alpha =  std::log(a/b)/(a-b);

	const double fu    = -std::log(alpha);

	const double ru    = -b*alpha+log_b+1.0;

	const double zeta  =  (fu+ru)/2.0-1.0;

	const double delta =  (fu-ru)/2.0;

	z.force_intersection(std::log(a), log_b);

	z.unary_op(x, alpha, zeta, delta);
}

void aa_sqr(affine& z, const affine& x) {

	dbg_consistency(z, x);

	const interval x_rng = x.range();

	ASSERT(!x_rng.degenerate());

	const double a = x_rng.inf();

	const double b = x_rng.sup();

	const double a_2_p_b_2 = std::pow(a, 2)+std::pow(b, 2);

	const double a_b = a*b;

	const double alpha = a + b;

	const double zeta  =  -(a_2_p_b_2+6.0*a_b)/8.0;

	const double delta =   (a_2_p_b_2-2.0*a_b)/8.0;

	z.force_intersection(sqr(x_rng));

	z.unary_op(x, alpha, zeta, delta);
}

//----------------------------------------------------------------------
//  Reciprocal function using minimax / Chebysev approximation
//  The implementation is based on the paper of L. V. Kolev
//  "An improved interval linearization for solving non-linear problems"
//  Numerical Algorithms, 37, pp.213-224 (2004)

void aa_reciprocal(affine& z, const affine& y) {

	y.dbg_consistency();

	const interval y_rng = y.range();

	ASSERT (!y_rng.degenerate());
	ASSERT2(!y_rng.contains(0),"division by zero: "<<y_rng);

	const double y_inf = y_rng.inf();
	const double y_sup = y_rng.sup();

	const double s  = -1.0/(y_inf*y_sup);

	const double y1 = -std::sqrt(-1.0/s);
	const double y2 = -y1;

	const double ys = (y_inf > 0.0) ? y2 : y1;

	const double f_inf = 1.0/ys    - s*ys   ;
	const double f_sup = 1.0/y_inf - s*y_inf;

	const double f0 = (f_inf + f_sup) / 2.0;

	const double rf = f_sup-f0;

	z.force_intersection(1.0/y_sup, 1.0/y_inf);

	z.unary_op(y, s, f0, rf);
}

void affine::dbg_consistency() const {

	ASSERT(get_range().valid());
    ASSERT(!noise_vars.empty());
	ASSERT(noise_vars.at(0).index==0);
}

void dbg_consistency(const affine& z, const affine& x) {

	ASSERT(z.range_index>=0);

	x.dbg_consistency();
}

void dbg_consistency(const affine& z, const affine& x, const affine& y) {

	dbg_consistency(z, x);

	y.dbg_consistency();
}

affine::~affine() {
	// Out of line dtor to make the compiler shut up
}

void affine::release_all() {

	delete lp;
	lp = 0;

	lp_solver::free_environment();

	v = 0;
}

}
