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
#include <iostream>
#include "affine.hpp"
#include "affine_pair_iterator.hpp"
#include "diagnostics.hpp"

namespace asol {

int affine::max_used_index(0);

std::vector<interval>* affine::v(0);

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

	e_i.coeff = rng.diameter();
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

	binary_operation(affine& z) : z(z), rad(0.0) { this->z.noise_vars.clear(); }

	void init_z0(double x0, double y0);

	void set_zi(int index, double x_i, double y_i);

	void set_z_range(const interval& x, const interval& y);

private:

	void intersect_range() {

		const double mid = z.noise_vars.at(0).coeff;

		z.intersect_range(mid - rad, mid + rad);
	}

	void process_zi(int index, double tmp) {

		z.add_noise_var(index, tmp);

		rad += std::fabs(tmp);
	}

	affine& z;
	double rad;
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

	z.intersect_range(x+y);

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

	z.intersect_range(x-y);

	intersect_range();
}

void aa_addition(affine& z, const affine& x, const affine& y) {

	binary_op(x, y, binary_operation<Add>(z));
}

void aa_substraction(affine& z, const affine& x, const affine& y) {

	binary_op(x, y, binary_operation<Sub>(z));
}

void aa_multiplication(affine& z, const affine& x, const affine& y) {

}

void aa_division(affine& z, const affine& x, const affine& y) {

}

void affine::equals(double value) {

}

void affine::less_than_or_equal_to(affine& rhs) {

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

	noise_vars.at(0).coeff += zeta;

	// TODO Intersect with the range of the AA form?
}

void aa_exp(affine& z, const affine& x) {

	dbg_consistency(z, x);

	const interval x_rng = x.range();

	ASSERT(!x_rng.degenerate());

	const double a = x_rng.inf();

	const double b = x_rng.sup();

	const double e_a = std::exp(a);

	const double e_b = std::exp(b);

	const double alpha = (e_b-e_a)/(b-a);

	const double ln_alpha = std::log(alpha);

	const double zeta  = (-alpha*ln_alpha+e_b-alpha*b+alpha)/2.0;

	const double delta = ( alpha*ln_alpha+e_b-alpha*b-alpha)/2.0;

	z.intersect_range(e_a, e_b);

	z.unary_op(x, alpha, zeta, delta);
}

void aa_log(affine& z, const affine& x) {

	dbg_consistency(z, x);

	const interval x_rng = x.range();

	ASSERT(!x_rng.degenerate());

	const double a = x_rng.inf();

	const double b = x_rng.sup();

	ASSERT2( a > 0, "invalid argument in ln(): " << x_rng);

	const double log_b =  std::log(b);

	const double alpha =  std::log(a/b)/(a-b);

	const double fu    = -std::log(alpha);

	const double ru    = -b*alpha+log_b+1.0;

	const double zeta  =  (fu+ru)/2.0-1.0;

	const double delta =  (fu-ru)/2.0;

	z.intersect_range(std::log(a), log_b);

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

	z.intersect_range(sqr(x_rng));

	z.unary_op(x, alpha, zeta, delta);
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

}
